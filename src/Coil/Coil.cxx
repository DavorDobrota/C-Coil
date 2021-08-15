#include "Coil.h"
#include "LegendreMatrix.h"
#include "ctpl.h"
#include "PrecisionGlobalVars.h"
#include "CoilType.h"

#include <cmath>


namespace
{
    const double g_MiReduced = 0.0000001;

    const double g_defaultCurrent = 1.0;
    const double g_defaultResistivity = 1.63e-8;
    const double g_defaultSineFrequency = 50;
}


Coil::Coil() : Coil(0.0, 0.0, 0.0, 3600, 0) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
           double wireResistivity, double sineFrequency, const PrecisionArguments &precisionSettings,
           int threadCount, vec3::CoordVector3 coordinatePosition, double xAxisAngle, double zAxisAngle) :
        innerRadius(innerRadius),
        thickness(thickness / innerRadius < g_thinCoilApproximationRatio ? 0.1 * innerRadius * g_thinCoilApproximationRatio : thickness),
        length(length / innerRadius < g_thinCoilApproximationRatio ? 0.1 * innerRadius * g_thinCoilApproximationRatio : length),
        numOfTurns(numOfTurns),
        defaultPrecision(precisionSettings),
        useFastMethod(length / innerRadius >= g_thinCoilApproximationRatio),
        positionVector(coordinatePosition),
        xAxisAngle(xAxisAngle), zAxisAngle(zAxisAngle)
{
    calculateTransformationMatrices();
    calculateCoilType();
    setCurrent(current);
    calculateAverageWireThickness();
    setWireResistivity(wireResistivity);
    setSineFrequency(sineFrequency);
    setThreadCount(threadCount);
}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
           double wireResistivity, double sineFrequency, PrecisionFactor precisionFactor,
           int threadCount, vec3::CoordVector3 coordinatePosition, double xAxisAngle, double zAxisAngle) :
        innerRadius(innerRadius),
        thickness(thickness / innerRadius < g_thinCoilApproximationRatio ? 0.1 * innerRadius * g_thinCoilApproximationRatio : thickness),
        length(length / innerRadius < g_thinCoilApproximationRatio ? 0.1 * innerRadius * g_thinCoilApproximationRatio : length),
        numOfTurns(numOfTurns),
        useFastMethod(length / innerRadius >= g_thinCoilApproximationRatio),
        positionVector(coordinatePosition),
        xAxisAngle(xAxisAngle), zAxisAngle(zAxisAngle)
{
    calculateTransformationMatrices();
    calculateCoilType();
    setCurrent(current);
    calculateAverageWireThickness();
    setWireResistivity(wireResistivity);
    setSineFrequency(sineFrequency);
    setPrecisionSettings(PrecisionArguments::getCoilPrecisionArgumentsCPU(*this, precisionFactor));
    setThreadCount(threadCount);
}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current, double sineFrequency,
           PrecisionFactor precisionFactor, int threadCount) :
           Coil(innerRadius, thickness, length, numOfTurns, current, g_defaultResistivity, sineFrequency,
                precisionFactor, threadCount) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current, double sineFrequency,
           const PrecisionArguments &precisionSettings, int threadCount) :
           Coil(innerRadius, thickness, length, numOfTurns, current,
                g_defaultResistivity, sineFrequency, precisionSettings, threadCount) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
           PrecisionFactor precisionFactor, int threadCount)  :
           Coil(innerRadius, thickness, length, numOfTurns, current, g_defaultResistivity,
                g_defaultSineFrequency, precisionFactor, threadCount) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
           const PrecisionArguments &precisionSettings, int threadCount) :
           Coil(innerRadius, thickness, length, numOfTurns, current, g_defaultResistivity,
                g_defaultSineFrequency, precisionSettings, threadCount) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, PrecisionFactor precisionFactor,
           int threadCount) :
           Coil(innerRadius, thickness, length, numOfTurns, g_defaultCurrent, g_defaultResistivity,
                g_defaultSineFrequency, precisionFactor, threadCount){}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns,
           const PrecisionArguments &precisionSettings, int threadCount) :
           Coil(innerRadius, thickness, length, numOfTurns, g_defaultCurrent, g_defaultResistivity,
                g_defaultSineFrequency, precisionSettings, threadCount) {}


double Coil::getCurrentDensity() const { return currentDensity; }

double Coil::getCurrent() const { return current; }

int Coil::getNumOfTurns() const { return numOfTurns; }

double Coil::getInnerRadius() const { return innerRadius; }

double Coil::getThickness() const { return thickness; }

double Coil::getLength() const { return length; }

double Coil::getAverageWireThickness() const { return averageWireThickness; }

bool Coil::isSineDriven1() const { return isSineDriven; }

double Coil::getSineFrequency() const { return sineFrequency; }

double Coil::getMagneticMoment()
{
    calculateMagneticMoment();
    return magneticMoment;
}

double Coil::getWireResistivity() const { return wireResistivity; }

double Coil::getSelfInductance() const { return selfInductance; }

double Coil::getResistance()
{
    calculateResistance();
    return resistance;
}

double Coil::getReactance()
{
    calculateReactance();
    return reactance;
}

double Coil::getImpedance()
{
    calculateImpedance();
    return impedance;
}

const PrecisionArguments &Coil::getPrecisionSettings() const { return defaultPrecision; }

int Coil::getThreadCount() const { return threadCount; }

bool Coil::isUsingFastMethod() const { return useFastMethod; }

CoilType Coil::getCoilType() const { return coilType; }


void Coil::setCurrentDensity(double currentDensity)
{
    this->currentDensity = currentDensity;
    current = currentDensity * length * thickness / numOfTurns;
}

void Coil::setCurrent(double current)
{
    this->current = current;
    currentDensity = current * numOfTurns / (length * thickness);
}

void Coil::setWireResistivity(double wireResistivity)
{
    this->wireResistivity = wireResistivity;
}

void Coil::setSineFrequency(double sineFrequency)
{
    if (sineFrequency > 0.0)
    {
        isSineDriven = true;
        this->sineFrequency = sineFrequency;
    }
    else
    {
        isSineDriven = false;
        this->sineFrequency = 0.0;
    }
    calculateImpedance();
}

void Coil::setPrecisionSettings(const PrecisionArguments &precisionSettings)
{
    this->defaultPrecision = precisionSettings;
}

void Coil::setSelfInductance(double selfInductance)
{
    this->selfInductance = selfInductance;
}


void Coil::setPositionAndOrientation(vec3::CoordVector3 positionVector, double xAxisAngle, double zAxisAngle)
{
    this->positionVector = positionVector;
    this->xAxisAngle = xAxisAngle;
    this->zAxisAngle = zAxisAngle;

    calculateTransformationMatrices();
}

void Coil::calculateMagneticMoment()
{
    magneticMoment = M_PI * current * numOfTurns *
            (pow(innerRadius, 2) + pow(thickness, 2) + innerRadius * thickness / 3);
}

void Coil::calculateAverageWireThickness()
{
    averageWireThickness = std::sqrt(length * thickness / numOfTurns);
}

void Coil::calculateResistance()
{
    double wireRadius = averageWireThickness * 0.5;
    double ohmicResistance = wireResistivity * numOfTurns * 2*M_PI *
            (innerRadius + thickness * 0.5) / (wireRadius * wireRadius * M_PI);
    double skinDepth = std::sqrt(wireResistivity / (M_PI * sineFrequency * g_MiReduced));

    double ohmicSurface = M_PI * wireRadius * wireRadius;
    double sineSurface = 2*M_PI * (
            skinDepth * skinDepth * (exp(-wireRadius / skinDepth) - 1) +
            skinDepth * wireRadius);

    resistance = ohmicResistance * (ohmicSurface / sineSurface);
}

void Coil::calculateReactance()
{
    reactance = selfInductance * 2*M_PI * sineFrequency;
}

void Coil::calculateImpedance()
{
    calculateResistance();
    calculateReactance();
    impedance = std::sqrt(resistance * resistance + reactance * reactance);
}

void Coil::calculateCoilType()
{
    if (thickness / innerRadius < g_thinCoilApproximationRatio && length / innerRadius < g_thinCoilApproximationRatio)
        coilType = CoilType::FILAMENT;

    else if (thickness / innerRadius < g_thinCoilApproximationRatio)
        coilType = CoilType::THIN;

    else if (length / innerRadius < g_thinCoilApproximationRatio)
        coilType = CoilType::FLAT;

    else
        coilType = CoilType::RECTANGULAR;
}

void Coil::calculateTransformationMatrices()
{

    double normalVectorX = std::sin(zAxisAngle) * sin(xAxisAngle);
    double normalVectorY = (-1) * std::sin(xAxisAngle) * std::cos(zAxisAngle);

    double rotationX = std::acos(normalVectorX) - M_PI_2;
    double rotationY = std::acos(normalVectorY) - M_PI_2;
    double rotationZ = std::cos(xAxisAngle);

    double cosX = std::cos(rotationX); double sinX = std::sin(rotationX);
    double cosY = std::cos(rotationY); double sinY = std::sin(rotationY);
    double cosZ = std::cos(rotationZ); double sinZ = std::sin(rotationZ);



    transformationMatrix = vec3::Matrix3(cosZ * cosY, cosZ * sinY * sinX - sinZ * cosX, cosZ * sinY * cosX + sinZ * sinX,
                                         sinZ * cosY, sinZ * sinY * sinX + cosZ * cosX, sinZ * sinY * cosX - cosZ * sinX,
                                         (-1) * sinY, cosY * sinX, cosY * cosY);

    inverseTransformationMatrix = vec3::Matrix3(cosZ * cosY, cosZ * sinY * sinX + sinZ * cosX, (-1) * cosZ * sinY * cosX + sinZ * sinX,
                                                (-1) * sinZ * cosY, (-1) * sinZ * sinY * sinX + cosZ * cosX, sinZ * sinY * cosX + cosZ * sinX,
                                                sinY,  (-1) * cosY * sinX, cosY * cosY);
}

double Coil::computeAndSetSelfInductance(PrecisionFactor precisionFactor, ComputeMethod method)
{
    if (coilType == CoilType::FILAMENT)
    {
        fprintf(stderr, "ERROR: The integral of a filament is divergent, try a thin rectangular coil\n");
        throw "Coil loop calculation not supported";
    }

    double inductance = Coil::computeMutualInductance(*this, *this, 0.0, precisionFactor, method);
    setSelfInductance(inductance);

    return inductance;
}


std::vector<std::pair<vec3::FieldVector3, vec3::FieldVector3>>
Coil::calculateRingIncrementPosition(int angularBlocks, int angularIncrements,
                                     double alpha, double beta, double ringIntervalSize)
{
    int numElements = angularBlocks * angularIncrements;

    std::vector<std::pair<vec3::FieldVector3, vec3::FieldVector3>> unitRingVector;
    unitRingVector.reserve(numElements);

    vec3::FieldVector3 ringPosition;
    vec3::FieldVector3 ringTangent;

    double angularBlock = ringIntervalSize / angularBlocks;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    angularBlocks--;
    angularIncrements--;

    for (int phiBlock = 0; phiBlock <= angularBlocks; ++phiBlock)
    {
        double blockPositionPhi = angularBlock * (phiBlock + 0.5);

        for (int phiIndex = 0; phiIndex <= angularIncrements; ++phiIndex)
        {
            // M_PI/2 added to readjust to an even interval so a shortcut can be used
            double phi = M_PI/2 + blockPositionPhi +
                         (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][phiIndex];

            ringPosition = vec3::FieldVector3(cos(beta) * cos(phi) - sin(beta) * cos(alpha) * sin(phi),
                                              sin(beta) * cos(phi) + cos(beta) * cos(alpha) * sin(phi),
                                              sin(alpha) * sin(phi));

            ringTangent = vec3::FieldVector3((-1) * cos(beta) * sin(phi) - sin(beta) * cos(alpha) * cos(phi),
                                             (-1) * sin(beta) * sin(phi) + cos(beta) * cos(alpha) * cos(phi),
                                             sin(alpha) * cos(phi));

            unitRingVector.emplace_back(ringPosition, ringTangent);
        }
    }
    return unitRingVector;
}