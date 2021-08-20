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
           int threadCount, vec3::CoordVector3 coordinatePosition, double yAxisAngle, double zAxisAngle) :
        innerRadius(innerRadius),
        thickness(thickness / innerRadius < g_thinCoilApproximationRatio ? 0.1 * innerRadius * g_thinCoilApproximationRatio : thickness),
        length(length / innerRadius < g_thinCoilApproximationRatio ? 0.1 * innerRadius * g_thinCoilApproximationRatio : length),
        numOfTurns(numOfTurns),
        defaultPrecision(precisionSettings),
        useFastMethod(length / innerRadius >= g_thinCoilApproximationRatio),
        positionVector(coordinatePosition),
        yAxisAngle(yAxisAngle), zAxisAngle(zAxisAngle)
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
           int threadCount, vec3::CoordVector3 coordinatePosition, double yAxisAngle, double zAxisAngle) :
        innerRadius(innerRadius),
        thickness(thickness / innerRadius < g_thinCoilApproximationRatio ? 0.1 * innerRadius * g_thinCoilApproximationRatio : thickness),
        length(length / innerRadius < g_thinCoilApproximationRatio ? 0.1 * innerRadius * g_thinCoilApproximationRatio : length),
        numOfTurns(numOfTurns),
        useFastMethod(length / innerRadius >= g_thinCoilApproximationRatio),
        positionVector(coordinatePosition),
        yAxisAngle(yAxisAngle), zAxisAngle(zAxisAngle)
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

vec3::CoordVector3 Coil::getPositionVector() const { return positionVector; }

std::pair<double, double> Coil::getRotationAngles() const { return std::make_pair(yAxisAngle, zAxisAngle); }


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


void Coil::setPositionAndOrientation(vec3::CoordVector3 positionVector, double yAxisAngle, double zAxisAngle)
{
    this->positionVector = positionVector;
    this->yAxisAngle = yAxisAngle;
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

    double cosY = std::cos(yAxisAngle); double sinY = std::sin(yAxisAngle);
    double cosZ = std::cos(zAxisAngle); double sinZ = std::sin(zAxisAngle);

    transformationMatrix = vec3::Matrix3(cosZ * cosZ * cosY - sinZ * sinZ, -sinZ * cosZ * cosY - sinZ * cosZ, cosZ * sinY,
                                         sinZ * cosZ * cosY + sinZ * cosZ, cosZ * cosZ - sinZ * sinZ * cosY, sinZ * sinY,
                                         -sinY * cosZ, sinY * sinZ, cosY);

    inverseTransformationMatrix = vec3::Matrix3(cosZ * cosZ * cosY - sinZ * sinZ, sinZ * cosZ * cosY + sinZ * cosZ, -cosZ * sinY,
                                                -sinZ * cosZ * cosY - sinZ * cosZ, cosZ * cosZ - sinZ * sinZ * cosY, sinZ * sinY,
                                                sinY * cosZ, sinY * sinZ, cosY);
}

double Coil::computeAndSetSelfInductance(PrecisionFactor precisionFactor, ComputeMethod method)
{
    if (coilType == CoilType::FILAMENT)
    {
        fprintf(stderr, "ERROR: The integral of a filament is divergent, try a thin rectangular coil\n");
        throw "Coil loop calculation not supported";
    }
    // centering the coil at (0, 0, 0) and setting angles to 0 improves the accuracy by leveraging the z-axis formula
    vec3::CoordVector3 tempPosition = getPositionVector();
    std::pair tempAngles = getRotationAngles();
    setPositionAndOrientation();

    double inductance = Coil::computeMutualInductance(*this, *this, precisionFactor, method);
    setSelfInductance(inductance);
    setPositionAndOrientation(tempPosition, tempAngles.first, tempAngles.second);

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

    double sinA = std::sin(alpha); double cosA = std::cos(alpha);
    double sinB = std::sin(beta); double cosB = std::cos(beta);

    for (int phiBlock = 0; phiBlock <= angularBlocks; ++phiBlock)
    {
        double blockPositionPhi = angularBlock * (phiBlock + 0.5);

        for (int phiIndex = 0; phiIndex <= angularIncrements; ++phiIndex)
        {
            double phi = blockPositionPhi +
                         (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][phiIndex];

            double sinPhi = std::sin(phi); double cosPhi = std::cos(phi);

            ringPosition = vec3::FieldVector3(cosB * cosA * cosPhi - sinB * sinPhi,
                                              sinB * cosA * cosPhi + cosB * sinPhi,
                                              (-1) * sinA * cosPhi);

            ringTangent = vec3::FieldVector3((-1) * cosB * cosA * sinPhi - sinB * cosPhi,
                                             (-1) * sinB * cosA * sinPhi + cosB * cosPhi,
                                             sinA * sinPhi);

            unitRingVector.emplace_back(ringPosition, ringTangent);
        }
    }
    return unitRingVector;
}

bool Coil::isZAxisCase(const Coil &primary, const Coil &secondary)
{
    vec3::FieldVector3 primPositionVec = vec3::CoordVector3::convertToFieldVector(primary.getPositionVector());
    vec3::FieldVector3 secPositionVec = vec3::CoordVector3::convertToFieldVector(secondary.getPositionVector());

    if (primPositionVec.xComponent / primary.innerRadius < g_zAxisApproximationRatio &&
    primPositionVec.yComponent / primary.innerRadius < g_zAxisApproximationRatio &&
    secPositionVec.xComponent / primary.innerRadius < g_zAxisApproximationRatio &&
    secPositionVec.yComponent / primary.innerRadius < g_zAxisApproximationRatio &&
    primary.yAxisAngle / (2 * M_PI) < g_zAxisApproximationRatio &&
    secondary.yAxisAngle / (2 * M_PI) < g_zAxisApproximationRatio)
    {
        return true;
    }
    return false;
}