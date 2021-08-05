#include "Coil.h"
#include "LegendreMatrix.h"
#include "ctpl.h"

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
           int threadCount) :
        innerRadius(innerRadius), thickness(thickness), length(length), numOfTurns(numOfTurns),
        defaultPrecision(precisionSettings)
{
    setCurrent(current);
    calculateAverageWireThickness();
    setWireResistivity(wireResistivity);
    setSineFrequency(sineFrequency);
    setThreadCount(threadCount);
}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
           double wireResistivity, double sineFrequency, PrecisionFactor precisionFactor,
           int threadCount) :
           innerRadius(innerRadius), thickness(thickness), length(length), numOfTurns(numOfTurns)
{
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