#include <cstdio>
#include <cmath>
#include <vector>
#include <functional>
#include <chrono>

#include "Coil.h"

namespace
{
    const double PI = 3.14159265357989323;
    const double g_MiReduced = 0.0000001;

    const int g_maxLegendrePol = 20;

    const double g_defaultCurrent = 1.0;
    const double g_defaultResistivity = 1.63e-8;
    const double g_defaultSineFrequency = 50;
    const PrecisionArguments g_defaultPrecision = PrecisionArguments(3, 2, 2, 12, 8, 8);

}

namespace Precision
{
    const PrecisionArguments defaultPrecision_ULTRAFAST = PrecisionArguments(1, 1, 1, 12, 8, 8);
    const PrecisionArguments defaultPrecision_FAST = PrecisionArguments(2, 1, 1, 12, 8, 8);
    const PrecisionArguments defaultPrecision_NORMAL = PrecisionArguments(3, 1, 1, 12, 12, 12);
    const PrecisionArguments defaultPrecision_PRECISE = PrecisionArguments(4, 2, 2, 12, 8, 8);
}


PrecisionArguments::PrecisionArguments(
        int numOfAngularBlocks, int numOfThicknessBlocks, int numOfLengthBlocks,
        int numOfAngularIncrements, int numOfThicknessIncrements, int numOfLengthIncrements) :
        numOfAngularBlocks(numOfAngularBlocks), numOfThicknessBlocks(numOfThicknessBlocks),
        numOfLengthBlocks(numOfLengthBlocks), numOfAngularIncrements(numOfAngularIncrements),
        numOfThicknessIncrements(numOfThicknessIncrements), numOfLengthIncrements(numOfLengthIncrements)
{
    if (numOfAngularIncrements > g_maxLegendrePol)
    {
        PrecisionArguments::numOfAngularIncrements = 12;
    }
    if(numOfThicknessIncrements > g_maxLegendrePol)
    {
        PrecisionArguments::numOfThicknessIncrements = 12;
    }
    if(numOfLengthIncrements > g_maxLegendrePol)
    {
        PrecisionArguments::numOfLengthIncrements = 12;
    }

    precisionFactor = 0.0;
    genPrecisionVectors();
}

PrecisionArguments::PrecisionArguments(double precisionFactor)
{
    genParametersFromPrecision();
    genPrecisionVectors();
}

void PrecisionArguments::genPrecisionVectors()
{
    //TO-DO
}

void PrecisionArguments::genParametersFromPrecision()
{
    //TO-DO
}


Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
           double wireResistivity, double sineFrequency, const PrecisionArguments &precisionSettings) :
           innerRadius(innerRadius), thickness(thickness), length(length), numOfTurns(numOfTurns),
           precisionSettings(precisionSettings)
{
    setCurrent(current);
    setWireResistivity(wireResistivity);
    setSineFrequency(sineFrequency);
}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current, double sineFrequency) :
        Coil(innerRadius, thickness, length, numOfTurns, current,
             g_defaultResistivity, sineFrequency, g_defaultPrecision) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current) :
        Coil(innerRadius, thickness, length, numOfTurns, current,
             g_defaultResistivity, g_defaultSineFrequency, g_defaultPrecision) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns) :
        Coil(innerRadius, thickness, length, numOfTurns,
             g_defaultCurrent, g_defaultResistivity, g_defaultSineFrequency, g_defaultPrecision) {}


double Coil::getCurrentDensity() const
{
    return currentDensity;
}

double Coil::getCurrent() const
{
    return current;
}

int Coil::getNumOfTurns() const
{
    return numOfTurns;
}

double Coil::getInnerRadius() const
{
    return innerRadius;
}

double Coil::getThickness() const
{
    return thickness;
}

double Coil::getLength() const
{
    return length;
}

double Coil::getAverageWireThickness() const
{
    return averageWireThickness;
}

bool Coil::isSineDriven1() const
{
    return isSineDriven;
}

double Coil::getSineFrequency() const
{
    return sineFrequency;
}

double Coil::getSelfInductance() const
{
    return selfInductance;
}

double Coil::getMagneticMoment() const
{
    return magneticMoment;
}

double Coil::getWireResistivity() const
{
    return wireResistivity;
}

double Coil::getResistance() const
{
    return resistance;
}

double Coil::getReactance() const
{
    return reactance;
}

double Coil::getImpedance() const
{
    return impedance;
}

const PrecisionArguments &Coil::getPrecisionSettings() const
{
    return precisionSettings;
}


void Coil::setCurrentDensity(double currentDensity)
{
    Coil::currentDensity = currentDensity;
    current = currentDensity * length * thickness / numOfTurns;
    calculateMagneticMoment();
}

void Coil::setCurrent(double current)
{
    Coil::current = current;
    currentDensity = current * numOfTurns / (length * thickness);
    calculateMagneticMoment();
}

void Coil::setWireResistivity(double wireResistivity)
{
    Coil::wireResistivity = wireResistivity;
    calculateImpedance();
}

void Coil::setSineFrequency(double sineFrequency)
{
    if (sineFrequency > 0.0)
    {
        isSineDriven = true;
        Coil::sineFrequency = sineFrequency;
    }
    else
    {
        isSineDriven = false;
        Coil::sineFrequency = 0.0;
    }
    calculateImpedance();
}

void Coil::setPrecisionSettings(const PrecisionArguments &precisionSettings)
{
    Coil::precisionSettings = precisionSettings;
}

void Coil::calculateMagneticMoment()
{
    magneticMoment = PI * current * numOfTurns *
            (pow(innerRadius, 2) + pow(thickness, 2) + innerRadius * thickness);
}

void Coil::calculateAverageWireThickness()
{
    averageWireThickness = sqrt(length * thickness / numOfTurns);
}

void Coil::calculateResistance()
{
    double wireRadius = averageWireThickness / 2;
    double ohmicResistance = numOfTurns * (innerRadius - thickness / 2);
    double skinDepth = sqrt(wireResistivity / (PI * sineFrequency * g_MiReduced));

    double ohmicSurface = PI * pow(wireRadius, 2);
    double sineSurface = 2*PI * (
            pow(skinDepth, 2) * (exp(-wireRadius / skinDepth) - 1) +
            skinDepth * wireRadius);

    resistance = ohmicResistance * (ohmicSurface / sineSurface);
}

void Coil::calculateReactance()
{
    reactance = selfInductance * 2*PI * sineFrequency;
}

void Coil::calculateImpedance()
{
    calculateResistance();
    calculateReactance();
    impedance = sqrt(pow(resistance, 2) + pow(reactance, 2));
}

void Coil::calculateSelfInductance()
{
    //TODO - complicated task, yet unresolved
    selfInductance = 0.0;
}
