
#include <cstdio>
#include <cmath>
#include <vector>
#include <functional>
#include <chrono>

#include "Coil.h"
#include "Polynomial.h"
#include "ComputeMethod.h"
#include "hardware_acceleration.h"
#include "LegendreMatrix.h"

namespace
{
    const double PI = 3.14159265357989323;
    const double g_MiReduced = 0.0000001;

    const int g_defaultLegendreOrder = 12;
    const int g_defaultBlockCount = 1;

    const double g_defaultCurrent = 1.0;
    const double g_defaultResistivity = 1.63e-8;
    const double g_defaultSineFrequency = 50;
    const PrecisionArguments g_defaultPrecision = PrecisionArguments(4, 1, 1, 14, 14, 14);

    const double g_minPrecisionFactor = 1.0;
    const double g_maxPrecisionFactor = 7.0;

    const int g_minPrimLinearIncrements = 12;
    const int g_minPrimAngularIncrements = 12;
    const int g_minPrimAngularBlocks = 1;
    const int g_minSecLinearIncrements = 8;
    const int g_minSecAngularIncrements = 8;
    const int g_minSecAngularBlocks = 1;

}

namespace Precision
{
    const PrecisionArguments defaultPrecision_ULTRAFAST = PrecisionArguments(1, 1, 1, 12, 8, 8);
    const PrecisionArguments defaultPrecision_FAST = PrecisionArguments(1, 1, 1, 12, 12, 12);
    const PrecisionArguments defaultPrecision_NORMAL = PrecisionArguments(2, 1, 1, 12, 12, 12);
    const PrecisionArguments defaultPrecision_PRECISE = PrecisionArguments(1, 1, 1, 48, 16, 16);
}

PrecisionArguments::PrecisionArguments() : PrecisionArguments(g_defaultPrecision) {}

PrecisionArguments::PrecisionArguments(
        int numOfAngularBlocks, int numOfThicknessBlocks, int numOfLengthBlocks,
        int numOfAngularIncrements, int numOfThicknessIncrements, int numOfLengthIncrements) :
        angularBlockCount(numOfAngularBlocks), thicknessBlockCount(numOfThicknessBlocks),
        lengthBlockCount(numOfLengthBlocks), angularIncrementCount(numOfAngularIncrements),
        thicknessIncrementCount(numOfThicknessIncrements), lengthIncrementCount(numOfLengthIncrements)
{
    //TODO - fix constructor calls from main
    if (numOfAngularIncrements > Legendre::maxLegendreOrder || numOfAngularIncrements < 1)
    {
        PrecisionArguments::angularIncrementCount = g_defaultLegendreOrder;
    }
    if (numOfThicknessIncrements >= Legendre::maxLegendreOrder || numOfThicknessIncrements < 1)
    {
        PrecisionArguments::thicknessIncrementCount = g_defaultLegendreOrder;
    }
    if (numOfLengthIncrements >= Legendre::maxLegendreOrder || numOfLengthIncrements < 1)
    {
        PrecisionArguments::lengthIncrementCount = g_defaultLegendreOrder;
    }

    if (numOfAngularBlocks < 1)
    {
        PrecisionArguments::angularBlockCount = g_defaultBlockCount;
    }
    if (numOfThicknessBlocks < 1)
    {
        PrecisionArguments::thicknessBlockCount = g_defaultBlockCount;
    }
    if (numOfLengthIncrements < 1)
    {
        PrecisionArguments::lengthBlockCount = g_defaultBlockCount;
    }

    precisionFactor = 0.0;
}

const int PrecisionArguments::blockPrecisionCPUArray[precisionArraySize] =
        {
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 9,
        9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14,
        15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20, 21, 21, 21, 22, 22, 22, 23,
        23, 23, 24, 24, 24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 32, 32, 33, 33, 34, 34, 35, 35,
        36, 36, 37, 37, 38, 38, 39, 39, 40, 40, 41, 41, 42, 42, 43, 43, 44, 44, 45, 45, 46, 46, 47, 47, 48, 48, 49,
        49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75,
        76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99
        };

const int PrecisionArguments::incrementPrecisionCPUArray[precisionArraySize] =
        {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
        31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 26, 27, 28, 29, 30, 31, 32,
        33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 34, 35, 36, 37, 38, 39, 40, 41, 42,
        43, 44, 45, 46, 47, 48, 49, 50, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 41, 42, 43, 44, 45, 46,
        47, 48, 49, 50, 42, 43, 44, 45, 46, 47, 48, 49, 50, 43, 44, 45, 46, 47, 48, 49, 50, 44, 45, 46, 47, 48, 49,
        50, 45, 46, 47, 48, 49, 50, 46, 47, 48, 49, 50, 46, 47, 48, 49, 50, 46, 47, 48, 49, 50, 47, 48, 49, 50, 47,
        48, 49, 50, 47, 48, 49, 50, 47, 48, 49, 50, 48, 49, 50, 48, 49, 50, 48, 49, 50, 48, 49, 50, 48, 49, 50, 48,
        49, 50, 48, 49, 50, 48, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49,
        50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50,
        49, 50, 49, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
        50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50
        };

PrecisionArguments::PrecisionArguments(double precisionFactor)
{
    genParametersFromPrecision();
}

void PrecisionArguments::genParametersFromPrecision()
{
    //TODO - when further analysis is complete a method will be devised, until then
    angularBlockCount = 2;
    lengthBlockCount = 1;
    thicknessBlockCount = 1;

    angularIncrementCount = 12;
    thicknessIncrementCount = 12;
    lengthIncrementCount = 12;
}

void PrecisionArguments::getMutualInductancePrecisionSettingsZCPU(const Coil &primCoil,
                                                                  const Coil &secCoil,
                                                                  double precisionFactor,
                                                                  PrecisionArguments &fieldPrecision,
                                                                  int &linearIncrements)
{
    if (precisionFactor < g_minPrecisionFactor)
        precisionFactor = g_minPrecisionFactor;

    if (precisionFactor > g_maxPrecisionFactor)
        precisionFactor = g_maxPrecisionFactor;

    int primLinearIncrements = g_minPrimLinearIncrements;
    int secLinearIncrements = g_minSecLinearIncrements;

    int primAngularArrayIndex = g_minPrimAngularIncrements - 1;

    int totalIncrements = pow(2, 16 + precisionFactor);
    int currentIncrements;

    do
    {
        double primAngularStep =
                PI * (primCoil.getInnerRadius() + primCoil.getThickness() * 0.5) /
                (blockPrecisionCPUArray[primAngularArrayIndex] *
                 incrementPrecisionCPUArray[primAngularArrayIndex]);

        double primLengthThicknessStep = sqrt(2 * primCoil.getThickness() * primCoil.getLength()) /
                                         primLinearIncrements;

        double secLengthThicknessStep = sqrt(secCoil.getThickness() * secCoil.getLength()) /
                                        secLinearIncrements;

        if (primAngularStep / primLengthThicknessStep >= 1.0)
            primAngularArrayIndex++;
        else
        {
            if (primLengthThicknessStep / secLengthThicknessStep >= 1.0)
                primLinearIncrements++;
            else
                secLinearIncrements++;
        }
        currentIncrements =
                primLinearIncrements * primLinearIncrements * secLinearIncrements * secLinearIncrements *
                blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex];
    }
    while (currentIncrements < totalIncrements);

    fieldPrecision = PrecisionArguments(blockPrecisionCPUArray[primAngularArrayIndex], 1, 1,
                                        incrementPrecisionCPUArray[primAngularArrayIndex],
                                        primLinearIncrements, primLinearIncrements);
    linearIncrements = secLinearIncrements;

}


Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
           double wireResistivity, double sineFrequency, const PrecisionArguments &precisionSettings) :
           innerRadius(innerRadius), thickness(thickness), length(length), numOfTurns(numOfTurns),
           precisionSettings(precisionSettings)
{

    setCurrent(current);
    calculateAverageWireThickness();
    setWireResistivity(wireResistivity);
    setSineFrequency(sineFrequency);
    calculateSelfInductance();
}

Coil::Coil() : Coil(0.0, 0.0, 0.0, 3600, 0,
                    g_defaultResistivity, g_defaultSineFrequency, g_defaultPrecision){}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current, double sineFrequency) :
        Coil(innerRadius, thickness, length, numOfTurns, current,
             g_defaultResistivity, sineFrequency, g_defaultPrecision) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current, double sineFrequency,
           const PrecisionArguments &precisionSettings) :
           Coil(innerRadius, thickness, length, numOfTurns, current,
                g_defaultResistivity, sineFrequency, precisionSettings) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current) :
        Coil(innerRadius, thickness, length, numOfTurns, current,
             g_defaultResistivity, g_defaultSineFrequency, g_defaultPrecision) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
           const PrecisionArguments &precisionSettings) :
           Coil(innerRadius, thickness, length, numOfTurns, current,
                g_defaultResistivity, g_defaultSineFrequency, precisionSettings) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns) :
        Coil(innerRadius, thickness, length, numOfTurns,
             g_defaultCurrent, g_defaultResistivity, g_defaultSineFrequency, g_defaultPrecision) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns,
           const PrecisionArguments &precisionSettings) :
           Coil(innerRadius, thickness, length, numOfTurns, g_defaultCurrent,
                g_defaultResistivity, g_defaultSineFrequency, precisionSettings) {}


double Coil::getCurrentDensity() const { return currentDensity; }

double Coil::getCurrent() const { return current; }

int Coil::getNumOfTurns() const { return numOfTurns; }

double Coil::getInnerRadius() const { return innerRadius; }

double Coil::getThickness() const { return thickness; }

double Coil::getLength() const { return length; }

double Coil::getAverageWireThickness() const { return averageWireThickness; }

bool Coil::isSineDriven1() const { return isSineDriven; }

double Coil::getSineFrequency() const { return sineFrequency; }

double Coil::getSelfInductance() const { return selfInductance; }

double Coil::getMagneticMoment() const { return magneticMoment; }

double Coil::getWireResistivity() const { return wireResistivity; }

double Coil::getResistance() const { return resistance; }

double Coil::getReactance() const { return reactance; }

double Coil::getImpedance() const { return impedance; }

const PrecisionArguments &Coil::getPrecisionSettings() const { return precisionSettings; }


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
            (pow(innerRadius, 2) + pow(thickness, 2) + innerRadius * thickness / 3);
}

void Coil::calculateAverageWireThickness()
{
    averageWireThickness = sqrt(length * thickness / numOfTurns);
}

void Coil::calculateResistance()
{
    double wireRadius = averageWireThickness * 0.5;
    double ohmicResistance = wireResistivity * numOfTurns * 2*PI *
            (innerRadius + thickness * 0.5) / (wireRadius * wireRadius * PI);
    double skinDepth = sqrt(wireResistivity / (PI * sineFrequency * g_MiReduced));

    double ohmicSurface = PI * wireRadius * wireRadius;
    double sineSurface = 2*PI * (
            skinDepth * skinDepth * (exp(-wireRadius / skinDepth) - 1) +
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
    impedance = sqrt(resistance * resistance + reactance * reactance);
}

void Coil::calculateSelfInductance()
{
    //TODO - complicated task, yet unresolved
    selfInductance = 0.0;
}

std::pair<double, double> Coil::calculateBField(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    double magneticFieldZ = 0.0;
    double magneticFieldH = 0.0;
    
    double lengthBlock = length / usedPrecision.lengthBlockCount;
    double thicknessBlock = thickness / usedPrecision.thicknessBlockCount;
    double angularBlock = PI / usedPrecision.angularBlockCount;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int lengthIncrements = usedPrecision.lengthIncrementCount - 1;
    int thicknessIncrements = usedPrecision.thicknessIncrementCount - 1;
    int angularIncrements = usedPrecision.angularIncrementCount - 1;

    // multiplication by 2 because cosine is an even function and by 0.125 for a triple change of interval (3 times 1/2)
    double constant = g_MiReduced * currentDensity * lengthBlock * thicknessBlock * angularBlock * 2 * 0.125;

    for (int indBlockL = 0; indBlockL < usedPrecision.lengthBlockCount; ++indBlockL)
    {
        for (int indBlockT = 0; indBlockT < usedPrecision.thicknessBlockCount; ++indBlockT)
        {
            double blockPositionL = (-1) * (length * 0.5) + lengthBlock * (indBlockL + 0.5);
            double blockPositionT = innerRadius + thicknessBlock * (indBlockT + 0.5);

            for (int incL = 0; incL <= lengthIncrements; ++incL)
            {
                for (int incT = 0; incT <= thicknessIncrements; ++incT)
                {
                    double incrementPositionL = blockPositionL +
                            (lengthBlock * 0.5) * Legendre::positionMatrix[lengthIncrements][incL];
                    double incrementPositionT = blockPositionT +
                            (thicknessBlock * 0.5) * Legendre::positionMatrix[thicknessIncrements][incT];

                    double incrementWeightS =
                            Legendre::weightsMatrix[lengthIncrements][incL] *
                            Legendre::weightsMatrix[thicknessIncrements][incT];

                    double tempConstA = incrementPositionT * incrementPositionT;
                    double tempConstB = incrementPositionT * (incrementPositionL + zAxis);
                    double tempConstC = incrementPositionT * rPolar;
                    double tempConstD = tempConstA + rPolar * rPolar +
                                        (incrementPositionL + zAxis) * (incrementPositionL + zAxis);
                    double tempConstE = constant * incrementWeightS;

                    for (int indBlockFi = 0; indBlockFi < usedPrecision.angularBlockCount; ++indBlockFi)
                    {
                        double blockPositionFi = angularBlock * (indBlockFi + 0.5);

                        for (int incFi = 0; incFi <= angularIncrements; ++incFi)
                        {
                            double incrementPositionFi = blockPositionFi +
                                    (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][incFi];

                            double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incFi];

                            double cosinePhi = cos(incrementPositionFi);
                            double tempConstF = 2 * tempConstC * cosinePhi;
                            double tempConstH = (tempConstD - tempConstF) * sqrt(tempConstD - tempConstF);
                            double tempConstG = tempConstE * incrementWeightFi / tempConstH;

                            magneticFieldZ += tempConstG * (tempConstA - tempConstC * cosinePhi);
                            magneticFieldH += tempConstG * (tempConstB * cosinePhi);
                        }
                    }
                }
            }
        }
    }
    std::pair<double, double> output;
    output.first = magneticFieldH;
    output.second = magneticFieldZ;

    return output;
}

double Coil::calculateBFieldVertical(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    double magneticFieldZ = 0.0;
    
    double lengthBlock = length / usedPrecision.lengthBlockCount;
    double thicknessBlock = thickness / usedPrecision.thicknessBlockCount;
    double angularBlock = PI / usedPrecision.angularBlockCount;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int lengthIncrements = usedPrecision.lengthIncrementCount - 1;
    int thicknessIncrements = usedPrecision.thicknessIncrementCount - 1;
    int angularIncrements = usedPrecision.angularIncrementCount - 1;

    // multiplication by 2 because cosine is an even function and by 0.125 for a triple change of interval (3 times 1/2)
    double constant = g_MiReduced * currentDensity * lengthBlock * thicknessBlock * angularBlock * 2 * 0.125;

    for (int indBlockL = 0; indBlockL < usedPrecision.lengthBlockCount; ++indBlockL)
    {
        for (int indBlockT = 0; indBlockT < usedPrecision.thicknessBlockCount; ++indBlockT)
        {
            double blockPositionL = (-1) * (length * 0.5) + lengthBlock * (indBlockL + 0.5);
            double blockPositionT = innerRadius + thicknessBlock * (indBlockT + 0.5);

            for (int incL = 0; incL <= lengthIncrements; ++incL)
            {
                for (int incT = 0; incT <= thicknessIncrements; ++incT)
                {
                    double incrementPositionL = blockPositionL +
                            (lengthBlock * 0.5) * Legendre::positionMatrix[lengthIncrements][incL];
                    double incrementPositionT = blockPositionT +
                            (thicknessBlock * 0.5) * Legendre::positionMatrix[thicknessIncrements][incT];

                    double incrementWeightS =
                            Legendre::weightsMatrix[lengthIncrements][incL] *
                            Legendre::weightsMatrix[thicknessIncrements][incT];

                    double tempConstA = incrementPositionT * incrementPositionT;
                    double tempConstB = incrementPositionT * rPolar;
                    double tempConstC = tempConstA + rPolar * rPolar +
                            (incrementPositionL + zAxis) * (incrementPositionL + zAxis);

                    for (int indBlockFi = 0; indBlockFi < usedPrecision.angularBlockCount; ++indBlockFi)
                    {
                        double blockPositionFi = angularBlock * (indBlockFi + 0.5);

                        for (int incFi = 0; incFi <= angularIncrements; ++incFi)
                        {
                            double incrementPositionFi = blockPositionFi +
                                    (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][incFi];

                            double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incFi];

                            double cosinePhi = cos(incrementPositionFi);
                            double tempConstD = tempConstC - 2 * tempConstB * cosinePhi;

                            magneticFieldZ += constant * incrementWeightS * incrementWeightFi *
                                    (tempConstA - tempConstB * cosinePhi) / (tempConstD * sqrt(tempConstD));
                        }
                    }
                }
            }
        }
    }
    return magneticFieldZ;
}

double Coil::calculateBFieldHorizontal(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    double magneticFieldH = 0.0;
    
    double lengthBlock = length / usedPrecision.lengthBlockCount;
    double thicknessBlock = thickness / usedPrecision.thicknessBlockCount;
    double angularBlock = PI / usedPrecision.angularBlockCount;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int lengthIncrements = usedPrecision.lengthIncrementCount - 1;
    int thicknessIncrements = usedPrecision.thicknessIncrementCount - 1;
    int angularIncrements = usedPrecision.angularIncrementCount - 1;

    // multiplication by 2 because cosine is an even function and by 0.125 for a triple change of interval (3 times 1/2)
    double constant = g_MiReduced * currentDensity * lengthBlock * thicknessBlock * angularBlock * 2 * 0.125;

    for (int indBlockL = 0; indBlockL < usedPrecision.lengthBlockCount; ++indBlockL)
    {
        for (int indBlockT = 0; indBlockT < usedPrecision.thicknessBlockCount; ++indBlockT)
        {
            double blockPositionL = (-1) * (length * 0.5) + lengthBlock * (indBlockL + 0.5);
            double blockPositionT = innerRadius + thicknessBlock * (indBlockT + 0.5);

            for (int incL = 0; incL <= lengthIncrements; ++incL)
            {
                for (int incT = 0; incT <= thicknessIncrements; ++incT)
                {
                    double incrementPositionL = blockPositionL +
                                                (lengthBlock * 0.5) * Legendre::positionMatrix[lengthIncrements][incL];
                    double incrementPositionT = blockPositionT +
                                                (thicknessBlock * 0.5) * Legendre::positionMatrix[thicknessIncrements][incT];

                    double incrementWeightS =
                            Legendre::weightsMatrix[lengthIncrements][incL] *
                            Legendre::weightsMatrix[thicknessIncrements][incT];

                    double tempConstA = incrementPositionT * (incrementPositionL + zAxis);
                    double tempConstB = incrementPositionT * rPolar;
                    double tempConstC = incrementPositionT * incrementPositionT + rPolar * rPolar +
                                        (incrementPositionL + zAxis) * (incrementPositionL + zAxis);

                    for (int indBlockFi = 0; indBlockFi < usedPrecision.angularBlockCount; ++indBlockFi)
                    {
                        double blockPositionFi = angularBlock * (indBlockFi + 0.5);

                        for (int incFi = 0; incFi <= angularIncrements; ++incFi)
                        {
                            double incrementPositionFi = blockPositionFi +
                                                         (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][incFi];

                            double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incFi];

                            double cosinePhi = cos(incrementPositionFi);
                            double tempConstD = tempConstC - 2 * tempConstB * cosinePhi;

                            magneticFieldH += constant * incrementWeightS * incrementWeightFi *
                                    (tempConstA * cosinePhi) /(tempConstD * sqrt(tempConstD));
                        }
                    }
                }
            }
        }
    }
    return magneticFieldH;
}

double Coil::calculateAPotential(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    double magneticPotential = 0.0;

    double lengthBlock = length / usedPrecision.lengthBlockCount;
    double thicknessBlock = thickness / usedPrecision.thicknessBlockCount;
    double angularBlock = PI / usedPrecision.angularBlockCount;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int lengthIncrements = usedPrecision.lengthIncrementCount - 1;
    int thicknessIncrements = usedPrecision.thicknessIncrementCount - 1;
    int angularIncrements = usedPrecision.angularIncrementCount - 1;

    // multiplication by 2 because cosine is an even function and by 0.125 for a triple change of interval (3 times 1/2)
    double constant = g_MiReduced * currentDensity * lengthBlock * thicknessBlock * angularBlock * 2 * 0.125;

    for (int indBlockL = 0; indBlockL < usedPrecision.lengthBlockCount; ++indBlockL)
    {
        for (int indBlockT = 0; indBlockT < usedPrecision.thicknessBlockCount; ++indBlockT)
        {
            double blockPositionL = (-1) * (length * 0.5) + lengthBlock * (indBlockL + 0.5);
            double blockPositionT = innerRadius + thicknessBlock * (indBlockT + 0.5);

            for (int incL = 0; incL <= lengthIncrements; ++incL)
            {
                for (int incT = 0; incT <= thicknessIncrements; ++incT)
                {
                    double incrementPositionL = blockPositionL +
                                                (lengthBlock * 0.5) * Legendre::positionMatrix[lengthIncrements][incL];
                    double incrementPositionT = blockPositionT +
                                                (thicknessBlock * 0.5) * Legendre::positionMatrix[thicknessIncrements][incT];

                    double incrementWeightS =
                            Legendre::weightsMatrix[lengthIncrements][incL] *
                            Legendre::weightsMatrix[thicknessIncrements][incT];

                    double tempConstA = incrementPositionT;
                    double tempConstB = incrementPositionT * rPolar;
                    double tempConstC = incrementPositionT * incrementPositionT + rPolar * rPolar +
                                        (incrementPositionL + zAxis) * (incrementPositionL + zAxis);

                    for (int indBlockFi = 0; indBlockFi < usedPrecision.angularBlockCount; ++indBlockFi)
                    {
                        double blockPositionFi = angularBlock * (indBlockFi + 0.5);

                        for (int incFi = 0; incFi <= angularIncrements; ++incFi)
                        {
                            double incrementPositionFi = blockPositionFi +
                                                         (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][incFi];

                            double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incFi];

                            double cosinePhi = cos(incrementPositionFi);

                            magneticPotential += constant * incrementWeightS * incrementWeightFi *
                                    (tempConstA * cosinePhi) /sqrt(tempConstC - 2*tempConstB * cosinePhi);
                        }
                    }
                }
            }
        }
    }
    return magneticPotential;
}

double Coil::computeBFieldX(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return calculateBFieldHorizontal(cylindricalZ, cylindricalR, precisionSettings) * cos(cylindricalPhi);
}

double Coil::computeBFieldX(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                            const PrecisionArguments &usedPrecision) const
{
    return calculateBFieldHorizontal(cylindricalZ, cylindricalR, usedPrecision) * cos(cylindricalPhi);
}

double Coil::computeBFieldY(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return calculateBFieldHorizontal(cylindricalZ, cylindricalR, precisionSettings) * sin(cylindricalPhi);
}

double Coil::computeBFieldY(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                            const PrecisionArguments &usedPrecision) const
{
    return calculateBFieldHorizontal(cylindricalZ, cylindricalR, usedPrecision) * sin(cylindricalPhi);
}

double Coil::computeBFieldH(double cylindricalZ, double cylindricalR) const
{
    return calculateBFieldHorizontal(cylindricalZ, cylindricalR, precisionSettings);
}

double Coil::computeBFieldH(double cylindricalZ, double cylindricalR, const PrecisionArguments &usedPrecision) const
{
    return calculateBFieldHorizontal(cylindricalZ, cylindricalR, usedPrecision);
}

double Coil::computeBFieldZ(double cylindricalZ, double cylindricalR) const
{
    return calculateBFieldVertical(cylindricalZ, cylindricalR, precisionSettings);
}

double Coil::computeBFieldZ(double cylindricalZ, double cylindricalR, const PrecisionArguments &usedPrecision) const
{
    return calculateBFieldVertical(cylindricalZ, cylindricalR, usedPrecision);
}



std::vector<double> Coil::computeBFieldVector(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                              const PrecisionArguments &usedPrecision) const
{
    std::vector<double> fieldVector;
    std::pair<double, double> fields = calculateBField(cylindricalZ, cylindricalR, usedPrecision);

    fieldVector.push_back(fields.first * cos(cylindricalPhi));
    fieldVector.push_back(fields.first * sin(cylindricalPhi));
    fieldVector.push_back(fields.second);

    return fieldVector;
}
std::vector<double> Coil::computeBFieldVector(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeBFieldVector(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

double Coil::computeAPotentialX(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return (-1) * calculateAPotential(cylindricalZ, cylindricalR, precisionSettings) * sin(cylindricalPhi);
}

double Coil::computeAPotentialX(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                const PrecisionArguments &usedPrecision) const
{
    return (-1) * calculateAPotential(cylindricalZ, cylindricalR, usedPrecision) * sin(cylindricalPhi);
}

double Coil::computeAPotentialY(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return calculateAPotential(cylindricalZ, cylindricalR, precisionSettings) * cos(cylindricalPhi);
}

double Coil::computeAPotentialY(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                const PrecisionArguments &usedPrecision) const
{
    return calculateAPotential(cylindricalZ, cylindricalR, usedPrecision) * cos(cylindricalPhi);
}

double Coil::computeAPotentialAbs(double cylindricalZ, double cylindricalR) const
{
    return calculateAPotential(cylindricalZ, cylindricalR, precisionSettings);
}

double Coil::computeAPotentialAbs(double cylindricalZ, double cylindricalR, PrecisionArguments &usedPrecision) const
{
    return calculateAPotential(cylindricalZ, cylindricalR, usedPrecision);
}

std::vector<double> Coil::computeAPotentialVector(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                                  const PrecisionArguments &usedPrecision) const
{
    std::vector<double> potentialVector;
    double potential = calculateAPotential(cylindricalZ, cylindricalR, usedPrecision);

    potentialVector.push_back(potential * (-sin(cylindricalPhi)));
    potentialVector.push_back(potential * cos(cylindricalPhi));
    return potentialVector;
}

std::vector<double> Coil::computeAPotentialVector(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeAPotentialVector(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

void Coil::calculateAllBFieldSINGLE(const std::vector<double> &cylindricalZArr,
                                    const std::vector<double> &cylindricalRArr,
                                    std::vector<double> &computedFieldHArr,
                                    std::vector<double> &computedFieldZArr,
                                    const PrecisionArguments &usedPrecision) const
{
    computedFieldHArr.resize(0);
    computedFieldZArr.resize(0);

    for (int i = 0; i < cylindricalZArr.size(); ++i)
    {
        std::pair<double, double> values = calculateBField(cylindricalZArr[i], cylindricalRArr[i], usedPrecision);
        computedFieldHArr.push_back(values.first);
        computedFieldZArr.push_back(values.second);
    }
}

void Coil::calculateAllBFieldVerticalSINGLE(const std::vector<double> &cylindricalZArr,
                                            const std::vector<double> &cylindricalRArr,
                                            std::vector<double> &computedFieldZArr,
                                            const PrecisionArguments &usedPrecision) const
{
    computedFieldZArr.resize(0);

    for (int i = 0; i < cylindricalZArr.size(); ++i)
    {
        double field = calculateBFieldVertical(cylindricalZArr[i], cylindricalRArr[i], usedPrecision);
        computedFieldZArr.push_back(field);
    }
}

void Coil::calculateAllBFieldHorizontalSINGLE(const std::vector<double> &cylindricalZArr,
                                              const std::vector<double> &cylindricalRArr,
                                              std::vector<double> &computedFieldHArr,
                                              const PrecisionArguments &usedPrecision) const
{
    computedFieldHArr.resize(0);

    for (int i = 0; i < cylindricalZArr.size(); ++i)
    {
        double field = calculateBFieldHorizontal(cylindricalZArr[i], cylindricalRArr[i], usedPrecision);
        computedFieldHArr.push_back(field);
    }
}

void Coil::calculateAllAPotentialSINGLE(const std::vector<double> &cylindricalZArr,
                                        const std::vector<double> &cylindricalRArr,
                                        std::vector<double> &computedPotentialArr,
                                        const PrecisionArguments &usedPrecision) const
{
    computedPotentialArr.resize(0);

    for (int i = 0; i < cylindricalZArr.size(); ++i)
    {
        double field = calculateAPotential(cylindricalZArr[i], cylindricalRArr[i], usedPrecision);
        computedPotentialArr.push_back(field);
    }
}

void Coil::convertPolarToCylindrical(double polarR, double polarTheta, double polarPhi,
                                     double &cylindricalZ, double &cylindricalR, double &cylindricalPhi)
{
    cylindricalZ = polarR * cos(polarTheta);
    cylindricalR = polarR * sin(polarTheta);
    cylindricalPhi = cylindricalPhi;
}

void Coil::convertAllPolarToCylindrical(const std::vector<double> &polarRArr,
                                        const std::vector<double> &polarThetaArr,
                                        const std::vector<double> &polarPhiArr,
                                        std::vector<double> &cylindricalZArr,
                                        std::vector<double> &cylindricalRArr,
                                        std::vector<double> &cylindricalPhiArr)
{
    if (polarRArr.size() == polarThetaArr.size() == polarPhiArr.size())
    {
        cylindricalZArr.resize(0);
        cylindricalRArr.resize(0);
        cylindricalPhiArr.resize(0);

        for (int i = 0; i < polarRArr.size(); ++i)
        {
            cylindricalZArr.push_back(polarRArr[i] * cos(polarThetaArr[i]));
            cylindricalRArr.push_back(polarRArr[i] * sin(polarThetaArr[i]));
            cylindricalPhiArr.push_back(polarPhiArr[i]);
        }
    }
}

void Coil::convertCylindricalToPolar(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                     double &polarR, double &polarTheta, double &polarPhi)
{
    polarR = sqrt(cylindricalZ * cylindricalZ + cylindricalR * cylindricalR);
    polarTheta = atan2(cylindricalR, cylindricalZ);
    polarPhi = cylindricalPhi;
}

void Coil::convertAllCylindricalToPolar(const std::vector<double> &cylindricalZArr,
                                        const std::vector<double> &cylindricalRArr,
                                        const std::vector<double> &cylindricalPhiArr,
                                        std::vector<double> &polarRArr,
                                        std::vector<double> &polarThetaArr,
                                        std::vector<double> &polarPhiArr)
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        polarRArr.resize(0);
        polarThetaArr.resize(0);
        polarPhiArr.resize(0);

        for (int i = 0; i < cylindricalZArr.size(); ++i)
        {
            polarRArr.push_back(sqrt(cylindricalZArr[i] * cylindricalZArr[i] + cylindricalRArr[i] * cylindricalRArr[i]));
            polarThetaArr.push_back(atan2(cylindricalRArr[i], cylindricalZArr[i]));
            polarPhiArr.push_back(cylindricalPhiArr[i]);
        }
    }
}

void Coil::calculateAllBFieldACCELERATED(const std::vector<double> &cylindricalZArr,
                                         const std::vector<double> &cylindricalRArr,
                                         std::vector<float> &computedFieldHArr,
                                         std::vector<float> &computedFieldZArr,
                                         const PrecisionArguments &usedPrecision) const
{
    std::vector<float> polarR;
    std::vector<float> polarTheta;

    for (int i = 0; i < cylindricalZArr.size(); ++i)
    {
        polarR.push_back(sqrt(cylindricalZArr[i] * cylindricalZArr[i] + cylindricalRArr[i] * cylindricalRArr[i]));
        polarTheta.push_back(atan2(cylindricalRArr[i], cylindricalZArr[i]));
    }

    computedFieldHArr.resize(polarR.size());
    computedFieldZArr.resize(polarR.size());

    Calculate_hardware_accelerated_b(polarR.size(), &polarTheta[0], &polarR[0],
                                     currentDensity, innerRadius, length, thickness,
                                     thickness/16, length/16, PI/48,
                                     &computedFieldHArr[0], &computedFieldZArr[0]);
}

void Coil::calculateAllAPotentialACCELERATED(const std::vector<double> &cylindricalZArr,
                                             const std::vector<double> &cylindricalRArr,
                                             std::vector<float> &computedPotentialArr,
                                             const PrecisionArguments &usedPrecision) const
{
    std::vector<float> polarR;
    std::vector<float> polarTheta;

    for (int i = 0; i < cylindricalZArr.size(); ++i)
    {
        polarR.push_back(sqrt(cylindricalZArr[i] * cylindricalZArr[i] + cylindricalRArr[i] * cylindricalRArr[i]));
        polarTheta.push_back(atan2(cylindricalRArr[i], cylindricalZArr[i]));
    }
    computedPotentialArr.resize(polarR.size());

    Calculate_hardware_accelerated_a(polarR.size(), &polarTheta[0], &polarR[0],
                                     currentDensity, innerRadius, length, thickness,
                                     thickness / 16, length / 16, PI / 48,
                                     nullptr, nullptr, &computedPotentialArr[0]);

    // TODO - fix frequency in GPU potential calculation, current temporary fix
    for (int i = 0; i < polarR.size(); ++i)
        computedPotentialArr[i] /= 2*PI;
}

void Coil::computeAllBFieldX(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             const PrecisionArguments &usedPrecision,
                             ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedFieldArr.resize(0);

        if (method == CPU_ST)
        {
            calculateAllBFieldHorizontalSINGLE(
                    cylindricalZArr, cylindricalRArr, computedFieldArr, usedPrecision);

            for (int i = 0; i < computedFieldArr.size(); ++i)
            {
                computedFieldArr[i] *= cos(cylindricalPhiArr[i]);
            }
        }
        else if (method == GPU)
        {
            std::vector<float> fieldH;
            std::vector<float> fieldZ;

            calculateAllBFieldACCELERATED(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision);

            for (int i = 0; i < fieldH.size(); ++i)
            {
                computedFieldArr.push_back(fieldH[i] * cos(cylindricalPhiArr[i]));
            }
        }
        // TODO - CPU_MT method
    }
    else
    {
        throw "Number of elements in input data vectors is not the same!";
    }
}

void Coil::computeAllBFieldX(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             ComputeMethod method) const
{
    computeAllBFieldX(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedFieldArr, precisionSettings, method);
}

void Coil::computeAllBFieldY(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             const PrecisionArguments &usedPrecision,
                             ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        if (method == CPU_ST)
        {
            calculateAllBFieldHorizontalSINGLE(
                    cylindricalZArr, cylindricalRArr, computedFieldArr, usedPrecision);

            for (int i = 0; i < computedFieldArr.size(); ++i)
            {
                computedFieldArr[i] *= sin(cylindricalPhiArr[i]);
            }
        }
        else if (method == GPU)
        {
            std::vector<float> fieldH;
            std::vector<float> fieldZ;

            calculateAllBFieldACCELERATED(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision);

            for (int i = 0; i < fieldH.size(); ++i)
            {
                computedFieldArr.push_back(fieldH[i] * cos(cylindricalPhiArr[i]));
            }
        }
        // TODO - CPU_MT method
    }
    else
    {
        throw "Number of elements in input data vectors is not the same!";
    }
}

void Coil::computeAllBFieldY(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             ComputeMethod method) const
{
    computeAllBFieldY(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedFieldArr, precisionSettings, method);
}

void Coil::computeAllBFieldH(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             const PrecisionArguments &usedPrecision,
                             ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        if (method == CPU_ST)
        {
            calculateAllBFieldHorizontalSINGLE(
                    cylindricalZArr, cylindricalRArr, computedFieldArr, usedPrecision);
        }
        else if (method == GPU)
        {
            std::vector<float> fieldH;
            std::vector<float> fieldZ;

            calculateAllBFieldACCELERATED(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision);

            for (float i : fieldH)
                computedFieldArr.push_back(i);
        }
        // TODO - CPU_MT method
    }
    else
    {
        throw "Number of elements in input data vectors is not the same!";
    }
}

void Coil::computeAllBFieldH(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             ComputeMethod method) const
{
    computeAllBFieldH(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedFieldArr, precisionSettings, method);
}

void Coil::computeAllBFieldZ(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             const PrecisionArguments &usedPrecision,
                             ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        if (method == CPU_ST)
        {
            calculateAllBFieldVerticalSINGLE(
                    cylindricalZArr, cylindricalRArr, computedFieldArr, usedPrecision);
        }
        else if (method == GPU)
        {
            std::vector<float> fieldH;
            std::vector<float> fieldZ;

            calculateAllBFieldACCELERATED(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision);

            for (float i : fieldZ)
                computedFieldArr.push_back(i);
        }
        // TODO - CPU_MT method
    }
    else
    {
        throw "Number of elements in input data vectors is not the same!";
    }
}

void Coil::computeAllBFieldZ(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             ComputeMethod method) const
{
    computeAllBFieldZ(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedFieldArr, precisionSettings, method);
}

void
Coil::computeAllBFieldComponents(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedFieldXArr,
                                 std::vector<double> &computedFieldYArr,
                                 std::vector<double> &computedFieldZArr,
                                 const PrecisionArguments &usedPrecision,
                                 ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedFieldXArr.resize(0);
        computedFieldYArr.resize(0);
        computedFieldZArr.resize(0);

        if (method == CPU_ST)
        {
            std::vector<double> hFieldArray;
            std::vector<double> zFieldArray;

            calculateAllBFieldSINGLE(
                    cylindricalZArr, cylindricalRArr, hFieldArray, zFieldArray, usedPrecision);

            for (int i = 0; i < hFieldArray.size(); ++i)
            {
                computedFieldXArr.push_back(hFieldArray[i] * cos(cylindricalPhiArr[i]));
                computedFieldYArr.push_back(hFieldArray[i] * sin(cylindricalPhiArr[i]));
                computedFieldZArr.push_back(zFieldArray[i]);
            }
        }
        else if (method == GPU)
        {
            std::vector<float> fieldH;
            std::vector<float> fieldZ;

            calculateAllBFieldACCELERATED(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision);

            for (int i = 0; i < fieldH.size(); ++i)
            {
                computedFieldXArr.push_back(fieldH[i] * cos(cylindricalPhiArr[i]));
                computedFieldYArr.push_back(fieldH[i] * sin(cylindricalPhiArr[i]));
                computedFieldZArr.push_back(fieldZ[i]);
            }
        }
        // TODO - CPU_MT method
    }
    else
    {
        throw "Number of elements in input data vectors is not the same!";
    }
}

void
Coil::computeAllBFieldComponents(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedFieldXArr,
                                 std::vector<double> &computedFieldYArr,
                                 std::vector<double> &computedFieldZArr,
                                 ComputeMethod method) const
{
    computeAllBFieldComponents(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
            computedFieldXArr, computedFieldYArr, computedFieldZArr,
            precisionSettings, method);
}

void Coil::computeAllAPotentialX(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedPotentialArr,
                                 const PrecisionArguments &usedPrecision,
                                 ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedPotentialArr.resize(0);

        if (method == CPU_ST)
        {
            calculateAllAPotentialSINGLE(
                    cylindricalZArr, cylindricalRArr, computedPotentialArr, usedPrecision);

            for (int i = 0; i < computedPotentialArr.size(); ++i)
            {
                computedPotentialArr[i] *= (-1) * sin(cylindricalPhiArr[i]);
            }
        }
        else if (method == GPU)
        {
            std::vector<float> potentialA;

            calculateAllAPotentialACCELERATED(cylindricalZArr, cylindricalRArr, potentialA, usedPrecision);

            for (int i = 0; i < potentialA.size(); ++i)
            {
                computedPotentialArr.push_back(potentialA[i] * (-1) * sin(cylindricalPhiArr[i]));
            }
        }
        // TODO - CPU_MT method
    }
    else
    {
        throw "Number of elements in input data vectors is not the same!";
    }
}

void Coil::computeAllAPotentialX(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedPotentialArr,
                                 ComputeMethod method) const
{
    computeAllAPotentialX(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedPotentialArr, precisionSettings, method);
}

void Coil::computeAllAPotentialY(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedPotentialArr,
                                 const PrecisionArguments &usedPrecision,
                                 ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedPotentialArr.resize(0);

        if (method == CPU_ST)
        {
            calculateAllAPotentialSINGLE(
                    cylindricalZArr, cylindricalRArr, computedPotentialArr, usedPrecision);

            for (int i = 0; i < computedPotentialArr.size(); ++i)
            {
                computedPotentialArr[i] *= cos(cylindricalPhiArr[i]);
            }
        }
        else if (method == GPU)
        {
            std::vector<float> potentialA;

            calculateAllAPotentialACCELERATED(cylindricalZArr, cylindricalRArr, potentialA, usedPrecision);

            for (int i = 0; i < potentialA.size(); ++i)
            {
                computedPotentialArr.push_back(potentialA[i] * cos(cylindricalPhiArr[i]));
            }
        }
        // TODO - CPU_MT method
    }
    else
    {
        throw "Number of elements in input data vectors is not the same!";
    }
}

void Coil::computeAllAPotentialY(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedPotentialArr,
                                 ComputeMethod method) const
{
    computeAllAPotentialY(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedPotentialArr, precisionSettings, method);
}

void
Coil::computeAllAPotentialAbs(const std::vector<double> &cylindricalZArr,
                              const std::vector<double> &cylindricalRArr,
                              std::vector<double> &computedPotentialArr,
                              const PrecisionArguments &usedPrecision,
                              ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size())
    {
        computedPotentialArr.resize(0);

        if (method == CPU_ST)
        {
            calculateAllAPotentialSINGLE(
                    cylindricalZArr, cylindricalRArr, computedPotentialArr, usedPrecision);
        }
        else if (method == GPU)
        {
            std::vector<float> potentialA;

            calculateAllAPotentialACCELERATED(cylindricalZArr, cylindricalRArr, potentialA, usedPrecision);

            for (float i : potentialA){
                computedPotentialArr.push_back(i);
            }

        }
        // TODO - CPU_MT method
    }
    else
    {
        throw "Number of elements in input data vectors is not the same!";
    }
}

void
Coil::computeAllAPotentialAbs(const std::vector<double> &cylindricalZArr,
                              const std::vector<double> &cylindricalRArr,
                              std::vector<double> &computedPotentialArr,
                              ComputeMethod method) const
{
    computeAllAPotentialAbs(
            cylindricalZArr, cylindricalRArr, computedPotentialArr, precisionSettings, method);
}

void Coil::computeAllAPotentialComponents(const std::vector<double> &cylindricalZArr,
                                          const std::vector<double> &cylindricalRArr,
                                          const std::vector<double> &cylindricalPhiArr,
                                          std::vector<double> &computedPotentialXArr,
                                          std::vector<double> &computedPotentialYArr,
                                          std::vector<double> &computedPotentialZArr,
                                          const PrecisionArguments &usedPrecision,
                                          ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedPotentialXArr.resize(0);
        computedPotentialYArr.resize(0);
        computedPotentialZArr.resize(0);

        if (method == CPU_ST)
        {
            std::vector<double> potentialArray;

            calculateAllAPotentialSINGLE(
                    cylindricalZArr, cylindricalRArr, potentialArray, usedPrecision);

            for (int i = 0; i < potentialArray.size(); ++i)
            {
                computedPotentialXArr.push_back(potentialArray[i] * (-1) * sin(cylindricalPhiArr[i]));
                computedPotentialYArr.push_back(potentialArray[i] * cos(cylindricalPhiArr[i]));
                computedPotentialZArr.push_back(0.0);
            }
        }
        else if (method == GPU)
        {
            std::vector<float> potentialA;

            calculateAllAPotentialACCELERATED(cylindricalZArr, cylindricalRArr, potentialA, usedPrecision);

            for (int i = 0; i < potentialA.size(); ++i)
            {
                computedPotentialXArr.push_back(potentialA[i] * (-1) * sin(cylindricalPhiArr[i]));
                computedPotentialYArr.push_back(potentialA[i] * cos(cylindricalPhiArr[i]));
                computedPotentialZArr.push_back(0.0);
            }
        }
        //TODO - other computation methods
    }
    else
    {
        throw "Number of elements in input data vectors is not the same!";
    }
}

void Coil::computeAllAPotentialComponents(const std::vector<double> &cylindricalZArr,
                                          const std::vector<double> &cylindricalRArr,
                                          const std::vector<double> &cylindricalPhiArr,
                                          std::vector<double> &computedPotentialXArr,
                                          std::vector<double> &computedPotentialYArr,
                                          std::vector<double> &computedPotentialZArr,
                                          ComputeMethod method) const
{
    computeAllAPotentialComponents(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
            computedPotentialXArr, computedPotentialYArr, computedPotentialZArr, precisionSettings, method);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary, double zDisplacement,
                                     double precisionFactor, ComputeMethod method)
{
    PrecisionArguments precisionArguments;
    int linearIncrements;

    PrecisionArguments::getMutualInductancePrecisionSettingsZCPU(primary, secondary,
                                                                 precisionFactor,
                                                                 precisionArguments, linearIncrements);
    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    linearIncrements--; //

    std::vector<double> zPositions;
    std::vector<double> rPositions;

    std::vector<double> weights;

    for (int zIndex = 0; zIndex <= linearIncrements; ++zIndex)
    {
        for (int rIndex = 0; rIndex <= linearIncrements; ++rIndex)
        {
            zPositions.push_back(zDisplacement + (secondary.length * 0.5) *
            Legendre::positionMatrix[linearIncrements][zIndex]);

            rPositions.push_back(secondary.innerRadius + secondary.thickness * 0.5 +
            (secondary.thickness * 0.5) * Legendre::positionMatrix[linearIncrements][rIndex]);

            weights.push_back(0.25 *
            Legendre::weightsMatrix[linearIncrements][zIndex] *
            Legendre::weightsMatrix[linearIncrements][rIndex]);
        }
    }

    std::vector<double> potentialA;
    double mutualInductance = 0.0;

    primary.computeAllAPotentialAbs(zPositions, rPositions, potentialA, precisionArguments, method);

    for (int i = 0; i < potentialA.size(); ++i)
    {
        mutualInductance += 2*PI * rPositions[i] * potentialA[i] * weights[i];
    }

    return mutualInductance * secondary.numOfTurns;
}
