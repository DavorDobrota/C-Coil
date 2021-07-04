#include <cstdio>
#include <cmath>
#include <vector>
#include <functional>
#include <chrono>

#include "Coil.h"
#include "Polynomial.h"

namespace
{
    const double PI = 3.14159265357989323;
    const double g_MiReduced = 0.0000001;

    const int g_maxLegendrePol = 20;

    const double g_defaultCurrent = 1.0;
    const double g_defaultResistivity = 1.63e-8;
    const double g_defaultSineFrequency = 50;
    const PrecisionArguments g_defaultPrecision = PrecisionArguments(1, 1, 1, 16, 12, 12);

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
    Polynomial::getLegendreParametersForN(numOfAngularIncrements,
                                          angularIncrementPositions, angularIncrementWeights);

    Polynomial::getLegendreParametersForN(numOfThicknessIncrements,
                                          thicknessIncrementPositions, thicknessIncrementWeights);

    Polynomial::getLegendreParametersForN(numOfLengthIncrements,
                                          lengthIncrementPositions, lengthIncrementWeights);
}

void PrecisionArguments::genParametersFromPrecision()
{
    //TODO - when further analysis is complete a method will be divised, until then
    numOfAngularBlocks = 2;
    numOfLengthBlocks = 1;
    numOfThicknessBlocks = 1;

    numOfAngularIncrements = 12;
    numOfThicknessIncrements = 12;
    numOfLengthIncrements = 12;
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
}

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
            (pow(innerRadius, 2) + pow(thickness, 2) + innerRadius * thickness / 3);
}

void Coil::calculateAverageWireThickness()
{
    averageWireThickness = sqrt(length * thickness / numOfTurns);
}

void Coil::calculateResistance()
{
    double wireRadius = averageWireThickness / 2;
    double ohmicResistance = wireResistivity * numOfTurns * 2*PI * (innerRadius + thickness / 2) / (pow(wireRadius, 2) * PI);
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

double Coil::calculateBFieldVertical(double zAxis, double rPolar)
{
    double magneticFieldZ = 0.0;

    double lengthBlock = length / precisionSettings.numOfLengthBlocks;
    double thicknessBlock = thickness / precisionSettings.numOfThicknessBlocks;
    double angularBlock = PI / precisionSettings.numOfAngularBlocks;

    double constant = g_MiReduced * currentDensity * lengthBlock * thicknessBlock * angularBlock * 2;

    for (int indBlockL = 0; indBlockL < precisionSettings.numOfLengthBlocks; ++indBlockL)
    {
        for (int indBlockT = 0; indBlockT < precisionSettings.numOfThicknessBlocks; ++indBlockT)
        {
            double blockPositionL = (-1) * (length / 2) + lengthBlock * (indBlockL + 0.5);
            double blockPositionT = innerRadius + thicknessBlock * (indBlockT + 0.5);

            for (int incL = 0; incL < precisionSettings.numOfLengthIncrements; ++incL)
            {
                for (int incT = 0; incT < precisionSettings.numOfThicknessIncrements; ++incT)
                {
                    double incrementPositionL = blockPositionL +
                            (lengthBlock / 2) * precisionSettings.lengthIncrementPositions[incL];
                    double incrementPositionT = blockPositionT +
                            (thicknessBlock / 2) * precisionSettings.thicknessIncrementPositions[incT];

                    double incrementWeightL = precisionSettings.lengthIncrementWeights[incL] / 2;
                    double incrementWeightT = precisionSettings.thicknessIncrementWeights[incT] / 2;

                    double incrementWeightS = incrementWeightL * incrementWeightT;

                    double tempConstA = pow(incrementPositionT, 2);
                    double tempConstB = incrementPositionT * rPolar;
                    double tempConstC = tempConstA + pow(rPolar, 2) + pow((incrementPositionL + zAxis), 2);

                    for (int indBlockFi = 0; indBlockFi < precisionSettings.numOfAngularBlocks; ++indBlockFi)
                    {
                        double blockPositionFi = angularBlock * (indBlockFi + 0.5);

                        for (int incFi = 0; incFi < precisionSettings.numOfAngularIncrements; ++incFi)
                        {
                            double incrementPositionFi = blockPositionFi +
                                    (angularBlock / 2) * precisionSettings.angularIncrementPositions[incFi];

                            double incrementWeightFi = precisionSettings.angularIncrementWeights[incFi] / 2;

                            double tempConstD = tempConstB * cos(incrementPositionFi);

//                            printf("%d %d %d %2d %2d %2d : ", indBlockL, indBlockT, indBlockFi, incL, incT, incFi);
//                            printf("%.10f %.10f %.10f - ", incrementWeightL, incrementWeightT, incrementWeightFi);
//                            printf("%.10f %.10f %.10f\n", incrementPositionL, incrementPositionT, incrementPositionFi);

                            magneticFieldZ += constant * incrementWeightS * incrementWeightFi *
                                    (tempConstA - tempConstD)/pow((tempConstC - 2*tempConstD), 1.5);
                        }
                    }
                }
            }
        }
    }
    return magneticFieldZ;
}

double Coil::calculateBFieldHorizontal(double zAxis, double rPolar)
{
    double magneticFieldH = 0.0;

    double lengthBlock = length / precisionSettings.numOfLengthBlocks;
    double thicknessBlock = thickness / precisionSettings.numOfThicknessBlocks;
    double angularBlock = PI / precisionSettings.numOfAngularBlocks;

    double constant = g_MiReduced * currentDensity * lengthBlock * thicknessBlock * angularBlock * 2;

    for (int indBlockL = 0; indBlockL < precisionSettings.numOfLengthBlocks; ++indBlockL)
    {
        for (int indBlockT = 0; indBlockT < precisionSettings.numOfThicknessBlocks; ++indBlockT)
        {
            double blockPositionL = (-1) * (length / 2) + lengthBlock * (indBlockL + 0.5);
            double blockPositionT = innerRadius + thicknessBlock * (indBlockT + 0.5);

            for (int incL = 0; incL < precisionSettings.numOfLengthIncrements; ++incL)
            {
                for (int incT = 0; incT < precisionSettings.numOfThicknessIncrements; ++incT)
                {
                    double incrementPositionL = blockPositionL +
                                                (lengthBlock / 2) * precisionSettings.lengthIncrementPositions[incL];
                    double incrementPositionT = blockPositionT +
                                                (thicknessBlock / 2) * precisionSettings.thicknessIncrementPositions[incT];

                    double incrementWeightL = precisionSettings.lengthIncrementWeights[incL] / 2;
                    double incrementWeightT = precisionSettings.thicknessIncrementWeights[incT] / 2;

                    double incrementWeightS = incrementWeightL * incrementWeightT;

                    double tempConstA = incrementPositionT * (incrementPositionL + zAxis);
                    double tempConstB = incrementPositionT * rPolar;
                    double tempConstC = pow(incrementPositionT, 2) +
                            pow(rPolar, 2) + pow((incrementPositionL + zAxis), 2);

                    for (int indBlockFi = 0; indBlockFi < precisionSettings.numOfAngularBlocks; ++indBlockFi)
                    {
                        double blockPositionFi = angularBlock * (indBlockFi + 0.5);

                        for (int incFi = 0; incFi < precisionSettings.numOfAngularIncrements; ++incFi)
                        {
                            double incrementPositionFi = blockPositionFi +
                                                         (angularBlock / 2) * precisionSettings.angularIncrementPositions[incFi];

                            double incrementWeightFi = precisionSettings.angularIncrementWeights[incFi] / 2;

                            double cosinePhi = cos(incrementPositionFi);

                            magneticFieldH += constant * incrementWeightS * incrementWeightFi *
                                    (tempConstA * cosinePhi) /pow((tempConstC - 2*tempConstB * cosinePhi), 1.5);
                        }
                    }
                }
            }
        }
    }
    return magneticFieldH;
}

double Coil::calculateAPotential(double zAxis, double rPolar)
{
    double magneticPotential = 0.0;

    double lengthBlock = length / precisionSettings.numOfLengthBlocks;
    double thicknessBlock = thickness / precisionSettings.numOfThicknessBlocks;
    double angularBlock = PI / precisionSettings.numOfAngularBlocks;

    double constant = g_MiReduced * currentDensity * lengthBlock * thicknessBlock * angularBlock * 2;

    for (int indBlockL = 0; indBlockL < precisionSettings.numOfLengthBlocks; ++indBlockL)
    {
        for (int indBlockT = 0; indBlockT < precisionSettings.numOfThicknessBlocks; ++indBlockT)
        {
            double blockPositionL = (-1) * (length / 2) + lengthBlock * (indBlockL + 0.5);
            double blockPositionT = innerRadius + thicknessBlock * (indBlockT + 0.5);

            for (int incL = 0; incL < precisionSettings.numOfLengthIncrements; ++incL)
            {
                for (int incT = 0; incT < precisionSettings.numOfThicknessIncrements; ++incT)
                {
                    double incrementPositionL = blockPositionL +
                                                (lengthBlock / 2) * precisionSettings.lengthIncrementPositions[incL];
                    double incrementPositionT = blockPositionT +
                                                (thicknessBlock / 2) * precisionSettings.thicknessIncrementPositions[incT];

                    double incrementWeightL = precisionSettings.lengthIncrementWeights[incL] / 2;
                    double incrementWeightT = precisionSettings.thicknessIncrementWeights[incT] / 2;

                    double incrementWeightS = incrementWeightL * incrementWeightT;

                    double tempConstA = incrementPositionT;
                    double tempConstB = incrementPositionT * rPolar;
                    double tempConstC = pow(incrementPositionT, 2) +
                                        pow(rPolar, 2) + pow((incrementPositionL + zAxis), 2);

                    for (int indBlockFi = 0; indBlockFi < precisionSettings.numOfAngularBlocks; ++indBlockFi)
                    {
                        double blockPositionFi = angularBlock * (indBlockFi + 0.5);

                        for (int incFi = 0; incFi < precisionSettings.numOfAngularIncrements; ++incFi)
                        {
                            double incrementPositionFi = blockPositionFi +
                                                         (angularBlock / 2) * precisionSettings.angularIncrementPositions[incFi];

                            double incrementWeightFi = precisionSettings.angularIncrementWeights[incFi] / 2;

                            double cosinePhi = cos(incrementPositionFi);

                            magneticPotential += constant * incrementWeightS * incrementWeightFi *
                                              (tempConstA * cosinePhi) /pow((tempConstC - 2*tempConstB * cosinePhi), 1.5);
                        }
                    }
                }
            }
        }
    }
    return magneticPotential;
}

double Coil::computeBFieldX(double cylindricalZ, double cylindricalR, double cylindricalPhi)
{
    return calculateBFieldHorizontal(cylindricalZ, cylindricalR) * cos(cylindricalPhi);
}

double Coil::computeBFieldY(double cylindricalZ, double cylindricalR, double cylindricalPhi)
{
    return calculateBFieldHorizontal(cylindricalZ, cylindricalR) * sin(cylindricalPhi);
}

double Coil::computeBFieldH(double cylindricalZ, double cylindricalR)
{
    return calculateBFieldHorizontal(cylindricalZ, cylindricalR);
}

double Coil::computeBFieldZ(double cylindricalZ, double cylindricalR)
{
    return calculateBFieldVertical(cylindricalZ, cylindricalR);
}

std::vector<double> Coil::computeBFieldVector(double cylindricalZ, double cylindricalR, double cylindricalPhi)
{
    std::vector<double> fieldVector;
    double FieldZ = calculateBFieldVertical(cylindricalZ, cylindricalR);
    double FieldH = calculateBFieldHorizontal(cylindricalZ, cylindricalR);

    fieldVector.push_back(FieldH * cos(cylindricalPhi));
    fieldVector.push_back(FieldH * sin(cylindricalPhi));
    fieldVector.push_back(FieldZ);

    return fieldVector;
}
