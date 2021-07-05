#include <cstdio>
#include <cmath>
#include <vector>
#include <functional>
#include <chrono>

#include "Coil.h"
#include "Polynomial.h"
#include "ComputeMethod.h"

namespace
{
    const double PI = 3.14159265357989323;
    const double g_MiReduced = 0.0000001;

    const int g_maxLegendrePol = 20;

    const double g_defaultCurrent = 1.0;
    const double g_defaultResistivity = 1.63e-8;
    const double g_defaultSineFrequency = 50;
    const PrecisionArguments g_defaultPrecision = PrecisionArguments(1, 1, 1, 12, 12, 12);

}

namespace Precision
{
    const PrecisionArguments defaultPrecision_ULTRAFAST = PrecisionArguments(1, 1, 1, 12, 8, 8);
    const PrecisionArguments defaultPrecision_FAST = PrecisionArguments(1, 1, 1, 12, 12, 12);
    const PrecisionArguments defaultPrecision_NORMAL = PrecisionArguments(2, 1, 1, 12, 12, 12);
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
    calculateSelfInductance();
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

std::pair<double, double> Coil::calculateBField(double zAxis, double rPolar, const PrecisionArguments &usedPrecision)
{
    double magneticFieldZ = 0.0;
    double magneticFieldH = 0.0;
    
    double lengthBlock = length / usedPrecision.numOfLengthBlocks;
    double thicknessBlock = thickness / usedPrecision.numOfThicknessBlocks;
    double angularBlock = PI / usedPrecision.numOfAngularBlocks;

    double constant = g_MiReduced * currentDensity * lengthBlock * thicknessBlock * angularBlock * 2;

    for (int indBlockL = 0; indBlockL < usedPrecision.numOfLengthBlocks; ++indBlockL)
    {
        for (int indBlockT = 0; indBlockT < usedPrecision.numOfThicknessBlocks; ++indBlockT)
        {
            double blockPositionL = (-1) * (length / 2) + lengthBlock * (indBlockL + 0.5);
            double blockPositionT = innerRadius + thicknessBlock * (indBlockT + 0.5);

            for (int incL = 0; incL < usedPrecision.numOfLengthIncrements; ++incL)
            {
                for (int incT = 0; incT < usedPrecision.numOfThicknessIncrements; ++incT)
                {
                    double incrementPositionL = blockPositionL +
                                                (lengthBlock / 2) * usedPrecision.lengthIncrementPositions[incL];
                    double incrementPositionT = blockPositionT +
                                                (thicknessBlock / 2) * usedPrecision.thicknessIncrementPositions[incT];

                    double incrementWeightL = usedPrecision.lengthIncrementWeights[incL] / 2;
                    double incrementWeightT = usedPrecision.thicknessIncrementWeights[incT] / 2;

                    double incrementWeightS = incrementWeightL * incrementWeightT;

                    double tempConstA = pow(incrementPositionT, 2);
                    double tempConstB = incrementPositionT * (incrementPositionL + zAxis);
                    double tempConstC = incrementPositionT * rPolar;
                    double tempConstD = tempConstA + pow(rPolar, 2) + pow((incrementPositionL + zAxis), 2);
                    double tempConstE = constant * incrementWeightS;

                    for (int indBlockFi = 0; indBlockFi < usedPrecision.numOfAngularBlocks; ++indBlockFi)
                    {
                        double blockPositionFi = angularBlock * (indBlockFi + 0.5);

                        for (int incFi = 0; incFi < usedPrecision.numOfAngularIncrements; ++incFi)
                        {
                            double incrementPositionFi = blockPositionFi +
                                                         (angularBlock / 2) * usedPrecision.angularIncrementPositions[incFi];

                            double incrementWeightFi = usedPrecision.angularIncrementWeights[incFi] / 2;

                            double cosinePhi = cos(incrementPositionFi);
                            double tempConstF = 2 * tempConstC * cosinePhi;
                            double tempConstG = tempConstE * incrementWeightFi;

                            magneticFieldZ += tempConstG *
                                              (tempConstA - tempConstC * cosinePhi)/pow((tempConstD - tempConstF), 1.5);
                            magneticFieldH += tempConstG *
                                              (tempConstB * cosinePhi) /pow((tempConstD - tempConstF), 1.5);
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

double Coil::calculateBFieldVertical(double zAxis, double rPolar, const PrecisionArguments &usedPrecision)
{
    double magneticFieldZ = 0.0;
    
    double lengthBlock = length / usedPrecision.numOfLengthBlocks;
    double thicknessBlock = thickness / usedPrecision.numOfThicknessBlocks;
    double angularBlock = PI / usedPrecision.numOfAngularBlocks;

    double constant = g_MiReduced * currentDensity * lengthBlock * thicknessBlock * angularBlock * 2;

    for (int indBlockL = 0; indBlockL < usedPrecision.numOfLengthBlocks; ++indBlockL)
    {
        for (int indBlockT = 0; indBlockT < usedPrecision.numOfThicknessBlocks; ++indBlockT)
        {
            double blockPositionL = (-1) * (length / 2) + lengthBlock * (indBlockL + 0.5);
            double blockPositionT = innerRadius + thicknessBlock * (indBlockT + 0.5);

            for (int incL = 0; incL < usedPrecision.numOfLengthIncrements; ++incL)
            {
                for (int incT = 0; incT < usedPrecision.numOfThicknessIncrements; ++incT)
                {
                    double incrementPositionL = blockPositionL +
                                                (lengthBlock / 2) * usedPrecision.lengthIncrementPositions[incL];
                    double incrementPositionT = blockPositionT +
                                                (thicknessBlock / 2) * usedPrecision.thicknessIncrementPositions[incT];

                    double incrementWeightL = usedPrecision.lengthIncrementWeights[incL] / 2;
                    double incrementWeightT = usedPrecision.thicknessIncrementWeights[incT] / 2;

                    double incrementWeightS = incrementWeightL * incrementWeightT;

                    double tempConstA = pow(incrementPositionT, 2);
                    double tempConstB = incrementPositionT * rPolar;
                    double tempConstC = tempConstA + pow(rPolar, 2) + pow((incrementPositionL + zAxis), 2);

                    for (int indBlockFi = 0; indBlockFi < usedPrecision.numOfAngularBlocks; ++indBlockFi)
                    {
                        double blockPositionFi = angularBlock * (indBlockFi + 0.5);

                        for (int incFi = 0; incFi < usedPrecision.numOfAngularIncrements; ++incFi)
                        {
                            double incrementPositionFi = blockPositionFi +
                                                         (angularBlock / 2) * usedPrecision.angularIncrementPositions[incFi];

                            double incrementWeightFi = usedPrecision.angularIncrementWeights[incFi] / 2;

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

double Coil::calculateBFieldHorizontal(double zAxis, double rPolar, const PrecisionArguments &usedPrecision)
{
    double magneticFieldH = 0.0;
    
    double lengthBlock = length / usedPrecision.numOfLengthBlocks;
    double thicknessBlock = thickness / usedPrecision.numOfThicknessBlocks;
    double angularBlock = PI / usedPrecision.numOfAngularBlocks;

    double constant = g_MiReduced * currentDensity * lengthBlock * thicknessBlock * angularBlock * 2;

    for (int indBlockL = 0; indBlockL < usedPrecision.numOfLengthBlocks; ++indBlockL)
    {
        for (int indBlockT = 0; indBlockT < usedPrecision.numOfThicknessBlocks; ++indBlockT)
        {
            double blockPositionL = (-1) * (length / 2) + lengthBlock * (indBlockL + 0.5);
            double blockPositionT = innerRadius + thicknessBlock * (indBlockT + 0.5);

            for (int incL = 0; incL < usedPrecision.numOfLengthIncrements; ++incL)
            {
                for (int incT = 0; incT < usedPrecision.numOfThicknessIncrements; ++incT)
                {
                    double incrementPositionL = blockPositionL +
                                                (lengthBlock / 2) * usedPrecision.lengthIncrementPositions[incL];
                    double incrementPositionT = blockPositionT +
                                                (thicknessBlock / 2) * usedPrecision.thicknessIncrementPositions[incT];

                    double incrementWeightL = usedPrecision.lengthIncrementWeights[incL] / 2;
                    double incrementWeightT = usedPrecision.thicknessIncrementWeights[incT] / 2;

                    double incrementWeightS = incrementWeightL * incrementWeightT;

                    double tempConstA = incrementPositionT * (incrementPositionL + zAxis);
                    double tempConstB = incrementPositionT * rPolar;
                    double tempConstC = pow(incrementPositionT, 2) +
                            pow(rPolar, 2) + pow((incrementPositionL + zAxis), 2);

                    for (int indBlockFi = 0; indBlockFi < usedPrecision.numOfAngularBlocks; ++indBlockFi)
                    {
                        double blockPositionFi = angularBlock * (indBlockFi + 0.5);

                        for (int incFi = 0; incFi < usedPrecision.numOfAngularIncrements; ++incFi)
                        {
                            double incrementPositionFi = blockPositionFi +
                                                         (angularBlock / 2) * usedPrecision.angularIncrementPositions[incFi];

                            double incrementWeightFi = usedPrecision.angularIncrementWeights[incFi] / 2;

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

double Coil::calculateAPotential(double zAxis, double rPolar, const PrecisionArguments &usedPrecision)
{
    double magneticPotential = 0.0;

    double lengthBlock = length / usedPrecision.numOfLengthBlocks;
    double thicknessBlock = thickness / usedPrecision.numOfThicknessBlocks;
    double angularBlock = PI / usedPrecision.numOfAngularBlocks;

    double constant = g_MiReduced * currentDensity * lengthBlock * thicknessBlock * angularBlock * 2;

    for (int indBlockL = 0; indBlockL < usedPrecision.numOfLengthBlocks; ++indBlockL)
    {
        for (int indBlockT = 0; indBlockT < usedPrecision.numOfThicknessBlocks; ++indBlockT)
        {
            double blockPositionL = (-1) * (length / 2) + lengthBlock * (indBlockL + 0.5);
            double blockPositionT = innerRadius + thicknessBlock * (indBlockT + 0.5);

            for (int incL = 0; incL < usedPrecision.numOfLengthIncrements; ++incL)
            {
                for (int incT = 0; incT < usedPrecision.numOfThicknessIncrements; ++incT)
                {
                    double incrementPositionL = blockPositionL +
                                                (lengthBlock / 2) * usedPrecision.lengthIncrementPositions[incL];
                    double incrementPositionT = blockPositionT +
                                                (thicknessBlock / 2) * usedPrecision.thicknessIncrementPositions[incT];

                    double incrementWeightL = usedPrecision.lengthIncrementWeights[incL] / 2;
                    double incrementWeightT = usedPrecision.thicknessIncrementWeights[incT] / 2;

                    double incrementWeightS = incrementWeightL * incrementWeightT;

                    double tempConstA = incrementPositionT;
                    double tempConstB = incrementPositionT * rPolar;
                    double tempConstC = pow(incrementPositionT, 2) +
                                        pow(rPolar, 2) + pow((incrementPositionL + zAxis), 2);

                    for (int indBlockFi = 0; indBlockFi < usedPrecision.numOfAngularBlocks; ++indBlockFi)
                    {
                        double blockPositionFi = angularBlock * (indBlockFi + 0.5);

                        for (int incFi = 0; incFi < usedPrecision.numOfAngularIncrements; ++incFi)
                        {
                            double incrementPositionFi = blockPositionFi +
                                                         (angularBlock / 2) * usedPrecision.angularIncrementPositions[incFi];

                            double incrementWeightFi = usedPrecision.angularIncrementWeights[incFi] / 2;

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
    return calculateBFieldHorizontal(cylindricalZ, cylindricalR, precisionSettings) * cos(cylindricalPhi);
}

double Coil::computeBFieldX(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                            const PrecisionArguments &usedPrecision)
{
    return calculateBFieldHorizontal(cylindricalZ, cylindricalR, usedPrecision) * cos(cylindricalPhi);
}

double Coil::computeBFieldY(double cylindricalZ, double cylindricalR, double cylindricalPhi)
{
    return calculateBFieldHorizontal(cylindricalZ, cylindricalR, precisionSettings) * sin(cylindricalPhi);
}

double Coil::computeBFieldY(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                            const PrecisionArguments &usedPrecision)
{
    return calculateBFieldHorizontal(cylindricalZ, cylindricalR, usedPrecision) * sin(cylindricalPhi);
}

double Coil::computeBFieldH(double cylindricalZ, double cylindricalR)
{
    return calculateBFieldHorizontal(cylindricalZ, cylindricalR, precisionSettings);
}

double Coil::computeBFieldH(double cylindricalZ, double cylindricalR, const PrecisionArguments &usedPrecision)
{
    return calculateBFieldHorizontal(cylindricalZ, cylindricalR, usedPrecision);
}

double Coil::computeBFieldZ(double cylindricalZ, double cylindricalR)
{
    return calculateBFieldVertical(cylindricalZ, cylindricalR, precisionSettings);
}

double Coil::computeBFieldZ(double cylindricalZ, double cylindricalR, const PrecisionArguments &usedPrecision)
{
    return calculateBFieldVertical(cylindricalZ, cylindricalR, usedPrecision);
}



std::vector<double> Coil::computeBFieldVector(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                              const PrecisionArguments &usedPrecision)
{
    std::vector<double> fieldVector;
    std::pair<double, double> fields = calculateBField(cylindricalZ, cylindricalR, usedPrecision);

    fieldVector.push_back(fields.first * cos(cylindricalPhi));
    fieldVector.push_back(fields.first * sin(cylindricalPhi));
    fieldVector.push_back(fields.second);

    return fieldVector;
}
std::vector<double> Coil::computeBFieldVector(double cylindricalZ, double cylindricalR, double cylindricalPhi)
{
    return computeBFieldVector(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

double Coil::computeAPotentialX(double cylindricalZ, double cylindricalR, double cylindricalPhi)
{
    return (-1) * calculateAPotential(cylindricalZ, cylindricalR, precisionSettings) * sin(cylindricalPhi);
}

double Coil::computeAPotentialX(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                const PrecisionArguments &usedPrecision)
{
    return (-1) * calculateAPotential(cylindricalZ, cylindricalR, usedPrecision) * sin(cylindricalPhi);
}

double Coil::computeAPotentialY(double cylindricalZ, double cylindricalR, double cylindricalPhi)
{
    return calculateAPotential(cylindricalZ, cylindricalR, precisionSettings) * cos(cylindricalPhi);
}

double Coil::computeAPotentialY(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                const PrecisionArguments &usedPrecision)
{
    return calculateAPotential(cylindricalZ, cylindricalR, usedPrecision) * cos(cylindricalPhi);
}

double Coil::computeAPotentialAbs(double cylindricalZ, double cylindricalR)
{
    return calculateAPotential(cylindricalZ, cylindricalR, precisionSettings);
}

double Coil::computeAPotentialAbs(double cylindricalZ, double cylindricalR, PrecisionArguments &usedPrecision)
{
    return calculateAPotential(cylindricalZ, cylindricalR, usedPrecision);
}

std::vector<double> Coil::computeAPotentialVector(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                                  const PrecisionArguments &usedPrecision)
{
    std::vector<double> potentialVector;
    double potential = calculateAPotential(cylindricalZ, cylindricalR, usedPrecision);

    potentialVector.push_back(potential * (-sin(cylindricalPhi)));
    potentialVector.push_back(potential * cos(cylindricalPhi));
    return potentialVector;
}

std::vector<double> Coil::computeAPotentialVector(double cylindricalZ, double cylindricalR, double cylindricalPhi)
{
    return computeAPotentialVector(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

void Coil::calculateAllBFieldSINGLE(const std::vector<double> &cylindricalZArr,
                                    const std::vector<double> &cylindricalRArr,
                                    std::vector<double> &computedFieldHArr,
                                    std::vector<double> &computedFieldZArr,
                                    const PrecisionArguments &usedPrecision) {
    if (cylindricalZArr.size() == cylindricalRArr.size())
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
    else
    {
        throw "Input arrays must be the same size!";
    }
}

void Coil::calculateAllBFieldVerticalSINGLE(const std::vector<double> &cylindricalZArr,
                                            const std::vector<double> &cylindricalRArr,
                                            std::vector<double> &computedFieldZArr,
                                            const PrecisionArguments &usedPrecision)
{
    if (cylindricalZArr.size() == cylindricalRArr.size())
    {
        computedFieldZArr.resize(0);

        for (int i = 0; i < cylindricalZArr.size(); ++i)
        {
            double field = calculateBFieldVertical(cylindricalZArr[i], cylindricalRArr[i], usedPrecision);
            computedFieldZArr.push_back(field);
        }
    }
    else
    {
        throw "Input arrays must be the same size!";
    }
}

void Coil::calculateAllBFieldHorizontalSINGLE(const std::vector<double> &cylindricalZArr,
                                              const std::vector<double> &cylindricalRArr,
                                              std::vector<double> &computedFieldHArr,
                                              const PrecisionArguments &usedPrecision)
{
    if (cylindricalZArr.size() == cylindricalRArr.size())
    {
        computedFieldHArr.resize(0);

        for (int i = 0; i < cylindricalZArr.size(); ++i)
        {
            double field = calculateBFieldHorizontal(cylindricalZArr[i], cylindricalRArr[i], usedPrecision);
            computedFieldHArr.push_back(field);
        }
    }
    else
    {
        throw "Input arrays must be the same size!";
    }
}

void Coil::calculateAllAPotentialSINGLE(const std::vector<double> &cylindricalZArr,
                                        const std::vector<double> &cylindricalRArr,
                                        std::vector<double> &computedPotentialArr,
                                        const PrecisionArguments &usedPrecision)
{
    if (cylindricalZArr.size() == cylindricalRArr.size())
    {
        computedPotentialArr.resize(0);

        for (int i = 0; i < cylindricalZArr.size(); ++i)
        {
            double field = calculateAPotential(cylindricalZArr[i], cylindricalRArr[i], usedPrecision);
            computedPotentialArr.push_back(field);
        }
    }
    else
    {
        throw "Input arrays must be the same size!";
    }
}

