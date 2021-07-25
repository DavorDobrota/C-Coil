#include "Coil.h"
#include "LegendreMatrix.h"

#include <cmath>

namespace
{
    const double g_MiReduced = 0.0000001;
}

void Coil::calculateSelfInductance(PrecisionFactor precisionFactor)
{
    // TODO - new solution applied: still convergence is slow and not entirely tested, an approximation, future improvement
    selfInductance = 0.0;

    auto precisionArguments = CoilPairArguments::getSelfInductanceArguments(*this, precisionFactor);
    PrecisionArguments combinedArguments = precisionArguments.primaryPrecision;

    double lengthBlock = length / combinedArguments.lengthBlockCount;
    double thicknessBlock = thickness / combinedArguments.thicknessBlockCount;
    double angularBlock = M_PI / combinedArguments.angularBlockCount;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int lengthIncrements = combinedArguments.lengthIncrementCount - 1;
    int thicknessIncrements = combinedArguments.thicknessIncrementCount - 1;
    int angularIncrements = combinedArguments.angularIncrementCount - 1;

    // multiplication by 2 because cosine is an even function and by 0.125 for a triple change of interval (3 times 1/2)
    double constant = g_MiReduced * currentDensity * lengthBlock * thicknessBlock * angularBlock * 2 * 0.125;

    for (int zBlock = 0; zBlock < combinedArguments.lengthBlockCount; zBlock++)
    {
        for (int rBlock = 0; rBlock < combinedArguments.thicknessBlockCount; rBlock++)
        {
            double zBlockPosition = (-1) * (length * 0.5) + lengthBlock * (zBlock + 0.5);
            double rBlockPosition = innerRadius + thicknessBlock * (rBlock + 0.5);

            for (int zIndex = 0; zIndex <= lengthIncrements; ++zIndex)
            {
                for (int rIndex = 0; rIndex <= thicknessIncrements; ++rIndex)
                {
                    double incrementPositionZ = zBlockPosition +
                                                (lengthBlock * 0.5) * Legendre::positionMatrix[lengthIncrements][zIndex];
                    double incrementPositionR = rBlockPosition +
                                                (thicknessBlock * 0.5) * Legendre::positionMatrix[thicknessIncrements][rIndex];

                    double potential = 0.0;

                    for (int indBlockL = 0; indBlockL < combinedArguments.lengthBlockCount; ++indBlockL)
                    {
                        for (int indBlockT = 0; indBlockT < combinedArguments.thicknessBlockCount; ++indBlockT)
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

                                    if (!(std::fabs(incrementPositionL - incrementPositionZ) / length < 1e-14 &&
                                          std::fabs(incrementPositionT - incrementPositionR) / thickness < 1e-14))
                                    {
                                        double incrementWeightS =
                                                Legendre::weightsMatrix[lengthIncrements][incL] *
                                                Legendre::weightsMatrix[thicknessIncrements][incT];

                                        double tempConstA = incrementPositionT;
                                        double tempConstB = incrementPositionT * incrementPositionR;
                                        double tempConstC =
                                                incrementPositionT * incrementPositionT +
                                                incrementPositionR * incrementPositionR +
                                                (incrementPositionL + incrementPositionZ) *
                                                (incrementPositionL + incrementPositionZ);

                                        for (int indBlockFi = 0; indBlockFi < combinedArguments.angularBlockCount; ++indBlockFi)
                                        {
                                            double blockPositionFi = angularBlock * (indBlockFi + 0.5);

                                            for (int incFi = 0; incFi <= angularIncrements; ++incFi)
                                            {
                                                double incrementPositionFi = blockPositionFi +
                                                                             (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][incFi];

                                                double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incFi];

                                                double cosinePhi = cos(incrementPositionFi);

                                                potential += constant * incrementWeightS * incrementWeightFi *
                                                             (tempConstA * cosinePhi) /sqrt(tempConstC - 2*tempConstB * cosinePhi);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    selfInductance += 2 * M_PI * rBlockPosition * potential * 0.25 *
                                      Legendre::weightsMatrix[lengthIncrements][zIndex] *
                                      Legendre::weightsMatrix[thicknessIncrements][rIndex];
                }
            }
        }
    }
    selfInductance *= (numOfTurns / current);
    calculateImpedance();
}

void Coil::calculateApproximateSelfInductance(PrecisionFactor precisionFactor, ComputeMethod method)
{
    // this is a simplified solution but produces a significant error because of the inherent divergent integral
    // however it can be accelerated beyond single thread and thus may prove useful in some instances
    auto precisionArguments = CoilPairArguments::getSelfInductanceArguments(*this, precisionFactor);
    selfInductance = computeMutualInductance(*this, *this, 0.0, precisionArguments, method);
    calculateImpedance();
}

double Coil::computeAndSetSelfInductance(PrecisionFactor precisionFactor)
{
    calculateSelfInductance(precisionFactor);
    return selfInductance;
}

double Coil::computeAndSetApproximateSelfInductance(PrecisionFactor precisionFactor, ComputeMethod method)
{
    calculateApproximateSelfInductance(precisionFactor, method);
    return selfInductance;
}
