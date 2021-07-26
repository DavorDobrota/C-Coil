#include "Coil.h"
#include "LegendreMatrix.h"

#include <cmath>

namespace
{
    const double g_MiReduced = 0.0000001;
}

std::pair<double, double> Coil::calculateBField(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    double magneticFieldZ = 0.0;
    double magneticFieldH = 0.0;

    double lengthBlock = length / usedPrecision.lengthBlockCount;
    double thicknessBlock = thickness / usedPrecision.thicknessBlockCount;
    double angularBlock = M_PI / usedPrecision.angularBlockCount;

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

double Coil::calculateAPotential(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    double magneticPotential = 0.0;

    double lengthBlock = length / usedPrecision.lengthBlockCount;
    double thicknessBlock = thickness / usedPrecision.thicknessBlockCount;
    double angularBlock = M_PI / usedPrecision.angularBlockCount;

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
