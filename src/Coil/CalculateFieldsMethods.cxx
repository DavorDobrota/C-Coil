#include "Coil.h"
#include "LegendreMatrix.h"

#include <cmath>
#include <cstdio>

namespace
{
    const double g_MiReduced = 0.0000001;
}

double Coil::calculateAPotential(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    double magneticPotential = 0.0;

    double thicknessBlock = thickness / usedPrecision.thicknessBlockCount;
    double angularBlock = M_PI / usedPrecision.angularBlockCount;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int lengthIncrements = usedPrecision.lengthIncrementCount - 1;
    int thicknessIncrements = usedPrecision.thicknessIncrementCount - 1;
    int angularIncrements = usedPrecision.angularIncrementCount - 1;

    // multiplication by 2 because cosine is an even function and by 0.125 for a double change of interval (2 times 1/2)
    double constant = g_MiReduced * currentDensity * thicknessBlock * angularBlock * 2 * 0.25;

    double topEdge = zAxis + length * 0.5;
    double bottomEdge = zAxis - length * 0.5;

    for (int indBlockT = 0; indBlockT < usedPrecision.thicknessBlockCount; ++indBlockT)
    {
        double blockPositionT = innerRadius + thicknessBlock * (indBlockT + 0.5);

        for (int incT = 0; incT <= thicknessIncrements; ++incT)
        {
            double incrementPositionT = blockPositionT +
                    (thicknessBlock * 0.5) * Legendre::positionMatrix[thicknessIncrements][incT];

            double incrementWeightT = Legendre::weightsMatrix[thicknessIncrements][incT];

            double tempConstA = 2 * incrementPositionT * rPolar;
            double tempConstB = incrementPositionT * incrementPositionT + rPolar * rPolar;

            for (int indBlockFi = 0; indBlockFi < usedPrecision.angularBlockCount; ++indBlockFi)
            {
                double blockPositionFi = angularBlock * (indBlockFi + 0.5);

                for (int incFi = 0; incFi <= angularIncrements; ++incFi)
                {
                    double incrementPositionFi = blockPositionFi +
                            (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][incFi];

                    double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incFi];

                    double cosinePhi = std::cos(incrementPositionFi);
                    double tempConstC = 1 / std::sqrt(tempConstB - tempConstA * cosinePhi);

                    magneticPotential += constant * incrementWeightT * incrementWeightFi *incrementPositionT *
                            cosinePhi * (std::asinh(topEdge * tempConstC) - std::asinh(bottomEdge * tempConstC));
                }
            }
        }
    }
    return magneticPotential;
}

std::pair<double, double> Coil::calculateBField(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    double magneticFieldZ = 0.0;
    double magneticFieldH = 0.0;

    double thicknessBlock = thickness / usedPrecision.thicknessBlockCount;
    double angularBlock = M_PI / usedPrecision.angularBlockCount;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int thicknessIncrements = usedPrecision.thicknessIncrementCount - 1;
    int angularIncrements = usedPrecision.angularIncrementCount - 1;

    // multiplication by 2 because cosine is an even function and by 0.125 for a double change of interval (2 times 1/2)
    double constant = g_MiReduced * currentDensity * thicknessBlock * angularBlock * 2 * 0.25;

    double topEdge = zAxis + length * 0.5;
    double bottomEdge = zAxis - length * 0.5;

    for (int indBlockT = 0; indBlockT < usedPrecision.thicknessBlockCount; ++indBlockT)
    {
        double blockPositionT = innerRadius + thicknessBlock * (indBlockT + 0.5);

        for (int incT = 0; incT <= thicknessIncrements; ++incT)
        {
            double incrementPositionT = blockPositionT +
                                        (thicknessBlock * 0.5) * Legendre::positionMatrix[thicknessIncrements][incT];

            double incrementWeightT = Legendre::weightsMatrix[thicknessIncrements][incT];

            double tempConstA = incrementPositionT * incrementPositionT;
            double tempConstB = incrementPositionT * rPolar;
            double tempConstC = tempConstA + rPolar * rPolar;

            double tempConstD1 = topEdge * topEdge + tempConstC;
            double tempConstD2 = bottomEdge * bottomEdge + tempConstC;

            double tempConstE = constant * incrementWeightT;

            for (int indBlockFi = 0; indBlockFi < usedPrecision.angularBlockCount; ++indBlockFi)
            {
                double blockPositionFi = angularBlock * (indBlockFi + 0.5);

                for (int incFi = 0; incFi <= angularIncrements; ++incFi)
                {
                    double incrementPositionFi = blockPositionFi +
                                                 (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][incFi];

                    double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incFi];

                    double cosinePhi = std::cos(incrementPositionFi);

                    double tempConstF = 2 * tempConstB * cosinePhi;
                    double tempConstG1 = 1 / std::sqrt(tempConstD1 - tempConstF);
                    double tempConstG2 = 1 / std::sqrt(tempConstD2 - tempConstF);

                    double tempConstH = tempConstE * incrementWeightFi;

                    magneticFieldH += tempConstH * incrementPositionT * cosinePhi * (tempConstG2 - tempConstG1);
                    magneticFieldZ += tempConstH *
                            ((tempConstA - 0.5 * tempConstF) / (tempConstC - tempConstF)) *
                            (topEdge * tempConstG1 - bottomEdge * tempConstG2);
                }
            }
        }
    }
    return std::make_pair(magneticFieldH, magneticFieldZ);
}

std::vector<double> Coil::calculateBGradient(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    double bufferValueRPhi = 0.0;
    double bufferValueRR = 0.0;
    double bufferValueZZ = 0.0;
    double bufferValueRZ = 0.0;

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
                    double tempConstB = rPolar * rPolar;
                    double tempConstC = (incrementPositionL + zAxis) * (incrementPositionL + zAxis);
                    double tempConstD = rPolar * incrementPositionT;
                    double tempConstE = incrementPositionT * (incrementPositionL + zAxis);
                    double tempConstF = 1 / rPolar;

                    double tempConstG = 2 * tempConstA + 2 * tempConstB - tempConstC;
                    double tempConstH = tempConstA + tempConstB + tempConstC;
                    double tempConstI = constant * incrementWeightS;

                    for (int indBlockFi = 0; indBlockFi < usedPrecision.angularBlockCount; ++indBlockFi)
                    {
                        double blockPositionFi = angularBlock * (indBlockFi + 0.5);

                        for (int incFi = 0; incFi <= angularIncrements; ++incFi)
                        {
                            double incrementPositionFi = blockPositionFi +
                                    (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][incFi];

                            double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incFi];

                            double cosinePhi = std::cos(incrementPositionFi);

                            double tempConstJ = tempConstH - 2 * tempConstD * cosinePhi;
                            double tempConstK = tempConstJ * tempConstJ * std::sqrt(tempConstJ);

                            double tempConstX = tempConstI * incrementWeightFi / (tempConstJ * std::sqrt(tempConstJ));
                            double tempConstY = tempConstI * incrementWeightFi / tempConstK;

                            bufferValueRPhi += tempConstX * (tempConstF * tempConstE * cosinePhi);
                            bufferValueRR += tempConstY * (-3 * tempConstE * (rPolar - incrementPositionT * cosinePhi)) * cosinePhi;
                            bufferValueZZ += tempConstY * (-3 * tempConstE * (incrementPositionT - rPolar * cosinePhi));
                            bufferValueRZ += tempConstY *
                                    (incrementPositionT * (tempConstG - tempConstD * cosinePhi) * cosinePhi -3 * tempConstA * rPolar);
                        }
                    }
                }
            }
        }
    }
    std::vector<double> bufferValues(4);

    bufferValues[0] = bufferValueRPhi;
    bufferValues[1] = bufferValueRR;
    bufferValues[2] = bufferValueRZ;
    bufferValues[3] = bufferValueZZ;

    return bufferValues;
}
