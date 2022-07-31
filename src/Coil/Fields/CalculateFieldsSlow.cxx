#include "Coil.h"
#include "LegendreMatrix.h"

#define _USE_MATH_DEFINES
#include <math.h>


namespace
{
    const double g_MiReduced = 0.0000001;
}

double Coil::calculateAPotentialSlow(double zCoord, double rCoord, const PrecisionArguments &usedPrecision) const
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

    std::vector<std::vector<double>> cosPhiPrecomputeMat(usedPrecision.angularBlockCount);

    for (int indBlockPhi = 0; indBlockPhi < usedPrecision.angularBlockCount; ++indBlockPhi)
    {
        double blockPositionPhi = angularBlock * (indBlockPhi + 0.5);
        cosPhiPrecomputeMat[indBlockPhi].resize(usedPrecision.angularIncrementCount);

        for (int incPhi = 0; incPhi <= angularIncrements; ++incPhi)
        {
            double incrementPositionFi = blockPositionPhi +
                                         (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][incPhi];
            cosPhiPrecomputeMat[indBlockPhi][incPhi] = std::cos(incrementPositionFi);
        }
    }

    for (int indBlockL = 0; indBlockL < usedPrecision.lengthBlockCount; ++indBlockL)
    {
        for (int indBlockT = 0; indBlockT < usedPrecision.thicknessBlockCount; ++indBlockT)
        {
            double blockPositionL = (-1.0) * (length * 0.5) + lengthBlock * (indBlockL + 0.5);
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

                    double tempConstA = 2.0 * incrementPositionT * rCoord;
                    double tempConstB = incrementPositionT * incrementPositionT + rCoord * rCoord +
                                        (incrementPositionL + zCoord) * (incrementPositionL + zCoord);

                    for (int indBlockFi = 0; indBlockFi < usedPrecision.angularBlockCount; ++indBlockFi)
                    {
                        for (int incFi = 0; incFi <= angularIncrements; ++incFi)
                        {
                            double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incFi];
                            double cosinePhi = cosPhiPrecomputeMat[indBlockFi][incFi];

                            magneticPotential += constant * incrementWeightS * incrementWeightFi *
                                    (incrementPositionT * cosinePhi) / std::sqrt(tempConstB - tempConstA * cosinePhi);
                        }
                    }
                }
            }
        }
    }
    return magneticPotential;
}

std::pair<double, double> Coil::calculateBFieldSlow(double zCoord, double rCoord, const PrecisionArguments &usedPrecision) const
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

    std::vector<std::vector<double>> cosPhiPrecomputeMat(usedPrecision.angularBlockCount);

    for (int indBlockPhi = 0; indBlockPhi < usedPrecision.angularBlockCount; ++indBlockPhi)
    {
        double blockPositionPhi = angularBlock * (indBlockPhi + 0.5);
        cosPhiPrecomputeMat[indBlockPhi].resize(usedPrecision.angularIncrementCount);

        for (int incPhi = 0; incPhi <= angularIncrements; ++incPhi)
        {
            double incrementPositionFi = blockPositionPhi +
                                         (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][incPhi];
            cosPhiPrecomputeMat[indBlockPhi][incPhi] = std::cos(incrementPositionFi);
        }
    }

    for (int indBlockL = 0; indBlockL < usedPrecision.lengthBlockCount; ++indBlockL)
    {
        for (int indBlockT = 0; indBlockT < usedPrecision.thicknessBlockCount; ++indBlockT)
        {
            double blockPositionL = (-1.0) * (length * 0.5) + lengthBlock * (indBlockL + 0.5);
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
                    double tempConstB = incrementPositionT * (incrementPositionL + zCoord);
                    double tempConstC = incrementPositionT * rCoord;
                    double tempConstD = tempConstA + rCoord * rCoord +
                                        (incrementPositionL + zCoord) * (incrementPositionL + zCoord);
                    double tempConstE = constant * incrementWeightS;

                    for (int indBlockPhi = 0; indBlockPhi < usedPrecision.angularBlockCount; ++indBlockPhi)
                    {
                        for (int incPhi = 0; incPhi <= angularIncrements; ++incPhi)
                        {
                            double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incPhi];
                            double cosinePhi = cosPhiPrecomputeMat[indBlockPhi][incPhi];

                            double tempConstF = 2.0 * tempConstC * cosinePhi;
                            double tempConstH = (tempConstD - tempConstF) * std::sqrt(tempConstD - tempConstF);
                            double tempConstG = tempConstE * incrementWeightFi / tempConstH;

                            magneticFieldZ += tempConstG * (tempConstA - tempConstC * cosinePhi);
                            magneticFieldH += tempConstG * (tempConstB * cosinePhi);
                        }
                    }
                }
            }
        }
    }
    return {magneticFieldH, magneticFieldZ};
}

std::vector<double> Coil::calculateBGradientSlow(double zCoord, double rCoord, const PrecisionArguments &usedPrecision) const
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

    std::vector<std::vector<double>> cosPhiPrecomputeMat(usedPrecision.angularBlockCount);

    for (int indBlockPhi = 0; indBlockPhi < usedPrecision.angularBlockCount; ++indBlockPhi)
    {
        double blockPositionPhi = angularBlock * (indBlockPhi + 0.5);
        cosPhiPrecomputeMat[indBlockPhi].resize(usedPrecision.angularIncrementCount);

        for (int incPhi = 0; incPhi <= angularIncrements; ++incPhi)
        {
            double incrementPositionFi = blockPositionPhi +
                                         (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][incPhi];
            cosPhiPrecomputeMat[indBlockPhi][incPhi] = std::cos(incrementPositionFi);
        }
    }

    for (int indBlockL = 0; indBlockL < usedPrecision.lengthBlockCount; ++indBlockL)
    {
        for (int indBlockT = 0; indBlockT < usedPrecision.thicknessBlockCount; ++indBlockT)
        {
            double blockPositionL = (-1.0) * (length * 0.5) + lengthBlock * (indBlockL + 0.5);
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
                    double tempConstB = rCoord * rCoord;
                    double tempConstC = (incrementPositionL + zCoord) * (incrementPositionL + zCoord);
                    double tempConstD = rCoord * incrementPositionT;
                    double tempConstE = incrementPositionT * (incrementPositionL + zCoord);
                    double tempConstF = 1.0 / rCoord;

                    double tempConstG = 2.0 * tempConstA + 2.0 * tempConstB - tempConstC;
                    double tempConstH = tempConstA + tempConstB + tempConstC;
                    double tempConstI = constant * incrementWeightS;

                    for (int indBlockPhi = 0; indBlockPhi < usedPrecision.angularBlockCount; ++indBlockPhi)
                    {
                        for (int incPhi = 0; incPhi <= angularIncrements; ++incPhi)
                        {
                            double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incPhi];

                            double cosinePhi = cosPhiPrecomputeMat[indBlockPhi][incPhi];

                            double tempConstJ = tempConstH - 2.0 * tempConstD * cosinePhi;
                            double tempConstK = tempConstJ * tempConstJ * std::sqrt(tempConstJ);

                            double tempConstX = tempConstI * incrementWeightFi / (tempConstJ * std::sqrt(tempConstJ));
                            double tempConstY = tempConstI * incrementWeightFi / tempConstK;

                            bufferValueRPhi += tempConstX * (tempConstF * tempConstE * cosinePhi);
                            bufferValueRR += tempConstY * (-3.0 * tempConstE * (rCoord - incrementPositionT * cosinePhi)) * cosinePhi;
                            bufferValueZZ += tempConstY * (-3.0 * tempConstE * (incrementPositionT - rCoord * cosinePhi));
                            bufferValueRZ += tempConstY *
                                    (incrementPositionT * (tempConstG - tempConstD * cosinePhi) * cosinePhi - 3.0 * tempConstA * rCoord);
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

