#include "Coil.h"
#include "LegendreMatrix.h"

#define _USE_MATH_DEFINES
#include <math.h>


namespace
{
    const double g_MiReduced = 0.0000001;
}

double Coil::calculateAPotentialSlow(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
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

                    double tempConstA = incrementPositionT;
                    double tempConstB = incrementPositionT * rPolar;
                    double tempConstC = incrementPositionT * incrementPositionT + rPolar * rPolar +
                            (incrementPositionL + zAxis) * (incrementPositionL + zAxis);

                    for (int indBlockFi = 0; indBlockFi < usedPrecision.angularBlockCount; ++indBlockFi)
                    {
                        for (int incFi = 0; incFi <= angularIncrements; ++incFi)
                        {
                            double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incFi];
                            double cosinePhi = cosPhiPrecomputeMat[indBlockFi][incFi];

                            magneticPotential += constant * incrementWeightS * incrementWeightFi *
                                    (tempConstA * cosinePhi) / std::sqrt(tempConstC - 2*tempConstB * cosinePhi);
                        }
                    }
                }
            }
        }
    }
    return magneticPotential;
}

std::pair<double, double> Coil::calculateBFieldSlow(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
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
                    double tempConstB = incrementPositionT * (incrementPositionL + zAxis);
                    double tempConstC = incrementPositionT * rPolar;
                    double tempConstD = tempConstA + rPolar * rPolar +
                            (incrementPositionL + zAxis) * (incrementPositionL + zAxis);
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

std::vector<double> Coil::calculateBGradientSlow(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
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
                    double tempConstB = rPolar * rPolar;
                    double tempConstC = (incrementPositionL + zAxis) * (incrementPositionL + zAxis);
                    double tempConstD = rPolar * incrementPositionT;
                    double tempConstE = incrementPositionT * (incrementPositionL + zAxis);
                    double tempConstF = 1.0 / rPolar;

                    double tempConstG = 2.0 * tempConstA + 2.0 * tempConstB - tempConstC;
                    double tempConstH = tempConstA + tempConstB + tempConstC;
                    double tempConstI = constant * incrementWeightS;

                    for (int indBlockPhi = 0; indBlockPhi < usedPrecision.angularBlockCount; ++indBlockPhi)
                    {
                        for (int incPhi = 0; incPhi <= angularIncrements; ++incPhi)
                        {
                            double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incPhi];

                            double cosinePhi = cosPhiPrecomputeMat[indBlockPhi][incPhi];

                            double tempConstJ = tempConstH - 2 * tempConstD * cosinePhi;
                            double tempConstK = tempConstJ * tempConstJ * std::sqrt(tempConstJ);

                            double tempConstX = tempConstI * incrementWeightFi / (tempConstJ * std::sqrt(tempConstJ));
                            double tempConstY = tempConstI * incrementWeightFi / tempConstK;

                            bufferValueRPhi += tempConstX * (tempConstF * tempConstE * cosinePhi);
                            bufferValueRR += tempConstY * (-3.0 * tempConstE * (rPolar - incrementPositionT * cosinePhi)) * cosinePhi;
                            bufferValueZZ += tempConstY * (-3.0 * tempConstE * (incrementPositionT - rPolar * cosinePhi));
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

