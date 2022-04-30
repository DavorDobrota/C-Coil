#include "Coil.h"
#include "LegendreMatrix.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include "immintrin.h"


namespace
{
    const double g_MiReduced = 0.0000001;
}

double Coil::calculateAPotentialFast(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    double magneticPotential = 0.0;

    double thicknessBlock = thickness / usedPrecision.thicknessBlockCount;
    double angularBlock = M_PI / usedPrecision.angularBlockCount;

    // initialising precompute array
    const int numPhiIncrements = usedPrecision.angularBlockCount * usedPrecision.angularIncrementCount;

    #if defined(__GNUC__)
        if(numPhiIncrements > 20480)
            throw "Number of increments too great (the maximum is 20480)";
        double cosPhiPrecomputeArr[numPhiIncrements];
        precomputeCosPhi(usedPrecision.angularBlockCount, usedPrecision.angularIncrementCount, cosPhiPrecomputeArr);
    #else
        std::vector<double> cosPhiPrecomputeArr(numPhiIncrements);
        precomputeCosPhi(usedPrecision.angularBlockCount, usedPrecision.angularIncrementCount, &cosPhiPrecomputeArr.front());
    #endif

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int thicknessIncrements = usedPrecision.thicknessIncrementCount - 1;
    int angularIncrements = usedPrecision.angularIncrementCount - 1;

    // multiplication by 2 because cosine is an even function and by 0.25 for a double change of interval (2 times 1/2)
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
                for (int incFi = 0; incFi <= angularIncrements; ++incFi)
                {
                    double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incFi];

                    int arrPos = indBlockFi * (angularIncrements + 1) + incFi;
                    double cosinePhi = cosPhiPrecomputeArr[arrPos];
                    double tempConstC = 1.0 / std::sqrt(tempConstB - tempConstA * cosinePhi);

                    double tempConstD1 = topEdge * tempConstC;
                    double tempConstD2 = bottomEdge * tempConstC;

                    double tempConstE1 = std::sqrt(tempConstD1 * tempConstD1 + 1.0);
                    double tempConstE2 = std::sqrt(tempConstD2 * tempConstD2 + 1.0);

                    double tempConstF = std::log((tempConstE1 + tempConstD1) / (tempConstE2 + tempConstD2));

                    magneticPotential += constant * incrementWeightT * incrementWeightFi * incrementPositionT *
                                         cosinePhi * tempConstF;
                }
            }
        }
    }
    return magneticPotential;
}

std::pair<double, double> Coil::calculateBFieldFast(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    double magneticFieldZ = 0.0;
    double magneticFieldH = 0.0;

    double thicknessBlock = thickness / usedPrecision.thicknessBlockCount;
    double angularBlock = M_PI / usedPrecision.angularBlockCount;

    // initialising precompute array
    const int numPhiIncrements = usedPrecision.angularBlockCount * usedPrecision.angularIncrementCount;

    #if defined(__GNUC__)
        if(numPhiIncrements > 20480)
            throw "Number of increments too great (the maximum is 20480)";
        double cosPhiPrecomputeArr[numPhiIncrements];
        precomputeCosPhi(usedPrecision.angularBlockCount, usedPrecision.angularIncrementCount, cosPhiPrecomputeArr);
    #else
        std::vector<double> cosPhiPrecomputeArr(numPhiIncrements);
        precomputeCosPhi(usedPrecision.angularBlockCount, usedPrecision.angularIncrementCount, &cosPhiPrecomputeArr.front());
    #endif

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int thicknessIncrements = usedPrecision.thicknessIncrementCount - 1;
    int angularIncrements = usedPrecision.angularIncrementCount - 1;

    // multiplication by 2 because cosine is an even function and by 0.25 for a double change of interval (2 times 1/2)
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
                for (int incFi = 0; incFi <= angularIncrements; ++incFi)
                {
                    double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incFi];
                    int arrPos = indBlockFi * (angularIncrements + 1) + incFi;
                    double cosinePhi = cosPhiPrecomputeArr[arrPos];

                    double tempConstF = 2 * tempConstB * cosinePhi;
                    double tempConstG1 = 1.0 / std::sqrt(tempConstD1 - tempConstF);
                    double tempConstG2 = 1.0 / std::sqrt(tempConstD2 - tempConstF);

                    double tempConstH = tempConstE * incrementWeightFi;

                    magneticFieldH += tempConstH * incrementPositionT * cosinePhi * (tempConstG2 - tempConstG1);
                    magneticFieldZ += tempConstH *
                            ((tempConstA - 0.5 * tempConstF) / (tempConstC - tempConstF)) *
                            (topEdge * tempConstG1 - bottomEdge * tempConstG2);
                }
            }
        }
    }
    return {magneticFieldH, magneticFieldZ};
}

std::vector<double> Coil::calculateBGradientFast(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    double bufferValueRPhi = 0.0;
    double bufferValueRR = 0.0;
    double bufferValueZZ = 0.0;
    double bufferValueRZ = 0.0;

    double thicknessBlock = thickness / usedPrecision.thicknessBlockCount;
    double angularBlock = M_PI / usedPrecision.angularBlockCount;

    // initialising precompute array
    const int numPhiIncrements = usedPrecision.angularBlockCount * usedPrecision.angularIncrementCount;

    #if defined(__GNUC__)
        if(numPhiIncrements > 20480)
            throw "Number of increments too great (the maximum is 20480)";
        double cosPhiPrecomputeArr[numPhiIncrements];
        precomputeCosPhi(usedPrecision.angularBlockCount, usedPrecision.angularIncrementCount, cosPhiPrecomputeArr);
    #else
        std::vector<double> cosPhiPrecomputeArr(numPhiIncrements);
        precomputeCosPhi(usedPrecision.angularBlockCount, usedPrecision.angularIncrementCount, &cosPhiPrecomputeArr.front());
    #endif

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int thicknessIncrements = usedPrecision.thicknessIncrementCount - 1;
    int angularIncrements = usedPrecision.angularIncrementCount - 1;

    // multiplication by 2 because cosine is an even function and by 0.25 for a double change of interval (2 times 1/2)
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
            double tempConstB = rPolar * rPolar;
            double tempConstC = incrementPositionT * rPolar;

            double tempConstD = tempConstA + tempConstB;
            double tempConstE = tempConstA * tempConstA + tempConstB * tempConstB;
            double tempConstF = tempConstC * tempConstC;

            double tempConstG1 = tempConstD + topEdge * topEdge;
            double tempConstG2 = tempConstD + bottomEdge * bottomEdge;

            double tempConstI = constant * incrementWeightT;

            for (int indBlockFi = 0; indBlockFi < usedPrecision.angularBlockCount; ++indBlockFi)
            {
                for (int incFi = 0; incFi <= angularIncrements; ++incFi)
                {
                    double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incFi];

                    int arrPos = indBlockFi * (angularIncrements + 1) + incFi;
                    double cosinePhi = cosPhiPrecomputeArr[arrPos];

                    double cosinePhi2 = cosinePhi * cosinePhi;
                    double phiExpression = 2 * tempConstC * cosinePhi;

                    double tempConstJ1 = tempConstG1 - phiExpression;
                    double tempConstJ2 = tempConstG2 - phiExpression;

                    double tempConstK1 = std::sqrt(tempConstJ1);
                    double tempConstK2 = std::sqrt(tempConstJ2);

                    double tempConstL1 = 1 / (tempConstJ1 * tempConstK1);
                    double tempConstL2 = 1 / (tempConstJ2 * tempConstK2);

                    double tempConstM = tempConstD - phiExpression;
                    double tempConstN =
                            2 * tempConstF * cosinePhi * (cosinePhi2 + 2.0) -
                            tempConstC * (3 * cosinePhi2 + 1) * tempConstD + cosinePhi * tempConstE;
                    double tempConstO = cosinePhi * tempConstD - 2*tempConstC;

                    double tempConstZ = tempConstI * incrementWeightFi;

                    bufferValueRPhi += tempConstZ * (incrementPositionT * cosinePhi / rPolar) * (1 / tempConstK2 - 1 / tempConstK1);
                    bufferValueRR += tempConstZ * (tempConstC - tempConstA * cosinePhi) * cosinePhi * (tempConstL1 - tempConstL2);
                    bufferValueZZ += tempConstZ * (tempConstA - tempConstC * cosinePhi) * (tempConstL1 - tempConstL2);
                    bufferValueRZ += tempConstZ * incrementPositionT / (tempConstM * tempConstM) *
                            (topEdge * tempConstL1 * (tempConstO * tempConstJ1 + tempConstN) -
                            bottomEdge * tempConstL2 * (tempConstO * tempConstJ2 + tempConstN));

                    //printf("%.15g: %.15g %.15g %.15g %.15g\n", incrementPositionFi, bufferValueRPhi, bufferValueRR, bufferValueZZ, bufferValueRZ);
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

