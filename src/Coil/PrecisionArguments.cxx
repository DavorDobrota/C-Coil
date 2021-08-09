#include "Coil.h"
#include "LegendreMatrix.h"
#include "CoilData.h"
#include "PrecisionGlobalVars.h"

#include <cstdio>
#include <cmath>

namespace
{
    const int g_defaultLegendreOrder = 12;
    const int g_defaultBlockCount = 1;
}

PrecisionArguments::PrecisionArguments() :
        PrecisionArguments(g_defaultBlockCount, g_defaultBlockCount, g_defaultBlockCount,
                           g_defaultLegendreOrder, g_defaultLegendreOrder, g_defaultLegendreOrder) {}

PrecisionArguments::PrecisionArguments(
        int angularBlocks, int thicknessBlocks, int lengthBlocks,
        int angularIncrements, int thicknessIncrements, int lengthIncrements) :
        angularBlockCount(angularBlocks), thicknessBlockCount(thicknessBlocks),
        lengthBlockCount(lengthBlocks), angularIncrementCount(angularIncrements),
        thicknessIncrementCount(thicknessIncrements), lengthIncrementCount(lengthIncrements)
{
    if (angularIncrements > Legendre::maxLegendreOrder || angularIncrements < 1)
        PrecisionArguments::angularIncrementCount = g_defaultLegendreOrder;

    if (thicknessIncrements >= Legendre::maxLegendreOrder || thicknessIncrements < 1)
        PrecisionArguments::thicknessIncrementCount = g_defaultLegendreOrder;

    if (lengthIncrements >= Legendre::maxLegendreOrder || lengthIncrements < 1)
        PrecisionArguments::lengthIncrementCount = g_defaultLegendreOrder;


    if (angularBlocks < 1)
        PrecisionArguments::angularBlockCount = g_defaultBlockCount;

    if (thicknessBlocks < 1)
        PrecisionArguments::thicknessBlockCount = g_defaultBlockCount;

    if (lengthIncrements < 1)
        PrecisionArguments::lengthBlockCount = g_defaultBlockCount;
}

PrecisionArguments PrecisionArguments::getCoilPrecisionArgumentsCPU(const Coil &coil, PrecisionFactor precisionFactor)
{
    int lengthArrayIndex = g_minPrimLengthIncrements - 1;
    int thicknessArrayIndex = g_minPrimThicknessIncrements - 1;
    int angularArrayIndex = g_minPrimAngularIncrements - 1;

    double precisionMultiplier = std::pow(2, precisionFactor.relativePrecision - 1.0);

    int totalIncrements = 0;
    int currentIncrements;
    int caseIndex;

    if (coil.getThickness() / coil.getInnerRadius() < g_thinCoilApproximationRatio &&
        coil.getLength() / coil.getInnerRadius() < g_thinCoilApproximationRatio)
    {
        caseIndex = 1;
        totalIncrements = (int) (g_baseLayerIncrements * precisionMultiplier);
    }
    else if (coil.getThickness() / coil.getLength() < g_thinCoilApproximationRatio)
    {
        caseIndex = 2;
        totalIncrements = (int) (g_baseLayerIncrements * precisionMultiplier);
    }
    else if (coil.getLength() / coil.getThickness() < g_thinCoilApproximationRatio)
    {
        caseIndex = 3;
        totalIncrements = (int) (g_baseLayerIncrements * g_baseLayerIncrements * precisionMultiplier);
    }
    else
    {
        caseIndex = 4;
        totalIncrements = (int) (g_baseLayerIncrements * g_baseLayerIncrements * precisionMultiplier);
    }

    do
    {
        double angularStep = M_PI * g_angularWeightModifier * (coil.getInnerRadius() + coil.getThickness() * 0.5) /
                             (blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex]);

        double lengthStep = coil.getLength() /
                            (blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex]);

        double thicknessStep = coil.getThickness() /
                               (blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex]);

        switch (caseIndex)
        {
            case (1):
                thicknessArrayIndex = 0; lengthArrayIndex = 0;
                angularArrayIndex++;
                break;
            case (2):
                thicknessArrayIndex = 0;
                if (angularStep / lengthStep >= 1.0)
                    angularArrayIndex++;
                else
                    lengthArrayIndex++;
                break;
            case (3):
                lengthArrayIndex = 0;
                if (angularStep / thicknessStep >= 1.0)
                    angularArrayIndex++;
                else
                    thicknessArrayIndex++;
                break;
            default:
                if (angularStep / thicknessStep >= 1.0)
                    angularArrayIndex++;
                else
                { lengthArrayIndex++; thicknessArrayIndex++;}
        }

        currentIncrements =
                blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex] *
                blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex];
    }
    while (currentIncrements <= totalIncrements);

    #if PRINT_ENABLED
        printf("case %d - %d : %d %d %d\n", caseIndex, currentIncrements,
               blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex],
               blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex],
               blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex]);
    #endif // PRINT_ENABLED

    return PrecisionArguments(blockPrecisionCPUArray[angularArrayIndex],
                              blockPrecisionCPUArray[thicknessArrayIndex],
                              blockPrecisionCPUArray[lengthArrayIndex],
                              incrementPrecisionCPUArray[angularArrayIndex],
                              incrementPrecisionCPUArray[thicknessArrayIndex],
                              incrementPrecisionCPUArray[lengthArrayIndex]);
}

PrecisionArguments PrecisionArguments::getCoilPrecisionArgumentsGPU(const Coil &coil, PrecisionFactor precisionFactor)
{
    const int linearIncrements = arrSize;
    const int angularIncrements = arrSize;

    int linearBlocks = 1;
    int angularBlocks = 1;

    int totalIncrements = (int) pow(2, 9 + precisionFactor.relativePrecision);
    int currentIncrements;

    do
    {
        double linearStep = std::sqrt(2 * coil.getLength() * coil.getThickness()) / (linearIncrements * linearBlocks);
        double angularStep = M_PI * (coil.getInnerRadius() + coil.getThickness() * 0.5) / (angularIncrements * angularBlocks);

        if (angularStep / linearStep >= 1.0)
            angularBlocks++;
        else
            linearBlocks++;

        currentIncrements =
                linearBlocks * linearIncrements * linearIncrements * linearBlocks * angularBlocks * angularIncrements;
    }
    while(currentIncrements < totalIncrements);

    #if PRINT_ENABLED
        printf("%d : %d %d %d\n", currentIncrements,
               linearBlocks * linearIncrements, linearBlocks * linearIncrements, angularBlocks * angularIncrements);
    #endif //PRINT_ENABLED
    return PrecisionArguments(angularBlocks, linearBlocks, linearBlocks, angularIncrements, linearIncrements, linearIncrements);
}