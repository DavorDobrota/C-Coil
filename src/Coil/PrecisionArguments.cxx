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
    int lengthArrayIndex, thicknessArrayIndex, angularArrayIndex;

    double angularRoot = std::sqrt(M_PI * (coil.getInnerRadius() + 0.5 * coil.getThickness()));
    double thicknessRoot = std::sqrt(coil.getThickness());

    int totalIncrements = 1;
    int currentIncrements;

    totalIncrements *= g_baseLayerIncrements;

    switch (coil.getCoilType())
    {
        case CoilType::RECTANGULAR:
        {
            totalIncrements *= g_baseLayerIncrements;
            thicknessArrayIndex = g_minPrimThicknessIncrements - 1;
            break;
        }
        case CoilType::THIN:
        {
            thicknessArrayIndex = 0;
            break;
        }
        case CoilType::FLAT:
        {
            thicknessArrayIndex = g_minPrimThicknessIncrements - 1;
            totalIncrements *= g_baseLayerIncrements;
            break;
        }
        case CoilType::FILAMENT:
        {
            thicknessArrayIndex = 0;
            break;
        }
    }
    angularArrayIndex = g_minPrimAngularIncrements - 1;
    lengthArrayIndex = 0;

    do
    {
        double angularStep = angularRoot /
                             (blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex]);

        double thicknessStep = thicknessRoot /
                (blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex]);

        if (angularStep >= thicknessStep)
            angularArrayIndex++;
        else
            thicknessArrayIndex++;

        currentIncrements =
                blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex] *
                blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex];
    }
    while (currentIncrements <= totalIncrements);

    #if PRINT_ENABLED
        printf("%d : %d %d %d\n", currentIncrements,
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