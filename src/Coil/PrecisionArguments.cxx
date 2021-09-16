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
    double lengthRoot = std::sqrt(coil.getLength());

    int totalIncrements = g_baseLayerIncrements;
    int currentIncrements;

    switch (coil.getCoilType())
    {
        case CoilType::RECTANGULAR:
        {
            totalIncrements *= g_baseLayerIncrements;
            thicknessArrayIndex = g_minPrimThicknessIncrements - 1;
            lengthArrayIndex = g_minPrimLengthIncrements - 1;
            break;
        }
        case CoilType::THIN:
        {
            thicknessArrayIndex = 0;
            thicknessArrayIndex = g_minPrimLengthIncrements - 1;
            break;
        }
        case CoilType::FLAT:
        {
            totalIncrements *= g_baseLayerIncrements;
            thicknessArrayIndex = g_minPrimThicknessIncrements - 1;
            lengthArrayIndex = 0;
            break;
        }
        case CoilType::FILAMENT:
        {
            thicknessArrayIndex = 0;
            lengthArrayIndex = 0;
            break;
        }
    }
    angularArrayIndex = g_minPrimAngularIncrements - 1;

    totalIncrements *= std::pow(2, precisionFactor.relativePrecision - 1.0);

    double angularStep, lengthStep, thicknessStep;

    do
    {
        angularStep = angularRoot /
                (blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex]);

        lengthStep = lengthRoot /
                (blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex]);

        thicknessStep = thicknessRoot /
                (blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex]);


        if (angularStep >= 0.5 * (thicknessStep + lengthStep))
            angularArrayIndex++;
        else
        {
            if (thicknessStep >= lengthStep)
                thicknessArrayIndex++;
            else
                lengthArrayIndex++;
        }

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
    const int lengthIncrements = arrSize;
    const int thicknessIncrements = arrSize;
    const int angularIncrements = arrSize;

    int lengthBlocks = 1;
    int thicknessBlocks = 1;
    int angularBlocks = 1;

    double angularRoot = std::sqrt(M_PI * (coil.getInnerRadius() + 0.5 * coil.getThickness()));
    double thicknessRoot = std::sqrt(coil.getThickness());
    double lengthRoot = std::sqrt(coil.getLength());

    int totalIncrements = (int) (arrSize * arrSize * arrSize * pow(2, precisionFactor.relativePrecision));
    int currentIncrements;

    double angularStep, thicknessStep, lengthStep;

    do
    {
        lengthStep = lengthRoot / (lengthIncrements * lengthBlocks);
        thicknessStep = thicknessRoot / (thicknessIncrements * thicknessBlocks);
        angularStep = angularRoot / (angularIncrements * angularBlocks);

        if (angularStep >= 0.5 * (thicknessStep + lengthStep))
            angularBlocks++;
        else
        {
            if (thicknessStep >= lengthStep)
                lengthBlocks++;
            else
                thicknessBlocks++;
        }

        currentIncrements = lengthBlocks*lengthIncrements * thicknessIncrements*thicknessBlocks * angularBlocks*angularIncrements;
    }
    while(currentIncrements < totalIncrements);

    #if PRINT_ENABLED
        printf("%d : %d %d %d\n", currentIncrements,
               lengthBlocks * lengthIncrements, thicknessIncrements * thicknessBlocks, angularBlocks*angularIncrements);
    #endif //PRINT_ENABLED
    return PrecisionArguments(angularBlocks, thicknessBlocks, lengthBlocks, angularIncrements, thicknessIncrements, lengthIncrements);
}