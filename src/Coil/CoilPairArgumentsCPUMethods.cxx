#include "Coil.h"

#include <cmath>
#include <cstdio>


CoilPairArguments CoilPairArguments::calculateCoilPairArgumentsCPU(const Coil &primary, const Coil &secondary,
                                                                   PrecisionFactor precisionFactor, bool zAxisCase)
{
    int primLengthArrayIndex, primThicknessArrayIndex, primAngularArrayIndex;
    int secLengthArrayIndex, secThicknessArrayIndex, secAngularArrayIndex;

    int totalIncrements = 1;
    int currentIncrements;

    double primAngularRoot = std::sqrt(M_PI * (primary.getInnerRadius() + 0.5 * primary.getThickness()));
    double primThicknessRoot = std::sqrt(primary.getThickness());
    double primLengthRoot = std::sqrt(primary.getLength());

    double secAngularRoot = std::sqrt(2 * M_PI * (secondary.getInnerRadius() + 0.5 * secondary.getThickness()));
    double secThicknessRoot = std::sqrt(secondary.getThickness());
    double secLengthRoot = std::sqrt(secondary.getLength());


    // multiplying for secondary angular increments, if present
    if (!zAxisCase)
        totalIncrements *= g_baseLayerIncrements;

    switch (primary.getCoilType())
    {
        case CoilType::RECTANGULAR:
        {
            totalIncrements *= g_baseLayerIncrements;
            primThicknessArrayIndex = g_minPrimThicknessIncrements - 1;
            break;
        }
        case CoilType::THIN:
        {
            primThicknessArrayIndex = 0;
            break;
        }
        case CoilType::FLAT:
        {
            primThicknessArrayIndex = g_minPrimThicknessIncrements - 1;
            totalIncrements *= g_baseLayerIncrements;
            break;
        }
        case CoilType::FILAMENT:
        {
            primThicknessArrayIndex = 0;
            break;
        }
    }
    primAngularArrayIndex = g_minPrimAngularIncrements - 1;

    double secDimensionRatio = secLengthRoot / secThicknessRoot;
    switch (secondary.getCoilType())
    {
        case CoilType::RECTANGULAR:
        {
            totalIncrements *= g_baseLayerIncrements * g_baseLayerIncrements;
            secLengthArrayIndex = g_minSecLengthIncrements - 1;
            secThicknessArrayIndex = g_minSecThicknessIncrements - 1;
            break;
        }
        case CoilType::THIN:
        {
            totalIncrements *= g_baseLayerIncrements;
            secThicknessArrayIndex = 0;
            secLengthArrayIndex = g_minSecLengthIncrements - 1;
            break;
        }
        case CoilType::FLAT:
        {
            totalIncrements *= g_baseLayerIncrements;
            secLengthArrayIndex = 0;
            secThicknessArrayIndex = g_minSecAngularIncrements - 1;
            break;
        }
        case CoilType::FILAMENT:
        {
            secLengthArrayIndex = 0;
            secThicknessArrayIndex = 0;
            break;
        }
    }
    secAngularArrayIndex = g_minSecAngularIncrements - 1;

    totalIncrements = (int) (totalIncrements * std::pow(2, precisionFactor.relativePrecision - 1));

    double primAngularStep, primThicknessStep;
    double secAngularStep, secLengthStep, secThicknessStep, secLinearStep;

    do
    {
        primAngularStep = g_primAngularWeightModifier * primAngularRoot /
                (blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex]);

        primThicknessStep = g_primLinearWeightModifier * primThicknessRoot /
                (blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex]);

        secAngularStep = g_secAngularWeightModifier * secAngularRoot /
                (blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);

        secLengthStep = g_secLinearWeightModifier * secLengthRoot /
                (blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex]);

        secThicknessStep = g_secLinearWeightModifier * secThicknessRoot /
                (blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

        secLinearStep = 0.5 * (secLengthStep + secThicknessStep);

        if (primAngularStep + primThicknessStep >= secAngularStep + secLinearStep)
        {
            if (primAngularStep >= primThicknessStep)
                primAngularArrayIndex++;
            else
                primThicknessArrayIndex++;
        }
        else
        {
            if (secAngularStep >= secLinearStep)
                secAngularArrayIndex++;
            else
            {
                if (secLengthStep >= secThicknessStep)
                    secLengthArrayIndex++;
                else
                    secThicknessArrayIndex++;
            }
        }

        currentIncrements =
                blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex] *
                blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex] *
                blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex] *
                blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex];

        if (!zAxisCase)
            currentIncrements *= blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex];
    }
    while (currentIncrements < totalIncrements);

    // length index chosen only for purpose of data representation
    primLengthArrayIndex = 0;

    auto primaryPrecision = PrecisionArguments(blockPrecisionCPUArray[primAngularArrayIndex],
                                               blockPrecisionCPUArray[primThicknessArrayIndex],
                                               blockPrecisionCPUArray[primLengthArrayIndex],
                                               incrementPrecisionCPUArray[primAngularArrayIndex],
                                               incrementPrecisionCPUArray[primThicknessArrayIndex],
                                               incrementPrecisionCPUArray[primLengthArrayIndex]);

    auto secondaryPrecision = PrecisionArguments(blockPrecisionCPUArray[secAngularArrayIndex],
                                                 blockPrecisionCPUArray[secThicknessArrayIndex],
                                                 blockPrecisionCPUArray[secLengthArrayIndex],
                                                 incrementPrecisionCPUArray[secAngularArrayIndex],
                                                 incrementPrecisionCPUArray[secThicknessArrayIndex],
                                                 incrementPrecisionCPUArray[secLengthArrayIndex]);

    #if PRINT_ENABLED
        printf("%d : %d %d %d | %d %d %d\n", currentIncrements,
               blockPrecisionCPUArray[primLengthArrayIndex] * incrementPrecisionCPUArray[primLengthArrayIndex],
               blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex],
               blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex],
               blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex],
               blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex],
               blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);
    #endif // PRINT_ENABLED

    return CoilPairArguments(primaryPrecision, secondaryPrecision);
}


