#include "Coil.h"

#include <cmath>
#include <cstdio>


CoilPairArguments CoilPairArguments::calculateCoilPairArgumentsCPU(const Coil &primary, const Coil &secondary,
                                                                   PrecisionFactor precisionFactor, bool zAxisCase)
{
    int primLengthArrayIndex, primThicknessArrayIndex, primAngularArrayIndex;
    int secLengthArrayIndex, secThicknessArrayIndex, secAngularArrayIndex;

    int totalIncrements = (int) std::round(std::pow(2, precisionFactor.relativePrecision - 1));
    int currentIncrements;

    double primRadius = primary.getInnerRadius();
    double primThickness = primary.getThickness();

    double secRadius = secondary.getInnerRadius();
    double secThickness = secondary.getThickness();
    double secLength = secondary.getLength();

    // multiplying for primary angular increments
    totalIncrements *= g_baseLayerIncrements;
    // multiplying for secondary angular increments, if present
    if (!zAxisCase)
    {
        totalIncrements *= g_baseLayerIncrements;
    }

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

    double secDimensionRatio = std::sqrt(secLength / secThickness);
    switch (secondary.getCoilType())
    {
        case CoilType::RECTANGULAR:
        {
            totalIncrements *= g_baseLayerIncrements * g_baseLayerIncrements;
            // ensuring the minimum distributions are 2 * nB or nA * 2, otherwise it is just a thin coil
            if (secDimensionRatio >= (double) g_minSecAreaIncrements / 2.0)
            {
                secLengthArrayIndex = g_minSecAreaIncrements / 2 - 1;
                secThicknessArrayIndex = 1;
            }
            else if (1.0 / secDimensionRatio >= (double) g_minSecAreaIncrements / 2.0)
            {
                secLengthArrayIndex = 1;
                secThicknessArrayIndex = g_minSecAreaIncrements / 2 - 1;
            }
            else if (secDimensionRatio >= 1.0)
            {
                secThicknessArrayIndex = (int) std::ceil(std::sqrt(g_minSecAreaIncrements / secDimensionRatio)) - 1;
                secLengthArrayIndex = g_minSecAreaIncrements / (secThicknessArrayIndex + 1) - 1;
            }
            else
            {
                secLengthArrayIndex = (int) std::ceil(std::sqrt(g_minSecAreaIncrements * secDimensionRatio)) - 1;
                secThicknessArrayIndex = g_minSecAreaIncrements / (secLengthArrayIndex + 1) - 1;
            }
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

    double primAngularStep, primThicknessStep;
    double secAngularStep, secLengthStep, secThicknessStep, secLinearStep;
    double secLengthRootStep, secThicknessRootStep;

    do
    {
        primAngularStep = M_PI * g_primAngularWeightModifier * (primRadius + 0.5 * primThickness) /
                (blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex]);

        primThicknessStep = g_primLinearWeightModifier * primThickness /
                (blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex]);

        secAngularStep = 2 * M_PI * g_secAngularWeightModifier * (secRadius + 0.5 * secThickness) /
                (blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);

        secLengthStep = g_secLinearWeightModifier * secLength /
                (blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex]);

        secThicknessStep = g_secAngularWeightModifier * secThickness /
                (blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

        secLinearStep = 0.5 * (secLengthStep + secThicknessStep);

        secLengthRootStep = std::sqrt(g_secLinearWeightModifier * secLength) /
                (blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex]);

        secThicknessRootStep = std::sqrt(g_secLinearWeightModifier * secThickness) /
                (blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

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
                if (secLengthRootStep >= secThicknessRootStep)
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


