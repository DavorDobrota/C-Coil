#include "Coil.h"
#include "CoilData.h"

#include <cmath>
#include <cstdio>


CoilPairArguments::CoilPairArguments() : CoilPairArguments(PrecisionArguments(), PrecisionArguments()) {}

CoilPairArguments::CoilPairArguments(const PrecisionArguments &primaryPrecision,
                                     const PrecisionArguments &secondaryPrecision)
{
    this->primaryPrecision = primaryPrecision;
    this->secondaryPrecision = secondaryPrecision;
}


CoilPairArguments CoilPairArguments::getAppropriateCoilPairArguments(const Coil &primary, const Coil &secondary,
                                                                     PrecisionFactor precisionFactor,
                                                                     ComputeMethod method, bool zAxisCase)
{
    if (method == GPU)
        return CoilPairArguments::calculateCoilPairArgumentsGPU(primary, secondary, precisionFactor, zAxisCase);
    else
        return CoilPairArguments::calculateCoilPairArgumentsCPU(primary, secondary, precisionFactor, zAxisCase);
}

CoilPairArguments CoilPairArguments::calculateCoilPairArgumentsCPU(const Coil &primary, const Coil &secondary,
                                                                   PrecisionFactor precisionFactor, bool zAxisCase)
{
    int primLengthArrayIndex, primThicknessArrayIndex, primAngularArrayIndex;
    int secLengthArrayIndex, secThicknessArrayIndex, secAngularArrayIndex;

    int totalIncrements = g_baseLayerIncrements;
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

    totalIncrements *= std::pow(2, precisionFactor.relativePrecision - 1.0);

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

CoilPairArguments CoilPairArguments::calculateCoilPairArgumentsGPU(const Coil &primary, const Coil &secondary,
                                                                   PrecisionFactor precisionFactor, bool zAxisCase)
{
    const int primLengthIncrements = arrSize;
    const int primThicknessIncrements = arrSize;
    const int primAngularIncrements = arrSize;

    int primLengthBlocks = 1;
    int primThicknessBlocks = 1;
    int primAngularBlocks = 1;

    int secLengthArrayIndex, secThicknessArrayIndex, secAngularArrayIndex;

    int totalIncrements = arrSize * arrSize * arrSize;
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

    totalIncrements *= std::pow(2, precisionFactor.relativePrecision);

    double primAngularStep, primThicknessStep, primLengthStep, primLinearStep;
    double secAngularStep, secLengthStep, secThicknessStep, secLinearStep;

    do
    {
        primLengthStep = g_primAngularWeightModifier * primLengthRoot / (primLengthIncrements * primLengthBlocks);
        primThicknessStep = primThicknessRoot / (primThicknessIncrements * primThicknessBlocks);
        primAngularStep = primAngularRoot / (primAngularIncrements * primAngularBlocks);

        secLengthStep = secLengthRoot /
                (blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex]);

        secThicknessStep = secThicknessRoot /
                (blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

        secAngularStep = g_secAngularWeightModifier * secAngularRoot /
                (blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);

        primLinearStep = 0.5 * (primLengthStep + primThicknessStep);
        secLinearStep = 0.5 * (secLengthStep + secThicknessStep);

        if (primAngularStep + primLinearStep >= secAngularStep + secLinearStep)
        {
            if (primAngularStep >= primLinearStep)
                primAngularBlocks++;
            else
            {
                if (primLengthStep >= primThicknessStep)
                    primLengthBlocks++;
                else
                    primThicknessBlocks++;
            }
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

        currentIncrements = primAngularBlocks * primAngularIncrements *
                primLengthBlocks * primLengthIncrements * primThicknessBlocks * primThicknessIncrements *
                blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex] *
                blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex];

        if (!zAxisCase)
            currentIncrements *= blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex];
    }
    while (currentIncrements < totalIncrements);

    auto primaryPrecision = PrecisionArguments(primAngularBlocks, primThicknessBlocks, primLengthBlocks,
                                               primAngularIncrements, primThicknessIncrements, primLengthIncrements);

    auto secondaryPrecision = PrecisionArguments(blockPrecisionCPUArray[secAngularArrayIndex],
                                                 blockPrecisionCPUArray[secThicknessArrayIndex],
                                                 blockPrecisionCPUArray[secLengthArrayIndex],
                                                 incrementPrecisionCPUArray[secAngularArrayIndex],
                                                 incrementPrecisionCPUArray[secThicknessArrayIndex],
                                                 incrementPrecisionCPUArray[secLengthArrayIndex]);
    #if PRINT_ENABLED
    printf("%d : %d %d %d | %d %d %d\n", currentIncrements,
           primLengthBlocks * primLengthIncrements,
           primThicknessBlocks * primThicknessIncrements,
           primAngularBlocks * primAngularIncrements,
           blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex],
           blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex],
           blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);
    #endif // PRINT_ENABLED

    return CoilPairArguments(primaryPrecision, secondaryPrecision);
}
