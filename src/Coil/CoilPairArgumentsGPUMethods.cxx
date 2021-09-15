#include "Coil.h"
#include "PrecisionGlobalVars.h"
#include "CoilData.h"

#include <cmath>
#include <cstdio>

CoilPairArguments CoilPairArguments::calculateCoilPairArgumentsZAxisGPU(const Coil &primary, const Coil &secondary,
                                                                        PrecisionFactor precisionFactor)
{
    const int primLinearIncrements = arrSize;
    const int primAngularIncrements = arrSize;

    int primLinearBlocks = 1;
    int primAngularBlocks = 1;

    int secLengthArrayIndex = g_minSecLengthIncrements - 1;
    int secThicknessArrayIndex = g_minSecThicknessIncrements - 1;

    int totalIncrements = 0;
    int currentIncrements;
    int caseIndex;

    getGeometryCaseAndIncrementsSingleCoil(secondary, caseIndex, totalIncrements);

    totalIncrements *= g_baseLayerIncrements * g_baseLayerIncrements * g_baseLayerIncrements *
            std::pow(2, precisionFactor.relativePrecision);


    do
    {
        double primAngularStep = M_PI * (primary.getInnerRadius() + primary.getThickness() * 0.5) /
                (primAngularBlocks * primAngularIncrements);

        double primLinearStep = std::sqrt(primary.getThickness() * primary.getLength()) /
                (primLinearBlocks * primLinearIncrements);

        double secLengthStep = secondary.getLength() /
                (blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex]);

        double secThicknessStep = secondary.getThickness() /
                (blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

        double secLinearStep = std::sqrt(secLengthStep * secThicknessStep);

        switch (caseIndex)
        {
            case (1):
                secLengthArrayIndex = 0; secThicknessArrayIndex = 0;

                if (primAngularStep / primLinearStep <= 1.0)
                    primLinearBlocks++;
                else
                    primAngularBlocks++;
                break;
            case (2):
                secThicknessArrayIndex = 0;

                if (primAngularStep / std::sqrt(secLengthStep * primLinearStep) <= 1.0)
                {
                    if (secLengthStep / primLinearStep >= 1.0)
                        secLengthArrayIndex++;
                    else
                        primLinearBlocks++;
                }
                else
                {
                    primAngularBlocks++;
                }
                break;
            case (3):
                secLengthArrayIndex = 0;

                if (primAngularStep / std::sqrt(secThicknessStep * primLinearStep) <= 1.0)
                {
                    if (secThicknessStep / primLinearStep >= 1.0)
                        secThicknessArrayIndex++;
                    else
                        primLinearBlocks++;
                }
                else
                    primAngularBlocks++;
                break;
            default:
                if (primAngularStep / std::sqrt(secLinearStep * primLinearStep) <= 1.0)
                {
                    if (secLinearStep / primLinearStep >= 1.0)
                        { secLengthArrayIndex++; secThicknessArrayIndex++; }
                    else
                        primLinearBlocks++;
                }
                else
                    primAngularBlocks++;
        }
        currentIncrements = primAngularBlocks * primAngularIncrements *
                primLinearBlocks * primLinearIncrements * primLinearBlocks * primLinearIncrements *
                blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex] *
                blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex];
    }
    while (currentIncrements < totalIncrements);

    auto primaryPrecision = PrecisionArguments(primAngularBlocks, primLinearBlocks, primLinearBlocks,
                                               primAngularIncrements, primLinearIncrements, primLinearIncrements);

    auto secondaryPrecision = PrecisionArguments(0,
                                                 blockPrecisionCPUArray[secThicknessArrayIndex],
                                                 blockPrecisionCPUArray[secLengthArrayIndex],
                                                 0,
                                                 incrementPrecisionCPUArray[secThicknessArrayIndex],
                                                 incrementPrecisionCPUArray[secLengthArrayIndex]);
    #if PRINT_ENABLED
        printf("case %d - %d : %d %d %d | %d %d\n", caseIndex, currentIncrements,
               primLinearBlocks * primLinearIncrements,
               primLinearBlocks * primLinearIncrements,
               primAngularBlocks * primAngularIncrements,
               blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex],
               blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);
    #endif // PRINT_ENABLED

    return CoilPairArguments(primaryPrecision, secondaryPrecision);
}

CoilPairArguments CoilPairArguments::calculateCoilPairArgumentsGeneralGPU(const Coil &primary, const Coil &secondary,
                                                                          PrecisionFactor precisionFactor)
{
    const int primLinearIncrements = arrSize;
    const int primAngularIncrements = arrSize;

    int primLinearBlocks = 1;
    int primAngularBlocks = 1;

    int secLengthArrayIndex = g_minSecLengthIncrements - 1;
    int secThicknessArrayIndex = g_minSecThicknessIncrements - 1;
    int secAngularArrayIndex = g_minSecAngularIncrements - 1;

    int totalIncrements = 0;
    int currentIncrements;
    int caseIndex;

    getGeometryCaseAndIncrementsSingleCoil(secondary, caseIndex, currentIncrements);

    totalIncrements *= g_baseLayerIncrements * g_baseLayerIncrements * g_baseLayerIncrements * g_baseLayerIncrements
            * std::pow(2, precisionFactor.relativePrecision);

    do
    {
        double primAngularStep = M_PI * (primary.getInnerRadius() + primary.getThickness() * 0.5) /
                (primAngularBlocks * primAngularIncrements);

        double primLinearStep = std::sqrt(primary.getThickness() * primary.getLength()) /
                (primLinearBlocks * primLinearIncrements);

        double secLengthStep = secondary.getLength() /
                (blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex]);

        double secThicknessStep = secondary.getThickness() /
                (blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

        double secAngularStep = M_PI * (secondary.getInnerRadius() + secondary.getThickness() * 0.5) /
                (blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);

        double secLinearStep = std::sqrt(secLengthStep * secThicknessStep);

        switch (caseIndex)
        {
            case (1):
                secLengthArrayIndex = 0; secThicknessArrayIndex = 0;

                if ((primLinearStep * primLinearStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularBlocks++;
                    else
                        secAngularArrayIndex++;
                }
                else
                    primLinearBlocks++;
                break;
            case (2):
                secThicknessArrayIndex = 0;

                if ((primLinearStep * secLengthStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularBlocks++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secLengthStep / primLinearStep <= 1.0)
                        primLinearBlocks++;
                    else
                        secLengthArrayIndex++;
                }
                break;
            case (3):
                secLengthArrayIndex = 0;

                if ((primLinearStep * secThicknessStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularBlocks++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secThicknessStep / primLinearStep <= 1.0)
                        primLinearBlocks++;
                    else
                        secThicknessArrayIndex++;
                }
                break;
            default:
                if ((primLinearStep * secLinearStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularBlocks++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secLinearStep / primLinearStep <= 1.0)
                        primLinearBlocks++;
                    else
                        { secLengthArrayIndex++; secThicknessArrayIndex++; }
                }
        }
        currentIncrements = primAngularBlocks * primAngularIncrements *
                primLinearBlocks * primLinearIncrements * primLinearBlocks * primLinearIncrements *
                blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex] *
                blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex] *
                blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex];
    }
    while (currentIncrements < totalIncrements);

    auto primaryPrecision = PrecisionArguments(primAngularBlocks, primLinearBlocks, primLinearBlocks,
                                               primAngularIncrements, primLinearIncrements, primLinearIncrements);

    auto secondaryPrecision = PrecisionArguments(blockPrecisionCPUArray[secAngularArrayIndex],
                                                 blockPrecisionCPUArray[secThicknessArrayIndex],
                                                 blockPrecisionCPUArray[secLengthArrayIndex],
                                                 incrementPrecisionCPUArray[secAngularArrayIndex],
                                                 incrementPrecisionCPUArray[secThicknessArrayIndex],
                                                 incrementPrecisionCPUArray[secLengthArrayIndex]);
    #if PRINT_ENABLED
        printf("case %d - %d : %d %d %d | %d %d %d\n", caseIndex, currentIncrements,
               primLinearBlocks * primLinearIncrements,
               primLinearBlocks * primLinearIncrements,
               primAngularBlocks * primAngularIncrements,
               blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex],
               blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex],
               blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);
    #endif // PRINT_ENABLED

    return CoilPairArguments(primaryPrecision, secondaryPrecision);
}
