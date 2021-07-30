#include "Coil.h"
#include "PrecisionGlobalVars.h"

#include <cmath>
#include <cstdio>

CoilPairArguments CoilPairArguments::calculateCoilPairArgumentsZAxisCPU(const Coil &primary, const Coil &secondary,
                                                                        PrecisionFactor precisionFactor)
                                                                        {
    int primLengthArrayIndex = g_minPrimLengthIncrements - 1;
    int primThicknessArrayIndex = g_minPrimThicknessIncrements -1;
    int primAngularArrayIndex = g_minPrimAngularIncrements - 1;

    int secLengthArrayIndex = g_minSecLengthIncrements - 1;
    int secThicknessArrayIndex = g_minSecThicknessIncrements - 1;

    int totalIncrements = 0;
    int currentIncrements;
    int caseIndex;

    getGeometryCaseAndIncrementsCoilPair(primary, secondary, precisionFactor, caseIndex, totalIncrements);

    do
    {
        double primAngularStep = M_PI * (primary.getInnerRadius() + primary.getThickness() * 0.5) /
                (blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex]);

        double primLengthStep = sqrt(2) * primary.getLength() /
                (blockPrecisionCPUArray[primLengthArrayIndex] * incrementPrecisionCPUArray[primLengthArrayIndex]);

        double primThicknessStep = sqrt(2) * primary.getThickness() /
                (blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex]);

        double secLengthStep = secondary.getLength() /
                (blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex]);

        double secThicknessStep = secondary.getThickness() /
                (blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

        double primLinearStep = sqrt(primLengthStep * primThicknessStep);
        double secLinearStep = sqrt(secLengthStep * secThicknessStep);

        switch (caseIndex)
        {
            case (1):
                primThicknessArrayIndex = 0;
                secThicknessArrayIndex = 0; secLengthArrayIndex = 0;

                if (primAngularStep / primLengthStep >= 1.0)
                    primAngularArrayIndex++;
                else
                    primLengthArrayIndex++;
                break;
            case (2):
                primThicknessArrayIndex = 0;
                secThicknessArrayIndex = 0;

                if (primAngularStep / primLengthStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primLengthStep / secLengthStep >= 1.0)
                        primLengthArrayIndex++;
                    else
                        secLengthArrayIndex++;
                }
                break;
            case (3):
                primThicknessArrayIndex = 0;
                secLengthArrayIndex = 0;

                if (primAngularStep / primLengthStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primLengthStep / secThicknessStep >= 1.0)
                        primLengthArrayIndex++;
                    else
                        secThicknessArrayIndex++;
                }
                break;
            case (4):
                primThicknessArrayIndex = 0;

                if (primAngularStep / primLengthStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primLengthStep / secLinearStep >= 1.0)
                        primLengthArrayIndex++;
                    else
                        { secLengthArrayIndex++; secThicknessArrayIndex++; }
                }
                break;
            case (5):
                primLengthArrayIndex = 0;
                secThicknessArrayIndex = 0; secLengthArrayIndex = 0;

                if (primAngularStep / primThicknessStep >= 1.0)
                    primAngularArrayIndex++;
                else
                    primThicknessArrayIndex++;
                break;
            case (6):
                primLengthArrayIndex = 0;
                secThicknessArrayIndex = 0;

                if (primAngularStep / primThicknessStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primThicknessStep / secLengthStep >= 1.0)
                        primThicknessArrayIndex++;
                    else
                        secLengthArrayIndex++;
                }
                break;
            case (7):
                primLengthArrayIndex = 0;
                secLengthArrayIndex = 0;

                if (primAngularStep / primThicknessStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primThicknessStep / secThicknessStep >= 1.0)
                        primThicknessArrayIndex++;
                    else
                        secThicknessArrayIndex++;
                }
                break;
            case (8):
                primLengthArrayIndex = 0;

                if (primAngularStep / primThicknessStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primThicknessStep / secLinearStep >= 1.0)
                        primThicknessArrayIndex++;
                    else
                        { secLengthArrayIndex++; secThicknessArrayIndex++; }
                }
                break;
            case (9):
                primLengthArrayIndex = 0; primThicknessArrayIndex = 0;
                secLengthArrayIndex = 0; secThicknessArrayIndex = 0;

                primAngularArrayIndex++;
                break;
            case (10):
                primLengthArrayIndex = 0; primThicknessArrayIndex = 0;
                secThicknessArrayIndex = 0;

                if (primAngularStep / secLengthStep >= 1.0)
                    primAngularArrayIndex++;
                else
                    secLengthArrayIndex++;
                break;
            case (11):
                primLengthArrayIndex = 0; primThicknessArrayIndex = 0;
                secLengthArrayIndex = 0;

                if (primAngularStep / secThicknessStep >= 1.0)
                    primAngularArrayIndex++;
                else
                    secThicknessArrayIndex++;
                break;
            case (12):
                primLengthArrayIndex = 0; primThicknessArrayIndex = 0;

                if (primAngularStep / secLinearStep >= 1.0)
                    primAngularArrayIndex++;
                else
                    { secLengthArrayIndex++; secThicknessArrayIndex++; }
                break;
            case (13):
                secLengthArrayIndex = 0; secThicknessArrayIndex = 0;

                if (primAngularStep / primLinearStep >= 1.0)
                    primAngularArrayIndex++;
                else
                    { primLengthArrayIndex++; primThicknessArrayIndex++; }
                break;
            case (14):
                secThicknessArrayIndex = 0;

                if (primAngularStep / primLinearStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primLinearStep / secLengthStep >= 1.0)
                        { primLengthArrayIndex++; primThicknessArrayIndex++; }
                    else
                        secLengthArrayIndex++;
                }
                break;
            case (15):
                secLengthArrayIndex = 0;

                if (primAngularStep / primLinearStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primLinearStep / secThicknessStep >= 1.0)
                        { primLengthArrayIndex++; primThicknessArrayIndex++; }
                    else
                        secThicknessArrayIndex++;
                }
                break;
            default:
                if (primAngularStep / primLinearStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primLinearStep / secLinearStep >= 1.0)
                        { primLengthArrayIndex++; primThicknessArrayIndex++; }
                    else
                        { secLengthArrayIndex++; secThicknessArrayIndex++; }
                }
        }

        currentIncrements =
                blockPrecisionCPUArray[primLengthArrayIndex] * incrementPrecisionCPUArray[primLengthArrayIndex] *
                blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex] *
                blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex] *
                blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex] *
                blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex];
    }
    while (currentIncrements < totalIncrements);

    auto primaryPrecision = PrecisionArguments(blockPrecisionCPUArray[primAngularArrayIndex],
                                               blockPrecisionCPUArray[primThicknessArrayIndex],
                                               blockPrecisionCPUArray[primLengthArrayIndex],
                                               incrementPrecisionCPUArray[primAngularArrayIndex],
                                               incrementPrecisionCPUArray[primThicknessArrayIndex],
                                               incrementPrecisionCPUArray[primLengthArrayIndex]);

    auto secondaryPrecision = PrecisionArguments(0,
                                                 blockPrecisionCPUArray[secThicknessArrayIndex],
                                                 blockPrecisionCPUArray[secLengthArrayIndex],
                                                 0,
                                                 incrementPrecisionCPUArray[secThicknessArrayIndex],
                                                 incrementPrecisionCPUArray[secLengthArrayIndex]);

    #if PRINT_ENABLED
        printf("case %d - %d : %d %d %d | %d %d\n", caseIndex, currentIncrements,
               blockPrecisionCPUArray[primLengthArrayIndex] * incrementPrecisionCPUArray[primLengthArrayIndex],
               blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex],
               blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex],
               blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex],
               blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);
    #endif // PRINT_ENABLED

    return CoilPairArguments(primaryPrecision, secondaryPrecision);
}

CoilPairArguments CoilPairArguments::calculateCoilPairArgumentsGeneralCPU(const Coil &primary, const Coil &secondary,
                                                                          PrecisionFactor precisionFactor)
{
    int primLengthArrayIndex = g_minPrimLengthIncrements - 1;
    int primThicknessArrayIndex = g_minPrimThicknessIncrements - 1;
    int primAngularArrayIndex = g_minPrimAngularIncrements - 1;

    int secLengthArrayIndex = g_minSecLengthIncrements - 1;
    int secThicknessArrayIndex = g_minSecThicknessIncrements - 1;
    int secAngularArrayIndex = g_minSecAngularIncrements - 1;

    int totalIncrements = 0;
    int currentIncrements;
    int caseIndex;

    getGeometryCaseAndIncrementsCoilPair(primary, secondary, precisionFactor, caseIndex, totalIncrements);
    // multiplying the number of increments by 16 because of one extra dimension of integration compared to z-axis
    totalIncrements *= 16;

    do
    {
        double primAngularStep = M_PI * (primary.getInnerRadius() + primary.getThickness() * 0.5) /
                (blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex]);

        double primLengthStep = sqrt(2) * primary.getLength() /
                (blockPrecisionCPUArray[primLengthArrayIndex] * incrementPrecisionCPUArray[primLengthArrayIndex]);

        double primThicknessStep = sqrt(2) * primary.getThickness() /
                (blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex]);

        double secLengthStep = secondary.getLength() /
                (blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex]);

        double secThicknessStep = secondary.getThickness() /
                (blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

        double secAngularStep = M_PI * (secondary.getInnerRadius() + secondary.getThickness() * 0.5) /
                (blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);

        double primLinearStep = sqrt(primLengthStep * primThicknessStep);
        double secLinearStep = sqrt(secLengthStep * secThicknessStep);

        switch (caseIndex)
        {
            case (1):
                primThicknessArrayIndex = 0;
                secThicknessArrayIndex = 0; secLengthArrayIndex = 0;

                if ((primLengthStep * primLengthStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                    primLengthArrayIndex++;
                break;
            case (2):
                primThicknessArrayIndex = 0;
                secThicknessArrayIndex = 0;

                if ((primLengthStep * secLengthStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secLengthStep / primLengthStep <= 1.0)
                        primLengthArrayIndex++;
                    else
                        secLengthArrayIndex++;
                }
                break;
            case (3):
                primThicknessArrayIndex = 0;
                secLengthArrayIndex = 0;

                if ((primLengthStep * secThicknessStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secThicknessStep / primLengthStep <= 1.0)
                        primLengthArrayIndex++;
                    else
                        secThicknessArrayIndex++;
                }
                break;
            case (4):
                primThicknessArrayIndex = 0;

                if ((primLengthStep * secLinearStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secLinearStep / primLengthStep <= 1.0)
                        primLengthArrayIndex++;
                    else
                        { secLengthArrayIndex++; secThicknessArrayIndex++; }
                }
                break;
            case (5):
                primLengthArrayIndex = 0;
                secThicknessArrayIndex = 0; secLengthArrayIndex = 0;

                if ((primThicknessStep * primThicknessStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                    primThicknessArrayIndex++;
                break;
            case (6):
                primLengthArrayIndex = 0;
                secThicknessArrayIndex = 0;

                if ((primThicknessStep * secLengthStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secLengthStep / primThicknessStep <= 1.0)
                        primThicknessArrayIndex++;
                    else
                        secLengthArrayIndex++;
                }
                break;
            case (7):
                primLengthArrayIndex = 0;
                secLengthArrayIndex = 0;

                if ((primThicknessStep * secThicknessStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secThicknessStep / primThicknessStep <= 1.0)
                        primThicknessArrayIndex++;
                    else
                        secThicknessArrayIndex++;
                }
                break;
            case (8):
                primLengthArrayIndex = 0;

                if ((primThicknessStep * secLinearStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secLinearStep / primThicknessStep <= 1.0)
                        primThicknessArrayIndex++;
                    else
                        { secLengthArrayIndex++; secThicknessArrayIndex++; }
                }
                break;
            case (9):
                primLengthArrayIndex = 0; primThicknessArrayIndex = 0;
                secLengthArrayIndex = 0; secThicknessArrayIndex = 0;

                if (secAngularStep / primAngularStep <= 1.0)
                    primAngularArrayIndex++;
                else
                    secAngularArrayIndex++;
                break;
            case (10):
                primLengthArrayIndex = 0; primThicknessArrayIndex = 0;
                secThicknessArrayIndex = 0;

                if ((secLengthStep * secLengthStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                    secLengthArrayIndex++;
                break;
            case (11):
                primLengthArrayIndex = 0; primThicknessArrayIndex = 0;
                secLengthArrayIndex = 0;

                if ((secThicknessStep * secThicknessStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                    secThicknessArrayIndex++;
                break;
            case (12):
                primLengthArrayIndex = 0; primThicknessArrayIndex = 0;

                if ((secLengthStep * secThicknessStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                    { secThicknessArrayIndex++; secThicknessArrayIndex++; }
                break;
            case (13):
                secLengthArrayIndex = 0; secThicknessArrayIndex = 0;

                if ((primLengthStep * primThicknessStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                    { primLengthArrayIndex++; primThicknessArrayIndex++; }
                break;
            case (14):
                secThicknessArrayIndex = 0;

                if ((primLinearStep * secLengthStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secLengthStep / primLinearStep <= 1.0)
                        { primLengthArrayIndex++; primThicknessArrayIndex++; }
                    else
                        secLengthArrayIndex++;
                }
                break;
            case (15):
                secLengthArrayIndex = 0;

                if ((primLinearStep * secThicknessStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secThicknessStep / primLinearStep <= 1.0)
                        { primLengthArrayIndex++; primThicknessArrayIndex++; }
                    else
                        secThicknessArrayIndex++;
                }
                break;
            default:
                if ((primLinearStep * secLinearStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secLinearStep / primLinearStep <= 1.0)
                        { primLengthArrayIndex++; primThicknessArrayIndex++; }
                    else
                        { secLengthArrayIndex++; secThicknessArrayIndex++; }
                }
        }

        currentIncrements =
                blockPrecisionCPUArray[primLengthArrayIndex] * incrementPrecisionCPUArray[primLengthArrayIndex] *
                blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex] *
                blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex] *
                blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex] *
                blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex] *
                blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex];
    }
    while (currentIncrements < totalIncrements);

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
        printf("case %d - %d : %d %d %d | %d %d %d\n", caseIndex, currentIncrements,
               blockPrecisionCPUArray[primLengthArrayIndex] * incrementPrecisionCPUArray[primLengthArrayIndex],
               blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex],
               blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex],
               blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex],
               blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex],
               blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);
    #endif // PRINT_ENABLED

    return CoilPairArguments(primaryPrecision, secondaryPrecision);
}

