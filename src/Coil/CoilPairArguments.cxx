#include "Coil.h"
#include "PrecisionGlobalVars.h"
#include "CoilData.h"

#include <cmath>
#include <cstdio>


CoilPairArguments::CoilPairArguments() : CoilPairArguments(PrecisionArguments(), PrecisionArguments()) {}

CoilPairArguments::CoilPairArguments(const PrecisionArguments &primaryPrecision,
                                     const PrecisionArguments &secondaryPrecision)
{
    CoilPairArguments::primaryPrecision = primaryPrecision;
    CoilPairArguments::secondaryPrecision = secondaryPrecision;
}

void CoilPairArguments::getGeometryCaseAndIncrementsSingleCoil(const Coil &coil, PrecisionFactor precisionFactor,
                                                               int &caseIndex, int &totalIncrements)
{
    double radius = coil.getInnerRadius();
    double thickness = coil.getThickness();
    double length = coil.getLength();

    if (thickness / radius < g_thinCoilApproximationRatio && length / radius < g_thinCoilApproximationRatio)
    { caseIndex = 1; totalIncrements = pow(2, 3 + precisionFactor.relativePrecision); }

    else if (thickness / length < g_thinCoilApproximationRatio)
    { caseIndex = 2; totalIncrements = pow(2, 6 + precisionFactor.relativePrecision); }

    else if (length / thickness < g_thinCoilApproximationRatio)
    { caseIndex = 3; totalIncrements = pow(2, 6 + precisionFactor.relativePrecision); }

    else
    { caseIndex = 4; totalIncrements = pow(2, 9 + precisionFactor.relativePrecision); }
}

void CoilPairArguments::getGeometryCaseAndIncrementsCoilPair(const Coil &primary, const Coil &secondary,
                                                             PrecisionFactor precisionFactor,
                                                             int &caseIndex, int &totalIncrements)
{
    double primRadius = primary.getInnerRadius();
    double primThickness = primary.getThickness();
    double primLength = primary.getLength();

    double secRadius = secondary.getInnerRadius();
    double secThickness = secondary.getThickness();
    double secLength = secondary.getLength();

    if (primThickness / primLength < g_thinCoilApproximationRatio)
    {
        getGeometryCaseAndIncrementsSingleCoil(secondary, precisionFactor, caseIndex, totalIncrements);

        totalIncrements *= 8;
    }
    else if (primLength / primThickness < g_thinCoilApproximationRatio)
    {
        getGeometryCaseAndIncrementsSingleCoil(secondary, precisionFactor, caseIndex, totalIncrements);

        caseIndex += 4;
        totalIncrements *= 8;
    }
    else if (primThickness / primRadius < g_thinCoilApproximationRatio && primLength / primRadius < g_thinCoilApproximationRatio)
    {
        getGeometryCaseAndIncrementsSingleCoil(secondary, precisionFactor, caseIndex, totalIncrements);

        caseIndex += 8;
    }
    else
    {
        getGeometryCaseAndIncrementsSingleCoil(secondary, precisionFactor, caseIndex, totalIncrements);

        caseIndex += 12;
        totalIncrements *= 8 * 8;
    }
}

CoilPairArguments CoilPairArguments::calculateMInductanceArgumentsZCPU(const Coil &primary, const Coil &secondary,
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

    if (PRINT_ENABLED)
        printf("case %d - %d : %d %d %d | %d %d\n", caseIndex, currentIncrements,
               blockPrecisionCPUArray[primLengthArrayIndex] * incrementPrecisionCPUArray[primLengthArrayIndex],
               blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex],
               blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex],
               blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex],
               blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

    return CoilPairArguments(primaryPrecision, secondaryPrecision);
}

CoilPairArguments CoilPairArguments::calculateMInductanceArgumentsGeneralCPU(const Coil &primary, const Coil &secondary,
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
    if (PRINT_ENABLED)
        printf("case %d - %d : %d %d %d | %d %d %d\n", caseIndex, currentIncrements,
               blockPrecisionCPUArray[primLengthArrayIndex] * incrementPrecisionCPUArray[primLengthArrayIndex],
               blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex],
               blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex],
               blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex],
               blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex],
               blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);

    return CoilPairArguments(primaryPrecision, secondaryPrecision);
}

CoilPairArguments CoilPairArguments::calculateMInductanceArgumentsZGPU(const Coil &primary, const Coil &secondary,
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

    getGeometryCaseAndIncrementsSingleCoil(secondary, precisionFactor, caseIndex, totalIncrements);
    //multiplying total increments by 2^7 to satisfy this specific workload
    totalIncrements *= 128;

    do
    {
        double primAngularStep = M_PI * (primary.getInnerRadius() + primary.getThickness() * 0.5) /
                                 (primAngularBlocks * primAngularIncrements);

        double primLinearStep = sqrt(primary.getThickness() * primary.getLength()) /
                                (primLinearBlocks * primLinearIncrements);

        double secLengthStep = secondary.getLength() /
                               (blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex]);

        double secThicknessStep = secondary.getThickness() /
                                  (blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

        double secLinearStep = sqrt(secLengthStep * secThicknessStep);

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
                if (primAngularStep / sqrt(secLengthStep * primLinearStep) <= 1.0)
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
                if (primAngularStep / sqrt(secThicknessStep * primLinearStep) <= 1.0)
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
                if (primAngularStep / sqrt(secLinearStep * primLinearStep) <= 1.0)
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
    if (PRINT_ENABLED)
        printf("case %d - %d : %d %d %d | %d %d\n", caseIndex, currentIncrements,
               primLinearBlocks * primLinearIncrements,
               primLinearBlocks * primLinearIncrements,
               primAngularBlocks * primAngularIncrements,
               blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex],
               blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

    return CoilPairArguments(primaryPrecision, secondaryPrecision);
}

CoilPairArguments CoilPairArguments::calculateMInductanceArgumentsGeneralGPU(const Coil &primary, const Coil &secondary,
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

    getGeometryCaseAndIncrementsSingleCoil(secondary, precisionFactor, caseIndex, currentIncrements);
    //multiplying total increments by 2^10 to satisfy this specific workload
    totalIncrements *= 1024;

    do
    {
        double primAngularStep = M_PI * (primary.getInnerRadius() + primary.getThickness() * 0.5) /
                                 (primAngularBlocks * primAngularIncrements);

        double primLinearStep = sqrt(primary.getThickness() * primary.getLength()) /
                                (primLinearBlocks * primLinearIncrements);

        double secLengthStep = secondary.getLength() /
                               (blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex]);

        double secThicknessStep = secondary.getThickness() /
                                  (blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

        double secAngularStep = M_PI * (secondary.getInnerRadius() + secondary.getThickness() * 0.5) /
                                (blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);

        double secLinearStep = sqrt(secLengthStep * secThicknessStep);

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
    if (PRINT_ENABLED)
        printf("case %d - %d : %d %d %d | %d %d %d\n", caseIndex, currentIncrements,
               primLinearBlocks * primLinearIncrements,
               primLinearBlocks * primLinearIncrements,
               primAngularBlocks * primAngularIncrements,
               blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex],
               blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex],
               blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);

    return CoilPairArguments(primaryPrecision, secondaryPrecision);
}

CoilPairArguments CoilPairArguments::getSelfInductanceArguments(const Coil &coil, PrecisionFactor precisionFactor)
{
    int lengthArrayIndex = g_minPrimLengthIncrements - 1;
    int thicknessArrayIndex = g_minPrimThicknessIncrements -1;
    int angularArrayIndex = g_minPrimAngularIncrements - 1;

    int totalIncrements = 0;
    int currentIncrements;
    int caseIndex;

    if (coil.getThickness() / coil.getInnerRadius() < g_thinCoilApproximationRatio &&
        coil.getLength() / coil.getInnerRadius() < g_thinCoilApproximationRatio)
    {
        caseIndex = 1; totalIncrements = pow(2, 5 + precisionFactor.relativePrecision);
    }
    else if (coil.getThickness() / coil.getLength() < g_thinCoilApproximationRatio)
    {
        caseIndex = 2; totalIncrements = pow(2, 11 + precisionFactor.relativePrecision);
    }
    else if (coil.getLength() / coil.getThickness() < g_thinCoilApproximationRatio)
    {
        caseIndex = 3; totalIncrements = pow(2, 11 + precisionFactor.relativePrecision);
    }
    else
    {
        caseIndex = 4; totalIncrements = pow(2, 17 + precisionFactor.relativePrecision);
    }

    do
    {
        double angularStep = M_PI * (coil.getInnerRadius() + coil.getThickness() * 0.5) /
                             (blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex]);

        double lengthStep = sqrt(2) * coil.getLength() /
                            (blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex]);

        double thicknessStep = sqrt(2) * coil.getThickness() /
                               (blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex]);

        double linearStep = sqrt(lengthStep * thicknessStep);

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
                blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex] *
                blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex] *
                blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex] *
                blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex] *
                blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex];
    }
    while (currentIncrements < totalIncrements);

    auto primaryPrecision = PrecisionArguments(blockPrecisionCPUArray[angularArrayIndex],
                                               blockPrecisionCPUArray[thicknessArrayIndex],
                                               blockPrecisionCPUArray[lengthArrayIndex],
                                               incrementPrecisionCPUArray[angularArrayIndex],
                                               incrementPrecisionCPUArray[thicknessArrayIndex],
                                               incrementPrecisionCPUArray[lengthArrayIndex]);

    auto secondaryPrecision = PrecisionArguments(0,
                                                 blockPrecisionCPUArray[thicknessArrayIndex],
                                                 blockPrecisionCPUArray[lengthArrayIndex],
                                                 0,
                                                 incrementPrecisionCPUArray[thicknessArrayIndex],
                                                 incrementPrecisionCPUArray[lengthArrayIndex]);
    if (PRINT_ENABLED)
        printf("case %d - %d : %d %d %d | %d %d\n", caseIndex, currentIncrements,
               blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex],
               blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex],
               blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex],
               blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex],
               blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex]);

    return CoilPairArguments(primaryPrecision, secondaryPrecision);
}

CoilPairArguments CoilPairArguments::getAppropriateMInductanceArguments(const Coil &primary, const Coil &secondary,
                                                                        PrecisionFactor precisionFactor,
                                                                        ComputeMethod method, bool isGeneral)
{
    if (!isGeneral)
    {
        if (method == GPU)
            return CoilPairArguments::calculateMInductanceArgumentsZGPU(primary, secondary, precisionFactor);
        else
            return CoilPairArguments::calculateMInductanceArgumentsZCPU(primary, secondary, precisionFactor);
    }
    else
    {
        if (method == GPU)
            return CoilPairArguments::calculateMInductanceArgumentsGeneralGPU(primary, secondary, precisionFactor);
        else
            return CoilPairArguments::calculateMInductanceArgumentsGeneralCPU(primary, secondary, precisionFactor);
    }
}