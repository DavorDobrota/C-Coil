#include "Coil.h"
#include "PrecisionGlobalVars.h"

#include <cmath>
#include <cstdio>


CoilPairArguments::CoilPairArguments() : CoilPairArguments(PrecisionArguments(), PrecisionArguments()) {}

CoilPairArguments::CoilPairArguments(const PrecisionArguments &primaryPrecision,
                                     const PrecisionArguments &secondaryPrecision)
{
    this->primaryPrecision = primaryPrecision;
    this->secondaryPrecision = secondaryPrecision;
}

void CoilPairArguments::getGeometryCaseAndIncrementsSingleCoil(const Coil &coil, PrecisionFactor precisionFactor,
                                                               int &caseIndex, int &totalIncrements)
{
    double radius = coil.getInnerRadius();
    double thickness = coil.getThickness();
    double length = coil.getLength();

    double precisionMultiplier = std::pow(2, precisionFactor.relativePrecision - 1);

    if (thickness / radius < g_thinCoilApproximationRatio && length / radius < g_thinCoilApproximationRatio)
    {
        caseIndex = 1;
        totalIncrements = (int) precisionMultiplier;
    }
    else if (thickness / length < g_thinCoilApproximationRatio)
    {
        caseIndex = 2;
        totalIncrements = (int) (g_baseLayerIncrements * precisionMultiplier);
    }
    else if (length / thickness < g_thinCoilApproximationRatio)
    {
        caseIndex = 3;
        totalIncrements = (int) (g_baseLayerIncrements * precisionMultiplier);
    }
    else
    {
        caseIndex = 4;
        totalIncrements = (int) (g_baseLayerIncrements * g_baseLayerIncrements * precisionMultiplier);
    }
}

void CoilPairArguments::getGeometryCaseAndIncrementsCoilPair(const Coil &primary, const Coil &secondary,
                                                             PrecisionFactor precisionFactor,
                                                             int &caseIndex, int &totalIncrements)
{
    double primRadius = primary.getInnerRadius();
    double primThickness = primary.getThickness();
    double primLength = primary.getLength();

    if (primThickness / primLength < g_thinCoilApproximationRatio)
    {
        getGeometryCaseAndIncrementsSingleCoil(secondary, precisionFactor, caseIndex, totalIncrements);

        totalIncrements *= g_baseLayerIncrements;
    }
    else if (primLength / primThickness < g_thinCoilApproximationRatio)
    {
        getGeometryCaseAndIncrementsSingleCoil(secondary, precisionFactor, caseIndex, totalIncrements);

        caseIndex += 4;
        totalIncrements *= g_baseLayerIncrements * g_baseLayerIncrements;
    }
    else if (primThickness / primRadius < g_thinCoilApproximationRatio && primLength / primRadius < g_thinCoilApproximationRatio)
    {
        getGeometryCaseAndIncrementsSingleCoil(secondary, precisionFactor, caseIndex, totalIncrements);

        caseIndex += 8;
        totalIncrements *= g_baseLayerIncrements;
    }
    else
    {
        getGeometryCaseAndIncrementsSingleCoil(secondary, precisionFactor, caseIndex, totalIncrements);

        caseIndex += 12;
        totalIncrements *= g_baseLayerIncrements * g_baseLayerIncrements;
    }
}

CoilPairArguments CoilPairArguments::getAppropriateCoilPairArguments(const Coil &primary, const Coil &secondary,
                                                                     PrecisionFactor precisionFactor,
                                                                     ComputeMethod method, bool isGeneral)
{
    if (!isGeneral)
    {
        if (method == GPU)
            return CoilPairArguments::calculateCoilPairArgumentsZAxisGPU(primary, secondary, precisionFactor);
        else
            return CoilPairArguments::calculateCoilPairArgumentsZAxisCPU(primary, secondary, precisionFactor);
    }
    else
    {
        if (method == GPU)
            return CoilPairArguments::calculateCoilPairArgumentsGeneralGPU(primary, secondary, precisionFactor);
        else
            return CoilPairArguments::calculateCoilPairArgumentsGeneralCPU(primary, secondary, precisionFactor);
    }
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
        caseIndex = 1; totalIncrements = std::pow(2, 5 + precisionFactor.relativePrecision);
    }
    else if (coil.getThickness() / coil.getLength() < g_thinCoilApproximationRatio)
    {
        caseIndex = 2; totalIncrements = std::pow(2, 11 + precisionFactor.relativePrecision);
    }
    else if (coil.getLength() / coil.getThickness() < g_thinCoilApproximationRatio)
    {
        caseIndex = 3; totalIncrements = std::pow(2, 11 + precisionFactor.relativePrecision);
    }
    else
    {
        caseIndex = 4; totalIncrements = std::pow(2, 17 + precisionFactor.relativePrecision);
    }

    do
    {
        double angularStep = M_PI * (coil.getInnerRadius() + coil.getThickness() * 0.5) /
                             (blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex]);

        double lengthStep = std::sqrt(2) * coil.getLength() /
                            (blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex]);

        double thicknessStep = std::sqrt(2) * coil.getThickness() /
                               (blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex]);

        double linearStep = std::sqrt(lengthStep * thicknessStep);

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
    #if PRINT_ENABLED
        printf("case %d - %d : %d %d %d | %d %d\n", caseIndex, currentIncrements,
               blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex],
               blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex],
               blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex],
               blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex],
               blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex]);
    #endif // PRINT_ENABLED

    return CoilPairArguments(primaryPrecision, secondaryPrecision);
}