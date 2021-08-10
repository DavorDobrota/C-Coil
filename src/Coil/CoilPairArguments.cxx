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

void CoilPairArguments::getGeometryCaseAndIncrementsSingleCoil(const Coil &coil, int &caseIndex, int &totalIncrements)
{
    switch (coil.getCoilType())
    {
        case CoilType::FILAMENT:
            caseIndex = 1; totalIncrements = 1;
            break;
        case CoilType::THIN:
            caseIndex = 2; totalIncrements = g_baseLayerIncrements;
            break;
        case CoilType::FLAT:
            caseIndex = 3; totalIncrements = g_baseLayerIncrements;
            break;
        default:
            caseIndex = 4;
            totalIncrements = g_baseLayerIncrements * g_baseLayerIncrements;
    }
}

void CoilPairArguments::getGeometryCaseAndIncrementsCoilPair(const Coil &primary, const Coil &secondary,
                                                             PrecisionFactor precisionFactor,
                                                             int &caseIndex, int &totalIncrements)
{
    double precisionMultiplier = std::pow(2, precisionFactor.relativePrecision - 1.0);

    getGeometryCaseAndIncrementsSingleCoil(secondary, caseIndex, totalIncrements);

    switch (primary.getCoilType())
    {
        case CoilType::THIN:
            totalIncrements *= g_baseLayerIncrements * precisionMultiplier;
            break;
        case CoilType::FLAT:
            caseIndex += 4;
            totalIncrements *= g_baseLayerIncrements * g_baseLayerIncrements * precisionMultiplier;
            break;
        case CoilType::FILAMENT:
            caseIndex += 8;
            totalIncrements *= g_baseLayerIncrements * precisionMultiplier;
            break;
        default:
            caseIndex += 12;
            totalIncrements *= g_baseLayerIncrements * g_baseLayerIncrements * precisionMultiplier;
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
