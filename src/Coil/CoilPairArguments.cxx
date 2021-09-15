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

CoilPairArguments CoilPairArguments::getAppropriateCoilPairArguments(const Coil &primary, const Coil &secondary,
                                                                     PrecisionFactor precisionFactor,
                                                                     ComputeMethod method, bool isGeneral)
{
    if (!isGeneral)
    {
        if (method == GPU)
            return CoilPairArguments::calculateCoilPairArgumentsZAxisGPU(primary, secondary, precisionFactor);
        else
            return CoilPairArguments::calculateCoilPairArgumentsCPU(primary, secondary, precisionFactor, !isGeneral);
    }
    else
    {
        if (method == GPU)
            return CoilPairArguments::calculateCoilPairArgumentsGeneralGPU(primary, secondary, precisionFactor);
        else
            return CoilPairArguments::calculateCoilPairArgumentsCPU(primary, secondary, precisionFactor, !isGeneral);
    }
}
