#include "Coil.h"

#define _USE_MATH_DEFINES
#include <math.h>


double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     CoilPairArguments inductanceArguments, ComputeMethod computeMethod)
{
    if (isZAxisCase(primary, secondary))
    {
        vec3::FieldVector3 secPositionVec = vec3::CoordVector3::convertToFieldVector(secondary.getPositionVector());
        return calculateMutualInductanceZAxis(primary, secondary, secPositionVec.z, inductanceArguments, computeMethod);
    }
    else
        return calculateMutualInductanceGeneral(primary, secondary, inductanceArguments, computeMethod);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     PrecisionFactor precisionFactor, ComputeMethod computeMethod)
{
    bool zAxisCase = isZAxisCase(primary, secondary);
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, computeMethod, zAxisCase);

    return computeMutualInductance(primary, secondary, args, computeMethod);
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, CoilPairArguments inductanceArguments,
                                            ComputeMethod computeMethod) const
{
    return computeMutualInductance(*this, secondary, inductanceArguments, computeMethod) * 2*M_PI * sineFrequency;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, PrecisionFactor precisionFactor,
                                            ComputeMethod computeMethod) const
{
    return computeMutualInductance(*this, secondary, precisionFactor, computeMethod) * 2*M_PI * sineFrequency;
}