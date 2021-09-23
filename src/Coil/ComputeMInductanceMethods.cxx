#include "Coil.h"

#include <cmath>


double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     CoilPairArguments inductanceArguments, ComputeMethod method)
{
    if (isZAxisCase(primary, secondary))
    {
        vec3::FieldVector3 secPositionVec = vec3::CoordVector3::convertToFieldVector(secondary.getPositionVector());
        return calculateMutualInductanceZAxis(primary, secondary, secPositionVec.zComponent, inductanceArguments, method);
    }
    else
        return calculateMutualInductanceGeneral(primary, secondary, inductanceArguments, method);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     PrecisionFactor precisionFactor, ComputeMethod method)
{
    bool zAxisCase = isZAxisCase(primary, secondary);
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, method, zAxisCase);

    return computeMutualInductance(primary, secondary, args, method);
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, CoilPairArguments inductanceArguments,
                                            ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, inductanceArguments, method) * 2*M_PI * sineFrequency;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, PrecisionFactor precisionFactor,
                                            ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, precisionFactor, method) * 2*M_PI * sineFrequency;
}