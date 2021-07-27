#include "Coil.h"

#include <cmath>

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary, double zDisplacement,
                                     CoilPairArguments inductanceArguments, ComputeMethod method)
{
    return calculateMutualInductanceGeneral(primary, secondary, zDisplacement,
                                            0.0, 0.0, 0.0, inductanceArguments, method);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary, double zDisplacement,
                                     PrecisionFactor precisionFactor, ComputeMethod method)
{
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, method,
                                                                   false);
    return computeMutualInductance(primary, secondary, zDisplacement, args, method);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     double zDisplacement, double rDisplacement,
                                     CoilPairArguments inductanceArguments, ComputeMethod method)
{
    return calculateMutualInductanceGeneral(primary, secondary, zDisplacement, rDisplacement,
                                            0.0, 0.0, inductanceArguments, method);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     double zDisplacement, double rDisplacement,
                                     PrecisionFactor precisionFactor, ComputeMethod method)
{
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, method);
    return computeMutualInductance(primary, secondary, zDisplacement, rDisplacement, args, method);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     double zDisplacement, double rDisplacement, double alphaAngle,
                                     CoilPairArguments inductanceArguments, ComputeMethod method)
{
    return calculateMutualInductanceGeneral(primary, secondary, zDisplacement, rDisplacement, alphaAngle,
                                            0.0, inductanceArguments, method);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     double zDisplacement, double rDisplacement, double alphaAngle,
                                     PrecisionFactor precisionFactor, ComputeMethod method)
{
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, method);
    return computeMutualInductance(primary, secondary, zDisplacement, rDisplacement, alphaAngle, args, method);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     double zDisplacement, double rDisplacement, double alphaAngle, double betaAngle,
                                     CoilPairArguments inductanceArguments, ComputeMethod method)
{
    return calculateMutualInductanceGeneral(primary, secondary, zDisplacement, rDisplacement, alphaAngle, betaAngle,
                                            inductanceArguments, method);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     double zDisplacement, double rDisplacement, double alphaAngle, double betaAngle,
                                     PrecisionFactor precisionFactor, ComputeMethod method)
{
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, method);
    return computeMutualInductance(primary, secondary, zDisplacement, rDisplacement, alphaAngle, betaAngle, args, method);
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement,
                                            PrecisionFactor precisionFactor, ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, precisionFactor, method) *
           2 * M_PI * sineFrequency;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement,
                                            CoilPairArguments inductanceArguments, ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, inductanceArguments, method) *
           2 * M_PI * sineFrequency;;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement, double rDisplacement,
                                            PrecisionFactor precisionFactor, ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, rDisplacement, precisionFactor, method) *
           2 * M_PI * sineFrequency;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement, double rDisplacement,
                                            CoilPairArguments inductanceArguments, ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, rDisplacement, inductanceArguments, method) *
           2 * M_PI * sineFrequency;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement, double rDisplacement,
                                            double alphaAngle, PrecisionFactor precisionFactor, ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, rDisplacement, alphaAngle,
                                   precisionFactor, method) * 2 * M_PI * sineFrequency;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement, double rDisplacement,
                                            double alphaAngle, CoilPairArguments inductanceArguments,
                                            ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, rDisplacement, alphaAngle,
                                   inductanceArguments, method) * 2 * M_PI * sineFrequency;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement, double rDisplacement,
                                            double alphaAngle, double betaAngle, PrecisionFactor precisionFactor,
                                            ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, rDisplacement, alphaAngle, betaAngle,
                                   precisionFactor, method) * 2 * M_PI * sineFrequency;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement, double rDisplacement,
                                            double alphaAngle, double betaAngle, CoilPairArguments inductanceArguments,
                                            ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, rDisplacement, alphaAngle, betaAngle,
                                   inductanceArguments, method) * 2 * M_PI * sineFrequency;
}