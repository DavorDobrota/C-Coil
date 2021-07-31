#include "Coil.h"

double Coil::computeAmpereForceZAxis(const Coil &primary, const Coil &secondary, double zDisplacement,
                                     PrecisionFactor precisionFactor, ComputeMethod method)
{
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, method,
                                                                   false);
    return computeAmpereForceZAxis(primary, secondary, zDisplacement, args, method);
}

double Coil::computeAmpereForceZAxis(const Coil &primary, const Coil &secondary, double zDisplacement,
                                     CoilPairArguments forceArguments, ComputeMethod method)
{
    return calculateAmpereForceZAxis(primary, secondary, zDisplacement, forceArguments, method);
}

std::vector<double> Coil::computeAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                                    double zDisplacement, double rDisplacement,
                                                    PrecisionFactor precisionFactor, ComputeMethod method)
{
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, method);
    return computeAmpereForceGeneral(primary, secondary, zDisplacement,rDisplacement, args, method);
}

std::vector<double> Coil::computeAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                                    double zDisplacement, double rDisplacement,
                                                    CoilPairArguments forceArguments, ComputeMethod method)
{
    return calculateAmpereForceGeneral(primary, secondary, zDisplacement, rDisplacement,
                                       0.0, 0.0, forceArguments, method);
}

std::vector<double> Coil::computeAmpereForceGeneral(const Coil &primary, const Coil &secondary, double zDisplacement,
                                                    double rDisplacement, double alphaAngle,
                                                    PrecisionFactor precisionFactor, ComputeMethod method)
{
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, method);
    return computeAmpereForceGeneral(primary, secondary, zDisplacement, rDisplacement, alphaAngle, args, method);
}

std::vector<double> Coil::computeAmpereForceGeneral(const Coil &primary, const Coil &secondary, double zDisplacement,
                                                    double rDisplacement, double alphaAngle,
                                                    CoilPairArguments forceArguments, ComputeMethod method)
{
    return calculateAmpereForceGeneral(primary, secondary, zDisplacement, rDisplacement,
                                       alphaAngle, 0.0, forceArguments, method);
}

std::vector<double> Coil::computeAmpereForceGeneral(const Coil &primary, const Coil &secondary, double zDisplacement,
                                                    double rDisplacement, double alphaAngle, double betaAngle,
                                                    PrecisionFactor precisionFactor, ComputeMethod method)
{
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, method);
    return computeAmpereForceGeneral(primary, secondary, zDisplacement, rDisplacement, alphaAngle, betaAngle, args, method);
}

std::vector<double> Coil::computeAmpereForceGeneral(const Coil &primary, const Coil &secondary, double zDisplacement,
                                                    double rDisplacement, double alphaAngle, double betaAngle,
                                                    CoilPairArguments forceArguments, ComputeMethod method)
{
    return calculateAmpereForceGeneral(primary, secondary, zDisplacement, rDisplacement,
                                       alphaAngle, betaAngle, forceArguments, method);
}