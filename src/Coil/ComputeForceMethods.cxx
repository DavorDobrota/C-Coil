#include "Coil.h"

double Coil::computeAmpereForce(const Coil &primary, const Coil &secondary, double zDisplacement,
                                PrecisionFactor precisionFactor, ComputeMethod method)
{
    auto args = CoilPairArguments::getAppropriateMInductanceArguments(primary, secondary, precisionFactor, method,
                                                                      false);
    return calculateAmpereForceZAxis(primary, secondary, zDisplacement, args, method);
}