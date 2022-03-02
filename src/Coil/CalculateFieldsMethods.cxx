#include "Coil.h"


double Coil::calculateAPotential(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    if (useFastMethod)
        return calculateAPotentialFast(zAxis, rPolar, usedPrecision);

    return calculateAPotentialSlow(zAxis, rPolar, usedPrecision);
}

std::pair<double, double> Coil::calculateBField(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    if (useFastMethod)
        return calculateBFieldFast(zAxis, rPolar, usedPrecision);

    return calculateBFieldSlow(zAxis, rPolar, usedPrecision);
}

std::vector<double> Coil::calculateBGradient(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    if (useFastMethod)
        return calculateBGradientFast(zAxis, rPolar, usedPrecision);

    return calculateBGradientSlow(zAxis, rPolar, usedPrecision);
}
