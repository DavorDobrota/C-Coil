#include "Coil.h"


double Coil::calculateAPotential(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    // TODO - implement a solution which can completely eliminate overhead
    if (usedPrecision.lengthIncrementCount == 1)
        return calculateAPotentialFast(zAxis, rPolar, usedPrecision);

    return calculateAPotentialSlow(zAxis, rPolar, usedPrecision);
}

std::pair<double, double> Coil::calculateBField(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    // TODO - implement a solution which can completely eliminate overhead
    if (usedPrecision.lengthIncrementCount == 1)
        return calculateBFieldFast(zAxis, rPolar, usedPrecision);

    return calculateBFieldSlow(zAxis, rPolar, usedPrecision);
}

std::vector<double> Coil::calculateBGradient(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    // TODO - implement a solution which can completely eliminate overhead
    if (usedPrecision.lengthIncrementCount == 1)
        return calculateBGradientFast(zAxis, rPolar, usedPrecision);

    return calculateBGradientSlow(zAxis, rPolar, usedPrecision);
}
