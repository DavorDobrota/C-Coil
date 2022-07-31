#include "Coil.h"

#include <sstream>


namespace
{
    const double g_minRelativePrecision = 1.0;
    const double g_maxRelativePrecision = 15.0;
    const double g_defaultRelativePrecision = 5.0;
}

PrecisionFactor::PrecisionFactor() : PrecisionFactor(g_defaultRelativePrecision) {}

PrecisionFactor::PrecisionFactor(double relativePrecision)
{
    if (relativePrecision < g_minRelativePrecision || relativePrecision > g_maxRelativePrecision)
        PrecisionFactor::relativePrecision = g_defaultRelativePrecision;
    else
        PrecisionFactor::relativePrecision = relativePrecision;
}

PrecisionFactor::operator std::string() const
{
    std::stringstream output;

    output << "PrecisionFactor(relative_precision=" << relativePrecision << ")";

    return output.str();
}
