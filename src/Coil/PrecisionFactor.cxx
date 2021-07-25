#include "Coil.h"

namespace
{
    const double g_minRelativePrecision = 1.0;
    const double g_maxRelativePrecision = 8.0;
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