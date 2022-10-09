#include "Coil.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace
{
    const double g_MiReduced = 0.0000001;
}


void Coil::calculateMagneticMoment()
{
    magneticMoment = M_PI * current * numOfTurns *
                     (innerRadius * innerRadius + innerRadius * thickness + thickness * thickness / 3);
}

void Coil::calculateAverageWireThickness()
{
    averageWireThickness = std::sqrt(length * thickness / numOfTurns);
}

void Coil::calculateResistance()
{
    double wireRadius = averageWireThickness * 0.5;
    double ohmicResistance = wireResistivity * numOfTurns * 2*M_PI *
                             (innerRadius + thickness * 0.5) / (wireRadius * wireRadius * M_PI);
    double skinDepth = std::sqrt(wireResistivity / (M_PI * sineFrequency * g_MiReduced));

    double ohmicSurface = M_PI * wireRadius * wireRadius;
    double sineSurface = 2*M_PI * (
            skinDepth * skinDepth * (exp(-wireRadius / skinDepth) - 1) +
            skinDepth * wireRadius);

    resistance = ohmicResistance * (ohmicSurface / sineSurface);
}

void Coil::calculateReactance()
{
    reactance = selfInductance * 2*M_PI * sineFrequency;
}

void Coil::calculateImpedance()
{
    calculateResistance();
    calculateReactance();
    impedance = std::sqrt(resistance * resistance + reactance * reactance);
}

void Coil::calculateCoilType()
{
    if (thickness / innerRadius < g_thinCoilApproximationRatio && length / innerRadius < g_thinCoilApproximationRatio)
        coilType = CoilType::FILAMENT;

    else if (thickness / innerRadius < g_thinCoilApproximationRatio)
        coilType = CoilType::THIN;

    else if (length / innerRadius < g_thinCoilApproximationRatio)
        coilType = CoilType::FLAT;

    else
        coilType = CoilType::RECTANGULAR;
}

void Coil::calculateTransformationMatrices()
{

    double cosY = std::cos(yAxisAngle); double sinY = std::sin(yAxisAngle);
    double cosZ = std::cos(zAxisAngle); double sinZ = std::sin(zAxisAngle);

    transformationMatrix = vec3::Matrix3(cosZ * cosZ * cosY - sinZ * sinZ, -sinZ * cosZ * cosY - sinZ * cosZ, cosZ * sinY,
                                         sinZ * cosZ * cosY + sinZ * cosZ, cosZ * cosZ - sinZ * sinZ * cosY, sinZ * sinY,
                                         -sinY * cosZ, sinY * sinZ, cosY);

    inverseTransformationMatrix = vec3::Matrix3(cosZ * cosZ * cosY - sinZ * sinZ, sinZ * cosZ * cosY + sinZ * cosZ, -cosZ * sinY,
                                                -sinZ * cosZ * cosY - sinZ * cosZ, cosZ * cosZ - sinZ * sinZ * cosY, sinZ * sinY,
                                                sinY * cosZ, sinY * sinZ, cosY);
}
