#include "Coil.h"
#include "LegendreMatrix.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>

namespace
{
    const double g_MiReduced = 0.0000001;
}


std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
Coil::calculateRingIncrementPosition(int angularBlocks, int angularIncrements, double alpha, double beta)
{
    int numElements = angularBlocks * angularIncrements;

    std::vector<std::pair<vec3::Vector3, vec3::Vector3>> unitRingVector;
    unitRingVector.reserve(numElements);

    vec3::Vector3 ringPosition;
    vec3::Vector3 ringTangent;

    double angularBlock = 2*M_PI / angularBlocks;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    angularBlocks--;
    angularIncrements--;

    double sinA = std::sin(alpha); double cosA = std::cos(alpha);
    double sinB = std::sin(beta); double cosB = std::cos(beta);

    for (int phiBlock = 0; phiBlock <= angularBlocks; ++phiBlock)
    {
        double blockPositionPhi = angularBlock * (phiBlock + 0.5);

        for (int phiIndex = 0; phiIndex <= angularIncrements; ++phiIndex)
        {
            double phi = blockPositionPhi +
                         (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][phiIndex];

            double sinPhi = std::sin(phi); double cosPhi = std::cos(phi);

            ringPosition = vec3::Vector3(cosB * cosA * cosPhi - sinB * sinPhi,
                                         sinB * cosA * cosPhi + cosB * sinPhi,
                                         (-1) * sinA * cosPhi);

            ringTangent = vec3::Vector3((-1) * cosB * cosA * sinPhi - sinB * cosPhi,
                                        (-1) * sinB * cosA * sinPhi + cosB * cosPhi,
                                        sinA * sinPhi);

            unitRingVector.emplace_back(ringPosition, ringTangent);
        }
    }
    return unitRingVector;
}


bool Coil::isZAxisCase(const Coil &primary, const Coil &secondary)
{
    vec3::Vector3 primPositionVec = primary.getPositionVector();
    vec3::Vector3 secPositionVec = secondary.getPositionVector();

    if (std::abs(primPositionVec.x / primary.innerRadius) < g_zAxisApproximationRatio &&
        std::abs(primPositionVec.y / primary.innerRadius) < g_zAxisApproximationRatio &&
        std::abs(secPositionVec.x / primary.innerRadius) < g_zAxisApproximationRatio &&
        std::abs(secPositionVec.y / primary.innerRadius) < g_zAxisApproximationRatio &&
        std::abs(primary.yAxisAngle / (2 * M_PI)) < g_zAxisApproximationRatio &&
        std::abs(secondary.yAxisAngle / (2 * M_PI)) < g_zAxisApproximationRatio)
    {
        return true;
    }
    return false;
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"

void Coil::generateCoilData(CoilData &coilData, const PrecisionArguments &usedPrecision) const
{
    if (useFastMethod)
        coilData.constFactor = g_MiReduced * currentDensity * thickness * M_PI * 0.5;
    else
        coilData.constFactor = g_MiReduced * currentDensity * thickness * length * M_PI * 0.5;

    coilData.current = current;
    coilData.frequency = sineFrequency;
    coilData.useFastMethod = useFastMethod;

    coilData.innerRadius = innerRadius;
    coilData.thickness = thickness;
    coilData.length = length;

    coilData.lengthIncrements = usedPrecision.lengthIncrementCount;
    coilData.thicknessIncrements = usedPrecision.thicknessIncrementCount;
    coilData.angularIncrements = usedPrecision.angularIncrementCount;

    for (int i = 0; i < coilData.angularIncrements; ++i)
    {
        double phiPosition = M_PI_2 * (1.0 + Legendre::positionMatrix[coilData.angularIncrements - 1][i]);

        coilData.cosPrecomputeArray[i] = std::cos(phiPosition);
        coilData.angularWeightArray[i] = Legendre::weightsMatrix[coilData.angularIncrements - 1][i];
    }

    for (int i = 0; i < coilData.thicknessIncrements; ++i)
    {
        coilData.thicknessPositionArray[i] = Legendre::positionMatrix[coilData.thicknessIncrements - 1][i];
        coilData.thicknessWeightArray[i] = Legendre::weightsMatrix[coilData.thicknessIncrements - 1][i];
    }

    coilData.positionVector[0] = positionVector.x;
    coilData.positionVector[1] = positionVector.y;
    coilData.positionVector[2] = positionVector.z;

    coilData.transformArray[0] = transformationMatrix.xx;
    coilData.transformArray[1] = transformationMatrix.xy;
    coilData.transformArray[2] = transformationMatrix.xz;
    coilData.transformArray[3] = transformationMatrix.yx;
    coilData.transformArray[4] = transformationMatrix.yy;
    coilData.transformArray[5] = transformationMatrix.yz;
    coilData.transformArray[6] = transformationMatrix.zx;
    coilData.transformArray[7] = transformationMatrix.zy;
    coilData.transformArray[8] = transformationMatrix.zz;

    coilData.invTransformArray[0] = inverseTransformationMatrix.xx;
    coilData.invTransformArray[1] = inverseTransformationMatrix.xy;
    coilData.invTransformArray[2] = inverseTransformationMatrix.xz;
    coilData.invTransformArray[3] = inverseTransformationMatrix.yx;
    coilData.invTransformArray[4] = inverseTransformationMatrix.yy;
    coilData.invTransformArray[5] = inverseTransformationMatrix.yz;
    coilData.invTransformArray[6] = inverseTransformationMatrix.zx;
    coilData.invTransformArray[7] = inverseTransformationMatrix.zy;
    coilData.invTransformArray[8] = inverseTransformationMatrix.zz;
}


void Coil::generateCoilPairArgumentsData(const Coil &primary, const Coil &secondary,
                                         CoilPairArgumentsData &coilPairArgumentsData,
                                         const CoilPairArguments &inductanceArguments, bool forceCalculation)
{
    if (primary.useFastMethod)
        coilPairArgumentsData.constFactor = g_MiReduced * primary.currentDensity * primary.thickness * M_PI * 0.5;
    else
        coilPairArgumentsData.constFactor = g_MiReduced * M_PI * 0.5 *
                                            primary.currentDensity * primary.thickness * primary.length;

    coilPairArgumentsData.useFastMethod = primary.useFastMethod;

    if (!forceCalculation)
        coilPairArgumentsData.correctionFactor = 2*M_PI * secondary.numOfTurns / primary.current;
    else
        coilPairArgumentsData.correctionFactor = 2*M_PI * secondary.numOfTurns * secondary.current;

    coilPairArgumentsData.primInnerRadius = primary.innerRadius;
    coilPairArgumentsData.primThickness = primary.thickness;
    coilPairArgumentsData.primLength = primary.length;

    coilPairArgumentsData.primLengthIncrements = inductanceArguments.primaryPrecision.lengthIncrementCount;
    coilPairArgumentsData.primThicknessIncrements = inductanceArguments.primaryPrecision.thicknessIncrementCount;
    coilPairArgumentsData.primAngularIncrements = inductanceArguments.primaryPrecision.angularIncrementCount;

    coilPairArgumentsData.secInnerRadius = secondary.innerRadius;
    coilPairArgumentsData.secThickness = secondary.thickness;
    coilPairArgumentsData.secLength = secondary.length;

    coilPairArgumentsData.secLengthIncrements = inductanceArguments.secondaryPrecision.lengthIncrementCount;
    coilPairArgumentsData.secThicknessIncrements = inductanceArguments.secondaryPrecision.thicknessIncrementCount;
    coilPairArgumentsData.secAngularIncrements = inductanceArguments.secondaryPrecision.angularIncrementCount;

    for (int i = 0; i < inductanceArguments.primaryPrecision.angularIncrementCount; ++i)
    {
        double phiPosition =
                M_PI_2 * (1.0 + Legendre::positionMatrix[inductanceArguments.primaryPrecision.angularIncrementCount - 1][i]);

        coilPairArgumentsData.primCosPrecomputeArray[i] = std::cos(phiPosition);
        coilPairArgumentsData.primAngularWeightArray[i] =
                Legendre::weightsMatrix[inductanceArguments.primaryPrecision.angularIncrementCount - 1][i];
    }
    for (int i = 0; i < inductanceArguments.primaryPrecision.thicknessIncrementCount; ++i)
    {
        coilPairArgumentsData.primThicknessPositionArray[i] =
                Legendre::positionMatrix[inductanceArguments.primaryPrecision.thicknessIncrementCount - 1][i];
        coilPairArgumentsData.primThicknessWeightArray[i] =
                Legendre::weightsMatrix[inductanceArguments.primaryPrecision.thicknessIncrementCount - 1][i];
    }

    for (int i = 0; i < inductanceArguments.secondaryPrecision.angularIncrementCount; ++i)
    {
        coilPairArgumentsData.secAngularPositionArray[i] =
                Legendre::positionMatrix[inductanceArguments.secondaryPrecision.angularIncrementCount - 1][i];
        coilPairArgumentsData.secAngularWeightArray[i] =
                Legendre::weightsMatrix[inductanceArguments.secondaryPrecision.angularIncrementCount - 1][i];
    }
    for (int i = 0; i < inductanceArguments.secondaryPrecision.thicknessIncrementCount; ++i)
    {
        coilPairArgumentsData.secThicknessPositionArray[i] =
                Legendre::positionMatrix[inductanceArguments.secondaryPrecision.thicknessIncrementCount - 1][i];
        coilPairArgumentsData.secThicknessWeightArray[i] =
                Legendre::weightsMatrix[inductanceArguments.secondaryPrecision.thicknessIncrementCount - 1][i];
    }
    for (int i = 0; i < inductanceArguments.secondaryPrecision.lengthIncrementCount; ++i)
    {
        coilPairArgumentsData.secLengthPositionArray[i] =
                Legendre::positionMatrix[inductanceArguments.secondaryPrecision.lengthIncrementCount - 1][i];
        coilPairArgumentsData.secLengthWeightArray[i] =
                Legendre::weightsMatrix[inductanceArguments.secondaryPrecision.lengthIncrementCount - 1][i];
    }
}
#pragma clang diagnostic pop
