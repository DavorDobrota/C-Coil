#include "Coil.h"
#include "LegendreMatrix.h"

#define _USE_MATH_DEFINES
#include <math.h>


double Coil::calculateMutualInductanceGeneral(const Coil &primary, const Coil &secondary,
                                              const CoilPairArguments &inductanceArguments, ComputeMethod computeMethod)
{
    vec3::Vector3 displacementVec = secondary.getPositionVector();

    double xDisplacement = displacementVec.x;
    double yDisplacement = displacementVec.y;
    double zDisplacement = displacementVec.z;
    double alphaAngle = secondary.yAxisAngle;
    double betaAngle = secondary.zAxisAngle;

    PrecisionArguments primaryPrecisionArguments = inductanceArguments.primaryPrecision;

    int lengthBlocks = inductanceArguments.secondaryPrecision.lengthBlocks;
    int lengthIncrements = inductanceArguments.secondaryPrecision.lengthIncrements;

    int thicknessBlocks = inductanceArguments.secondaryPrecision.thicknessBlocks;
    int thicknessIncrements = inductanceArguments.secondaryPrecision.thicknessIncrements;

    int angularBlocks = inductanceArguments.secondaryPrecision.angularBlocks;
    int angularIncrements = inductanceArguments.secondaryPrecision.angularIncrements;

    int numElements = lengthBlocks * lengthIncrements * thicknessBlocks * thicknessIncrements * angularBlocks * angularIncrements;

    std::vector<std::pair<vec3::Vector3, vec3::Vector3>> unitRingValues =
            calculateRingIncrementPosition(angularBlocks, angularIncrements, alphaAngle, betaAngle);

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int maxLengthIndex = lengthIncrements - 1;
    int maxThicknessIndex = thicknessIncrements - 1;
    int maxAngularIncrementIndex = angularIncrements - 1;

    double lengthBlockSize = secondary.length / lengthBlocks;
    double thicknessBlockSize = secondary.thickness / thicknessBlocks;

    vec3::Vector3Array positionVectors;
    std::vector<double> weights;

    positionVectors.reserve(numElements);
    weights.reserve(numElements);

    for (int zBlock = 0; zBlock < lengthBlocks; ++zBlock)
    {
        for (int rBlock = 0; rBlock < thicknessBlocks; ++rBlock)
        {
            double zBlockPosition = (-1) * (secondary.length * 0.5) + lengthBlockSize * (zBlock + 0.5);
            double rBlockPosition = secondary.innerRadius + thicknessBlockSize * (rBlock + 0.5);

            for (int zIndex = 0; zIndex < lengthIncrements; ++zIndex)
            {
                for (int rIndex = 0; rIndex < thicknessIncrements; ++rIndex)
                {
                    double ringRadius = rBlockPosition +
                                        (thicknessBlockSize * 0.5) * Legendre::positionMatrix[maxThicknessIndex][rIndex];

                    double lengthDisplacement = zBlockPosition +
                                                (lengthBlockSize * 0.5) * Legendre::positionMatrix[maxLengthIndex][zIndex];

                    for (int phiBlock = 0; phiBlock < angularBlocks; ++phiBlock)
                    {
                        for (int phiIndex = 0; phiIndex < angularIncrements; ++phiIndex)
                        {
                            int phiPosition = phiBlock * angularIncrements + phiIndex;

                            double displacementX = xDisplacement +
                                                   lengthDisplacement * std::sin(alphaAngle) * std::cos(betaAngle) +
                                                   ringRadius * unitRingValues[phiPosition].first.x;

                            double displacementY = yDisplacement +
                                                   lengthDisplacement * std::sin(alphaAngle) * std::sin(betaAngle) +
                                                   ringRadius * unitRingValues[phiPosition].first.y;

                            double displacementZ = zDisplacement +
                                                   lengthDisplacement * std::cos(alphaAngle) +
                                                   ringRadius * unitRingValues[phiPosition].first.z;

                            positionVectors.append(displacementX, displacementY, displacementZ);

                            weights.emplace_back(
                                0.125 * ringRadius *
                                Legendre::weightsMatrix[maxLengthIndex][zIndex] *
                                Legendre::weightsMatrix[maxThicknessIndex][rIndex] *
                                Legendre::weightsMatrix[maxAngularIncrementIndex][phiIndex]
                            );
                        }
                    }
                }
            }
        }
    }
    vec3::Vector3Array potentialVectors = primary.computeAllAPotentialVectors(
        positionVectors, primaryPrecisionArguments, computeMethod
    );
    double mutualInductance = 0.0;

    for (int i = 0; i < numElements; ++i)
    {
        int p = i % (angularBlocks * angularIncrements);
        mutualInductance += vec3::Vector3::scalarProduct(potentialVectors[i], unitRingValues[p].second) * weights[i];
    }

    mutualInductance /= (lengthBlocks * thicknessBlocks * angularBlocks);
    return mutualInductance * 2*M_PI * secondary.numOfTurns / primary.current;
}
