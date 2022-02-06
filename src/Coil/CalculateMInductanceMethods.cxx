#include "Coil.h"
#include "LegendreMatrix.h"
#include "PrecisionGlobalVars.h"
#include "Math/CustomMath.h"

#include <cmath>


namespace
{
    const double g_MiReduced = 0.0000001;
}

double Coil::calculateMutualInductanceZAxis(const Coil &primary, const Coil &secondary, double zDisplacement,
                                            CoilPairArguments inductanceArguments, ComputeMethod method)
{
    PrecisionArguments primaryPrecisionArguments = inductanceArguments.primaryPrecision;

    int lengthBlocks = inductanceArguments.secondaryPrecision.lengthBlockCount;
    int lengthIncrements = inductanceArguments.secondaryPrecision.lengthIncrementCount;

    int thicknessBlocks = inductanceArguments.secondaryPrecision.thicknessBlockCount;
    int thicknessIncrements = inductanceArguments.secondaryPrecision.thicknessIncrementCount;

    int numElements = lengthBlocks * lengthIncrements * thicknessBlocks * thicknessIncrements;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int maxLengthIndex = lengthIncrements - 1;
    int maxThicknessIndex = thicknessIncrements - 1;

    double lengthBlockSize = secondary.length / lengthBlocks;
    double thicknessBlockSize = secondary.thickness / thicknessBlocks;

    std::vector<vec3::CoordVector3> positionVectors;
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
                    double incrementPositionZ = zDisplacement + zBlockPosition +
                                                (lengthBlockSize * 0.5) * Legendre::positionMatrix[maxLengthIndex][zIndex];
                    double incrementPositionR = rBlockPosition +
                                                (thicknessBlockSize * 0.5) * Legendre::positionMatrix[maxThicknessIndex][rIndex];

                    positionVectors.emplace_back(vec3::CYLINDRICAL, incrementPositionZ, incrementPositionR, 0.0);

                    weights.push_back(
                            0.25 * incrementPositionR *
                            Legendre::weightsMatrix[maxLengthIndex][zIndex] *
                            Legendre::weightsMatrix[maxThicknessIndex][rIndex]);
                }
            }
        }
    }
    std::vector<double> potentialA = primary.computeAllAPotentialAbs(positionVectors, primaryPrecisionArguments, method);
    double mutualInductance = 0.0;

    for (int i = 0; i < potentialA.size(); ++i)
    {
        mutualInductance += potentialA[i] * weights[i];
    }
    mutualInductance /= (lengthBlocks * thicknessBlocks);
    return mutualInductance * 2*M_PI * secondary.numOfTurns / primary.current;
}

double Coil::calculateMutualInductanceGeneral(const Coil &primary, const Coil &secondary,
                                              CoilPairArguments inductanceArguments, ComputeMethod method)
{
    vec3::FieldVector3 displacementVec = vec3::CoordVector3::convertToFieldVector(secondary.getPositionVector());
    vec3::FieldVector3 offsetVec = vec3::CoordVector3::convertToFieldVector(primary.getPositionVector());

    double xDisplacement = displacementVec.xComponent;
    double yDisplacement = displacementVec.yComponent;
    double zDisplacement = displacementVec.zComponent;
    double alphaAngle = secondary.yAxisAngle;
    double betaAngle = secondary.zAxisAngle;

    vec3::FieldVector3 relativeVec = primary.inverseTransformationMatrix * (displacementVec - offsetVec);
    double relativeAlpha = primary.yAxisAngle - secondary.yAxisAngle;
    double relativeBeta = primary.zAxisAngle - secondary.zAxisAngle;

    PrecisionArguments primaryPrecisionArguments = inductanceArguments.primaryPrecision;

    int lengthBlocks = inductanceArguments.secondaryPrecision.lengthBlockCount;
    int lengthIncrements = inductanceArguments.secondaryPrecision.lengthIncrementCount;

    int thicknessBlocks = inductanceArguments.secondaryPrecision.thicknessBlockCount;
    int thicknessIncrements = inductanceArguments.secondaryPrecision.thicknessIncrementCount;

    int angularBlocks = inductanceArguments.secondaryPrecision.angularBlockCount;
    int angularIncrements = inductanceArguments.secondaryPrecision.angularIncrementCount;

    int numElements = lengthBlocks * lengthIncrements * thicknessBlocks * thicknessIncrements * angularBlocks * angularIncrements;

    // sometimes the function is even so a shortcut can be used to improve performance and efficiency
    double ringIntervalSize;

    if (relativeVec.xComponent / primary.innerRadius < g_zAxisApproximationRatio &&
        relativeVec.yComponent / primary.innerRadius < g_zAxisApproximationRatio ||
        relativeAlpha < g_zAxisApproximationRatio || relativeBeta < g_zAxisApproximationRatio)
    {
        ringIntervalSize = M_PI;
    }
    else
        ringIntervalSize = 2 * M_PI;

    std::vector<std::pair<vec3::FieldVector3, vec3::FieldVector3>> unitRingValues =
            calculateRingIncrementPosition(angularBlocks, angularIncrements, alphaAngle, betaAngle, ringIntervalSize);

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int maxLengthIndex = lengthIncrements - 1;
    int maxThicknessIndex = thicknessIncrements - 1;
    int maxAngularIncrementIndex = angularIncrements - 1;

    double lengthBlockSize = secondary.length / lengthBlocks;
    double thicknessBlockSize = secondary.thickness / thicknessBlocks;

    std::vector<vec3::CoordVector3> positionVectors;
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

                            double displacementX = xDisplacement + lengthDisplacement * sin(alphaAngle) * cos(betaAngle) +
                                                   ringRadius * unitRingValues[phiPosition].first.xComponent;

                            double displacementY = yDisplacement + lengthDisplacement * sin(alphaAngle) * sin(betaAngle) +
                                                   ringRadius * unitRingValues[phiPosition].first.yComponent;

                            double displacementZ = zDisplacement + lengthDisplacement * cos(alphaAngle) +
                                                   ringRadius * unitRingValues[phiPosition].first.zComponent;

                            positionVectors.emplace_back(vec3::CARTESIAN, displacementX, displacementY, displacementZ);

                            weights.push_back(
                                    0.125 * ringRadius *
                                    Legendre::weightsMatrix[maxLengthIndex][zIndex] *
                                    Legendre::weightsMatrix[maxThicknessIndex][rIndex] *
                                    Legendre::weightsMatrix[maxAngularIncrementIndex][phiIndex]);
                        }
                    }
                }
            }
        }
    }
    std::vector<vec3::FieldVector3> potentialVectors =
            primary.computeAllAPotentialComponents(positionVectors, primaryPrecisionArguments, method);

    double mutualInductance = 0.0;

    for (int i = 0; i < numElements; ++i)
    {
        int p = i % (angularBlocks * angularIncrements);
        mutualInductance += vec3::FieldVector3::scalarProduct(potentialVectors[i], unitRingValues[p].second) * weights[i];
    }

    mutualInductance /= (lengthBlocks * thicknessBlocks * angularBlocks);
    return mutualInductance * 2*M_PI * secondary.numOfTurns / primary.current;
}

double Coil::calculateSelfInductance(CoilPairArguments inductanceArguments, ComputeMethod method) const
{
    double calculatedSelfInductance = 0.0;

    double thicknessBlock = thickness / inductanceArguments.primaryPrecision.thicknessBlockCount;
    double angularBlock = M_PI / inductanceArguments.primaryPrecision.angularBlockCount;
    double radialBlock = 1.0 / inductanceArguments.secondaryPrecision.thicknessBlockCount;

    // initialising precompute array
    const int numPhiIncrements =
            inductanceArguments.primaryPrecision.angularBlockCount * inductanceArguments.primaryPrecision.angularIncrementCount;
    double cosPhiPrecomputeArr[numPhiIncrements];
    precomputeCosPhi(inductanceArguments.primaryPrecision.angularBlockCount,
                     inductanceArguments.primaryPrecision.angularIncrementCount,
                     cosPhiPrecomputeArr);

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int thicknessIncrements = inductanceArguments.primaryPrecision.thicknessIncrementCount - 1;
    int angularIncrements = inductanceArguments.primaryPrecision.angularIncrementCount - 1;
    int radialIncrements = inductanceArguments.secondaryPrecision.thicknessIncrementCount - 1;

    // multiplication by 2 because cosine is an even function and by 0.125 for a triple change of interval (3 times 1/2)
    double constant = g_MiReduced * currentDensity * thicknessBlock * angularBlock * radialBlock * 2 * 0.125;
    double lengthSquared = length * length;

    for (int indBlockR = 0; indBlockR < inductanceArguments.secondaryPrecision.thicknessBlockCount; ++indBlockR)
    {
        double blockPositionR = innerRadius + radialBlock * thickness * (indBlockR + 0.5);

        for (int incR = 0; incR <= radialIncrements; ++incR)
        {
            double incrementPositionR = blockPositionR +
                                        (radialBlock * thickness * 0.5) * Legendre::positionMatrix[radialIncrements][incR];

            double incrementWeightR = Legendre::weightsMatrix[radialIncrements][incR];

            for (int indBlockT = 0; indBlockT < inductanceArguments.primaryPrecision.thicknessBlockCount; ++indBlockT)
            {
                double blockPositionT = innerRadius + thicknessBlock * (indBlockT + 0.5);

                for (int incT = 0; incT <= thicknessIncrements; ++incT)
                {
                    double incrementPositionT = blockPositionT +
                                                (thicknessBlock * 0.5) *
                                                Legendre::positionMatrix[thicknessIncrements][incT];

                    double incrementWeightT = Legendre::weightsMatrix[thicknessIncrements][incT];

                    double tempConstA = 2 * incrementPositionT * incrementPositionR;
                    double tempConstB = incrementPositionT * incrementPositionT + incrementPositionR * incrementPositionR;

                    for (int indBlockFi = 0; indBlockFi < inductanceArguments.primaryPrecision.angularBlockCount; ++indBlockFi)
                    {
                        for (int incFi = 0; incFi <= angularIncrements; ++incFi)
                        {
                            double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incFi];

                            int arrPos = indBlockFi * (angularIncrements + 1) + incFi;
                            double cosinePhi = cosPhiPrecomputeArr[arrPos];
                            double tempConstC = tempConstB - tempConstA * cosinePhi;

                            double tempConstD = std::sqrt(tempConstC);
                            double tempConstE = std::sqrt(tempConstC + lengthSquared);

                            double tempConstF = length / tempConstD;
                            double tempConstG = std::sqrt(tempConstF * tempConstF + 1.0);
                            double tempConstH = LN(tempConstF + tempConstG);

                            calculatedSelfInductance +=
                                    constant * incrementWeightR * incrementWeightT * incrementWeightFi *
                                    tempConstA * cosinePhi * (tempConstD - tempConstE + length * tempConstH);

                        }
                    }
                }
            }
        }
    }
    return calculatedSelfInductance * 2*M_PI * numOfTurns / (current * length);
}
