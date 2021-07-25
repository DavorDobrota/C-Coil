#include "Coil.h"

#include <cmath>

void Coil::calculateRingIncrementPosition(int angularBlocks, int angularIncrements,
                                          double alpha, double beta, double ringIntervalSize,
                                          std::vector<double> &ringXPosition,
                                          std::vector<double> &ringYPosition,
                                          std::vector<double> &ringZPosition,
                                          std::vector<double> &ringXTangent,
                                          std::vector<double> &ringYTangent,
                                          std::vector<double> &ringZTangent)
{
    ringXPosition.resize(0);
    ringYPosition.resize(0);
    ringZPosition.resize(0);

    ringXTangent.resize(0);
    ringYTangent.resize(0);
    ringZTangent.resize(0);

    double angularBlock = ringIntervalSize / angularBlocks;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    angularBlocks--;
    angularIncrements--;

    for (int phiBlock = 0; phiBlock <= angularBlocks; ++phiBlock)
    {
        double blockPositionPhi = angularBlock * (phiBlock + 0.5);

        for (int phiIndex = 0; phiIndex <= angularIncrements; ++phiIndex)
        {
            // M_PI/2 added to readjust to an even interval so a shortcut can be used
            double phi = M_PI/2 + blockPositionPhi +
                         (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][phiIndex];

            ringXPosition.push_back(cos(beta) * cos(phi) - sin(beta) * cos(alpha) * sin(phi));
            ringYPosition.push_back(sin(beta) * cos(phi) + cos(beta) * cos(alpha) * sin(phi));
            ringZPosition.push_back(sin(alpha) * sin(phi));

            ringXTangent.push_back((-1) * cos(beta) * sin(phi) - sin(beta) * cos(alpha) * cos(phi));
            ringYTangent.push_back((-1) * sin(beta) * sin(phi) + cos(beta) * cos(alpha) * cos(phi));
            ringZTangent.push_back(sin(alpha) * cos(phi));
        }
    }
}


double Coil::calculateMutualInductanceZAxis(const Coil &primary, const Coil &secondary, double zDisplacement,
                                            CoilPairArguments inductanceArguments, ComputeMethod method)
{
    PrecisionArguments primaryPrecisionArguments = inductanceArguments.primaryPrecision;

    int lengthBlocks = inductanceArguments.secondaryPrecision.lengthBlockCount;
    int lengthIncrements = inductanceArguments.secondaryPrecision.lengthIncrementCount;

    int thicknessBlocks = inductanceArguments.secondaryPrecision.thicknessBlockCount;
    int thicknessIncrements = inductanceArguments.secondaryPrecision.thicknessIncrementCount;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int maxLengthIndex = lengthIncrements - 1;
    int maxThicknessIndex = thicknessIncrements - 1;

    double lengthBlockSize = secondary.length / lengthBlocks;
    double thicknessBlockSize = secondary.thickness / thicknessBlocks;

    std::vector<double> zPositions;
    std::vector<double> rPositions;

    std::vector<double> weights;

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

                    zPositions.push_back(incrementPositionZ);
                    rPositions.push_back(incrementPositionR);

                    weights.push_back(
                            0.25 * Legendre::weightsMatrix[maxLengthIndex][zIndex] *
                            Legendre::weightsMatrix[maxThicknessIndex][rIndex]);
                }
            }
        }
    }

    std::vector<double> potentialA;
    double mutualInductance = 0.0;

    primary.computeAllAPotentialAbs(zPositions, rPositions, potentialA, primaryPrecisionArguments, method);

    for (int i = 0; i < potentialA.size(); ++i)
    {
        mutualInductance += 2*M_PI * rPositions[i] * potentialA[i] * weights[i];
    }
    mutualInductance /= (lengthBlocks * thicknessBlocks);
    return mutualInductance * secondary.numOfTurns / primary.current;
}

double Coil::calculateMutualInductanceGeneral(const Coil &primary, const Coil &secondary,
                                              double zDisplacement, double rDisplacement,
                                              double alphaAngle, double betaAngle,
                                              CoilPairArguments inductanceArguments, ComputeMethod method)
{
    if (rDisplacement == 0.0 && alphaAngle == 0.0)
    {
        return calculateMutualInductanceZAxis(primary, secondary, zDisplacement, inductanceArguments, method);
    }
    else {
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

        if (rDisplacement == 0.0 || alphaAngle == 0.0 || betaAngle == 0.0)
            ringIntervalSize = M_PI;
        else
            ringIntervalSize = 2 * M_PI;

        std::vector<double> unitRingPointsX, unitRingPointsY, unitRingPointsZ;
        std::vector<double> unitRingTangentsX, unitRingTangentsY, unitRingTangentsZ;

        calculateRingIncrementPosition(angularBlocks, angularIncrements, alphaAngle, betaAngle, ringIntervalSize,
                                       unitRingPointsX, unitRingPointsY, unitRingPointsZ,
                                       unitRingTangentsX, unitRingTangentsY, unitRingTangentsZ);

        // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
        int maxLengthIndex = lengthIncrements - 1;
        int maxThicknessIndex = thicknessIncrements - 1;
        int maxAngularIncrementIndex = angularIncrements - 1;

        double lengthBlockSize = secondary.length / lengthBlocks;
        double thicknessBlockSize = secondary.thickness / thicknessBlocks;

        std::vector<double> zPositions;
        std::vector<double> rPositions;

        std::vector<double> weights;

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

                                double displacementX = lengthDisplacement * sin(alphaAngle) * sin(betaAngle) +
                                                       ringRadius * unitRingPointsX[phiPosition];

                                double displacementY = rDisplacement - lengthDisplacement * sin(alphaAngle) * cos(betaAngle) +
                                                       ringRadius * unitRingPointsY[phiPosition];

                                double displacementZ = zDisplacement + lengthDisplacement * cos(alphaAngle) +
                                                       ringRadius * unitRingPointsZ[phiPosition];

                                zPositions.push_back(displacementZ);
                                rPositions.push_back(sqrt(displacementX * displacementX + displacementY * displacementY));

                                double rhoAngle = atan2(displacementY, displacementX);

                                double orientationFactor =
                                        - sin(rhoAngle) * unitRingTangentsX[phiPosition] +
                                        cos(rhoAngle) * unitRingTangentsY[phiPosition];

                                weights.push_back(
                                        0.125 * orientationFactor * 2 * M_PI * ringRadius *
                                        Legendre::weightsMatrix[maxLengthIndex][zIndex] *
                                        Legendre::weightsMatrix[maxThicknessIndex][rIndex] *
                                        Legendre::weightsMatrix[maxAngularIncrementIndex][phiIndex]);
                            }
                        }
                    }
                }
            }
        }
        std::vector<double> potentialArray;
        double mutualInductance = 0.0;

        primary.computeAllAPotentialAbs(zPositions, rPositions, potentialArray, primaryPrecisionArguments, method);

        for (int i = 0; i < numElements; ++i)
            mutualInductance += potentialArray[i] * weights[i];

        mutualInductance /= (lengthBlocks * thicknessBlocks * angularBlocks);
        return mutualInductance * secondary.numOfTurns / primary.current;
    }
}

CoilPairArguments Coil::calculateAppropriateMInductanceArguments(const Coil &primary, const Coil &secondary,
                                                                 PrecisionFactor precisionFactor,
                                                                 ComputeMethod method, bool isGeneral)
{
    if (!isGeneral)
    {
        if (method == GPU)
            return CoilPairArguments::getMInductanceArgumentsZGPU(primary, secondary, precisionFactor);
        else
            return CoilPairArguments::getMInductanceArgumentsZCPU(primary, secondary, precisionFactor);
    }
    else
    {
        if (method == GPU)
            return CoilPairArguments::getMInductanceArgumentsGeneralGPU(primary, secondary, precisionFactor);
        else
            return CoilPairArguments::getMInductanceArgumentsGeneralCPU(primary, secondary, precisionFactor);
    }
}


