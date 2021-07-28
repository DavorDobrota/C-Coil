#include "Coil.h"
#include "LegendreMatrix.h"

#include <cmath>
#include <cstdio>

double Coil::calculateAmpereForceZAxis(const Coil &primary, const Coil &secondary, double zDisplacement,
                                       CoilPairArguments forceArguments, ComputeMethod method)
{
    PrecisionArguments primaryPrecisionArguments = forceArguments.primaryPrecision;

    int lengthBlocks = forceArguments.secondaryPrecision.lengthBlockCount;
    int lengthIncrements = forceArguments.secondaryPrecision.lengthIncrementCount;

    int thicknessBlocks = forceArguments.secondaryPrecision.thicknessBlockCount;
    int thicknessIncrements = forceArguments.secondaryPrecision.thicknessIncrementCount;

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

    std::vector<double> fieldH;
    double ampereForce = 0.0;

    primary.computeAllBFieldH(zPositions, rPositions, fieldH, primaryPrecisionArguments, method);

    for (int i = 0; i < fieldH.size(); ++i)
    {
        ampereForce += 2 * M_PI * rPositions[i] * fieldH[i] * weights[i];
    }
    ampereForce /= (lengthBlocks * thicknessBlocks);
    return (-1) * ampereForce * secondary.numOfTurns * secondary.current;
}

std::vector<double> Coil::calculateAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                                      double zDisplacement, double rDisplacement,
                                                      double alphaAngle, double betaAngle,
                                                      CoilPairArguments forceArguments, ComputeMethod method) {
    std::vector<double> forceAndTorqueComponents;

    if (rDisplacement == 0.0 && alphaAngle == 0.0) {
        // to enforce standard 6 arguments returned, zeros are passed
        forceAndTorqueComponents.push_back(0);
        forceAndTorqueComponents.push_back(0);

        forceAndTorqueComponents.push_back(
                calculateAmpereForceZAxis(primary, secondary, zDisplacement, forceArguments, method));

        forceAndTorqueComponents.push_back(0);
        forceAndTorqueComponents.push_back(0);
        forceAndTorqueComponents.push_back(0);

        return forceAndTorqueComponents;
    }
    else
    {
        PrecisionArguments primaryPrecisionArguments = forceArguments.primaryPrecision;

        int lengthBlocks = forceArguments.secondaryPrecision.lengthBlockCount;
        int lengthIncrements = forceArguments.secondaryPrecision.lengthIncrementCount;

        int thicknessBlocks = forceArguments.secondaryPrecision.thicknessBlockCount;
        int thicknessIncrements = forceArguments.secondaryPrecision.thicknessIncrementCount;

        //multiplying by 2 because integration needs to be from 0 to 2PI
        int angularBlocks = forceArguments.secondaryPrecision.angularBlockCount;
        int angularIncrements = forceArguments.secondaryPrecision.angularIncrementCount;

        int numElements = lengthBlocks * lengthIncrements * thicknessBlocks * thicknessIncrements * angularBlocks *
                          angularIncrements;

        double ringIntervalSize = 2 * M_PI;

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
        std::vector<double> phiPositions;

        std::vector<double> radii;
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
                                            (thicknessBlockSize * 0.5) *
                                            Legendre::positionMatrix[maxThicknessIndex][rIndex];

                        double lengthDisplacement = zBlockPosition +
                                                    (lengthBlockSize * 0.5) *
                                                    Legendre::positionMatrix[maxLengthIndex][zIndex];

                        for (int phiBlock = 0; phiBlock < angularBlocks; ++phiBlock)
                        {
                            for (int phiIndex = 0; phiIndex < angularIncrements; ++phiIndex)
                            {
                                int phiPosition = phiBlock * angularIncrements + phiIndex;

                                double displacementX = lengthDisplacement * sin(alphaAngle) * sin(betaAngle) +
                                                       ringRadius * unitRingPointsX[phiPosition];

                                double displacementY =
                                        rDisplacement - lengthDisplacement * sin(alphaAngle) * cos(betaAngle) +
                                        ringRadius * unitRingPointsY[phiPosition];

                                double displacementZ = zDisplacement + lengthDisplacement * cos(alphaAngle) +
                                                       ringRadius * unitRingPointsZ[phiPosition];

                                zPositions.push_back(displacementZ);
                                rPositions.push_back(
                                        sqrt(displacementX * displacementX + displacementY * displacementY));
                                phiPositions.push_back(atan2(displacementY, displacementX));

                                radii.push_back(ringRadius);
                                weights.push_back(
                                        0.125 * 2 * M_PI * ringRadius *
                                        Legendre::weightsMatrix[maxLengthIndex][zIndex] *
                                        Legendre::weightsMatrix[maxThicknessIndex][rIndex] *
                                        Legendre::weightsMatrix[maxAngularIncrementIndex][phiIndex]);
                            }
                        }
                    }
                }
            }
        }
        std::vector<double> magneticFieldX;
        std::vector<double> magneticFieldY;
        std::vector<double> magneticFieldZ;

        primary.computeAllBFieldComponents(zPositions, rPositions, phiPositions,
                                           magneticFieldX, magneticFieldY, magneticFieldZ,
                                           primaryPrecisionArguments, method);

        double forceX = 0.0, forceY = 0.0, forceZ = 0.0;
        double torqueX = 0.0, torqueY = 0.0, torqueZ = 0.0;

        for (int i = 0; i < numElements; ++i)
        {
            int p = i % (angularBlocks * angularIncrements);

            double tempForceX = weights[i] *
                                (unitRingTangentsY[p] * magneticFieldZ[i] - unitRingTangentsZ[p] * magneticFieldY[i]);

            double tempForceY = weights[i] *
                                (unitRingTangentsZ[p] * magneticFieldX[i] - unitRingTangentsX[p] * magneticFieldZ[i]);

            double tempForceZ = weights[i] *
                                (unitRingTangentsX[p] * magneticFieldY[i] - unitRingTangentsY[p] * magneticFieldX[i]);

        //    printf("%.15f %.15f %.15f\n", zPositions[i], rPositions[i], phiPositions[i]);

            forceX += tempForceX;
            forceY += tempForceY;
            forceZ += tempForceZ;

            torqueX += radii[i] * (tempForceY * unitRingPointsZ[p] - tempForceZ * unitRingPointsY[p]);
            torqueY += radii[i] * (tempForceZ * unitRingPointsX[p] - tempForceX * unitRingPointsZ[p]);
            torqueZ += radii[i] * (tempForceX * unitRingPointsY[p] - tempForceY * unitRingPointsX[p]);
        }

        forceX /= (lengthBlocks * thicknessBlocks * angularBlocks);
        forceY /= (lengthBlocks * thicknessBlocks * angularBlocks);
        forceZ /= (lengthBlocks * thicknessBlocks * angularBlocks);

        torqueX /= (lengthBlocks * thicknessBlocks * angularBlocks);
        torqueY /= (lengthBlocks * thicknessBlocks * angularBlocks);
        torqueZ /= (lengthBlocks * thicknessBlocks * angularBlocks);

        forceAndTorqueComponents.push_back(forceX * secondary.current * secondary.numOfTurns);
        forceAndTorqueComponents.push_back(forceY * secondary.current * secondary.numOfTurns);
        forceAndTorqueComponents.push_back(forceZ * secondary.current * secondary.numOfTurns);

        forceAndTorqueComponents.push_back(torqueX * secondary.current * secondary.numOfTurns);
        forceAndTorqueComponents.push_back(torqueY * secondary.current * secondary.numOfTurns);
        forceAndTorqueComponents.push_back(torqueZ * secondary.current * secondary.numOfTurns);

        return forceAndTorqueComponents;
    }
}
