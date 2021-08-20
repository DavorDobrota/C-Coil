#include "Coil.h"
#include "LegendreMatrix.h"

#include <cmath>

double Coil::calculateAmpereForceZAxis(const Coil &primary, const Coil &secondary, double zDisplacement,
                                       CoilPairArguments forceArguments, ComputeMethod method)
{
    PrecisionArguments primaryPrecisionArguments = forceArguments.primaryPrecision;

    int lengthBlocks = forceArguments.secondaryPrecision.lengthBlockCount;
    int lengthIncrements = forceArguments.secondaryPrecision.lengthIncrementCount;

    int thicknessBlocks = forceArguments.secondaryPrecision.thicknessBlockCount;
    int thicknessIncrements = forceArguments.secondaryPrecision.thicknessIncrementCount;

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
                            incrementPositionR * 0.25 *
                            Legendre::weightsMatrix[maxLengthIndex][zIndex] *
                            Legendre::weightsMatrix[maxThicknessIndex][rIndex]);
                }
            }
        }
    }

    double ampereForce = 0.0;

    std::vector<double> fieldH = primary.computeAllBFieldX(positionVectors, primaryPrecisionArguments, method);

    for (int i = 0; i < fieldH.size(); ++i)
    {
        ampereForce += fieldH[i] * weights[i];
    }
    ampereForce /= (lengthBlocks * thicknessBlocks);
    return (-1) * ampereForce * 2*M_PI * secondary.numOfTurns * secondary.current;
}

std::pair<vec3::FieldVector3, vec3::FieldVector3>
Coil::calculateAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                  CoilPairArguments forceArguments, ComputeMethod method)
{
    std::vector<double> forceAndTorqueComponents(6);

    vec3::FieldVector3 displacementVec = vec3::CoordVector3::convertToFieldVector(secondary.getPositionVector());

    double xDisplacement = displacementVec.xComponent;
    double yDisplacement = displacementVec.yComponent;
    double zDisplacement = displacementVec.zComponent;
    double alphaAngle = secondary.yAxisAngle;
    double betaAngle = secondary.zAxisAngle;

    PrecisionArguments primaryPrecisionArguments = forceArguments.primaryPrecision;

    int lengthBlocks = forceArguments.secondaryPrecision.lengthBlockCount;
    int lengthIncrements = forceArguments.secondaryPrecision.lengthIncrementCount;

    int thicknessBlocks = forceArguments.secondaryPrecision.thicknessBlockCount;
    int thicknessIncrements = forceArguments.secondaryPrecision.thicknessIncrementCount;

    int angularBlocks = forceArguments.secondaryPrecision.angularBlockCount;
    int angularIncrements = forceArguments.secondaryPrecision.angularIncrementCount;

    int numElements = lengthBlocks * lengthIncrements * thicknessBlocks * thicknessIncrements * angularBlocks * angularIncrements;

    double ringIntervalSize = 2 * M_PI;

    std::vector<std::pair<vec3::FieldVector3, vec3::FieldVector3>> unitRingValues =
            calculateRingIncrementPosition(angularBlocks, angularIncrements, alphaAngle, betaAngle, ringIntervalSize);

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int maxLengthIndex = lengthIncrements - 1;
    int maxThicknessIndex = thicknessIncrements - 1;
    int maxAngularIncrementIndex = angularIncrements - 1;

    double lengthBlockSize = secondary.length / lengthBlocks;
    double thicknessBlockSize = secondary.thickness / thicknessBlocks;

    std::vector<vec3::CoordVector3> positionVectors;

    std::vector<double> radii;
    std::vector<double> weights;

    positionVectors.reserve(numElements);
    radii.reserve(numElements);
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

                            double displacementX = xDisplacement + lengthDisplacement * sin(alphaAngle) * cos(betaAngle) +
                                                   ringRadius * unitRingValues[phiPosition].first.xComponent;

                            double displacementY = yDisplacement - lengthDisplacement * sin(alphaAngle) * sin(betaAngle) +
                                                   ringRadius * unitRingValues[phiPosition].first.yComponent;

                            double displacementZ = zDisplacement + lengthDisplacement * cos(alphaAngle) +
                                                   ringRadius * unitRingValues[phiPosition].first.zComponent;

                            positionVectors.emplace_back(vec3::CARTESIAN, displacementX, displacementY, displacementZ);

                            radii.push_back(ringRadius);
                            weights.push_back(
                                    0.125 * 2*M_PI * ringRadius *
                                    Legendre::weightsMatrix[maxLengthIndex][zIndex] *
                                    Legendre::weightsMatrix[maxThicknessIndex][rIndex] *
                                    Legendre::weightsMatrix[maxAngularIncrementIndex][phiIndex]);
                        }
                    }
                }
            }
        }
    }
    std::vector<vec3::FieldVector3> magneticFields = primary.computeAllBFieldComponents(positionVectors,
                                                                                        primaryPrecisionArguments, method);

    vec3::FieldVector3 forceVector;
    vec3::FieldVector3 torqueVector;

    for (int i = 0; i < numElements; ++i)
    {
        int p = i % (angularBlocks * angularIncrements);

        vec3::FieldVector3 tempForce = vec3::FieldVector3::crossProduct(unitRingValues[p].second, magneticFields[i]);
        tempForce *= weights[i];
        forceVector += tempForce;

        vec3::FieldVector3 tempTorque = vec3::FieldVector3::crossProduct(unitRingValues[p].first, tempForce);
        torqueVector += tempTorque * radii[i];
    }
    double forceFactor = (secondary.current * secondary.numOfTurns) / (lengthBlocks * thicknessBlocks * angularBlocks);

    forceVector *= forceFactor;
    torqueVector *= forceFactor;

    return std::make_pair(forceVector, torqueVector);

}
