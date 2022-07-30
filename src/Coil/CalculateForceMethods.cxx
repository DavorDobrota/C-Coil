#include "Coil.h"
#include "LegendreMatrix.h"
#include "ThreadPool.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <numeric>

namespace
{
    const double g_MiReduced = 0.0000001;
    threadPool::ThreadPoolControl g_threadPool;
}


double Coil::calculateAmpereForceZAxisSlow(const Coil &primary, const Coil &secondary, double zDisplacement,
                                           CoilPairArguments forceArguments, ComputeMethod computeMethod)
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
                    double incrementPositionZ = zDisplacement + zBlockPosition +
                                                (lengthBlockSize * 0.5) * Legendre::positionMatrix[maxLengthIndex][zIndex];
                    double incrementPositionR = rBlockPosition +
                                                (thicknessBlockSize * 0.5) * Legendre::positionMatrix[maxThicknessIndex][rIndex];

                    positionVectors.append(incrementPositionR, 0.0, incrementPositionZ);

                    weights.emplace_back(
                            incrementPositionR * 0.25 *
                            Legendre::weightsMatrix[maxLengthIndex][zIndex] *
                            Legendre::weightsMatrix[maxThicknessIndex][rIndex]);
                }
            }
        }
    }

    double ampereForce = 0.0;
    vec3::Vector3Array fieldH = primary.computeAllBFieldVectors(positionVectors, primaryPrecisionArguments, computeMethod);

    for (int i = 0; i < fieldH.size(); ++i)
    {
        ampereForce += fieldH[i].x * weights[i];
    }
    ampereForce /= (lengthBlocks * thicknessBlocks);
    return (-1.0) * ampereForce * 2*M_PI * secondary.numOfTurns * secondary.current;
}

double Coil::calculateAmpereForceZAxisFast(const Coil &primary, const Coil &secondary, double zDisplacement,
                                           CoilPairArguments forceArguments, ComputeMethod computeMethod)
{

    double ampereForce = 0.0;

    double thicknessBlock = primary.thickness / forceArguments.primaryPrecision.thicknessBlockCount;
    double angularBlock = M_PI / forceArguments.primaryPrecision.angularBlockCount;
    double radialBlock = secondary.thickness / forceArguments.secondaryPrecision.thicknessBlockCount;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int thicknessIncrements = forceArguments.primaryPrecision.thicknessIncrementCount - 1;
    int angularIncrements = forceArguments.primaryPrecision.angularIncrementCount - 1;
    int radialIncrements = forceArguments.secondaryPrecision.thicknessIncrementCount - 1;

    // multiplication by 2 because cosine is an even function and by 0.125 for a triple change of interval (3 times 1/2)
    double constant = g_MiReduced * primary.currentDensity * secondary.currentDensity * thicknessBlock * angularBlock * radialBlock * 0.125;

    std::vector<std::vector<double>> cosPhiPrecomputeMat(forceArguments.primaryPrecision.angularBlockCount);

    for (int indBlockPhi = 0; indBlockPhi < forceArguments.primaryPrecision.angularBlockCount; ++indBlockPhi)
    {
        double blockPositionPhi = angularBlock * (indBlockPhi + 0.5);
        cosPhiPrecomputeMat[indBlockPhi].resize(angularIncrements + 1);

        for (int incPhi = 0; incPhi <= angularIncrements; ++incPhi)
        {
            double incrementPositionFi = blockPositionPhi +
                                         (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][incPhi];
            cosPhiPrecomputeMat[indBlockPhi][incPhi] = std::cos(incrementPositionFi);
        }
    }

    zDisplacement -= primary.getPositionVector().z;

    double constZ1 = zDisplacement + secondary.length * 0.5 + primary.length * 0.5;
    double constZ2 = zDisplacement + secondary.length * 0.5 - primary.length * 0.5;
    double constZ3 = zDisplacement - secondary.length * 0.5 - primary.length * 0.5;
    double constZ4 = zDisplacement - secondary.length * 0.5 + primary.length * 0.5;

    auto calculate = [&](int threadIndex, int startIndex, int endIndex, double &result)
    {
        for (int indBlockR = 0; indBlockR < forceArguments.secondaryPrecision.thicknessBlockCount; ++indBlockR)
        {
            double blockPositionR = secondary.innerRadius + radialBlock * (indBlockR + 0.5);

            for (int incR = startIndex; incR < endIndex; ++incR)
            {
                double incrementPositionR = blockPositionR +
                                            (radialBlock * 0.5) *
                                            Legendre::positionMatrix[radialIncrements][incR];

                double incrementWeightR = Legendre::weightsMatrix[radialIncrements][incR];

                for (int indBlockT = 0; indBlockT < forceArguments.primaryPrecision.thicknessBlockCount; ++indBlockT)
                {
                    double blockPositionT = primary.innerRadius + thicknessBlock * (indBlockT + 0.5);

                    for (int incT = 0; incT <= thicknessIncrements; ++incT)
                    {
                        double incrementPositionT = blockPositionT +
                                                    (thicknessBlock * 0.5) *
                                                    Legendre::positionMatrix[thicknessIncrements][incT];

                        double incrementWeightT = Legendre::weightsMatrix[thicknessIncrements][incT];

                        double tempConst = constant * incrementWeightR * incrementWeightT;
                        double tempConstA = 2.0 * incrementPositionT * incrementPositionR;
                        double tempConstB = incrementPositionT * incrementPositionT + incrementPositionR * incrementPositionR;

                        for (int indBlockPhi = 0; indBlockPhi < forceArguments.primaryPrecision.angularBlockCount; ++indBlockPhi)
                        {
                            for (int incPhi = 0; incPhi <= angularIncrements; ++incPhi)
                            {
                                double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incPhi];
                                double cosinePhi = cosPhiPrecomputeMat[indBlockPhi][incPhi];

                                double tempConstC = tempConstB - tempConstA * cosinePhi;
                                double tempConstD = std::sqrt(tempConstC);

                                double tempConstE1 = constZ1 / tempConstD;
                                double tempConstE2 = constZ2 / tempConstD;
                                double tempConstE3 = constZ3 / tempConstD;
                                double tempConstE4 = constZ4 / tempConstD;

                                double tempConstF1 = std::sqrt(tempConstE1 * tempConstE1 + 1.0);
                                double tempConstF2 = std::sqrt(tempConstE2 * tempConstE2 + 1.0);
                                double tempConstF3 = std::sqrt(tempConstE3 * tempConstE3 + 1.0);
                                double tempConstF4 = std::sqrt(tempConstE4 * tempConstE4 + 1.0);

                                double tempConstG1 = tempConstE1 + tempConstF1;
                                double tempConstG2 = tempConstE2 + tempConstF2;
                                double tempConstG3 = tempConstE3 + tempConstF3;
                                double tempConstG4 = tempConstE4 + tempConstF4;

                                double tempConstH = std::log(tempConstG2 * tempConstG4 / (tempConstG1 * tempConstG3));

                                result += tempConst * incrementWeightFi * tempConstA * cosinePhi * tempConstH;
                            }
                        }
                    }
                }
            }
        }
        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    if(computeMethod == CPU_ST)
    {
        calculate(0, 0, radialIncrements + 1, ampereForce);
        g_threadPool.getCompletedTasks().store(0ull);
    }
    else
    {
        int incrementCount = radialIncrements + 1;
        int threadCount = std::max(primary.getThreadCount(), secondary.getThreadCount());
        int mean = incrementCount / threadCount;

        g_threadPool.setTaskCount(threadCount);
        g_threadPool.getCompletedTasks().store(0ull);

        std::vector<double> results(threadCount);

        int delegatedIncrements = 0;
        for(int i = 0; i < threadCount; i++)
        {
            int remainingIncrements = incrementCount - delegatedIncrements;

            if(remainingIncrements % (threadCount - i) == 0)
            {
                g_threadPool.push(
                        calculate,
                        delegatedIncrements, delegatedIncrements + mean,
                        std::ref(results[i])
                );
                delegatedIncrements += mean;
            }
            else
            {
                g_threadPool.push(
                        calculate,
                        delegatedIncrements, delegatedIncrements + mean + 1,
                        std::ref(results[i])
                );
                delegatedIncrements += mean + 1;
            }
        }

        g_threadPool.synchronizeThreads();

        ampereForce = std::accumulate(results.begin(), results.end(), 0.0);
    }

    return (-1.0) * ampereForce * 2*M_PI;
}


std::pair<vec3::Vector3, vec3::Vector3>
Coil::calculateAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                  CoilPairArguments forceArguments, ComputeMethod computeMethod)
{
    std::vector<double> forceAndTorqueComponents(6);

    vec3::Vector3 displacementVec = secondary.getPositionVector();

    double xDisplacement = displacementVec.x;
    double yDisplacement = displacementVec.y;
    double zDisplacement = displacementVec.z;
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

    std::vector<std::pair<vec3::Vector3, vec3::Vector3>> unitRingValues =
            calculateRingIncrementPosition(angularBlocks, angularIncrements, alphaAngle, betaAngle);

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int maxLengthIndex = lengthIncrements - 1;
    int maxThicknessIndex = thicknessIncrements - 1;
    int maxAngularIncrementIndex = angularIncrements - 1;

    double lengthBlockSize = secondary.length / lengthBlocks;
    double thicknessBlockSize = secondary.thickness / thicknessBlocks;

    vec3::Vector3Array positionVectors;

    std::vector<double> radii;
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
                                                   ringRadius * unitRingValues[phiPosition].first.x;

                            double displacementY = yDisplacement + lengthDisplacement * sin(alphaAngle) * sin(betaAngle) +
                                                   ringRadius * unitRingValues[phiPosition].first.y;

                            double displacementZ = zDisplacement + lengthDisplacement * cos(alphaAngle) +
                                                   ringRadius * unitRingValues[phiPosition].first.z;

                            positionVectors.append(displacementX, displacementY, displacementZ);

                            weights.emplace_back(
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
    vec3::Vector3Array magneticFields = primary.computeAllBFieldVectors(positionVectors,
                                                                        primaryPrecisionArguments,
                                                                        computeMethod);

    vec3::Vector3 forceVector;
    vec3::Vector3 torqueVector;

    for (int i = 0; i < numElements; ++i)
    {
        int p = i % (angularBlocks * angularIncrements);

        vec3::Vector3 tempForce = vec3::Vector3::crossProduct(unitRingValues[p].second, magneticFields[i]);
        tempForce *= weights[i];
        forceVector += tempForce;

        vec3::Vector3 tempTorque = vec3::Vector3::crossProduct(positionVectors[i] - displacementVec, tempForce);
        torqueVector += tempTorque;
    }
    double forceFactor = 2*M_PI * (secondary.current * secondary.numOfTurns) / (lengthBlocks * thicknessBlocks * angularBlocks);

    forceVector *= forceFactor;
    torqueVector *= forceFactor;

    return {forceVector, torqueVector};

}
