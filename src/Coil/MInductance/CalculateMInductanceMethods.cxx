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


double Coil::calculateMutualInductanceZAxisSlow(const Coil &primary, const Coil &secondary, double zDisplacement,
                                                const CoilPairArguments &inductanceArguments, ComputeMethod computeMethod)
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
                            0.25 * incrementPositionR *
                            Legendre::weightsMatrix[maxLengthIndex][zIndex] *
                            Legendre::weightsMatrix[maxThicknessIndex][rIndex]);
                }
            }
        }
    }
    vec3::Vector3Array potentialA = primary.computeAllAPotentialVectors(positionVectors,
                                                                        primaryPrecisionArguments,
                                                                        computeMethod);
    double mutualInductance = 0.0;

    for (int i = 0; i < potentialA.size(); ++i)
    {
        mutualInductance += potentialA[i].abs() * weights[i];
    }
    mutualInductance /= (lengthBlocks * thicknessBlocks);
    return mutualInductance * 2*M_PI * secondary.numOfTurns / primary.current;
}

double Coil::calculateMutualInductanceZAxisFast(const Coil &primary, const Coil &secondary, double zDisplacement,
                                                const CoilPairArguments &inductanceArguments, ComputeMethod computeMethod)
{
    double mutualInductance = 0.0;

    double thicknessBlock = primary.thickness / inductanceArguments.primaryPrecision.thicknessBlockCount;
    double angularBlock = M_PI / inductanceArguments.primaryPrecision.angularBlockCount;
    double radialBlock = secondary.thickness / inductanceArguments.secondaryPrecision.thicknessBlockCount;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int thicknessIncrements = inductanceArguments.primaryPrecision.thicknessIncrementCount - 1;
    int angularIncrements = inductanceArguments.primaryPrecision.angularIncrementCount - 1;
    int radialIncrements = inductanceArguments.secondaryPrecision.thicknessIncrementCount - 1;

    // multiplication by 2 because cosine is an even function and by 0.125 for a triple change of interval (3 times 1/2)
    double constant = g_MiReduced * primary.currentDensity * secondary.currentDensity * thicknessBlock * angularBlock * radialBlock * 0.125;

    std::vector<std::vector<double>> cosPhiPrecomputeMat(inductanceArguments.primaryPrecision.angularBlockCount);

    for (int indBlockPhi = 0; indBlockPhi < inductanceArguments.primaryPrecision.angularBlockCount; ++indBlockPhi)
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

    double constZ1Squared = constZ1 * constZ1;
    double constZ2Squared = constZ2 * constZ2;
    double constZ3Squared = constZ3 * constZ3;
    double constZ4Squared = constZ4 * constZ4;

    auto calculate = [&](int threadIndex, int startIndex, int endIndex, double &result)
    {
        for (int indBlockR = 0; indBlockR < inductanceArguments.secondaryPrecision.thicknessBlockCount; ++indBlockR)
        {
            double blockPositionR = secondary.innerRadius + radialBlock * (indBlockR + 0.5);

            for (int incR = startIndex; incR < endIndex; ++incR)
            {
                double incrementPositionR = blockPositionR +
                                            (radialBlock * 0.5) * Legendre::positionMatrix[radialIncrements][incR];

                double incrementWeightR = Legendre::weightsMatrix[radialIncrements][incR];

                for (int indBlockT = 0; indBlockT < inductanceArguments.primaryPrecision.thicknessBlockCount; ++indBlockT)
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

                        for (int indBlockPhi = 0; indBlockPhi < inductanceArguments.primaryPrecision.angularBlockCount; ++indBlockPhi)
                        {
                            for (int incPhi = 0; incPhi <= angularIncrements; ++incPhi)
                            {
                                double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incPhi];
                                double cosinePhi = cosPhiPrecomputeMat[indBlockPhi][incPhi];

                                double tempConstC = tempConstB - tempConstA * cosinePhi;
                                double tempConstD = std::sqrt(tempConstC);

                                double tempConstE1 = std::sqrt(tempConstC + constZ1Squared);
                                double tempConstE2 = std::sqrt(tempConstC + constZ2Squared);
                                double tempConstE3 = std::sqrt(tempConstC + constZ3Squared);
                                double tempConstE4 = std::sqrt(tempConstC + constZ4Squared);

                                double tempConstF1 = constZ1 / tempConstD;
                                double tempConstF2 = constZ2 / tempConstD;
                                double tempConstF3 = constZ3 / tempConstD;
                                double tempConstF4 = constZ4 / tempConstD;

                                double tempConstG1 = std::sqrt(tempConstF1 * tempConstF1 + 1.0);
                                double tempConstG2 = std::sqrt(tempConstF2 * tempConstF2 + 1.0);
                                double tempConstG3 = std::sqrt(tempConstF3 * tempConstF3 + 1.0);
                                double tempConstG4 = std::sqrt(tempConstF4 * tempConstF4 + 1.0);

                                double tempConstH1 = constZ1 * std::log(tempConstF1 + tempConstG1);
                                double tempConstH2 = constZ2 * std::log(tempConstF2 + tempConstG2);
                                double tempConstH3 = constZ3 * std::log(tempConstF3 + tempConstG3);
                                double tempConstH4 = constZ4 * std::log(tempConstF4 + tempConstG4);

                                result +=
                                        tempConst * incrementWeightFi * tempConstA * cosinePhi *
                                        (tempConstH1 - tempConstH2 + tempConstH3 - tempConstH4 -
                                        tempConstE1 + tempConstE2 - tempConstE3 + tempConstE4);

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
        calculate(0, 0, radialIncrements + 1, mutualInductance);
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

        mutualInductance = std::accumulate(results.begin(), results.end(), 0.0);
    }

    return mutualInductance * 2*M_PI / primary.current;
}

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

    int lengthBlocks = inductanceArguments.secondaryPrecision.lengthBlockCount;
    int lengthIncrements = inductanceArguments.secondaryPrecision.lengthIncrementCount;

    int thicknessBlocks = inductanceArguments.secondaryPrecision.thicknessBlockCount;
    int thicknessIncrements = inductanceArguments.secondaryPrecision.thicknessIncrementCount;

    int angularBlocks = inductanceArguments.secondaryPrecision.angularBlockCount;
    int angularIncrements = inductanceArguments.secondaryPrecision.angularIncrementCount;

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
    vec3::Vector3Array potentialVectors = primary.computeAllAPotentialVectors(positionVectors,
                                                                              primaryPrecisionArguments,
                                                                              computeMethod);
    double mutualInductance = 0.0;

    for (int i = 0; i < numElements; ++i)
    {
        int p = i % (angularBlocks * angularIncrements);
        mutualInductance += vec3::Vector3::scalarProduct(potentialVectors[i], unitRingValues[p].second) * weights[i];
    }

    mutualInductance /= (lengthBlocks * thicknessBlocks * angularBlocks);
    return mutualInductance * 2*M_PI * secondary.numOfTurns / primary.current;
}

double Coil::calculateSelfInductance(CoilPairArguments inductanceArguments) const
{
    double calculatedSelfInductance = 0.0;

    double thicknessBlock = thickness / inductanceArguments.primaryPrecision.thicknessBlockCount;
    double angularBlock = M_PI / inductanceArguments.primaryPrecision.angularBlockCount;
    double radialBlock = 1.0 / inductanceArguments.secondaryPrecision.thicknessBlockCount;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int thicknessIncrements = inductanceArguments.primaryPrecision.thicknessIncrementCount - 1;
    int angularIncrements = inductanceArguments.primaryPrecision.angularIncrementCount - 1;
    int radialIncrements = inductanceArguments.secondaryPrecision.thicknessIncrementCount - 1;

    // multiplication by 2 because cosine is an even function and by 0.125 for a triple change of interval (3 times 1/2)
    double constant = g_MiReduced * currentDensity * thicknessBlock * angularBlock * radialBlock * 2 * 0.125;
    double lengthSquared = length * length;

    std::vector<std::vector<double>> cosPhiPrecomputeMat(inductanceArguments.primaryPrecision.angularBlockCount);

    for (int indBlockPhi = 0; indBlockPhi < inductanceArguments.primaryPrecision.angularBlockCount; ++indBlockPhi)
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

                    double tempConst = constant * incrementWeightR * incrementWeightT;
                    double tempConstA = 2 * incrementPositionT * incrementPositionR;
                    double tempConstB = incrementPositionT * incrementPositionT + incrementPositionR * incrementPositionR;

                    for (int indBlockFi = 0; indBlockFi < inductanceArguments.primaryPrecision.angularBlockCount; ++indBlockFi)
                    {
                        for (int incFi = 0; incFi <= angularIncrements; ++incFi)
                        {
                            double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incFi];
                            double cosinePhi = cosPhiPrecomputeMat[indBlockFi][incFi];

                            double tempConstC = tempConstB - tempConstA * cosinePhi;

                            double tempConstD = std::sqrt(tempConstC);
                            double tempConstE = std::sqrt(tempConstC + lengthSquared);

                            double tempConstF = length / tempConstD;
                            double tempConstG = std::sqrt(tempConstF * tempConstF + 1.0);
                            double tempConstH = std::log(tempConstF + tempConstG);

                            calculatedSelfInductance +=
                                    tempConst * incrementWeightFi *
                                    tempConstA * cosinePhi * (tempConstD - tempConstE + length * tempConstH);
                        }
                    }
                }
            }
        }
    }
    return calculatedSelfInductance * 2*M_PI * numOfTurns / (current * length);
}
