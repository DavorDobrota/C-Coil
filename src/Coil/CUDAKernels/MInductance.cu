#include "CoilAcceleration.h"

#include "Timing.h"
#include "CUDAUtils/ErrorCheck/CUDAErrorCheck.h"
#include "CUDAUtils/MemoryManagement/GPUMemoryManagement.h"

#include <cstdio>


__global__
void CalculateMutualInductanceConfigurations(long long configCount, long long pointCount,
                                             const CoilPairArgumentsData *pairData,
                                             const CoilPairPositionData *configArr, TYPE *inductanceArr)
{
    unsigned int index = threadIdx.x;
    long long global_index = (long long) blockIdx.x * blockDim.x + index;

    if(global_index >= configCount * pointCount)
        return;

    __shared__ CoilPairArgumentsData coilPair;
    coilPair = *pairData;

    int pairIndex = int(global_index % configCount);
    int lengthIndex = int((global_index / configCount) % coilPair.secLengthIncrements);
    int thicknessIndex = int(((global_index / configCount) / coilPair.secLengthIncrements) % coilPair.secThicknessIncrements);
    int angularIndex = int((((global_index / configCount) / coilPair.secLengthIncrements) / coilPair.secThicknessIncrements) % coilPair.secAngularIncrements);

    CoilPairPositionData position = configArr[pairIndex];

    TYPE lengthPosition = coilPair.secLength * 0.5f * coilPair.secLengthPositionArray[lengthIndex];
    TYPE thicknessPosition = coilPair.secInnerRadius +
                             coilPair.secThickness * 0.5f * (1.0f + coilPair.secThicknessPositionArray[thicknessIndex]);
    TYPE angularPosition = PI * (1.0f + coilPair.secAngularPositionArray[angularIndex]);

    TYPE sinAlpha = sin(position.secAlphaAngle);
    TYPE cosAlpha = cos(position.secAlphaAngle);
    TYPE sinBeta = sin(position.secBetaAngle);
    TYPE cosBeta = cos(position.secBetaAngle);
    TYPE sinPhi = sin(angularPosition);
    TYPE cosPhi = cos(angularPosition);

    TYPE ringPositionX = cosBeta * cosAlpha * cosPhi - sinBeta * sinPhi;
    TYPE ringPositionY = sinBeta * cosAlpha * cosPhi + cosBeta * sinPhi;
    TYPE ringPositionZ = (-1.0f) * sinAlpha * cosPhi;

    TYPE ringTangentX = (-1.0f) * cosBeta * cosAlpha * sinPhi - sinBeta * cosPhi;
    TYPE ringTangentY = (-1.0f) * sinBeta * cosAlpha * sinPhi + cosBeta * cosPhi;
    TYPE ringTangentZ = sinAlpha * sinPhi;

    TYPE posX = position.secPositionVector[0] - position.primPositionVector[0] +
                lengthPosition * sinAlpha * cosBeta + thicknessPosition * ringPositionX;
    TYPE posY = position.secPositionVector[1] - position.primPositionVector[1] +
                lengthPosition * sinAlpha * sinBeta + thicknessPosition * ringPositionY;
    TYPE posZ = position.secPositionVector[2] - position.primPositionVector[2] +
                lengthPosition * cosAlpha + thicknessPosition * ringPositionZ;

    TYPE sinY = sin(position.primAlphaAngle);
    TYPE cosY = cos(position.primAlphaAngle);
    TYPE sinZ = sin(position.primBetaAngle);
    TYPE cosZ = cos(position.primBetaAngle);

    TYPE x = posX * (cosZ * cosZ * cosY - sinZ * sinZ) +
             posY * (sinZ * cosZ * cosY + sinZ * cosZ) +
             posZ * (-1.0f * cosZ * sinY);
    TYPE y = posX * (-1.0f * sinZ * cosZ * cosY - sinZ * cosZ) +
             posY * (cosZ * cosZ - sinZ * sinZ * cosY) +
             posZ * (sinZ * sinY);
    TYPE z = posX * (sinY * cosZ) +
             posY * (sinY * sinZ) +
             posZ * cosY;

    TYPE zCoord = z;
    TYPE rCoord = sqrt(x * x + y * y);
    TYPE phiCord = atan2(y, x);

    TYPE potential = 0.0f;

    TYPE topEdge = zCoord + 0.5f * coilPair.primLength;
    TYPE bottomEdge = zCoord - 0.5f * coilPair.primLength;

    if (coilPair.useFastMethod)
    {
        for (int incT = 0; incT < coilPair.primThicknessIncrements; ++incT)
        {
            TYPE incrementPositionT = coilPair.primInnerRadius +
                                      0.5f * coilPair.primThickness * (1.0f + coilPair.primThicknessPositionArray[incT]);

            TYPE tempConstA = incrementPositionT * incrementPositionT + rCoord * rCoord;
            TYPE tempConstB = 2.0f * incrementPositionT * rCoord;

            for (int incF = 0; incF < coilPair.primAngularIncrements; ++incF)
            {
                TYPE cosinePhi = coilPair.primCosPrecomputeArray[incF];

                TYPE tempConstC = rsqrt(tempConstA - tempConstB * cosinePhi);

                TYPE tempConstD1 = topEdge * tempConstC;
                TYPE tempConstD2 = bottomEdge * tempConstC;

                TYPE tempConstE1 = sqrt(tempConstD1 * tempConstD1 + 1.0f);
                TYPE tempConstE2 = sqrt(tempConstD2 * tempConstD2 + 1.0f);

                TYPE tempConstF = log((tempConstE1 + tempConstD1) / (tempConstE2 + tempConstD2));

                potential += coilPair.constFactor *
                             coilPair.primThicknessWeightArray[incT] * coilPair.primAngularWeightArray[incF] *
                             incrementPositionT * cosinePhi * tempConstF;
            }
        }
    }
    else
    {
        for (int incT = 0; incT < coilPair.primThicknessIncrements; ++incT)
        {
            TYPE incrementPositionT = coilPair.primInnerRadius +
                                      0.5f * coilPair.primThickness * (1.0f + coilPair.primThicknessPositionArray[incT]);

            TYPE tempConstA = incrementPositionT * incrementPositionT + rCoord * rCoord + zCoord * zCoord;
            TYPE tempConstB = 2.0f * incrementPositionT * rCoord;

            for (int incF = 0; incF < coilPair.primAngularIncrements; ++incF)
            {
                TYPE cosinePhi = coilPair.primCosPrecomputeArray[incF];

                TYPE tempConstC = rsqrt(tempConstA - tempConstB * cosinePhi);

                potential += coilPair.constFactor *
                             coilPair.primThicknessWeightArray[incT] * coilPair.primAngularWeightArray[incF] *
                             incrementPositionT * cosinePhi * tempConstC;
            }
        }
    }

    TYPE xPot = (-1.0f) * sin(phiCord) * potential;
    TYPE yPot = potential * cos(phiCord);
    TYPE zPot = 0.0f;

    TYPE potentialX = xPot * (cosZ * cosZ * cosY - sinZ * sinZ) +
                      yPot * (-1.0f * sinZ * cosZ * cosY - sinZ * cosZ) +
                      zPot * (cosZ * sinY);
    TYPE potentialY = xPot * (sinZ * cosZ * cosY + sinZ * cosZ) +
                      yPot * (cosZ * cosZ - sinZ * sinZ * cosY) +
                      zPot * (sinZ * sinY);
    TYPE potentialZ = xPot * (-1.0f * sinY * cosZ) +
                      yPot * (sinY * sinZ) +
                      zPot * cosY;

    TYPE weight = 0.125f * thicknessPosition *
                  coilPair.secLengthWeightArray[lengthIndex] *
                  coilPair.secThicknessWeightArray[thicknessIndex] *
                  coilPair.secAngularWeightArray[angularIndex];

    TYPE mInductance = weight * coilPair.correctionFactor *
                       (potentialX * ringTangentX + potentialY * ringTangentY + potentialZ * ringTangentZ);

    atomicAdd(&inductanceArr[pairIndex], mInductance);
}

namespace
{
    CoilPairPositionData *g_configArr = nullptr;
    CoilPairArgumentsData *g_coilPair = nullptr;
    TYPE *g_inductanceArr = nullptr;

    void getBuffers(long long configs)
    {
        std::vector<void*> buffers = GPUMem::getBuffers(
                { configs * (long long)sizeof(CoilPairPositionData),
                  1 * (long long)sizeof(CoilPairArgumentsData),
                  configs * (long long)sizeof(TYPE)}
        );

        g_configArr = static_cast<CoilPairPositionData *>(buffers[0]);
        g_coilPair = static_cast<CoilPairArgumentsData *>(buffers[1]);
        g_inductanceArr = static_cast<TYPE *>(buffers[2]);
    }

    #if DEBUG_TIMINGS
        double g_duration;
    #endif
}

void Calculate_mutual_inductance_configurations(long long configCount, long long pointCount,
                                                const CoilPairArgumentsData *coilPair,
                                                const CoilPairPositionData *configArr,
                                                TYPE *inductanceArr)
{
    #if DEBUG_TIMINGS
        recordStartPoint();
        recordStartPoint();
    #endif

    int blocks = int(std::ceil(double((long long) configCount * pointCount) / NTHREADS));

    getBuffers(configCount);

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tResource startup:         %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    gpuErrchk(cudaMemcpy(g_configArr, configArr, configCount * sizeof(CoilPairPositionData), cudaMemcpyHostToDevice))
    gpuErrchk(cudaMemcpy(g_coilPair, coilPair, 1 * sizeof(CoilPairArgumentsData), cudaMemcpyHostToDevice))

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tMemory initialization:    %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    gpuErrchk(cudaMemset(g_inductanceArr, 0, configCount * sizeof(TYPE)))

    CalculateMutualInductanceConfigurations<<<blocks, NTHREADS>>>(
        configCount, pointCount, g_coilPair, g_configArr, g_inductanceArr
    );
    gpuErrchk(cudaDeviceSynchronize())

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tCalculations:             %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    if(inductanceArr != nullptr)
        gpuErrchk(cudaMemcpy(inductanceArr, g_inductanceArr, configCount * sizeof(TYPE), cudaMemcpyDeviceToHost))

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tWriting to output array:  %.9g s\n\n", g_duration);

        g_duration = getIntervalDuration();
        printf("\tDevice buffer size:       %.3lf MB\n", (configCount * double(sizeof(CoilPairPositionData) + sizeof(TYPE)) / 1.0e6));
        printf("\tTotal blocks:             %d\n", blocks);
        printf("\tThreads per calculation:  %i\n", NTHREADS);
        printf("\tTotal issued threads:     %i\n", NTHREADS * blocks);
        printf("\tTotal configurations:     %lli\n", configCount);
        printf("\tPoints per configuration: %lli\n", pointCount);
        printf("\tTotal calculations:       %lli\n", pointCount * configCount);
        printf("\n\tPerformance:              %.1lf kConfigs/s\n", double(0.001 * configCount / g_duration));
        printf("\tEffective Performance:    %.1lf kPoints/s\n", double(0.001 * pointCount * configCount / g_duration));
        printf("---------------------------------------------------\n\n");
    #endif
}