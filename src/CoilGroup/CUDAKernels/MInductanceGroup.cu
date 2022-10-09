#include "CoilGroupAcceleration.h"

#include "Timing.h"
#include "CUDAUtils/ErrorCheck/CUDAErrorCheck.h"
#include "CUDAUtils/MemoryManagement/GPUMemoryManagement.h"

#include <cstdio>


__global__
void CalculateMutualInductanceConfigurationsGroup(long long coilIndex, long long configCount, long long pointCount,
                                                  SecondaryCoilData secondaryCoil,
                                                  const CoilData *primaryCoils,
                                                  const SecondaryCoilPositionData *secondaryPositions,
                                                  TYPE *inductanceArr)
{
    unsigned int index = threadIdx.x;
    long long globalIndex = (long long) blockIdx.x * blockDim.x + index;

    if(globalIndex >= configCount * pointCount)
        return;

    int configIndex = int(globalIndex % configCount);
    int lengthIndex = int((globalIndex / configCount) % secondaryCoil.lengthIncrements);
    int thicknessIndex = int(((globalIndex / configCount) / secondaryCoil.lengthIncrements) % secondaryCoil.thicknessIncrements);
    int angularIndex = int((((globalIndex / configCount) / secondaryCoil.lengthIncrements) / secondaryCoil.thicknessIncrements) % secondaryCoil.angularIncrements);

    __shared__ CoilData primCoil;
    primCoil = primaryCoils[coilIndex];

    SecondaryCoilPositionData position = secondaryPositions[configIndex];

    TYPE lengthPosition = secondaryCoil.length * 0.5f * secondaryCoil.lengthPositionArray[lengthIndex];
    TYPE thicknessPosition = secondaryCoil.innerRadius +
                             secondaryCoil.thickness * 0.5f * (1.0f + secondaryCoil.thicknessPositionArray[thicknessIndex]);
    TYPE angularPosition = PI * (1.0f + secondaryCoil.angularPositionArray[angularIndex]);

    TYPE sinAlpha = sin(position.alphaAngle);
    TYPE cosAlpha = cos(position.alphaAngle);
    TYPE sinBeta = sin(position.betaAngle);
    TYPE cosBeta = cos(position.betaAngle);
    TYPE sinPhi = sin(angularPosition);
    TYPE cosPhi = cos(angularPosition);

    TYPE ringPositionX = cosBeta * cosAlpha * cosPhi - sinBeta * sinPhi;
    TYPE ringPositionY = sinBeta * cosAlpha * cosPhi + cosBeta * sinPhi;
    TYPE ringPositionZ = (-1.0f) * sinAlpha * cosPhi;

    TYPE ringTangentX = (-1.0f) * cosBeta * cosAlpha * sinPhi - sinBeta * cosPhi;
    TYPE ringTangentY = (-1.0f) * sinBeta * cosAlpha * sinPhi + cosBeta * cosPhi;
    TYPE ringTangentZ = sinAlpha * sinPhi;

    TYPE posX = position.positionVector[0] - primCoil.positionVector[0] +
                lengthPosition * sinAlpha * cosBeta + thicknessPosition * ringPositionX;
    TYPE posY = position.positionVector[1] - primCoil.positionVector[1] +
                lengthPosition * sinAlpha * sinBeta + thicknessPosition * ringPositionY;
    TYPE posZ = position.positionVector[2] - primCoil.positionVector[2] +
                lengthPosition * cosAlpha + thicknessPosition * ringPositionZ;

    TYPE x = posX * primCoil.invTransformArray[0] + posY * primCoil.invTransformArray[1] + posZ * primCoil.invTransformArray[2];
    TYPE y = posX * primCoil.invTransformArray[3] + posY * primCoil.invTransformArray[4] + posZ * primCoil.invTransformArray[5];
    TYPE z = posX * primCoil.invTransformArray[6] + posY * primCoil.invTransformArray[7] + posZ * primCoil.invTransformArray[8];

    TYPE zCoord = z;
    TYPE rCoord = sqrt(x * x + y * y);
    TYPE phiCord = atan2(y, x);

    TYPE potential = 0.0f;

    TYPE topEdge = zCoord + 0.5f * primCoil.length;
    TYPE bottomEdge = zCoord - 0.5f * primCoil.length;

    if (primCoil.useFastMethod)
    {
        for (int incT = 0; incT < primCoil.thicknessIncrements; ++incT)
        {
            TYPE incrementPositionT = primCoil.innerRadius +
                                      0.5f * primCoil.thickness * (1.0f + primCoil.thicknessPositionArray[incT]);

            TYPE tempConstA = incrementPositionT * incrementPositionT + rCoord * rCoord;
            TYPE tempConstB = 2.0f * incrementPositionT * rCoord;

            for (int incF = 0; incF < primCoil.angularIncrements; ++incF)
            {
                TYPE cosinePhi = primCoil.cosPrecomputeArray[incF];

                TYPE tempConstC = rsqrt(tempConstA - tempConstB * cosinePhi);

                TYPE tempConstD1 = topEdge * tempConstC;
                TYPE tempConstD2 = bottomEdge * tempConstC;

                TYPE tempConstE1 = sqrt(tempConstD1 * tempConstD1 + 1.0f);
                TYPE tempConstE2 = sqrt(tempConstD2 * tempConstD2 + 1.0f);

                TYPE tempConstF = log((tempConstE1 + tempConstD1) / (tempConstE2 + tempConstD2));

                potential += primCoil.constFactor *
                             primCoil.thicknessWeightArray[incT] * primCoil.angularWeightArray[incF] *
                             incrementPositionT * cosinePhi * tempConstF;
            }
        }
    }
    else
    {
        for (int incT = 0; incT < primCoil.thicknessIncrements; ++incT)
        {
            TYPE incrementPositionT = primCoil.innerRadius +
                                      0.5f * primCoil.thickness * (1.0f + primCoil.thicknessPositionArray[incT]);

            TYPE tempConstA = incrementPositionT * incrementPositionT + rCoord * rCoord + zCoord * zCoord;
            TYPE tempConstB = 2.0f * incrementPositionT * rCoord;

            for (int incF = 0; incF < primCoil.angularIncrements; ++incF)
            {
                TYPE cosinePhi = primCoil.cosPrecomputeArray[incF];

                TYPE tempConstC = rsqrt(tempConstA - tempConstB * cosinePhi);

                potential += primCoil.constFactor *
                             primCoil.thicknessWeightArray[incT] * primCoil.angularWeightArray[incF] *
                             incrementPositionT * cosinePhi * tempConstC;
            }
        }
    }

    TYPE xPot = (-1.0f) * sin(phiCord) * potential;
    TYPE yPot = potential * cos(phiCord);
    TYPE zPot = 0.0f;

    TYPE potentialX = xPot * primCoil.transformArray[0] + yPot * primCoil.transformArray[1] + zPot * primCoil.transformArray[2];
    TYPE potentialY = xPot * primCoil.transformArray[3] + yPot * primCoil.transformArray[4] + zPot * primCoil.transformArray[5];
    TYPE potentialZ = xPot * primCoil.transformArray[6] + yPot * primCoil.transformArray[7] + zPot * primCoil.transformArray[8];

    TYPE weight = 0.125f * thicknessPosition / primCoil.current *
                  secondaryCoil.lengthWeightArray[lengthIndex] *
                  secondaryCoil.thicknessWeightArray[thicknessIndex] *
                  secondaryCoil.angularWeightArray[angularIndex];

    TYPE mInductance = weight * secondaryCoil.correctionFactor *
                       (potentialX * ringTangentX + potentialY * ringTangentY + potentialZ * ringTangentZ);

    atomicAdd(&inductanceArr[configIndex], mInductance);
}

namespace
{
    CoilData *g_coilArr = nullptr;
    SecondaryCoilPositionData *g_secondaryPositionArr = nullptr;
    TYPE *g_inductanceArr = nullptr;

    void getBuffers(long long coilCount, long long configs)
    {
        std::vector<void*> buffers = GPUMem::getBuffers(
                { coilCount * (long long)sizeof(CoilData),
                  configs * (long long)sizeof(SecondaryCoilPositionData),
                  configs * (long long)sizeof(TYPE)}
        );

        g_coilArr = static_cast<CoilData *>(buffers[0]);
        g_secondaryPositionArr = static_cast<SecondaryCoilPositionData *>(buffers[1]);
        g_inductanceArr = static_cast<TYPE*>(buffers[2]);
    }

    #if DEBUG_TIMINGS
        double g_duration;
    #endif
}


void Calculate_mutual_inductance_configurations_group(long long coilCount, long long configCount, long long pointCount,
                                                      SecondaryCoilData secondaryCoil,
                                                      const CoilData *coils,
                                                      const SecondaryCoilPositionData *secondaryPositions,
                                                      TYPE *inductanceArr)
{
    #if DEBUG_TIMINGS
        recordStartPoint();
        recordStartPoint();
    #endif

    int blocks = int(std::ceil(double((long long) configCount * pointCount) / NTHREADS));

    getBuffers(coilCount, configCount);

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tResource startup:         %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    gpuErrchk(cudaMemcpy(g_coilArr, coils, coilCount * sizeof(CoilData), cudaMemcpyHostToDevice))
    gpuErrchk(cudaMemcpy(g_secondaryPositionArr, secondaryPositions, configCount * sizeof(SecondaryCoilPositionData), cudaMemcpyHostToDevice))

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tMemory initialization:    %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    gpuErrchk(cudaMemset(g_inductanceArr, 0, configCount * sizeof(TYPE)))

    for (int i = 0; i < coilCount; ++i)
        CalculateMutualInductanceConfigurationsGroup<<<blocks, NTHREADS>>> (
            i, configCount, pointCount, secondaryCoil, g_coilArr,
            g_secondaryPositionArr, g_inductanceArr
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
        printf("\tDevice buffer size:       %.3lf MB\n",
               (double(coilCount * sizeof(CoilData) + configCount * (sizeof(SecondaryCoilPositionData) + sizeof(TYPE))) / 1e6)
        );
        printf("\tTotal blocks:             %d * %lli\n", blocks, coilCount);
        printf("\tThreads per calculation:  %i\n", NTHREADS);
        printf("\tTotal issued threads:     %i * %lli\n", NTHREADS * blocks, coilCount);
        printf("\tTotal coils:              %lli\n", coilCount);
        printf("\tTotal configurations:     %lli\n", configCount);
        printf("\tPoints per configuration: %lli\n", pointCount);
        printf("\tTotal calculations:       %lli\n", coilCount * pointCount * configCount);
        printf("\n\tPerformance:              %.0f Configs/s\n", double(configCount / g_duration));
        printf("\tEffective Performance:    %.1f kConfigs/s\n", double(1e-3 * coilCount * configCount / g_duration));
        printf("\tEquivalent Performance:   %.1f MPoints/s\n", double(1e-6 * coilCount * pointCount * configCount / g_duration));
        printf("---------------------------------------------------\n\n");
    #endif
}
