#include "CoilAcceleration.h"

#include "Timing.h"
#include "CUDAUtils/ErrorCheck/CUDAErrorCheck.h"
#include "CUDAUtils/MemoryManagement/GPUMemoryManagement.h"

#include <cstdio>


__global__
void CalculateForceAndTorqueConfigurations(long long configCount, long long pointCount, CoilPairArgumentsData coilPair,
                                           const CoilPairPositionData *configArr,
                                           ForceTorqueData *forceTorqueArr)
{
    unsigned int index = threadIdx.x;
    long long global_index = blockIdx.x * blockDim.x + index;

    if(global_index >= configCount * pointCount)
        return;

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

    TYPE ringX = lengthPosition * sinAlpha * cosBeta + thicknessPosition * ringPositionX;
    TYPE ringY = lengthPosition * sinAlpha * sinBeta + thicknessPosition * ringPositionY;
    TYPE ringZ = lengthPosition * cosAlpha + thicknessPosition * ringPositionZ;

    TYPE posX = position.secPositionVector[0] - position.primPositionVector[0] + ringX;
    TYPE posY = position.secPositionVector[1] - position.primPositionVector[1] + ringY;
    TYPE posZ = position.secPositionVector[2] - position.primPositionVector[2] + ringZ;

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

    TYPE fieldH = 0.0f;
    TYPE fieldZ = 0.0f;

    TYPE topEdge = zCoord + 0.5f * coilPair.primLength;
    TYPE bottomEdge = zCoord - 0.5f * coilPair.primLength;

    if (coilPair.useFastMethod)
    {
        for (int incT = 0; incT < coilPair.primThicknessIncrements; ++incT)
        {
            TYPE incrementPositionT = coilPair.primInnerRadius +
                                      0.5f * coilPair.primThickness * (1.0f + coilPair.primThicknessPositionArray[incT]);

            TYPE tempConstA = incrementPositionT * incrementPositionT;
            TYPE tempConstB = 2.0f * incrementPositionT * rCoord;
            TYPE tempConstC = tempConstA + rCoord * rCoord;

            TYPE tempConstD1 = topEdge * topEdge + tempConstC;
            TYPE tempConstD2 = bottomEdge * bottomEdge + tempConstC;

            for (int incF = 0; incF < coilPair.primAngularIncrements; ++incF)
            {
                TYPE cosinePhi = coilPair.primCosPrecomputeArray[incF];

                TYPE tempConstE = tempConstB * cosinePhi;

                TYPE tempConstF1 = rsqrt(tempConstD1 - tempConstE);
                TYPE tempConstF2 = rsqrt(tempConstD2 - tempConstE);

                TYPE tempConstG = coilPair.constFactor *
                                  coilPair.primThicknessWeightArray[incT] * coilPair.primAngularWeightArray[incF];

                fieldH += tempConstG * incrementPositionT * cosinePhi * (tempConstF2 - tempConstF1);
                fieldZ += tempConstG *
                          ((tempConstA - 0.5f * tempConstE) / (tempConstC - tempConstE)) *
                          (topEdge * tempConstF1 - bottomEdge * tempConstF2);
            }
        }
    }
    else
    {
        for (int incT = 0; incT < coilPair.primThicknessIncrements; ++incT)
        {
            TYPE incrementPositionT = coilPair.primInnerRadius +
                                      0.5f * coilPair.primThickness * (1.0f + coilPair.primThicknessPositionArray[incT]);

            TYPE tempConstA = incrementPositionT * incrementPositionT;
            TYPE tempConstB = incrementPositionT * rCoord;
            TYPE tempConstC = tempConstA + rCoord * rCoord + zCoord * zCoord;
            TYPE tempConstD = incrementPositionT * zCoord;

            for (int incF = 0; incF < coilPair.primAngularIncrements; ++incF)
            {
                TYPE cosinePhi = coilPair.primCosPrecomputeArray[incF];

                TYPE tempConstE = tempConstC - 2.0f * tempConstB * cosinePhi;
                TYPE tempConstF = tempConstE * sqrt(tempConstE);
                TYPE tempConstG = (coilPair.constFactor / tempConstF) *
                                  coilPair.primThicknessWeightArray[incT] * coilPair.primAngularWeightArray[incF];

                fieldH += tempConstG * (tempConstD * cosinePhi);
                fieldZ += tempConstG * (tempConstA - tempConstB * cosinePhi);
            }
        }
    }

    TYPE xComp = fieldH * cos(phiCord);
    TYPE yComp = fieldH * sin(phiCord);
    TYPE zComp = fieldZ;

    TYPE xField = xComp * (cosZ * cosZ * cosY - sinZ * sinZ) +
                  yComp * (-1.0f * sinZ * cosZ * cosY - sinZ * cosZ) +
                  zComp * (cosZ * sinY);
    TYPE yField = xComp * (sinZ * cosZ * cosY + sinZ * cosZ) +
                  yComp * (cosZ * cosZ - sinZ * sinZ * cosY) +
                  zComp * (sinZ * sinY);
    TYPE zField = xComp * (-1.0f * sinY * cosZ) +
                  yComp * (sinY * sinZ) +
                  zComp * cosY;

    TYPE weight = 0.125f * thicknessPosition *
                  coilPair.secLengthWeightArray[lengthIndex] *
                  coilPair.secThicknessWeightArray[thicknessIndex] *
                  coilPair.secAngularWeightArray[angularIndex];

    TYPE forceX = coilPair.correctionFactor * weight * (ringTangentY * zField - ringTangentZ * yField);
    TYPE forceY = coilPair.correctionFactor * weight * (ringTangentZ * xField - ringTangentX * zField);
    TYPE forceZ = coilPair.correctionFactor * weight * (ringTangentX * yField - ringTangentY * xField);

    TYPE torqueX = ringY * forceZ - ringZ * forceY;
    TYPE torqueY = ringZ * forceX - ringX * forceZ;
    TYPE torqueZ = ringX * forceY - ringY * forceX;

    atomicAdd(&forceTorqueArr[pairIndex].forceX, forceX);
    atomicAdd(&forceTorqueArr[pairIndex].forceY, forceY);
    atomicAdd(&forceTorqueArr[pairIndex].forceZ, forceZ);

    atomicAdd(&forceTorqueArr[pairIndex].torqueX, torqueX);
    atomicAdd(&forceTorqueArr[pairIndex].torqueY, torqueY);
    atomicAdd(&forceTorqueArr[pairIndex].torqueZ, torqueZ);
}

namespace
{
    CoilPairPositionData *g_configArr = nullptr;

    ForceTorqueData *g_forceTorqueArr = nullptr;

    void getBuffers(long long configs)
    {
        std::vector<void*> buffers = GPUMem::getBuffers(
                { configs * (long long)sizeof(CoilPairPositionData),
                  configs * (long long)sizeof(ForceTorqueData)}
        );

        g_configArr = static_cast<CoilPairPositionData*>(buffers[0]);
        g_forceTorqueArr = static_cast<ForceTorqueData*>(buffers[1]);

    }

#if DEBUG_TIMINGS
    double g_duration;
#endif
}

void Calculate_force_and_torque_configurations(long long configCount, long long pointCount,
                                               CoilPairArgumentsData coilPair,
                                               const CoilPairPositionData *configArr,
                                               ForceTorqueData *forceTorqueArr)
{
    #if DEBUG_TIMINGS
        recordStartPoint();
        recordStartPoint();
    #endif

    int blocks = ceil(double(configCount * pointCount) / NTHREADS);

    getBuffers(configCount);

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tResource startup:         %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    gpuErrchk(cudaMemcpy(g_configArr, configArr, configCount * sizeof(CoilPairPositionData), cudaMemcpyHostToDevice));

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tMemory initialization:    %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    gpuErrchk(cudaMemset(g_forceTorqueArr, 0, configCount * sizeof(ForceTorqueData)))

    CalculateForceAndTorqueConfigurations<<<blocks, NTHREADS>>>(configCount, pointCount, coilPair, g_configArr, g_forceTorqueArr);
    gpuErrchk(cudaDeviceSynchronize())

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tCalculations:             %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    if(g_forceTorqueArr != nullptr)
    gpuErrchk(cudaMemcpy(forceTorqueArr, g_forceTorqueArr, configCount * sizeof(ForceTorqueData), cudaMemcpyDeviceToHost))

#if DEBUG_TIMINGS
    g_duration = getIntervalDuration();
    printf("\tWriting to output array:  %.9g s\n\n", g_duration);

    g_duration = getIntervalDuration();
    printf("\tDevice buffer size:       %.3lf MB\n", (configCount * double(sizeof(CoilPairPositionData) + sizeof(TYPE)) / 1.0e6));
    printf("\tTotal blocks:             %d\n", blocks);
    printf("\tThreads per calculation:  %i\n", NTHREADS);
    printf("\tTotal issued threads:     %i\n", NTHREADS * blocks);
    printf("\tTotal configurations:     %d\n", configCount);
    printf("\tPoints per configuration: %d\n", pointCount);
    printf("\tTotal calculations:       %d\n", pointCount * configCount);
    printf("\n\tPerformance:              %.1f kCoils/s\n", double(0.001 * configCount / g_duration));
    printf("\n\tEffective Performance:    %.1f kPoints/s\n", double(0.001 * pointCount * configCount / g_duration));
    printf("---------------------------------------------------\n\n");
#endif
}