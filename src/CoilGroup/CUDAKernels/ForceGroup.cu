#include "CoilGroupAcceleration.h"

#include "Timing.h"
#include "CUDAUtils/ErrorCheck/CUDAErrorCheck.h"
#include "CUDAUtils/MemoryManagement/GPUMemoryManagement.h"

#include <cstdio>


__global__
void CalculateForceAndTorqueConfigurationsGroup(long long coilIndex, long long configCount, long long pointCount,
                                                SecondaryCoilData secondaryCoil,
                                                const CoilData *primaryCoils,
                                                const SecondaryCoilPositionData *secondaryPositions,
                                                ForceTorqueData *forceTorqueArr)
{
    unsigned int index = threadIdx.x;
    long long global_index = blockIdx.x * blockDim.x + index;

    if(global_index >= configCount * pointCount)
        return;

    int configIndex = int(global_index % configCount);
    int lengthIndex = int((global_index / configCount) % secondaryCoil.lengthIncrements);
    int thicknessIndex = int(((global_index / configCount) / secondaryCoil.lengthIncrements) % secondaryCoil.thicknessIncrements);
    int angularIndex = int((((global_index / configCount) / secondaryCoil.lengthIncrements) / secondaryCoil.thicknessIncrements) % secondaryCoil.angularIncrements);

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

    TYPE ringX = lengthPosition * sinAlpha * cosBeta + thicknessPosition * ringPositionX;
    TYPE ringY = lengthPosition * sinAlpha * sinBeta + thicknessPosition * ringPositionY;
    TYPE ringZ = lengthPosition * cosAlpha + thicknessPosition * ringPositionZ;

    TYPE posX = position.positionVector[0] - primCoil.positionVector[0] + ringX;
    TYPE posY = position.positionVector[1] - primCoil.positionVector[1] + ringY;
    TYPE posZ = position.positionVector[2] - primCoil.positionVector[2] + ringZ;

    TYPE x = posX * primCoil.invTransformArray[0] + posY * primCoil.invTransformArray[1] + posZ * primCoil.invTransformArray[2];
    TYPE y = posX * primCoil.invTransformArray[3] + posY * primCoil.invTransformArray[4] + posZ * primCoil.invTransformArray[5];
    TYPE z = posX * primCoil.invTransformArray[6] + posY * primCoil.invTransformArray[7] + posZ * primCoil.invTransformArray[8];

    TYPE zCoord = z;
    TYPE rCoord = sqrt(x * x + y * y);
    TYPE phiCord = atan2(y, x);

    TYPE fieldH = 0.0f;
    TYPE fieldZ = 0.0f;

    TYPE topEdge = zCoord + 0.5f * primCoil.length;
    TYPE bottomEdge = zCoord - 0.5f * primCoil.length;

    if (primCoil.useFastMethod)
    {
        for (int incT = 0; incT < primCoil.thicknessIncrements; ++incT)
        {
            TYPE incrementPositionT = primCoil.innerRadius + 0.5f * primCoil.thickness * (1.0f + primCoil.thicknessPositionArray[incT]);

            TYPE tempConstA = incrementPositionT * incrementPositionT;
            TYPE tempConstB = 2.0f * incrementPositionT * rCoord;
            TYPE tempConstC = tempConstA + rCoord * rCoord;

            TYPE tempConstD1 = topEdge * topEdge + tempConstC;
            TYPE tempConstD2 = bottomEdge * bottomEdge + tempConstC;

            for (int incF = 0; incF < primCoil.angularIncrements; ++incF)
            {
                TYPE cosinePhi = primCoil.cosPrecomputeArray[incF];

                TYPE tempConstE = tempConstB * cosinePhi;

                TYPE tempConstF1 = rsqrt(tempConstD1 - tempConstE);
                TYPE tempConstF2 = rsqrt(tempConstD2 - tempConstE);

                TYPE tempConstG = primCoil.constFactor * primCoil.thicknessWeightArray[incT] * primCoil.angularWeightArray[incF];

                fieldH += tempConstG * incrementPositionT * cosinePhi * (tempConstF2 - tempConstF1);
                fieldZ += tempConstG *
                          ((tempConstA - 0.5f * tempConstE) / (tempConstC - tempConstE)) *
                          (topEdge * tempConstF1 - bottomEdge * tempConstF2);
            }
        }
    }
    else
    {
        for (int incT = 0; incT < primCoil.thicknessIncrements; ++incT)
        {
            TYPE incrementPositionT = primCoil.innerRadius + 0.5f * primCoil.thickness * (1.0f + primCoil.thicknessPositionArray[incT]);

            TYPE tempConstA = incrementPositionT * incrementPositionT;
            TYPE tempConstB = incrementPositionT * rCoord;
            TYPE tempConstC = tempConstA + rCoord * rCoord + zCoord * zCoord;
            TYPE tempConstD = incrementPositionT * zCoord;

            for (int incF = 0; incF < primCoil.angularIncrements; ++incF)
            {
                TYPE cosinePhi = primCoil.cosPrecomputeArray[incF];

                TYPE tempConstE = tempConstC - 2.0f * tempConstB * cosinePhi;
                TYPE tempConstF = tempConstE * sqrt(tempConstE);
                TYPE tempConstG = primCoil.constFactor * primCoil.thicknessWeightArray[incT] * primCoil.angularWeightArray[incF] / tempConstF;

                fieldH += tempConstG * (tempConstD * cosinePhi);
                fieldZ += tempConstG * (tempConstA - tempConstB * cosinePhi);
            }
        }
    }

    TYPE xComp = fieldH * cos(phiCord);
    TYPE yComp = fieldH * sin(phiCord);
    TYPE zComp = fieldZ;

    TYPE xField = xComp * primCoil.transformArray[0] + yComp * primCoil.transformArray[1] + zComp * primCoil.transformArray[2];
    TYPE yField = xComp * primCoil.transformArray[3] + yComp * primCoil.transformArray[4] + zComp * primCoil.transformArray[5];
    TYPE zField = xComp * primCoil.transformArray[6] + yComp * primCoil.transformArray[7] + zComp * primCoil.transformArray[8];

    TYPE weight = 0.125f * thicknessPosition *
                  secondaryCoil.lengthWeightArray[lengthIndex] *
                  secondaryCoil.thicknessWeightArray[thicknessIndex] *
                  secondaryCoil.angularWeightArray[angularIndex];

    TYPE forceX = secondaryCoil.correctionFactor * weight * (ringTangentY * zField - ringTangentZ * yField);
    TYPE forceY = secondaryCoil.correctionFactor * weight * (ringTangentZ * xField - ringTangentX * zField);
    TYPE forceZ = secondaryCoil.correctionFactor * weight * (ringTangentX * yField - ringTangentY * xField);

    TYPE torqueX = ringY * forceZ - ringZ * forceY;
    TYPE torqueY = ringZ * forceX - ringX * forceZ;
    TYPE torqueZ = ringX * forceY - ringY * forceX;

    atomicAdd(&forceTorqueArr[configIndex].forceX, forceX);
    atomicAdd(&forceTorqueArr[configIndex].forceY, forceY);
    atomicAdd(&forceTorqueArr[configIndex].forceZ, forceZ);

    atomicAdd(&forceTorqueArr[configIndex].torqueX, torqueX);
    atomicAdd(&forceTorqueArr[configIndex].torqueY, torqueY);
    atomicAdd(&forceTorqueArr[configIndex].torqueZ, torqueZ);
}


namespace
{
    CoilData *g_coilArr = nullptr;
    SecondaryCoilPositionData *g_secondaryPositionArr = nullptr;
    ForceTorqueData *g_forceTorqueArr = nullptr;

    void getBuffers(long long coilCount, long long configs)
    {
        std::vector<void*> buffers = GPUMem::getBuffers(
                { coilCount * (long long)sizeof(CoilData),
                  configs * (long long)sizeof(SecondaryCoilPositionData),
                  configs * (long long)sizeof(ForceTorqueData)}
        );

        g_coilArr = static_cast<CoilData *>(buffers[0]);
        g_secondaryPositionArr = static_cast<SecondaryCoilPositionData *>(buffers[1]);
        g_forceTorqueArr = static_cast<ForceTorqueData *>(buffers[2]);
    }

    #if DEBUG_TIMINGS
        double g_duration;
    #endif
}

void Calculate_force_and_torque_configurations_group(long long coilCount, long long configCount, long long pointCount,
                                                     SecondaryCoilData secondaryCoil,
                                                     const CoilData *coils,
                                                     const SecondaryCoilPositionData *secondaryPositions,
                                                     ForceTorqueData *forceTorqueArr)
{
    #if DEBUG_TIMINGS
        recordStartPoint();
        recordStartPoint();
    #endif

    int blocks = ceil(double(configCount * pointCount) / NTHREADS);

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

    gpuErrchk(cudaMemset(g_forceTorqueArr, 0, configCount * sizeof(ForceTorqueData)))

    for (int i = 0; i < coilCount; ++i)
        CalculateForceAndTorqueConfigurationsGroup<<<blocks, NTHREADS>>>(
            i, configCount, pointCount, secondaryCoil, g_coilArr,
            g_secondaryPositionArr, g_forceTorqueArr
        );
    gpuErrchk(cudaDeviceSynchronize())

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tCalculations:             %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    if(forceTorqueArr != nullptr)
        gpuErrchk(cudaMemcpy(forceTorqueArr, g_forceTorqueArr, configCount * sizeof(TYPE), cudaMemcpyDeviceToHost))

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tWriting to output array:  %.9g s\n\n", g_duration);

        g_duration = getIntervalDuration();
        printf("\tDevice buffer size:       %.3lf MB\n", (double(coilCount * sizeof(CoilData) + configCount * (sizeof(SecondaryCoilPositionData) + sizeof(ForceTorqueData)) / 1.0e6)));
        printf("\tTotal blocks:             %d\n", blocks);
        printf("\tThreads per calculation:  %i\n", NTHREADS);
        printf("\tTotal issued threads:     %i\n", NTHREADS * blocks);
        printf("\tTotal coils:              %d\n", coilCount);
        printf("\tTotal configurations:     %d\n", configCount);
        printf("\tPoints per configuration: %d\n", pointCount);
        printf("\tTotal calculations:       %d\n", coilCount * pointCount * configCount);
        printf("\n\tPerformance:              %.0f Configs/s\n", double(configCount / g_duration));
        printf("\n\tEffective Performance:    %.1f kPoints/s\n", double(0.001 * coilCount * pointCount * configCount / g_duration));
        printf("---------------------------------------------------\n\n");
    #endif
}