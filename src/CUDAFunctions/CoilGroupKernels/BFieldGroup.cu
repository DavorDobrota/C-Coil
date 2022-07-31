#include "CoilGroupAcceleration.h"

#include "Timing.h"
#include "CUDAFunctions/ErrorCheck/CUDAErrorCheck.h"
#include "CUDAFunctions/MemoryManagement/GPUMemoryManagement.h"

#include <cstdio>


__global__
void calculateFieldGroup(long long numOps, long long coilIndex,
                         const CoilData *coilArr,
                         const DataVector *posArr,
                         DataVector *resArr)
{
    unsigned int index = threadIdx.x;
    long long global_index = blockIdx.x * blockDim.x + index;

    if(global_index >= numOps)
        return;

    __shared__ CoilData coil;
    coil = coilArr[coilIndex];

    TYPE x1 = posArr[global_index].x;
    TYPE y1 = posArr[global_index].y;
    TYPE z1 = posArr[global_index].z;

    x1 -= coil.positionVector[0];
    y1 -= coil.positionVector[1];
    z1 -= coil.positionVector[2];

    TYPE x = x1 * coil.invTransformArray[0] + y1 * coil.invTransformArray[1] + z1 * coil.invTransformArray[2];
    TYPE y = x1 * coil.invTransformArray[3] + y1 * coil.invTransformArray[4] + z1 * coil.invTransformArray[5];
    TYPE z = x1 * coil.invTransformArray[6] + y1 * coil.invTransformArray[7] + z1 * coil.invTransformArray[8];

    TYPE zCoord = z;
    TYPE rCoord = sqrt(x * x + y * y);
    TYPE phiCord = atan2(y, x);

    TYPE fieldH = 0.0f;
    TYPE fieldZ = 0.0f;

    TYPE topEdge = zCoord + 0.5f * coil.length;
    TYPE bottomEdge = zCoord - 0.5f * coil.length;

    if (coil.useFastMethod)
    {
        for (int incT = 0; incT < coil.thicknessIncrements; ++incT)
        {
            TYPE incrementPositionT = coil.innerRadius + 0.5f * coil.thickness * (1.0f + coil.thicknessPositionArray[incT]);

            TYPE tempConstA = incrementPositionT * incrementPositionT;
            TYPE tempConstB = 2.0f * incrementPositionT * rCoord;
            TYPE tempConstC = tempConstA + rCoord * rCoord;

            TYPE tempConstD1 = topEdge * topEdge + tempConstC;
            TYPE tempConstD2 = bottomEdge * bottomEdge + tempConstC;

            for (int incF = 0; incF < coil.angularIncrements; ++incF)
            {
                TYPE cosinePhi = coil.cosPrecomputeArray[incF];

                TYPE tempConstE = tempConstB * cosinePhi;

                TYPE tempConstF1 = rsqrt(tempConstD1 - tempConstE);
                TYPE tempConstF2 = rsqrt(tempConstD2 - tempConstE);

                TYPE tempConstG = coil.constFactor * coil.thicknessWeightArray[incT] * coil.angularWeightArray[incF];

                fieldH += tempConstG * incrementPositionT * cosinePhi * (tempConstF2 - tempConstF1);
                fieldZ += tempConstG *
                          ((tempConstA - 0.5f * tempConstE) / (tempConstC - tempConstE)) *
                          (topEdge * tempConstF1 - bottomEdge * tempConstF2);
            }
        }
    }
    else
    {
        for (int incT = 0; incT < coil.thicknessIncrements; ++incT)
        {
            TYPE incrementPositionT = coil.innerRadius + 0.5f * coil.thickness * (1.0f + coil.thicknessPositionArray[incT]);

            TYPE tempConstA = incrementPositionT * incrementPositionT;
            TYPE tempConstB = incrementPositionT * rCoord;
            TYPE tempConstC = tempConstA + rCoord * rCoord + zCoord * zCoord;
            TYPE tempConstD = incrementPositionT * zCoord;

            for (int incF = 0; incF < coil.angularIncrements; ++incF)
            {
                TYPE cosinePhi = coil.cosPrecomputeArray[incF];

                TYPE tempConstE = tempConstC - 2.0f * tempConstB * cosinePhi;
                TYPE tempConstF = tempConstE * sqrt(tempConstE);
                TYPE tempConstG = coil.constFactor * coil.thicknessWeightArray[incT] * coil.angularWeightArray[incF] / tempConstF;

                fieldH += tempConstG * (tempConstD * cosinePhi);
                fieldZ += tempConstG * (tempConstA - tempConstB * cosinePhi);
            }
        }
    }

    TYPE xField = fieldH * cos(phiCord);
    TYPE yField = fieldH * sin(phiCord);
    TYPE zField = fieldZ;

    TYPE xRes = xField * coil.transformArray[0] + yField * coil.transformArray[1] + zField * coil.transformArray[2];
    TYPE yRes = xField * coil.transformArray[3] + yField * coil.transformArray[4] + zField * coil.transformArray[5];
    TYPE zRes = xField * coil.transformArray[6] + yField * coil.transformArray[7] + zField * coil.transformArray[8];

    resArr[global_index].x += xRes;
    resArr[global_index].y += yRes;
    resArr[global_index].z += zRes;
}


namespace
{
    CoilData *g_coilArr = nullptr;
    DataVector *g_posArr = nullptr;
    DataVector *g_resArr = nullptr;

    void getBuffers(long long numCoils, long long numOps)
    {
        std::vector<void*> buffers = GPUMem::getBuffers(
                { numCoils * (long long)sizeof(CoilData),
                  numOps * (long long)sizeof(DataVector),
                  numOps * (long long)sizeof(DataVector)}
        );

        g_coilArr = static_cast<CoilData*>(buffers[0]);
        g_posArr = static_cast<DataVector*>(buffers[1]);
        g_resArr = static_cast<DataVector*>(buffers[2]);
    }

    #if DEBUG_TIMINGS
        double g_duration;
    #endif
}


void Calculate_hardware_accelerated_b_group(long long numCoils, long long numOps,
                                            const CoilData *coilArr,
                                            const DataVector *posArr,
                                            DataVector *resArr)
{
    #if DEBUG_TIMINGS
        recordStartPoint();
        recordStartPoint();
    #endif

    long long blocks = ceil(double(numOps) / NTHREADS);

    getBuffers(numCoils, numOps);

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tResource startup:         %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    gpuErrchk(cudaMemcpy(g_coilArr, coilArr, numCoils * sizeof(CoilData), cudaMemcpyHostToDevice))
    gpuErrchk(cudaMemcpy(g_posArr, posArr, numOps * sizeof(DataVector), cudaMemcpyHostToDevice))

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tMemory initialization:    %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    gpuErrchk(cudaMemset(g_resArr, 0, numOps * sizeof(DataVector)))

    for (int i = 0; i < numCoils; ++i)
    {
        calculateFieldGroup<<<blocks, NTHREADS>>>(numOps, i, g_coilArr, g_posArr, g_resArr);
        gpuErrchk(cudaDeviceSynchronize())
    }

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tCalculations:             %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    if(resArr != nullptr)
        gpuErrchk(cudaMemcpy(resArr, g_resArr, numOps * sizeof(DataVector), cudaMemcpyDeviceToHost))

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tWriting to output array:  %.9g s\n\n", g_duration);

        g_duration = getIntervalDuration();
        printf("\tDevice buffer size:       %.3lf MB\n", (6.0 * double(numOps * sizeof(TYPE) + numCoils * sizeof(CoilData)) / 1.0e6));
        printf("\tTotal blocks:             %lli\n", blocks);
        printf("\tThreads per calculation:  %i\n", NTHREADS);
        printf("\tTotal coils:              %lli\n", numCoils);
        printf("\tTotal points:             %lli\n", numOps);
        printf("\tTotal calculations:       %lli\n", numOps * numCoils);
        printf("\n\tPerformance:              %.1f kPoints/s\n", double(0.001 * numOps / g_duration));
        printf("\n\tEffectivePerformance:     %.1f kPoints/s\n", double(0.001 * numOps * numCoils / g_duration));
        printf("---------------------------------------------------\n\n");
    #endif
}