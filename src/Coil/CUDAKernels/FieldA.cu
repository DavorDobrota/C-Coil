#include "CoilAcceleration.h"

#include "Timing.h"
#include "CUDAUtils/ErrorCheck/CUDAErrorCheck.h"
#include "CUDAUtils/MemoryManagement/GPUMemoryManagement.h"

#include <cstdio>


__global__
void calculatePotentialSlow(long long opCount, CoilData coil, const VectorData *posArr, VectorData *resArr)
{
    unsigned int index = threadIdx.x;
    long long global_index = (long long) blockIdx.x * blockDim.x + index;

    if(global_index >= opCount)
        return;

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

    TYPE potential = 0.0f;

    for (int incT = 0; incT < coil.thicknessIncrements; ++incT)
    {
        TYPE incrementPositionT = coil.innerRadius + 0.5f * coil.thickness * (1.0f + coil.thicknessPositionArray[incT]);

        TYPE tempConstA = incrementPositionT * incrementPositionT + rCoord * rCoord + zCoord * zCoord;
        TYPE tempConstB = 2.0f * incrementPositionT * rCoord;

        for (int incF = 0; incF < coil.angularIncrements; ++incF)
        {
            TYPE cosinePhi = coil.cosPrecomputeArray[incF];

            TYPE tempConstC = rsqrt(tempConstA - tempConstB * cosinePhi);

            potential += coil.constFactor *
                         coil.thicknessWeightArray[incT] * coil.angularWeightArray[incF] *
                         incrementPositionT * cosinePhi * tempConstC;
        }
    }

    TYPE xPot = (-1.0f) * sin(phiCord) * potential;
    TYPE yPot = potential * cos(phiCord);
    TYPE zPot = 0.0f;

    TYPE xRes = xPot * coil.transformArray[0] + yPot * coil.transformArray[1] + zPot * coil.transformArray[2];
    TYPE yRes = xPot * coil.transformArray[3] + yPot * coil.transformArray[4] + zPot * coil.transformArray[5];
    TYPE zRes = xPot * coil.transformArray[6] + yPot * coil.transformArray[7] + zPot * coil.transformArray[8];

    resArr[global_index].x += xRes;
    resArr[global_index].y += yRes;
    resArr[global_index].z += zRes;
}

__global__
void calculatePotentialFast(long long opCount, CoilData coil, const VectorData *posArr, VectorData *resArr)
{
    unsigned int index = threadIdx.x;
    long long global_index = blockIdx.x * blockDim.x + index;

    if(global_index >= opCount)
        return;

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

    TYPE potential = 0.0f;

    TYPE topEdge = zCoord + 0.5f * coil.length;
    TYPE bottomEdge = zCoord - 0.5f * coil.length;

    for (int incT = 0; incT < coil.thicknessIncrements; ++incT)
    {
        TYPE incrementPositionT = coil.innerRadius + 0.5f * coil.thickness * (1.0f + coil.thicknessPositionArray[incT]);

        TYPE tempConstA = incrementPositionT * incrementPositionT + rCoord * rCoord;
        TYPE tempConstB = (-2.0f) * incrementPositionT * rCoord;

        for (int incF = 0; incF < coil.angularIncrements; ++incF)
        {
            TYPE cosinePhi = coil.cosPrecomputeArray[incF];

            TYPE tempConstC = rsqrt(fma(tempConstB, cosinePhi, tempConstA));

            TYPE tempConstD1 = topEdge * tempConstC;
            TYPE tempConstD2 = bottomEdge * tempConstC;

            TYPE tempConstE1 = sqrt(tempConstD1 * tempConstD1 + 1.0f);
            TYPE tempConstE2 = sqrt(tempConstD2 * tempConstD2 + 1.0f);

            TYPE tempConstF = log((tempConstE1 + tempConstD1) / (tempConstE2 + tempConstD2));

            potential += coil.constFactor *
                    coil.thicknessWeightArray[incT] * coil.angularWeightArray[incF] *
                    incrementPositionT * cosinePhi * tempConstF;
        }
    }

    TYPE xPot = (-1.0f) * sin(phiCord) * potential;
    TYPE yPot = potential * cos(phiCord);
    TYPE zPot = 0.0f;

    TYPE xRes = xPot * coil.transformArray[0] + yPot * coil.transformArray[1] + zPot * coil.transformArray[2];
    TYPE yRes = xPot * coil.transformArray[3] + yPot * coil.transformArray[4] + zPot * coil.transformArray[5];
    TYPE zRes = xPot * coil.transformArray[6] + yPot * coil.transformArray[7] + zPot * coil.transformArray[8];

    resArr[global_index].x += xRes;
    resArr[global_index].y += yRes;
    resArr[global_index].z += zRes;
}

	
namespace 
{
    VectorData *g_posArr = nullptr;
    VectorData *g_resArr = nullptr;

    void getBuffers(long long opCount)
    {
        std::vector<VectorData*> buffers = GPUMem::getBuffers<VectorData>({opCount, opCount});

        g_posArr = buffers[0];
        g_resArr = buffers[1];
    }
    
    #if DEBUG_TIMINGS
        double g_duration;
    #endif
}


void Calculate_hardware_accelerated_a (long long opCount, CoilData coil,
                                       const VectorData *posArr,
                                       VectorData *resArr)
{
    #if DEBUG_TIMINGS
        recordStartPoint();
        recordStartPoint();
        printf("---------------------------------------------------\n");
    #endif

    int blocks = int(std::ceil(double(opCount) / NTHREADS));

    getBuffers(opCount);

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tResource startup:         %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    gpuErrchk(cudaMemcpy(g_posArr, posArr, opCount * sizeof(VectorData), cudaMemcpyHostToDevice));

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tMemory initialization:    %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    gpuErrchk(cudaMemset(g_resArr, 0, opCount * sizeof(VectorData)));

    if (coil.useFastMethod)
        calculatePotentialFast<<<blocks, NTHREADS>>>(opCount, coil, g_posArr, g_resArr);
    else
        calculatePotentialSlow<<<blocks, NTHREADS>>>(opCount, coil, g_posArr, g_resArr);

	gpuErrchk(cudaDeviceSynchronize());

	#if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tCalculations:             %.9g s\n", g_duration);

        recordStartPoint();
    #endif

	if(resArr != nullptr)
        gpuErrchk(cudaMemcpy(resArr, g_resArr, opCount * sizeof(VectorData), cudaMemcpyDeviceToHost));


    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tWriting to output array:  %.9g s\n\n", g_duration);

        g_duration = getIntervalDuration();
        printf("\tDevice buffer size:       %.3lf MB\n", (6.0 * double(opCount * sizeof(TYPE)) / 1.0e6));
        printf("\tTotal blocks:             %lli\n", blocks);
        printf("\tThreads per calculation:  %i\n", NTHREADS);
        printf("\tPrecision:                %dx%d\n", coil.thicknessIncrements, coil.angularIncrements);
        printf("\tTotal calculations:       %lli\n", opCount);
        printf("\tTotal MegaIncrements:     %.f\n", 1e-6 * double(opCount * coil.thicknessIncrements * coil.angularIncrements));
        printf("\n\tPerformance:              %.1f kPoints/s\n", double(0.001 * opCount / g_duration));
        printf("---------------------------------------------------\n\n");
    #endif
}
