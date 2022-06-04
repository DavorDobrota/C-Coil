#include "hardware_acceleration.h"

#include "CUDAConstants.h"
#include "Timing.h"
#include "CUDAErrorCheck.h"
#include "CoilData.h"

#include <cstdio>


__global__
void calculateA(long long numOps, CoilData coil,
                const TYPE *xPosArr, const TYPE *yPosArr, const TYPE *zPosArr,
                TYPE *xResArr, TYPE *yResArr, TYPE *zResArr)
{
    unsigned int index = threadIdx.x;
    long long global_index = blockIdx.x * blockDim.x + index;

    if(global_index >= numOps)
        return;

    TYPE x1 = xPosArr[global_index];
    TYPE y1 = yPosArr[global_index];
    TYPE z1 = zPosArr[global_index];

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
    TYPE constant = coil.constFactor;

    TYPE topEdge = zCoord + 0.5f * coil.length;
    TYPE bottomEdge = zCoord - 0.5f * coil.length;

    for (int incT = 0; incT < coil.thicknessIncrements; ++incT)
    {
        TYPE incrementPositionT = coil.innerRadius + 0.5f * coil.thickness * (1.0f + coil.positionArray[incT]);

        for (int incF = 0; incF < coil.angularIncrements; ++incF)
        {
            TYPE cosinePhi = coil.cosPrecomputeArray[incF];

            TYPE tempConstC = rsqrt(incrementPositionT * incrementPositionT + rCoord * rCoord -
                                     2.0f * incrementPositionT * rCoord * cosinePhi);

            TYPE tempConstD1 = topEdge * tempConstC;
            TYPE tempConstD2 = bottomEdge * tempConstC;

            TYPE tempConstE1 = sqrt(tempConstD1 * tempConstD1 + 1.0f);
            TYPE tempConstE2 = sqrt(tempConstD2 * tempConstD2 + 1.0f);

            TYPE tempConstF = log((tempConstE1 + tempConstD1) / (tempConstE2 + tempConstD2));

            potential += constant *
                    coil.weightArray[incT] * coil.weightArray[incF] *
                    incrementPositionT * cosinePhi * tempConstF;
        }
    }
    TYPE xPot = (-1.f) * sin(phiCord) * potential;
    TYPE yPot = potential * cos(phiCord);
    TYPE zPot = 0.f;

    TYPE xRes = xPot * coil.transformArray[0] + yPot * coil.transformArray[1] + zPot * coil.transformArray[2];
    TYPE yRes = xPot * coil.transformArray[3] + yPot * coil.transformArray[4] + zPot * coil.transformArray[5];
    TYPE zRes = xPot * coil.transformArray[6] + yPot * coil.transformArray[7] + zPot * coil.transformArray[8];

    xResArr[global_index] = xRes;
    yResArr[global_index] = yRes;
    zResArr[global_index] = zRes;
}
	
namespace 
{
    long long g_last_num_ops = 0;

    TYPE *g_xPosArr = nullptr;
    TYPE *g_yPosArr = nullptr;
    TYPE *g_zPosArr = nullptr;

    TYPE *g_xResArr = nullptr;
    TYPE *g_yResArr = nullptr;
    TYPE *g_zResArr = nullptr;
    
    #if DEBUG_TIMINGS
        double g_duration;
    #endif
}

void resourceCleanupA()
{
	gpuErrchk(cudaFree(g_xPosArr));
    gpuErrchk(cudaFree(g_yPosArr));
    gpuErrchk(cudaFree(g_zPosArr));

    gpuErrchk(cudaFree(g_xResArr));
    gpuErrchk(cudaFree(g_yResArr));
    gpuErrchk(cudaFree(g_zResArr));

    g_xPosArr = nullptr;
    g_yPosArr = nullptr;
    g_zPosArr = nullptr;

    g_xResArr = nullptr;
    g_yResArr = nullptr;
    g_zResArr = nullptr;
}

void resourceStartupA(long long numOps)
{
    resourceCleanupA();
    
	gpuErrchk(cudaMalloc(&g_xPosArr, numOps * sizeof(TYPE)));
    gpuErrchk(cudaMalloc(&g_yPosArr, numOps * sizeof(TYPE)));
    gpuErrchk(cudaMalloc(&g_zPosArr, numOps * sizeof(TYPE)));

    gpuErrchk(cudaMalloc(&g_xResArr, numOps * sizeof(TYPE)));
    gpuErrchk(cudaMalloc(&g_yResArr, numOps * sizeof(TYPE)));
    gpuErrchk(cudaMalloc(&g_zResArr, numOps * sizeof(TYPE)));
}


void Calculate_hardware_accelerated_a (long long numOps, CoilData coil,
                                       const TYPE *xPosArr, const TYPE *yPosArr, const TYPE *zPosArr,
                                       TYPE *xResArr, TYPE *yResArr, TYPE *zResArr)
{
    #if DEBUG_TIMINGS
        recordStartPoint();
    #endif

    long long blocks = ceil(double(numOps) / NTHREADS);

    if (numOps > g_last_num_ops)
    {
        resourceStartupA(numOps);
        g_last_num_ops = numOps;
    }
    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tResource startup:         %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    gpuErrchk(cudaMemcpy(g_xPosArr, xPosArr, numOps * sizeof(TYPE), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(g_yPosArr, yPosArr, numOps * sizeof(TYPE), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(g_zPosArr, zPosArr, numOps * sizeof(TYPE), cudaMemcpyHostToDevice));

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tMemory initialization:    %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    calculateA<<<blocks, NTHREADS>>>(numOps, coil,
                                     g_xPosArr, g_yPosArr, g_zPosArr,
                                     g_xResArr, g_yResArr, g_zResArr);
	gpuErrchk(cudaDeviceSynchronize());

	#if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\t\tCalculations:     %.9g s\n", g_duration);
        printf("\t\tEstimated TFLOPS: %.2f\n",
               1e-12 * double(100 * numOps * coil.thicknessIncrements * coil.angularIncrements) / g_duration);

        recordStartPoint();
    #endif

	if(xResArr != nullptr)
    {
        gpuErrchk(cudaMemcpy(xResArr, g_xResArr, numOps * sizeof(TYPE), cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(yResArr, g_yResArr, numOps * sizeof(TYPE), cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(zResArr, g_zResArr, numOps * sizeof(TYPE), cudaMemcpyDeviceToHost));
    }

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tWriting to output array:  %.9g s\n\n", g_duration);
    #endif

	#if DEBUG_TIMINGS
        printf("\tDevice buffer size:       %.3lf MB\n", (6.0 * double(numOps * sizeof(TYPE)) / 1.0e6));
        printf("\tTotal blocks:             %lli\n", blocks);
        printf("\tThreads per calculation:  %i\n", NTHREADS);
        printf("\tPrecision:                %dx%d\n", coil.thicknessIncrements, coil.angularIncrements);
        printf("\tTotal calculations:       %lli\n", numOps);
        printf("\tTotal MegaIncrements:     %.f\n", 1e-6 * double(numOps * coil.thicknessIncrements * coil.angularIncrements));
    #endif
}
