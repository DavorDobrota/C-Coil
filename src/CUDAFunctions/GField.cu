#include "hardware_acceleration.h"

#include "CUDAConstants.h"
#include "Timing.h"
#include "CUDAErrorCheck.h"
#include "CoilData.h"

#include <cstdio>


__global__
void calculateG(long long numOps, CoilData coil, const DataVector *posArr, DataMatrix *resArr)
{
    unsigned int index = threadIdx.x;
    long long global_index = blockIdx.x * blockDim.x + index;

    if(global_index >= numOps)
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

    TYPE bufferValueRP = 0.f;
    TYPE bufferValueRR = 0.f;
    TYPE bufferValueRZ = 0.f;
    TYPE bufferValueZZ = 0.f;

    TYPE topEdge = zCoord + 0.5f * coil.length;
    TYPE bottomEdge = zCoord - 0.5f * coil.length;

    for (int incT = 0; incT < coil.thicknessIncrements; ++incT)
    {
        TYPE incrementPositionT = coil.innerRadius + 0.5f * coil.thickness * (1.0f + coil.positionArray[incT]);

        TYPE tempConstA = incrementPositionT * incrementPositionT;
        TYPE tempConstB = rCoord * rCoord;
        TYPE tempConstC = incrementPositionT * rCoord;

        TYPE tempConstD = tempConstA + tempConstB;
        TYPE tempConstE = tempConstA * tempConstA + tempConstB * tempConstB;
        TYPE tempConstF = tempConstC * tempConstC;

        TYPE tempConstG1 = tempConstD + topEdge * topEdge;
        TYPE tempConstG2 = tempConstD + bottomEdge * bottomEdge;

        for (int incF = 0; incF < coil.angularIncrements; ++incF)
        {
            TYPE cosinePhi = coil.cosPrecomputeArray[incF];
            TYPE cosinePhi2 = cosinePhi * cosinePhi;
            TYPE phiExpression = 2.f * tempConstC * cosinePhi;

            TYPE tempConstI = coil.constFactor * coil.weightArray[incT] * coil.weightArray[incF];

            TYPE tempConstJ1 = tempConstG1 - phiExpression;
            TYPE tempConstJ2 = tempConstG2 - phiExpression;

            TYPE tempConstK1 = rsqrtf(tempConstJ1);
            TYPE tempConstK2 = rsqrtf(tempConstJ2);

            TYPE tempConstL1 = tempConstK1 / (tempConstJ1);
            TYPE tempConstL2 = tempConstK2 / (tempConstJ2);

            TYPE tempConstM = tempConstD - phiExpression;
            TYPE tempConstN =
                    2.f * tempConstF * cosinePhi * (cosinePhi2 + 2.f) -
                    tempConstC * (3.f * cosinePhi2 + 1.f) * tempConstD + cosinePhi * tempConstE;
            TYPE tempConstO = cosinePhi * tempConstD - 2.f * tempConstC;

            bufferValueRP += tempConstI * (incrementPositionT * cosinePhi / rCoord) * (tempConstK2 - tempConstK1);
            bufferValueRR += tempConstI * (tempConstC - tempConstA * cosinePhi) * cosinePhi * (tempConstL1 - tempConstL2);
            bufferValueZZ += tempConstI * (tempConstA - tempConstC * cosinePhi) * (tempConstL1 - tempConstL2);
            bufferValueRZ += tempConstI * incrementPositionT / (tempConstM * tempConstM) *
                             (topEdge * tempConstL1 * (tempConstO * tempConstJ1 + tempConstN) -
                              bottomEdge * tempConstL2 * (tempConstO * tempConstJ2 + tempConstN));
        }
    }

    TYPE xxGrad, xyGrad, xzGrad, yxGrad, yyGrad, yzGrad, zxGrad, zyGrad, zzGrad;

    if (rCoord / coil.innerRadius > 1e-5)
    {
        TYPE sinPhi = sin(phiCord);
        TYPE cosPhi = cos(phiCord);

        xxGrad = bufferValueRR * cosPhi * cosPhi + bufferValueRP * sinPhi * sinPhi;
        yyGrad = bufferValueRR * sinPhi * sinPhi + bufferValueRP * cosPhi * cosPhi;
        zzGrad = bufferValueZZ;

        xyGrad = 0.5 * sin(2 * phiCord) * (bufferValueRR - bufferValueRP);
        xzGrad = bufferValueRZ * cosPhi;
        yzGrad = bufferValueRZ * sinPhi;

        yxGrad = xyGrad;
        zxGrad = xzGrad;
        zyGrad = yzGrad;
    }
    else
    {
        xxGrad = bufferValueRR;
        yyGrad = bufferValueRR;
        zzGrad = bufferValueZZ;

        xyGrad = 0.f;
        xzGrad = 0.f;
        yxGrad = 0.f;
        yzGrad = 0.f;
        zxGrad = 0.f;
        zyGrad = 0.f;
    }

    TYPE xxRes = coil.transformArray[0] * xxGrad + coil.transformArray[1] * yxGrad + coil.transformArray[2] * zxGrad;
    TYPE xyRes = coil.transformArray[0] * xyGrad + coil.transformArray[1] * yyGrad + coil.transformArray[2] * zyGrad;
    TYPE xzRes = coil.transformArray[0] * xzGrad + coil.transformArray[1] * yzGrad + coil.transformArray[2] * zzGrad;
    TYPE yxRes = coil.transformArray[3] * xxGrad + coil.transformArray[4] * yxGrad + coil.transformArray[5] * zxGrad;
    TYPE yyRes = coil.transformArray[3] * xyGrad + coil.transformArray[4] * yyGrad + coil.transformArray[5] * zyGrad;
    TYPE yzRes = coil.transformArray[3] * xzGrad + coil.transformArray[4] * yzGrad + coil.transformArray[5] * zzGrad;
    TYPE zxRes = coil.transformArray[6] * xxGrad + coil.transformArray[7] * yxGrad + coil.transformArray[8] * zxGrad;
    TYPE zyRes = coil.transformArray[6] * xyGrad + coil.transformArray[7] * yyGrad + coil.transformArray[8] * zyGrad;
    TYPE zzRes = coil.transformArray[6] * xzGrad + coil.transformArray[7] * yzGrad + coil.transformArray[8] * zzGrad;

    resArr[global_index].xx = xxRes;
    resArr[global_index].xy = xyRes;
    resArr[global_index].xz = xzRes;
    resArr[global_index].yx = yxRes;
    resArr[global_index].yy = yyRes;
    resArr[global_index].yz = yzRes;
    resArr[global_index].zx = zxRes;
    resArr[global_index].zy = zyRes;
    resArr[global_index].zz = zzRes;
}
	
namespace 
{
    long long g_last_num_ops = 0;

    DataVector *g_posArr = nullptr;
    DataMatrix *g_resArr = nullptr;
    
    #if DEBUG_TIMINGS
        double g_duration;
    #endif
}

void resourceCleanupG()
{
    gpuErrchk(cudaFree(g_posArr));
    gpuErrchk(cudaFree(g_resArr));

    g_posArr = nullptr;
    g_resArr = nullptr;
}

void resourceStartupG(long long numOps)
{
    resourceCleanupG();

    gpuErrchk(cudaMalloc(&g_posArr, numOps * sizeof(DataVector)));
    gpuErrchk(cudaMalloc(&g_resArr, numOps * sizeof(DataMatrix)));
}


void Calculate_hardware_accelerated_g(long long numOps, CoilData coil, const DataVector *posArr, DataMatrix *resArr)
{
    #if DEBUG_TIMINGS
        recordStartPoint();
        recordStartPoint();
    #endif

    long long blocks = ceil(double(numOps) / NTHREADS);

    if (numOps > g_last_num_ops)
    {
        resourceStartupG(numOps);
        g_last_num_ops = numOps;
    }
    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tResource startup:         %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    gpuErrchk(cudaMemcpy(g_posArr, posArr, numOps * sizeof(DataVector), cudaMemcpyHostToDevice));

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tMemory initialization:    %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    calculateG<<<blocks, NTHREADS>>>(numOps, coil, g_posArr, g_resArr);
	gpuErrchk(cudaDeviceSynchronize());

	#if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tCalculations:             %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    if(resArr != nullptr)
        gpuErrchk(cudaMemcpy(resArr, g_resArr, numOps * sizeof(DataMatrix), cudaMemcpyDeviceToHost));

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tWriting to output array:  %.9g s\n\n", g_duration);

        g_duration = getIntervalDuration();
        printf("\tDevice buffer size:       %.3lf MB\n", (12.0 * double(numOps * sizeof(TYPE)) / 1.0e6));
        printf("\tTotal blocks:             %lli\n", blocks);
        printf("\tThreads per calculation:  %i\n", NTHREADS);
        printf("\tPrecision:                %dx%d\n", coil.thicknessIncrements, coil.angularIncrements);
        printf("\tTotal calculations:       %lli\n", numOps);
        printf("\tTotal MegaIncrements:     %.f\n", 1e-6 * double(numOps * coil.thicknessIncrements * coil.angularIncrements));
        printf("\n\tPerformance:              %.1f kPoints/s\n", double(0.001 * numOps / g_duration));
        printf("---------------------------------------------------\n\n");
    #endif
}
