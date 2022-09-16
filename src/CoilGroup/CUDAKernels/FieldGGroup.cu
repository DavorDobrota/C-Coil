#include "CoilGroupAcceleration.h"

#include "Timing.h"
#include "CUDAUtils/ErrorCheck/CUDAErrorCheck.h"
#include "CUDAUtils/MemoryManagement/GPUMemoryManagement.h"

#include <cstdio>


__global__
void calculateGradientGroup(long long opCount, long long coilIndex,
                            const CoilData *coilArr,
                            const VectorData *posArr,
                            MatrixData *resArr)
{
    unsigned int index = threadIdx.x;
    long long globalIndex = blockIdx.x * blockDim.x + index;

    if(globalIndex >= opCount)
        return;

    __shared__ CoilData coil;
    coil = coilArr[coilIndex];

    TYPE x1 = posArr[globalIndex].x;
    TYPE y1 = posArr[globalIndex].y;
    TYPE z1 = posArr[globalIndex].z;

    x1 -= coil.positionVector[0];
    y1 -= coil.positionVector[1];
    z1 -= coil.positionVector[2];

    TYPE x = x1 * coil.invTransformArray[0] + y1 * coil.invTransformArray[1] + z1 * coil.invTransformArray[2];
    TYPE y = x1 * coil.invTransformArray[3] + y1 * coil.invTransformArray[4] + z1 * coil.invTransformArray[5];
    TYPE z = x1 * coil.invTransformArray[6] + y1 * coil.invTransformArray[7] + z1 * coil.invTransformArray[8];

    TYPE zCoord = z;
    TYPE rCoord = sqrt(x * x + y * y);
    TYPE phiCord = atan2(y, x);

    TYPE bufferValueRP = 0.0f;
    TYPE bufferValueRR = 0.0f;
    TYPE bufferValueRZ = 0.0f;
    TYPE bufferValueZZ = 0.0f;

    TYPE topEdge = zCoord + 0.5f * coil.length;
    TYPE bottomEdge = zCoord - 0.5f * coil.length;

    if (coil.useFastMethod)
    {
        for (int incT = 0; incT < coil.thicknessIncrements; ++incT)
        {
            TYPE incrementPositionT = coil.innerRadius + 0.5f * coil.thickness * (1.0f + coil.thicknessPositionArray[incT]);

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
                TYPE phiExpression = 2.0f * tempConstC * cosinePhi;

                TYPE tempConstI = coil.constFactor * coil.thicknessWeightArray[incT] * coil.angularWeightArray[incF];

                TYPE tempConstJ1 = tempConstG1 - phiExpression;
                TYPE tempConstJ2 = tempConstG2 - phiExpression;

                TYPE tempConstK1 = rsqrt(tempConstJ1);
                TYPE tempConstK2 = rsqrt(tempConstJ2);

                TYPE tempConstL1 = tempConstK1 / (tempConstJ1);
                TYPE tempConstL2 = tempConstK2 / (tempConstJ2);

                TYPE tempConstM = tempConstD - phiExpression;
                TYPE tempConstN =
                        2.0f * tempConstF * cosinePhi * (cosinePhi2 + 2.0f) -
                        tempConstC * (3.0f * cosinePhi2 + 1.0f) * tempConstD + cosinePhi * tempConstE;
                TYPE tempConstO = cosinePhi * tempConstD - 2.0f * tempConstC;

                bufferValueRP += tempConstI * (incrementPositionT * cosinePhi / rCoord) * (tempConstK2 - tempConstK1);
                bufferValueRR += tempConstI * (tempConstC - tempConstA * cosinePhi) * cosinePhi * (tempConstL1 - tempConstL2);
                bufferValueZZ += tempConstI * (tempConstA - tempConstC * cosinePhi) * (tempConstL1 - tempConstL2);
                bufferValueRZ += tempConstI * incrementPositionT / (tempConstM * tempConstM) *
                                 (topEdge * tempConstL1 * (tempConstO * tempConstJ1 + tempConstN) -
                                  bottomEdge * tempConstL2 * (tempConstO * tempConstJ2 + tempConstN));
            }
        }
    }
    else
    {
        for (int incT = 0; incT < coil.thicknessIncrements; ++incT)
        {
            TYPE incrementPositionT = coil.innerRadius + 0.5f * coil.thickness * (1.0f + coil.thicknessPositionArray[incT]);

            TYPE tempConstA = incrementPositionT * incrementPositionT;
            TYPE tempConstB = rCoord * rCoord;
            TYPE tempConstC = zCoord * zCoord;
            TYPE tempConstD = rCoord * incrementPositionT;
            TYPE tempConstE = incrementPositionT * zCoord;

            TYPE tempConstF = 2.0f * tempConstA + 2.0f * tempConstB - tempConstC;
            TYPE tempConstG = tempConstA + tempConstB + tempConstC;

            for (int incF = 0; incF < coil.angularIncrements; ++incF)
            {
                TYPE cosinePhi = coil.cosPrecomputeArray[incF];
                TYPE tempConstH = coil.constFactor * coil.thicknessWeightArray[incT] * coil.angularWeightArray[incF];

                TYPE tempConstI = tempConstG - 2.0f * tempConstD * cosinePhi;
                TYPE tempConstJ = tempConstI * sqrt(tempConstI);

                TYPE tempConstX = tempConstH / (tempConstJ);
                TYPE tempConstY = tempConstH / (tempConstJ * tempConstI);

                bufferValueRP += tempConstX * tempConstE * cosinePhi / rCoord;
                bufferValueRR += tempConstY * (-3.0f * tempConstE * (rCoord - incrementPositionT * cosinePhi)) * cosinePhi;
                bufferValueZZ += tempConstY * (-3.0f * tempConstE * (incrementPositionT - rCoord * cosinePhi));
                bufferValueRZ += tempConstY *
                                 (incrementPositionT * (tempConstF - tempConstD * cosinePhi) * cosinePhi - 3.0f * tempConstA * rCoord);
            }
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

        xyGrad = 0.5f * sin(2.0f * phiCord) * (bufferValueRR - bufferValueRP);
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

        xyGrad = 0.0f;
        xzGrad = 0.0f;
        yxGrad = 0.0f;
        yzGrad = 0.0f;
        zxGrad = 0.0f;
        zyGrad = 0.0f;
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

    resArr[globalIndex].xx += xxRes;
    resArr[globalIndex].xy += xyRes;
    resArr[globalIndex].xz += xzRes;
    resArr[globalIndex].yx += yxRes;
    resArr[globalIndex].yy += yyRes;
    resArr[globalIndex].yz += yzRes;
    resArr[globalIndex].zx += zxRes;
    resArr[globalIndex].zy += zyRes;
    resArr[globalIndex].zz += zzRes;
}


namespace
{
    CoilData *g_coilArr = nullptr;
    VectorData *g_posArr = nullptr;
    MatrixData *g_resArr = nullptr;

    void getBuffers(long long coilCount, long long opCount)
    {
        std::vector<void*> buffers = GPUMem::getBuffers(
                { coilCount * (long long)sizeof(CoilData),
                  opCount * (long long)sizeof(VectorData),
                  opCount * (long long)sizeof(MatrixData)}
        );

        g_coilArr = static_cast<CoilData*>(buffers[0]);
        g_posArr = static_cast<VectorData*>(buffers[1]);
        g_resArr = static_cast<MatrixData*>(buffers[2]);
    }

    #if DEBUG_TIMINGS
        double g_duration;
    #endif
}

void Calculate_hardware_accelerated_g_group(long long coilCount, long long opCount,
                                            const CoilData *coilArr,
                                            const VectorData *posArr,
                                            MatrixData *resArr)
{
    #if DEBUG_TIMINGS
        recordStartPoint();
        recordStartPoint();
    #endif

    long long blocks = ceil(double(opCount) / NTHREADS);

    getBuffers(coilCount, opCount);

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tResource startup:         %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    gpuErrchk(cudaMemcpy(g_coilArr, coilArr, coilCount * sizeof(CoilData), cudaMemcpyHostToDevice))
    gpuErrchk(cudaMemcpy(g_posArr, posArr, opCount * sizeof(VectorData), cudaMemcpyHostToDevice))

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tMemory initialization:    %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    gpuErrchk(cudaMemset(g_resArr, 0, opCount * sizeof(MatrixData)))

    for (int i = 0; i < coilCount; ++i)
    {
        calculateGradientGroup<<<blocks, NTHREADS>>>(
            opCount, i, g_coilArr, g_posArr, g_resArr
        );
        gpuErrchk(cudaDeviceSynchronize())
    }

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tCalculations:             %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    if(resArr != nullptr)
        gpuErrchk(cudaMemcpy(resArr, g_resArr, opCount * sizeof(MatrixData), cudaMemcpyDeviceToHost))

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tWriting to output array:  %.9g s\n\n", g_duration);

        g_duration = getIntervalDuration();
        printf("\tDevice buffer size:       %.3lf MB\n", (12.0 * double(opCount * sizeof(TYPE) + coilCount * sizeof(CoilData)) / 1.0e6));
        printf("\tTotal blocks:             %lli\n", blocks);
        printf("\tThreads per calculation:  %i\n", NTHREADS);
        printf("\tTotal coils:              %lli\n", coilCount);
        printf("\tTotal points:             %lli\n", opCount);
        printf("\tTotal calculations:       %lli\n", opCount * coilCount);
        printf("\n\tPerformance:              %.1f kPoints/s\n", double(0.001 * opCount / g_duration));
        printf("\n\tEffectivePerformance:     %.1f kPoints/s\n", double(0.001 * opCount * coilCount / g_duration));
        printf("---------------------------------------------------\n\n");
    #endif
}
