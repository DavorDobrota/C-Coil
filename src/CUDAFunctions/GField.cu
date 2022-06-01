#include "hardware_acceleration.h"

#include <cstdio>

#include "CUDAConstants.h"
#include "Timing.h"
#include "CUDAErrorCheck.h"


struct ParamG
{
    ParamG() :

        lengthIncrements{0},
        thicknessIncrements{0},
        angularIncrements{0},

        constantTerm(0.f),

        length(0.f),
        thickness(0.f),
        innerRadius(0.f),
        z{nullptr},
        r{nullptr},

        positions{-0.98940093499164993259615417, -0.94457502307323257607798842, -0.8656312023878317438804679,
              -0.7554044083550030338951012, -0.61787624440264374844667176, -0.4580167776572273863424194,
              -0.2816035507792589132304605, -0.09501250983763744018531934, 0.09501250983763744018531934,
              0.2816035507792589132304605, 0.4580167776572273863424194, 0.61787624440264374844667176,
              0.7554044083550030338951012, 0.8656312023878317438804679, 0.94457502307323257607798842,
              0.98940093499164993259615417},

        weights{0.02715245941175409485178057, 0.062253523938647892862843837, 0.09515851168249278480992511,
            0.12462897125553387205247628, 0.1495959888165767320815017, 0.16915651939500253818931208,
            0.18260341504492358886676367, 0.18945061045506849628539672, 0.1894506104550684962853967,
            0.1826034150449235888667637, 0.1691565193950025381893121, 0.14959598881657673208150173,
            0.1246289712555338720524763, 0.0951585116824927848099251, 0.06225352393864789286284384,
            0.027152459411754094851780572},

        cosPrecompute{0.999861409060662, 0.996212553862336, 0.977808137941974, 0.927094888788312,
                  0.825200872426836, 0.658971885399085, 0.428057063588869, 0.148691865890404,
                  -0.148691865890404, -0.428057063588869, -0.658971885399085, -0.825200872426836,
                  -0.927094888788312, -0.977808137941974, -0.996212553862336, -0.999861409060662},

        resG1{nullptr},
        resG2{nullptr},
        resG3{nullptr},
        resG4{nullptr}
    {}

    int lengthIncrements;
    int thicknessIncrements;
    int angularIncrements;

    TYPE constantTerm;

    TYPE length;
    TYPE thickness;
    TYPE innerRadius;

    TYPE *z;
    TYPE *r;

    TYPE positions[16];
    TYPE weights[16];
    TYPE cosPrecompute[16];

    TYPE *resG1;
    TYPE *resG2;
    TYPE *resG3;
    TYPE *resG4;
};


__global__
void calculateG(long long numOps, ParamG par)
{
    unsigned int index = threadIdx.x;
    long long global_index = blockIdx.x * blockDim.x + index;

    if(global_index >= numOps)
        return;

    TYPE zCoord = par.z[global_index];
    TYPE rCoord = par.r[global_index];

    TYPE bufferValueRP = 0.f;
    TYPE bufferValueRR = 0.f;
    TYPE bufferValueRZ = 0.f;
    TYPE bufferValueZZ = 0.f;

    TYPE constant = par.constantTerm;

    TYPE topEdge = zCoord + 0.5f * par.length;
    TYPE bottomEdge = zCoord - 0.5f * par.length;

    TYPE incrementPositionT, cosinePhi, cosinePhi2, phiExpression;
    TYPE tempConstA, tempConstB, tempConstC, tempConstD, tempConstE, tempConstF;
    TYPE tempConstG1, tempConstG2, tempConstI, tempConstJ1, tempConstJ2, tempConstK1, tempConstK2;
    TYPE tempConstL1, tempConstL2, tempConstM, tempConstN, tempConstO;

    for (int incT = 0; incT < par.thicknessIncrements; ++incT)
    {
        incrementPositionT = par.innerRadius + 0.5f * par.thickness * (1.0f + par.positions[incT]);

        tempConstA = incrementPositionT * incrementPositionT;
        tempConstB = rCoord * rCoord;
        tempConstC = incrementPositionT * rCoord;

        tempConstD = tempConstA + tempConstB;
        tempConstE = tempConstA * tempConstA + tempConstB * tempConstB;
        tempConstF = tempConstC * tempConstC;

        tempConstG1 = tempConstD + topEdge * topEdge;
        tempConstG2 = tempConstD + bottomEdge * bottomEdge;

        for (int incF = 0; incF < par.angularIncrements; ++incF)
        {
            cosinePhi = par.cosPrecompute[incF];
            cosinePhi2 = cosinePhi * cosinePhi;
            phiExpression = 2.f * tempConstC * cosinePhi;

            tempConstI = constant * par.weights[incT] * par.weights[incF];

            tempConstJ1 = tempConstG1 - phiExpression;
            tempConstJ2 = tempConstG2 - phiExpression;

            tempConstK1 = rsqrtf(tempConstJ1);
            tempConstK2 = rsqrtf(tempConstJ2);

            tempConstL1 = tempConstK1 / (tempConstJ1);
            tempConstL2 = tempConstK2 / (tempConstJ2);

            tempConstM = tempConstD - phiExpression;
            tempConstN =
                    2.f * tempConstF * cosinePhi * (cosinePhi2 + 2.f) -
                    tempConstC * (3.f * cosinePhi2 + 1.f) * tempConstD + cosinePhi * tempConstE;
            tempConstO = cosinePhi * tempConstD - 2.f * tempConstC;

            bufferValueRP += tempConstI * (incrementPositionT * cosinePhi / rCoord) * (tempConstK2 - tempConstK1);
            bufferValueRR += tempConstI * (tempConstC - tempConstA * cosinePhi) * cosinePhi * (tempConstL1 - tempConstL2);
            bufferValueZZ += tempConstI * (tempConstA - tempConstC * cosinePhi) * (tempConstL1 - tempConstL2);
            bufferValueRZ += tempConstI * incrementPositionT / (tempConstM * tempConstM) *
                             (topEdge * tempConstL1 * (tempConstO * tempConstJ1 + tempConstN) -
                              bottomEdge * tempConstL2 * (tempConstO * tempConstJ2 + tempConstN));
        }
    }
    par.resG1[global_index] = bufferValueRP;
    par.resG2[global_index] = bufferValueRR;
    par.resG3[global_index] = bufferValueRZ;
    par.resG4[global_index] = bufferValueZZ;
}
	
namespace 
{
    long long g_last_num_ops = 0;
//    long long g_last_blocks = 0;

    TYPE *g_zCoordArr = nullptr;
    TYPE *g_rCoordArr = nullptr;

    TYPE *g_tensorG1Arr = nullptr;
    TYPE *g_tensorG2Arr = nullptr;
    TYPE *g_tensorG3Arr = nullptr;
    TYPE *g_tensorG4Arr = nullptr;
    
    #if DEBUG_TIMINGS
        double g_duration;
    #endif
}

void resourceCleanupG()
{
	gpuErrchk(cudaFree(g_zCoordArr));
    gpuErrchk(cudaFree(g_rCoordArr));

    gpuErrchk(cudaFree(g_tensorG1Arr));
    gpuErrchk(cudaFree(g_tensorG2Arr));
    gpuErrchk(cudaFree(g_tensorG3Arr));
    gpuErrchk(cudaFree(g_tensorG4Arr));

    g_zCoordArr = nullptr;
    g_rCoordArr = nullptr;

    g_tensorG1Arr = nullptr;
    g_tensorG2Arr = nullptr;
    g_tensorG3Arr = nullptr;
    g_tensorG4Arr = nullptr;
}

void resourceStartupG(long long num_ops)
{
    resourceCleanupG();
    
	gpuErrchk(cudaMalloc(&g_zCoordArr, num_ops * sizeof(TYPE)));
    gpuErrchk(cudaMalloc(&g_rCoordArr, num_ops * sizeof(TYPE)));

    gpuErrchk(cudaMalloc(&g_tensorG1Arr, num_ops * sizeof(TYPE)));
    gpuErrchk(cudaMalloc(&g_tensorG2Arr, num_ops * sizeof(TYPE)));
    gpuErrchk(cudaMalloc(&g_tensorG3Arr, num_ops * sizeof(TYPE)));
    gpuErrchk(cudaMalloc(&g_tensorG4Arr, num_ops * sizeof(TYPE)));
}


void Calculate_hardware_accelerated_g
        (long long num_ops, const TYPE *z_ar, const TYPE *r_ar, TYPE current_density,
         TYPE innerRadius, TYPE length, TYPE thickness, int lengthIncrements,
         int thicknessIncrements, int angularIncrements, TYPE *gradientRPArr,
         TYPE *gradientRRArr, TYPE *gradientRZArr, TYPE *gradientZZArr)
{
    #if DEBUG_TIMINGS
        recordStartPoint();
    #endif

	ParamG par;

    par.lengthIncrements = INCREMENTCOUNT;
    par.thicknessIncrements = INCREMENTCOUNT;
    par.angularIncrements = INCREMENTCOUNT;

    par.length = length;
    par.thickness = thickness;
    par.innerRadius = innerRadius;

    par.constantTerm = MI * current_density * par.thickness * PI * 2 * 0.25;

    long long blocks = ceil(double(num_ops) / NTHREADS);

    if (num_ops > g_last_num_ops)
    {
        resourceStartupG(num_ops);
//        g_last_blocks = blocks;
    }

    par.z = g_zCoordArr;
    par.r = g_rCoordArr;

    par.resG1 = g_tensorG1Arr;
    par.resG2 = g_tensorG2Arr;
    par.resG3 = g_tensorG3Arr;
    par.resG4 = g_tensorG4Arr;

    gpuErrchk(cudaMemcpy(g_zCoordArr,z_ar,num_ops * sizeof(TYPE),cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(g_rCoordArr,r_ar,num_ops * sizeof(TYPE),cudaMemcpyHostToDevice));

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tMemory initialization:    %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    calculateG<<<blocks, NTHREADS>>>(num_ops, par);
	gpuErrchk(cudaDeviceSynchronize());

	#if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\t\tCalculations:     %.9g s\n", g_duration);
        printf("\t\tEstimated TFLOPS: %.2f\n",
               1e-12 * double(120 * num_ops * par.thicknessIncrements * par.angularIncrements) / g_duration);

        recordStartPoint();
    #endif

	if(gradientRPArr != nullptr)
    {
        gpuErrchk(cudaMemcpy(gradientRPArr, g_tensorG1Arr, num_ops * sizeof(TYPE), cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(gradientRRArr, g_tensorG2Arr, num_ops * sizeof(TYPE), cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(gradientRZArr, g_tensorG3Arr, num_ops * sizeof(TYPE), cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(gradientZZArr, g_tensorG4Arr, num_ops * sizeof(TYPE), cudaMemcpyDeviceToHost));
    }


    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tWriting to output array:  %.9g s\n\n", g_duration);
    #endif

	#if DEBUG_TIMINGS
        printf("\tDevice buffer size:       %.3lf MB\n", (6.0 * double(num_ops * sizeof(TYPE)) / 1.0e6));
        printf("\tTotal blocks:             %lli\n", blocks);
        printf("\tThreads per calculation:  %i\n", NTHREADS);
        printf("\tPrecision:                %dx%d\n", par.thicknessIncrements, par.angularIncrements);
        printf("\tTotal calculations:       %lli\n", num_ops);
        printf("\tTotal MegaIncrements:     %.f\n", 1e-6 * double(num_ops * par.thicknessIncrements * par.angularIncrements));
    #endif
}
