#include "hardware_acceleration.h"

#include <cstdio>

#include "CUDAConstants.h"
#include "Timing.h"
#include "CUDAErrorCheck.h"

#define DEBUG_TIMINGS 1


struct ParamA
{
    ParamA() :
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

        res{nullptr}
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

    TYPE *res;
};


__global__
void calculateA(long long numOps, ParamA par)
{
    unsigned int index = threadIdx.x;
    long long global_index = blockIdx.x * blockDim.x + index;

    if(global_index >= numOps)
        return;

    TYPE zCoord = par.z[global_index];
    TYPE rCoord = par.r[global_index];

    TYPE potential = 0.0f;
    TYPE constant = par.constantTerm;

    TYPE topEdge = zCoord + 0.5f * par.length;
    TYPE bottomEdge = zCoord - 0.5f * par.length;

    TYPE incrementPositionT, cosinePhi, tempConstC, tempConstD1, tempConstD2, tempConstE1, tempConstE2, tempConstF;

    for (int incT = 0; incT < par.thicknessIncrements; ++incT)
    {
        incrementPositionT = par.innerRadius + 0.5f * par.thickness * (1.0f + par.positions[incT]);

        for (int incF = 0; incF < par.angularIncrements; ++incF)
        {
            cosinePhi = par.cosPrecompute[incF];

            tempConstC = rsqrtf(incrementPositionT * incrementPositionT + rCoord * rCoord -
                                     2.0f * incrementPositionT * rCoord * cosinePhi);

            tempConstD1 = topEdge * tempConstC;
            tempConstD2 = bottomEdge * tempConstC;

            tempConstE1 = sqrtf(tempConstD1 * tempConstD1 + 1.0f);
            tempConstE2 = sqrtf(tempConstD2 * tempConstD2 + 1.0f);

            tempConstF = logf((tempConstE1 + tempConstD1) / (tempConstE2 + tempConstD2));

            potential += constant * par.weights[incT] * par.weights[incF] * incrementPositionT * cosinePhi * tempConstF;
        }
    }
    par.res[global_index] = potential;
}
	
namespace 
{
    long long g_last_num_ops = 0;
//    long long g_last_blocks = 0;
    TYPE *g_resultArr = nullptr;
    TYPE *g_zCoordArr = nullptr;
    TYPE *g_rCoordArr = nullptr;
    
    #if DEBUG_TIMINGS
        double g_duration;
    #endif
}

void resourceCleanupA()
{
	gpuErrchk(cudaFree(g_zCoordArr));
    gpuErrchk(cudaFree(g_rCoordArr));
    gpuErrchk(cudaFree(g_resultArr));

    g_zCoordArr = nullptr;
    g_rCoordArr = nullptr;
    g_resultArr = nullptr;
}

void resourceStartupA(long long num_ops)
{
    resourceCleanupA();
    
	gpuErrchk(cudaMalloc(&g_zCoordArr, num_ops * sizeof(TYPE)));
    gpuErrchk(cudaMalloc(&g_rCoordArr, num_ops * sizeof(TYPE)));
    gpuErrchk(cudaMalloc(&g_resultArr, num_ops * sizeof(TYPE)));
}


void Calculate_hardware_accelerated_a
        (long long num_ops, const TYPE *z_ar, const TYPE *r_ar, TYPE current_density,
         TYPE innerRadius, TYPE length, TYPE thickness, int lengthIncrements,
         int thicknessIncrements, int angularIncrements, TYPE *potentialArray) {
#if DEBUG_TIMINGS
    recordStartPoint();
#endif

    ParamA par;

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
        resourceStartupA(num_ops);
        g_last_num_ops = num_ops;
//        g_last_blocks = blocks;
    }

    par.z = g_zCoordArr;
    par.r = g_rCoordArr;
    par.res = g_resultArr;

    gpuErrchk(cudaMemcpy(g_zCoordArr,z_ar,num_ops * sizeof(TYPE),cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(g_rCoordArr,r_ar,num_ops * sizeof(TYPE),cudaMemcpyHostToDevice));

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tMemory initialization:    %.9g s\n", g_duration);

        recordStartPoint();
    #endif

    calculateA<<<blocks, NTHREADS>>>(num_ops, par);
	gpuErrchk(cudaDeviceSynchronize());

	#if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\t\tCalculations:     %.9g s\n", g_duration);
        printf("\t\tEstimated TFLOPS: %.2f\n",
               1e-12 * double(80 * num_ops * par.thicknessIncrements * par.angularIncrements) / g_duration);

        recordStartPoint();
    #endif

	if(potentialArray != nullptr)
        gpuErrchk(cudaMemcpy(potentialArray, g_resultArr, num_ops * sizeof(TYPE), cudaMemcpyDeviceToHost));

    #if DEBUG_TIMINGS
        g_duration = getIntervalDuration();
        printf("\tWriting to output array:  %.9g s\n\n", g_duration);
    #endif

	#if DEBUG_TIMINGS
        printf("\tDevice buffer size:       %.3lf MB\n", (3.0 * double(num_ops * sizeof(TYPE)) / 1.0e6));
        printf("\tTotal blocks:             %lli\n", blocks);
        printf("\tThreads per calculation:  %i\n", NTHREADS);
        printf("\tPrecision:                %dx%d\n", par.thicknessIncrements, par.angularIncrements);
        printf("\tTotal calculations:       %lli\n", num_ops);
        printf("\tTotal MegaIncrements:     %.f\n", 1e-6 * double(num_ops * par.thicknessIncrements * par.angularIncrements));
    #endif
}
