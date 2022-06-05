#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>

#include "Benchmark.h"
#include "Coil.h"
#include "CoilGroup.h"
#include "Test.h"
#include "Compare.h"


#pragma clang diagnostic push
#pragma ide diagnostic ignored "Simplify"
int main()
{
//    benchComputeAllFieldsWorkloadScalingMT(PrecisionFactor(3.0), 8, 22);
//    benchComputeAllFieldsWorkloadScalingGPU(PrecisionFactor(3.0), 26);
//    benchComputeAllFields();

//    benchComputeAllFieldsEveryCoilType(200003, 8);

//    benchMInductanceZAxisMTScaling(16);
    benchMInductanceGeneralMTScaling(16);
//    benchSelfInductance();
//    benchForceZAxisMTScaling(16);
//    benchForceGeneralMTScaling(16);

//    compMInductanceGeneralCase();

    compPrecisionCPUvsGPU();


    return 0;
}
#pragma clang diagnostic pop
