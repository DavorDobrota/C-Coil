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
    benchComputeAllFieldsWorkloadScalingMT(PrecisionFactor(1.0), 12, 26);
    benchComputeAllFieldsWorkloadScalingGPU(PrecisionFactor(3.0), 26);
    testRawGPUPerformance();

    return 0;
}
#pragma clang diagnostic pop
