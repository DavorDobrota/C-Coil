#include "Benchmark.h"
#include "Coil.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>
#include <chrono>


void benchMInductanceGeneral(ComputeMethod computeMethod, int threadCount)
{
    using namespace std::chrono;

    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);

    primary.setThreadCount(threadCount);
    primary.setPositionAndOrientation();
    secondary.setPositionAndOrientation(vec3::Vector3(0.1, 0.0, 0.2));

    int opCount = 1024;
    double temp;

    printf("Expected execution time for one MInductance general case calculation of specified precision\n");

    for (int i = 1; i <= 9; ++i)
    {
        int currentOperations = opCount / (int) pow(2, i);

        high_resolution_clock::time_point begin_time = high_resolution_clock::now();
        for (int j = 0; j < currentOperations; ++j)
            temp = Coil::computeMutualInductance(primary, secondary, PrecisionFactor(i), computeMethod);
        double interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("precisionFactor(%.1f) : %6.3f ms/op\n", (double) i, 1'000.0 * interval / currentOperations);

    }
}

void benchMInductanceGeneralMTScaling(int maxThreadCount)
{
    printf("Performance comparison between different numbers of threads:\n");

    printf(" -> single thread:\n");
    benchMInductanceGeneral(CPU_ST);
    printf("\n");

    for (int i = 2; i <= maxThreadCount; ++i)
    {
        printf(" -> %2d threads:\n", i);
        benchMInductanceGeneral(CPU_MT, i);
        printf("\n");
    }
}
