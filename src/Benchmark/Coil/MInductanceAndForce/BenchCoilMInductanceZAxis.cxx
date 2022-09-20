#include "Benchmark.h"
#include "Coil.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>
#include <chrono>


void Benchmark::mInductanceZAxis(ComputeMethod computeMethod, int threadCount)
{
    using namespace std::chrono;

    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);

    primary.setThreadCount(threadCount);
    primary.setPositionAndOrientation();
    secondary.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.2));

    int opCount = 8192;
    double temp;

    printf("Expected execution time for one MInductance z-axis calculation of specified precision\n");

    for (int i = 1; i <= 9; ++i)
    {
        int currentOperations = opCount / (int) pow(2, i);

        high_resolution_clock::time_point begin_time = high_resolution_clock::now();
        for (int j = 0; j < currentOperations; ++j)
            temp = Coil::computeMutualInductance(primary, secondary, PrecisionFactor(i), computeMethod);
        double interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("precisionFactor(%.1f) : %6.4f ms/op\n", (double) i, 1'000.0 * interval / currentOperations);
    }
}

void Benchmark::mInductanceZAxisMTScaling(int maxThreadCount)
{
    printf("Performance comparison between different numbers of threads:\n");

    printf(" -> single thread:\n");
    Benchmark::mInductanceZAxis(CPU_ST);
    printf("\n");

    for (int i = 2; i <= maxThreadCount; ++i)
    {
        printf(" -> %2d threads:\n", i);
        Benchmark::mInductanceZAxis(CPU_MT, i);
        printf("\n");
    }
}

void Benchmark::selfInductance()
{
    using namespace std::chrono;

    printf("Expected execution time for one self inductance calculation of specified precision\n");

    Coil coil = Coil(0.1, 0.1, 0.1, 10000);
    double temp;

    int opCount = 131'072;

    for (int i = 1; i <= 15; ++i)
    {
        int currentOperations = opCount / (int) pow(2, i);

        high_resolution_clock::time_point begin_time = high_resolution_clock::now();
        for (int j = 0; j < currentOperations; ++j)
            temp = coil.computeAndSetSelfInductance(PrecisionFactor(i));
        double interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();

        printf("precisionFactor(%.1f) : %6.4f ms/op\n", (double) i,  1'000.0 * interval / currentOperations);
    }
    printf("\n");
}
