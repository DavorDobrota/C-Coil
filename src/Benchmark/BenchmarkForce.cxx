#include "Benchmark.h"
#include "Coil.h"
#include "Tensor.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>
#include <chrono>


void benchmarkAmpereForceZAxis(ComputeMethod computeMethod, int nThreads)
{
    using namespace std::chrono;

    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);
    primary.setThreadCount(nThreads);
    secondary.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.2));

    int nOps = 8192;
    double temp;

    printf("Expected execution time for one Ampere force z-axis calculation of specified precision\n");

    for (int i = 1; i <= 9; ++i)
    {
        int currentOperations = nOps / (int) pow(2, i);

        high_resolution_clock::time_point begin_time = high_resolution_clock::now();
        for (int j = 0; j < currentOperations; ++j)
            temp = Coil::computeAmpereForce(primary, secondary, PrecisionFactor(i), computeMethod).first.z;
        double interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("precisionFactor(%.1f) : %6.4f ms/op\n", (double) i, 1'000.0 * interval / currentOperations);
    }
}

void benchmarkAmpereForceGeneral(ComputeMethod computeMethod, int nThreads)
{
    using namespace std::chrono;

    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);
    primary.setThreadCount(nThreads);

    primary.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.1, 0.0, 0.0));
    secondary.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.1, 0.0, 0.2));

    int nOps = 1024;
    std::pair<vec3::FieldVector3, vec3::FieldVector3> temp;

    printf("Expected execution time for one Ampere force general case calculation of specified precision\n");

    for (int i = 1; i <= 9; ++i)
    {
        int currentOperations = nOps / (int) pow(2, i);

        high_resolution_clock::time_point begin_time = high_resolution_clock::now();
        for (int j = 0; j < currentOperations; ++j)
            temp = Coil::computeAmpereForce(primary, secondary, PrecisionFactor(i), computeMethod);
        double interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("precisionFactor(%.1f) : %6.3f ms/op\n", (double) i, 1'000.0 * interval / currentOperations);
    }
}

void benchmarkAmpereForceZAxisMTScaling(int maxThreads)
{
    printf("Performance comparison between different numbers of threads:\n");

    printf(" -> single thread:\n");
    benchmarkAmpereForceZAxis(CPU_ST);
    printf("\n");

    for (int i = 2; i <= maxThreads; ++i)
    {
        printf(" -> %2d threads:\n", i);
        benchmarkAmpereForceZAxis(CPU_MT, i);
        printf("\n");
    }
}

void benchmarkAmpereForceGeneralMTScaling(int maxThreads)
{
    printf("Performance comparison between different numbers of threads:\n");

    printf(" -> single thread:\n");
    benchmarkAmpereForceGeneral(CPU_ST);
    printf("\n");

    for (int i = 2; i <= maxThreads; ++i)
    {
        printf(" -> %2d threads:\n", i);
        benchmarkAmpereForceGeneral(CPU_MT, i);
        printf("\n");
    }
}
