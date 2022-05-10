#include "Benchmark.h"
#include "Coil.h"
#include "Tensor.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>
#include <chrono>


void benchForceZAxis(ComputeMethod computeMethod, int threadCount)
{
    using namespace std::chrono;

    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);
    primary.setThreadCount(threadCount);
    secondary.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.2));

    int opCount = 8192;
    double temp;

    printf("Expected execution time for one Ampere force z-axis calculation of specified precision\n");

    for (int i = 1; i <= 9; ++i)
    {
        int currentOperations = opCount / (int) pow(2, i);

        high_resolution_clock::time_point begin_time = high_resolution_clock::now();
        for (int j = 0; j < currentOperations; ++j)
            temp = Coil::computeAmpereForce(primary, secondary, PrecisionFactor(i), computeMethod).first.z;
        double interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("precisionFactor(%.1f) : %6.4f ms/op\n", (double) i, 1'000.0 * interval / currentOperations);
    }
}

void benchForceGeneral(ComputeMethod computeMethod, int threadCount)
{
    using namespace std::chrono;

    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);
    primary.setThreadCount(threadCount);

    primary.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.1, 0.0, 0.0));
    secondary.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.1, 0.0, 0.2));

    int opCount = 1024;
    std::pair<vec3::FieldVector3, vec3::FieldVector3> temp;

    printf("Expected execution time for one Ampere force general case calculation of specified precision\n");

    for (int i = 1; i <= 9; ++i)
    {
        int currentOperations = opCount / (int) pow(2, i);

        high_resolution_clock::time_point begin_time = high_resolution_clock::now();
        for (int j = 0; j < currentOperations; ++j)
            temp = Coil::computeAmpereForce(primary, secondary, PrecisionFactor(i), computeMethod);
        double interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("precisionFactor(%.1f) : %6.3f ms/op\n", (double) i, 1'000.0 * interval / currentOperations);
    }
}

void benchForceZAxisMTScaling(int maxThreadCount)
{
    printf("Performance comparison between different numbers of threads:\n");

    printf(" -> single thread:\n");
    benchForceZAxis(CPU_ST);
    printf("\n");

    for (int i = 2; i <= maxThreadCount; ++i)
    {
        printf(" -> %2d threads:\n", i);
        benchForceZAxis(CPU_MT, i);
        printf("\n");
    }
}

void benchForceGeneralMTScaling(int maxThreadCount)
{
    printf("Performance comparison between different numbers of threads:\n");

    printf(" -> single thread:\n");
    benchForceGeneral(CPU_ST);
    printf("\n");

    for (int i = 2; i <= maxThreadCount; ++i)
    {
        printf(" -> %2d threads:\n", i);
        benchForceGeneral(CPU_MT, i);
        printf("\n");
    }
}
