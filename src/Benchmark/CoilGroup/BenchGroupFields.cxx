#include "Benchmark.h"
#include "Compare.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <cstdio>
#include <chrono>


void Benchmark::coilGroupComputeAllFields(PrecisionFactor precisionFactor, int coilCount,
                                          int opCount, int threadCount)
{
    using namespace std::chrono;

    high_resolution_clock::time_point begin_time;
    double interval;

    printf("Benchmarking MT and GPU field compute performance for %d coils in %d points and precision factor %.1f\n\n",
           coilCount, opCount, precisionFactor.relativePrecision);

    double torusRadius = 1.0;

    CoilGroup torusGroupThick = CoilGroup();
    CoilGroup torusGroupFlat = CoilGroup();

    vec3::Vector3Array fieldPoints(opCount);
    vec3::Vector3Array computedAPotential;
    vec3::Vector3Array computedBField;
    vec3::Matrix3Array computedGradient;

    for (int i = 0; i < opCount; ++i)
        fieldPoints[i] = vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / opCount);

    for (int i = 0; i < coilCount; ++i)
    {
        torusGroupThick.addCoil(
                torusRadius / 10.0, torusRadius / 100.0, torusRadius / 100.0, 10000,
                10, PrecisionFactor(), 8,
                vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / coilCount),
                M_PI_2, 2*M_PI * i / coilCount + M_PI_2
        );
    }
    torusGroupThick.setThreadCount(threadCount);
    torusGroupThick.setDefaultPrecisionFactor(precisionFactor);

    for (int i = 0; i < coilCount; ++i)
    {
        torusGroupFlat.addCoil(
                torusRadius / 10.0, torusRadius / 100.0, 0.0, 10000,
                10, PrecisionFactor(), 8,
                vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / coilCount),
                M_PI_2, 2*M_PI * i / coilCount + M_PI_2
        );
    }
    torusGroupFlat.setThreadCount(threadCount);
    torusGroupFlat.setDefaultPrecisionFactor(precisionFactor);

    // MT slow tests
    begin_time = high_resolution_clock::now();
    computedAPotential = torusGroupFlat.computeAllAPotentialVectors(fieldPoints, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_MT slow : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    begin_time = high_resolution_clock::now();
    computedBField = torusGroupFlat.computeAllBFieldVectors(fieldPoints, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Field     B CPU_MT slow : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    begin_time = high_resolution_clock::now();
    computedGradient = torusGroupFlat.computeAllBGradientMatrices(fieldPoints, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G CPU_MT slow : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    printf("\n");

    // MT fast tests
    begin_time = high_resolution_clock::now();
    computedAPotential = torusGroupThick.computeAllAPotentialVectors(fieldPoints, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_MT fast : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    begin_time = high_resolution_clock::now();
    computedBField = torusGroupThick.computeAllBFieldVectors(fieldPoints, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Field     B CPU_MT fast : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    begin_time = high_resolution_clock::now();
    computedGradient = torusGroupThick.computeAllBGradientMatrices(fieldPoints, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G CPU_MT fast : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    printf("\n");

    // GPU slow tests
    computedAPotential = torusGroupThick.computeAllAPotentialVectors(fieldPoints, GPU); // warmup

    begin_time = high_resolution_clock::now();
    computedAPotential = torusGroupFlat.computeAllAPotentialVectors(fieldPoints, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A GPU    slow : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    begin_time = high_resolution_clock::now();
    computedBField = torusGroupFlat.computeAllBFieldVectors(fieldPoints, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Field     B GPU    slow : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    begin_time = high_resolution_clock::now();
    computedGradient = torusGroupFlat.computeAllBGradientMatrices(fieldPoints, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G GPU    slow : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    printf("\n");

    // GPU fast tests

    begin_time = high_resolution_clock::now();
    computedAPotential = torusGroupThick.computeAllAPotentialVectors(fieldPoints, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A GPU    fast : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    begin_time = high_resolution_clock::now();
    computedBField = torusGroupThick.computeAllBFieldVectors(fieldPoints, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Field     B GPU    fast : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    begin_time = high_resolution_clock::now();
    computedGradient = torusGroupThick.computeAllBGradientMatrices(fieldPoints, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G GPU    fast : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    printf("\n");
}

void Benchmark::coilGroupComputeAllFieldsMTvsMTD(int threadCount, int pointCount)
{
    using namespace std::chrono;

    high_resolution_clock::time_point begin_time;
    double interval;

    int coilCount1 = threadCount;
    int coilCount2 = 3 * threadCount;

    printf("Quick performance benchmark for %d coils and %d points\n\n", threadCount, pointCount);

    begin_time = high_resolution_clock::now();
    Compare::fieldsCoilGroupMTD(coilCount1, pointCount, 1, false);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("ST  perf : %.0f kPoints/s\n", 1e-3 * coilCount1 * pointCount / interval);

    begin_time = high_resolution_clock::now();
    Compare::fieldsCoilGroupMTD(coilCount1, pointCount, threadCount, false);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MT  perf : %.0f kPoints/s\n", 1e-3 * coilCount1 * pointCount / interval);

    begin_time = high_resolution_clock::now();
    Compare::fieldsCoilGroupMTD(coilCount2, pointCount, threadCount, false);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MTD perf : %.0f kPoints/s\n", 1e-3 * coilCount2 * pointCount / interval);

    printf("\n");
}

void Benchmark::coilGroupComputeAllFieldsGPU(int coilCount, int opCount)
{
    using namespace std::chrono;

    high_resolution_clock::time_point begin_time;
    double interval;

    printf("Benchmarking GPU field compute performance for %d coils in %d points\n\n", coilCount, opCount);

    double torusRadius = 1.0;

    CoilGroup torusGroupThick = CoilGroup();
    CoilGroup torusGroupFlat = CoilGroup();

    vec3::Vector3Array fieldPoints(opCount);
    vec3::Vector3Array computedAPotential;
    vec3::Vector3Array computedBField;
    vec3::Matrix3Array computedGradient;

    for (int i = 0; i < opCount; ++i)
        fieldPoints[i] = vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / opCount);

    for (int i = 0; i < coilCount; ++i)
    {
        torusGroupThick.addCoil(
                torusRadius / 10.0, torusRadius / 100.0, torusRadius / 100.0, 10000,
                10, PrecisionFactor(), 8,
                vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / coilCount),
                M_PI_2, 2*M_PI * i / coilCount + M_PI_2
        );
    }

    for (int i = 0; i < coilCount; ++i)
    {
        torusGroupFlat.addCoil(
                torusRadius / 10.0, torusRadius / 100.0, 0.0, 10000,
                10, PrecisionFactor(), 8,
                vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / coilCount),
                M_PI_2, 2*M_PI * i / coilCount + M_PI_2
        );
    }

    computedAPotential = torusGroupFlat.computeAllAPotentialVectors(fieldPoints, GPU); // warmup

    for (int i = 1; i <= 7; ++i)
    {
        torusGroupFlat.setDefaultPrecisionFactor(PrecisionFactor(double(i)));
        torusGroupThick.setDefaultPrecisionFactor(PrecisionFactor(double(i)));

        printf("Precision factor %.1f\n", double(i));

        begin_time = high_resolution_clock::now();
        computedAPotential = torusGroupFlat.computeAllAPotentialVectors(fieldPoints, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("Potential A  slow : %.1f kPoints/s | eff %.1f MPoints/s\n",
               0.001 * opCount / interval, 1e-6 * opCount * coilCount / interval);

        begin_time = high_resolution_clock::now();
        computedBField = torusGroupFlat.computeAllBFieldVectors(fieldPoints, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("Field     B  slow : %.1f kPoints/s | eff %.1f MPoints/s\n",
               0.001 * opCount / interval, 1e-6 * opCount * coilCount / interval);

        begin_time = high_resolution_clock::now();
        computedGradient = torusGroupFlat.computeAllBGradientMatrices(fieldPoints, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("Gradient  G  slow : %.1f kPoints/s | eff %.1f MPoints/s\n",
               0.001 * opCount / interval, 1e-6 * opCount * coilCount / interval);

        begin_time = high_resolution_clock::now();
        computedAPotential = torusGroupThick.computeAllAPotentialVectors(fieldPoints, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("Potential A  fast : %.1f kPoints/s | eff %.1f MPoints/s\n",
               0.001 * opCount / interval, 1e-6 * opCount * coilCount / interval);

        begin_time = high_resolution_clock::now();
        computedBField = torusGroupThick.computeAllBFieldVectors(fieldPoints, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("Field     B  fast : %.1f kPoints/s | eff %.1f MPoints/s\n",
               0.001 * opCount / interval, 1e-6 * opCount * coilCount / interval);

        begin_time = high_resolution_clock::now();
        computedGradient = torusGroupThick.computeAllBGradientMatrices(fieldPoints, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("Gradient  G  fast : %.1f kPoints/s | eff %.1f MPoints/s\n",
               0.001 * opCount / interval, 1e-6 * opCount * coilCount / interval);

        printf("\n");
    }
}
