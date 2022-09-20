#include "Benchmark.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <cstdio>
#include <chrono>


void Benchmark::coilGroupComputeAllFieldsMTScaling(PrecisionFactor precisionFactor,
                                                   int threadCount, int coilCount, int maxPointsLog2)
{
    using namespace std::chrono;

    high_resolution_clock::time_point beginTime;
    double interval;
    double pointsPerSec;

    FILE *output = fopen("output.txt", "w");

    printf("Benchmarking expected CPU CoilGroup performance for a given number of points\n\n");

    double torusRadius = 1.0;

    CoilGroup torusGroupThick = CoilGroup();
    CoilGroup torusGroupFlat = CoilGroup();

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

    printf("Vector potential performance slow for precision factor %.1f and %d threads\n",
           precisionFactor.relativePrecision, threadCount);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int pointCount = int(std::pow(2, i));

        vec3::Vector3Array positions(pointCount);
        vec3::Vector3Array results;
        for (int j = 0; j < pointCount; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        results = torusGroupFlat.computeAllAPotentialVectors(positions, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = coilCount * pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic field performance slow for precision factor %.1f and %d threads\n",
           precisionFactor.relativePrecision, threadCount);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int pointCount = int(std::pow(2, i));

        vec3::Vector3Array positions(pointCount);
        vec3::Vector3Array results;
        for (int j = 0; j < pointCount; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        results = torusGroupFlat.computeAllBFieldVectors(positions, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = coilCount * pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic gradient performance slow for precision factor %.1f and %d threads\n",
           precisionFactor.relativePrecision, threadCount);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int pointCount = int(std::pow(2, i));

        vec3::Vector3Array positions(pointCount);
        vec3::Matrix3Array results;
        for (int j = 0; j < pointCount; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        results = torusGroupFlat.computeAllBGradientMatrices(positions, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = coilCount * pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Vector potential performance fast for precision factor %.1f and %d threads\n",
           precisionFactor.relativePrecision, threadCount);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int pointCount = int(std::pow(2, i));

        vec3::Vector3Array positions(pointCount);
        vec3::Vector3Array results;
        for (int j = 0; j < pointCount; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        results = torusGroupThick.computeAllAPotentialVectors(positions, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = coilCount * pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic field performance fast for precision factor %.1f and %d threads\n",
           precisionFactor.relativePrecision, threadCount);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int pointCount = int(std::pow(2, i));

        vec3::Vector3Array positions(pointCount);
        vec3::Vector3Array results;
        for (int j = 0; j < pointCount; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        results = torusGroupThick.computeAllBFieldVectors(positions, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = coilCount * pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic gradient performance fast for precision factor %.1f and %d threads\n",
           precisionFactor.relativePrecision, threadCount);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int pointCount = int(std::pow(2, i));

        vec3::Vector3Array positions(pointCount);
        vec3::Matrix3Array results;
        for (int j = 0; j < pointCount; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        results = torusGroupThick.computeAllBGradientMatrices(positions, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = coilCount * pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");

    fclose(output);
}

void Benchmark::coilGroupComputeAllFieldsGPUScaling(PrecisionFactor precisionFactor,
                                                    int coilCount, int maxPointsLog2)
{
    using namespace std::chrono;

    high_resolution_clock::time_point beginTime;
    double interval;
    double pointsPerSec;

    FILE *output = fopen("output.txt", "w");

    printf("Benchmarking expected GPU CoilGroup performance for a given number of points\n\n");

    double torusRadius = 1.0;

    CoilGroup torusGroupThick = CoilGroup();
    CoilGroup torusGroupFlat = CoilGroup();

    for (int i = 0; i < coilCount; ++i)
    {
        torusGroupThick.addCoil(
                torusRadius / 10.0, torusRadius / 100.0, torusRadius / 100.0, 10000,
                10, PrecisionFactor(), 8,
                vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / coilCount),
                M_PI_2, 2*M_PI * i / coilCount + M_PI_2
        );
    }
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
    torusGroupFlat.setDefaultPrecisionFactor(precisionFactor);

    printf("Vector potential performance slow for precision factor %.1f\n", precisionFactor.relativePrecision);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int pointCount = int(std::pow(2, i));

        vec3::Vector3Array positions(pointCount);
        vec3::Vector3Array results;
        for (int j = 0; j < pointCount; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        results = torusGroupFlat.computeAllAPotentialVectors(positions, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = coilCount * pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic field performance slow for precision factor %.1f\n", precisionFactor.relativePrecision);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int pointCount = int(std::pow(2, i));

        vec3::Vector3Array positions(pointCount);
        vec3::Vector3Array results;
        for (int j = 0; j < pointCount; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        results = torusGroupFlat.computeAllBFieldVectors(positions, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = coilCount * pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic gradient performance slow for precision factor %.1f\n", precisionFactor.relativePrecision);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int pointCount = int(std::pow(2, i));

        vec3::Vector3Array positions(pointCount);
        vec3::Matrix3Array results;
        for (int j = 0; j < pointCount; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        results = torusGroupFlat.computeAllBGradientMatrices(positions, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = coilCount * pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Vector potential performance fast for precision factor %.1f\n", precisionFactor.relativePrecision);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int pointCount = int(std::pow(2, i));

        vec3::Vector3Array positions(pointCount);
        vec3::Vector3Array results;
        for (int j = 0; j < pointCount; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        results = torusGroupThick.computeAllAPotentialVectors(positions, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = coilCount * pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic field performance fast for precision factor %.1f\n", precisionFactor.relativePrecision);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int pointCount = int(std::pow(2, i));

        vec3::Vector3Array positions(pointCount);
        vec3::Vector3Array results;
        for (int j = 0; j < pointCount; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        results = torusGroupThick.computeAllBFieldVectors(positions, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = coilCount * pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic gradient performance fast for precision factor %.1f\n", precisionFactor.relativePrecision);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int pointCount = int(std::pow(2, i));

        vec3::Vector3Array positions(pointCount);
        vec3::Matrix3Array results;
        for (int j = 0; j < pointCount; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        results = torusGroupThick.computeAllBGradientMatrices(positions, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = coilCount * pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");

    fclose(output);
}