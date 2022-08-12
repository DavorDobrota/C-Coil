#include "Benchmark.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <chrono>


void benchComputeAllFieldsWorkloadScalingMT(PrecisionFactor precisionFactor, int threadCount, int maxPointsLog2)
{
    using namespace std::chrono;

    FILE *output = fopen("output.txt", "w");

    printf("Benchmarking expected CPU performance for a given number of points\n\n");

    Coil coil = Coil(0.1, 0.1, 0.1, 10000);
    coil.setThreadCount(threadCount);
    coil.setDefaultPrecision(precisionFactor);

    Coil flat = Coil(0.1, 0.1, 0, 100);
    flat.setThreadCount(threadCount);
    flat.setDefaultPrecision(precisionFactor);

    high_resolution_clock::time_point beginTime;
    double interval;
    double pointsPerSec;

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
        results = flat.computeAllAPotentialVectors(positions, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = pointCount / interval;

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
        results = flat.computeAllBFieldVectors(positions, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = pointCount / interval;

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
        results = flat.computeAllBGradientMatrices(positions, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = pointCount / interval;

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
        results = coil.computeAllAPotentialVectors(positions, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = pointCount / interval;

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
        results = coil.computeAllBFieldVectors(positions, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = pointCount / interval;

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
        results = coil.computeAllBGradientMatrices(positions, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");

    fclose(output);
}

void benchComputeAllFieldsWorkloadScalingGPU(PrecisionFactor precisionFactor, int maxPointsLog2)
{
    using namespace std::chrono;

    FILE *output = fopen("output.txt", "w");

    printf("Benchmarking expected GPU performance for a given number of points\n\n");

    Coil coil = Coil(0.1, 0.1, 0.1, 10000);
    coil.setDefaultPrecision(precisionFactor);

    Coil flat = Coil(0.1, 0.1, 0, 100);
    flat.setDefaultPrecision(precisionFactor);

    high_resolution_clock::time_point beginTime;
    double interval;
    double pointsPerSec;

    printf("Vector potential performance slow for precision factor %.1f\n",
           precisionFactor.relativePrecision);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int pointCount = int(std::pow(2, i));

        vec3::Vector3Array positions(pointCount);
        vec3::Vector3Array potentialArr;
        for (int j = 0; j < pointCount; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        potentialArr = flat.computeAllAPotentialVectors(positions, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic field performance slow for precision factor %.1f\n",
           precisionFactor.relativePrecision);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int pointCount = int(std::pow(2, i));

        vec3::Vector3Array positions(pointCount);
        vec3::Vector3Array fieldArr;
        for (int j = 0; j < pointCount; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        fieldArr = flat.computeAllBFieldVectors(positions, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic gradient performance slow for precision factor %.1f\n",
           precisionFactor.relativePrecision);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int pointCount = int(std::pow(2, i));

        vec3::Vector3Array positions(pointCount);
        vec3::Matrix3Array gradientArr;
        for (int j = 0; j < pointCount; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        gradientArr = flat.computeAllBGradientMatrices(positions, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Vector potential performance fast for precision factor %.1f\n",
           precisionFactor.relativePrecision);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int pointCount = int(std::pow(2, i));

        vec3::Vector3Array positions(pointCount);
        vec3::Vector3Array potentialArr;
        for (int j = 0; j < pointCount; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        potentialArr = coil.computeAllAPotentialVectors(positions, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic field performance fast for precision factor %.1f\n",
           precisionFactor.relativePrecision);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int pointCount = int(std::pow(2, i));

        vec3::Vector3Array positions(pointCount);
        vec3::Vector3Array fieldArr;
        for (int j = 0; j < pointCount; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        fieldArr = coil.computeAllBFieldVectors(positions, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic gradient performance fast for precision factor %.1f\n",
           precisionFactor.relativePrecision);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int pointCount = int(std::pow(2, i));

        vec3::Vector3Array positions(pointCount);
        vec3::Matrix3Array gradientArr;
        for (int j = 0; j < pointCount; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        gradientArr = coil.computeAllBGradientMatrices(positions, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = pointCount / interval;

        printf("%8d : %.1f\n", pointCount, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", pointCount, 0.001 * pointsPerSec);
    }
    printf("\n");

    fclose(output);
}
