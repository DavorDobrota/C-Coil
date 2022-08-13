#include "Benchmark.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <chrono>


void benchComputeAllFieldsEveryCoilType(int opCount, int threadCount)
{
    using namespace std::chrono;

    FILE *output = fopen("output.txt", "w");

    Coil loop = Coil(0.1, 0.0, 0.0, 1);
    Coil pancake = Coil(0.1, 0.4, 0.0, 40);
    Coil thin = Coil(0.1, 0.0, 0.1, 100);
    Coil thick = Coil(0.1, 0.4, 0.4, 1600);


    high_resolution_clock::time_point beginTime;
    double interval, incrementsPerSec, pointsPerSec;
    vec3::Vector3Array positions(opCount);

    for (int i = 0; i < opCount; ++i)
        positions[i] = vec3::Vector3::getFromSphericalCoords(1.0, M_PI * i / opCount, 0.0);

    vec3::Vector3Array potentialArr;
    vec3::Vector3Array fieldArr;
    vec3::Matrix3Array gradientArr;

    printf("This test is created for the purpose of generating performance charts\n\n");

    printf("Loop potential performance:\n");
    for(int i = 1; i <= 8; ++i)
    {
        auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(loop, PrecisionFactor(i));
        int numIterations = precision.angularBlocks * precision.angularIncrements;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        potentialArr = loop.computeAllAPotentialVectors(positions, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = double(totalIterations) / interval;
        pointsPerSec = opCount / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= threadCount; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision =
                PrecisionArguments::getCoilPrecisionArgumentsCPU(loop, PrecisionFactor(i));
            int numIterations = precision.angularBlocks * precision.angularIncrements;
            long long totalIterations = (long long) opCount * numIterations;
            loop.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            potentialArr = loop.computeAllAPotentialVectors(positions, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = double(totalIterations) / interval;
            pointsPerSec = opCount / interval;
            printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
            fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
        }
        printf("\n");
        fprintf(output, "\n");
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Loop field performance:\n");
    for (int i = 1; i <= 8; ++i)
    {
        auto precision =
            PrecisionArguments::getCoilPrecisionArgumentsCPU(loop, PrecisionFactor(i));
        int numIterations = precision.angularBlocks * precision.angularIncrements;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        fieldArr = loop.computeAllBFieldVectors(positions, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = double(totalIterations) / interval;
        pointsPerSec = opCount / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= threadCount; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(loop, PrecisionFactor(i));
            int numIterations = precision.angularBlocks * precision.angularIncrements;
            long long totalIterations = (long long) opCount * numIterations;
            loop.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            fieldArr = loop.computeAllBFieldVectors(positions, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = double(totalIterations) / interval;
            pointsPerSec = opCount / interval;
            printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
            fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
        }
        printf("\n");
        fprintf(output, "\n");
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Loop gradient performance:\n");
    for (int i = 1; i <= 8; ++i)
    {
        auto precision =
            PrecisionArguments::getCoilPrecisionArgumentsCPU(loop, PrecisionFactor(i));
        int numIterations = precision.angularBlocks * precision.angularIncrements;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        gradientArr = loop.computeAllBGradientMatrices(positions, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = double(totalIterations) / interval;
        pointsPerSec = opCount / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= threadCount; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision =
                PrecisionArguments::getCoilPrecisionArgumentsCPU(loop, PrecisionFactor(i));
            int numIterations = precision.angularBlocks * precision.angularIncrements;
            long long totalIterations = (long long) opCount * numIterations;
            loop.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            gradientArr = loop.computeAllBGradientMatrices(positions, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = double(totalIterations) / interval;
            pointsPerSec = opCount / interval;
            printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
            fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
        }
        printf("\n");
        fprintf(output, "\n");
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Pancake potential performance:\n");
    for(int i = 1; i <= 8; ++i)
    {
        auto precision =
            PrecisionArguments::getCoilPrecisionArgumentsCPU(pancake, PrecisionFactor(i));
        int numIterations =
                precision.angularBlocks * precision.angularIncrements
                * precision.thicknessBlocks * precision.thicknessIncrements;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        potentialArr = pancake.computeAllAPotentialVectors(positions, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = double(totalIterations) / interval;
        pointsPerSec = opCount / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= threadCount; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision =
                PrecisionArguments::getCoilPrecisionArgumentsCPU(pancake, PrecisionFactor(i));
            int numIterations = precision.angularBlocks * precision.angularIncrements *
                                precision.thicknessBlocks * precision.thicknessIncrements;
            long long totalIterations = (long long) opCount * numIterations;
            pancake.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            potentialArr = pancake.computeAllAPotentialVectors(positions, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = double(totalIterations) / interval;
            pointsPerSec = opCount / interval;
            printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
            fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
        }
        printf("\n");
        fprintf(output, "\n");
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Pancake field performance:\n");
    for (int i = 1; i <= 8; ++i)
    {
        auto precision =
            PrecisionArguments::getCoilPrecisionArgumentsCPU(pancake, PrecisionFactor(i));
        int numIterations = precision.angularBlocks * precision.angularIncrements *
                            precision.thicknessBlocks * precision.thicknessIncrements;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        fieldArr = pancake.computeAllBFieldVectors(positions, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = double(totalIterations) / interval;
        pointsPerSec = opCount / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= threadCount; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision =
                PrecisionArguments::getCoilPrecisionArgumentsCPU(pancake, PrecisionFactor(i));
            int numIterations = precision.angularBlocks * precision.angularIncrements *
                                precision.thicknessBlocks * precision.thicknessIncrements;
            long long totalIterations = (long long) opCount * numIterations;
            pancake.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            fieldArr = pancake.computeAllBFieldVectors(positions, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = double(totalIterations) / interval;
            pointsPerSec = opCount / interval;
            printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
            fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
        }
        printf("\n");
        fprintf(output, "\n");
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Pancake gradient performance:\n");
    for (int i = 1; i <= 8; ++i)
    {
        auto precision =
            PrecisionArguments::getCoilPrecisionArgumentsCPU(pancake, PrecisionFactor(i));
        int numIterations = precision.angularBlocks * precision.angularIncrements *
                            precision.thicknessBlocks * precision.thicknessIncrements;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        gradientArr = pancake.computeAllBGradientMatrices(positions, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = double(totalIterations) / interval;
        pointsPerSec = opCount / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= threadCount; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision =
                PrecisionArguments::getCoilPrecisionArgumentsCPU(pancake, PrecisionFactor(i));
            int numIterations = precision.angularBlocks * precision.angularIncrements *
                                precision.thicknessBlocks * precision.thicknessIncrements;
            long long totalIterations = (long long) opCount * numIterations;
            pancake.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            gradientArr = pancake.computeAllBGradientMatrices(positions, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = double(totalIterations) / interval;
            pointsPerSec = opCount / interval;
            printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
            fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
        }
        printf("\n");
        fprintf(output, "\n");
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Thin potential performance:\n");
    for(int i = 1; i <= 8; ++i)
    {
        auto precision =
            PrecisionArguments::getCoilPrecisionArgumentsCPU(thin, PrecisionFactor(i));
        int numIterations = precision.angularBlocks * precision.angularIncrements *
                            precision.thicknessBlocks * precision.thicknessIncrements;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        potentialArr = thin.computeAllAPotentialVectors(positions, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = double(totalIterations) / interval;
        pointsPerSec = opCount / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= threadCount; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision =
                PrecisionArguments::getCoilPrecisionArgumentsCPU(thin, PrecisionFactor(i));
            int numIterations = precision.angularBlocks * precision.angularIncrements *
                                precision.thicknessBlocks * precision.thicknessIncrements;
            long long totalIterations = (long long) opCount * numIterations;
            thin.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            potentialArr = thin.computeAllAPotentialVectors(positions, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = double(totalIterations)/ interval;
            pointsPerSec = opCount / interval;
            printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
            fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
        }
        printf("\n");
        fprintf(output, "\n");
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Thin field performance:\n");
    for (int i = 1; i <= 8; ++i)
    {
        auto precision =
            PrecisionArguments::getCoilPrecisionArgumentsCPU(thin, PrecisionFactor(i));
        int numIterations = precision.angularBlocks * precision.angularIncrements *
                            precision.thicknessBlocks * precision.thicknessIncrements;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        fieldArr = thin.computeAllBFieldVectors(positions, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = double(totalIterations) / interval;
        pointsPerSec = opCount / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= threadCount; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision =
                PrecisionArguments::getCoilPrecisionArgumentsCPU(thin, PrecisionFactor(i));
            int numIterations = precision.angularBlocks * precision.angularIncrements *
                                precision.thicknessBlocks * precision.thicknessIncrements;
            long long totalIterations = (long long) opCount * numIterations;
            thin.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            fieldArr = thin.computeAllBFieldVectors(positions, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = double(totalIterations) / interval;
            pointsPerSec = opCount / interval;
            printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
            fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
        }
        printf("\n");
        fprintf(output, "\n");
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Thin gradient performance:\n");
    for (int i = 1; i <= 8; ++i)
    {
        auto precision =
            PrecisionArguments::getCoilPrecisionArgumentsCPU(thin, PrecisionFactor(i));
        int numIterations = precision.angularBlocks * precision.angularIncrements *
                            precision.thicknessBlocks * precision.thicknessIncrements;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        gradientArr = thin.computeAllBGradientMatrices(positions, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = double(totalIterations) / interval;
        pointsPerSec = opCount / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= threadCount; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision =
                PrecisionArguments::getCoilPrecisionArgumentsCPU(thin, PrecisionFactor(i));
            int numIterations = precision.angularBlocks * precision.angularIncrements *
                                precision.thicknessBlocks * precision.thicknessIncrements;
            long long totalIterations = (long long) opCount * numIterations;
            thin.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            gradientArr = thin.computeAllBGradientMatrices(positions, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = double(totalIterations) / interval;
            pointsPerSec = opCount / interval;
            printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
            fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
        }
        printf("\n");
        fprintf(output, "\n");
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Thick potential performance:\n");
    for(int i = 1; i <= 8; ++i)
    {
        auto precision =
            PrecisionArguments::getCoilPrecisionArgumentsCPU(thick, PrecisionFactor(i));
        int numIterations = precision.angularBlocks * precision.angularIncrements *
                            precision.thicknessBlocks * precision.thicknessIncrements;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        potentialArr = thick.computeAllAPotentialVectors(positions, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = double(totalIterations) / interval;
        pointsPerSec = opCount / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= threadCount; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision =
                PrecisionArguments::getCoilPrecisionArgumentsCPU(thick, PrecisionFactor(i));
            int numIterations = precision.angularBlocks * precision.angularIncrements *
                                precision.thicknessBlocks * precision.thicknessIncrements;
            long long totalIterations = (long long) opCount * numIterations;
            thick.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            potentialArr = thick.computeAllAPotentialVectors(positions, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = double(totalIterations)/ interval;
            pointsPerSec = opCount / interval;
            printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
            fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
        }
        printf("\n");
        fprintf(output, "\n");
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Thick field performance:\n");
    for (int i = 1; i <= 8; ++i)
    {
        auto precision =
            PrecisionArguments::getCoilPrecisionArgumentsCPU(thick, PrecisionFactor(i));
        int numIterations = precision.angularBlocks * precision.angularIncrements *
                            precision.thicknessBlocks * precision.thicknessIncrements;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        fieldArr = thick.computeAllBFieldVectors(positions, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = double(totalIterations) / interval;
        pointsPerSec = opCount / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= threadCount; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision =
                PrecisionArguments::getCoilPrecisionArgumentsCPU(thick, PrecisionFactor(i));
            int numIterations = precision.angularBlocks * precision.angularIncrements *
                                precision.thicknessBlocks * precision.thicknessIncrements;
            long long totalIterations = (long long) opCount * numIterations;
            thick.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            fieldArr = thick.computeAllBFieldVectors(positions, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = double(totalIterations) / interval;
            pointsPerSec = opCount / interval;
            printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
            fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
        }
        printf("\n");
        fprintf(output, "\n");
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Thick gradient performance:\n");
    for (int i = 1; i <= 8; ++i)
    {
        auto precision =
            PrecisionArguments::getCoilPrecisionArgumentsCPU(thick, PrecisionFactor(i));
        int numIterations = precision.angularBlocks * precision.angularIncrements *
                            precision.thicknessBlocks * precision.thicknessIncrements;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        gradientArr = thick.computeAllBGradientMatrices(positions, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = double(totalIterations) / interval;
        pointsPerSec = opCount / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= threadCount; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision =
                PrecisionArguments::getCoilPrecisionArgumentsCPU(thick, PrecisionFactor(i));
            int numIterations = precision.angularBlocks * precision.angularIncrements *
                                precision.thicknessBlocks * precision.thicknessIncrements;
            long long totalIterations = (long long) opCount * numIterations;
            thick.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            gradientArr = thick.computeAllBGradientMatrices(positions, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = double(totalIterations) / interval;
            pointsPerSec = opCount / interval;
            printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
            fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
        }
        printf("\n");
        fprintf(output, "\n");
    }
    printf("\n");

    fclose(output);
}
