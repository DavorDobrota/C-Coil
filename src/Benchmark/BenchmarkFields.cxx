#include "Benchmark.h"
#include "Coil.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <chrono>


void benchComputeFieldsST(int opCount)
{
    using namespace std::chrono;

    vec3::Vector3 temp1;
    vec3::Vector3 temp2;
    vec3::Matrix3 temp3;

    high_resolution_clock::time_point begin_time;
    double interval;

    Coil testCoil = Coil(0.03, 0.03, 0.12, 3600);

    PrecisionArguments precision = testCoil.getPrecisionSettingsCPU();

    int numOperations = opCount *
                        precision.thicknessBlockCount * precision.thicknessIncrementCount *
                        precision.lengthBlockCount * precision.lengthIncrementCount *
                        precision.angularBlockCount * precision.angularIncrementCount;

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < opCount; ++i){
        temp1 = testCoil.computeAPotentialVector(vec3::Vector3(0.0, 0.0, i*0.000001));
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("potential A : %.1f MInc/s\n", 1e-6 * numOperations / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < opCount; ++i){
        temp2 = testCoil.computeBFieldVector(vec3::Vector3(0.0, 0.0, i*0.000001));
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("combined B  : %.1f MInc/s\n", 1e-6 * numOperations / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < opCount; ++i){
        temp3 = testCoil.computeBGradientMatrix(vec3::Vector3(0.0, 0.0, i * 0.000001));
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("gradient G  : %.1f MInc/s\n", 1e-6 * numOperations / interval);
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"
void benchComputeAllFields(PrecisionFactor precisionFactor, int opCount, int repeatCount, int threadCount)
{
    using namespace std::chrono;

    printf("Testing performance of compute all methods for field calculation\n");

    Coil testCoilFast = Coil(0.03, 0.03, 0.12, 3600, precisionFactor, threadCount);
    Coil testCoilSlow = Coil(0.03, 0.03, 0.0, 30, precisionFactor, threadCount);

    const double radius = 0.1;

    PrecisionArguments precisionSlowCPU = testCoilSlow.getPrecisionSettingsCPU();
    PrecisionArguments precisionFastCPU = testCoilFast.getPrecisionSettingsCPU();
    PrecisionArguments precisionSlowGPU = testCoilSlow.getPrecisionSettingsGPU();
    PrecisionArguments precisionFastGPU = testCoilFast.getPrecisionSettingsGPU();

    const long long numOperationsSlowCPU = (long long) opCount *
                                           precisionSlowCPU.thicknessBlockCount * precisionSlowCPU.thicknessIncrementCount *
                                           precisionSlowCPU.angularBlockCount * precisionSlowCPU.angularIncrementCount;

    const long long numOperationsFastCPU1 = (long long) opCount *
                                            precisionFastCPU.thicknessBlockCount * precisionFastCPU.thicknessIncrementCount *
                                            precisionFastCPU.angularBlockCount * precisionFastCPU.angularIncrementCount;
    const long long numOperationsFastCPU2 = numOperationsFastCPU1 * precisionFastCPU.lengthBlockCount * precisionFastCPU.lengthIncrementCount;

    const long long numOperationsSlowGPU = (long long) opCount *
                                           precisionSlowGPU.thicknessIncrementCount * precisionSlowGPU.angularIncrementCount;

    const long long numOperationsFastGPU1 = opCount * precisionFastGPU.thicknessIncrementCount * precisionFastGPU.angularIncrementCount;
    const long long numOperationsFastGPU2 = numOperationsFastGPU1 * precisionFastGPU.lengthIncrementCount;

    vec3::Vector3Array positionValues(opCount);

    for (int i = 0; i < opCount; i++)
        positionValues[i] = vec3::Vector3::getFromCylindricalCoords(0.1, i * M_PI / opCount, 0.0);

    vec3::Vector3Array potentialVectors;

    vec3::Vector3Array fieldVectors;

    vec3::Matrix3Array gradientMatrices;

    high_resolution_clock::time_point begin_time;
    double interval;

    printf("\n");

    // ST methods slow - repetitions removed for single core, it is just too slow compared to the rest
    begin_time = high_resolution_clock::now();
    potentialVectors = testCoilSlow.computeAllAPotentialVectors(positionValues, CPU_ST);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_ST slow : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * numOperationsSlowCPU / interval, 1e-6 * numOperationsSlowCPU / interval, opCount / interval);

    begin_time = high_resolution_clock::now();
    fieldVectors = testCoilSlow.computeAllBFieldVectors(positionValues, CPU_ST);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Field     B CPU_ST slow : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * numOperationsSlowCPU / interval, 1e-6 * numOperationsSlowCPU / interval, opCount / interval);

    begin_time = high_resolution_clock::now();
    gradientMatrices = testCoilSlow.computeAllBGradientMatrices(positionValues, CPU_ST);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G CPU_ST slow : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * numOperationsSlowCPU / interval, 1e-6 * numOperationsSlowCPU / interval, opCount / interval);

    printf("\n");
    // ST methods fast - repetitions removed for single core, it is just too slow compared to the rest
    begin_time = high_resolution_clock::now();
    potentialVectors = testCoilFast.computeAllAPotentialVectors(positionValues, CPU_ST);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_ST fast : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * numOperationsFastCPU1 / interval, 1e-6 * numOperationsFastCPU2 / interval, opCount / interval);

    begin_time = high_resolution_clock::now();
    fieldVectors = testCoilFast.computeAllBFieldVectors(positionValues, CPU_ST);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Field     B CPU_ST fast : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * numOperationsFastCPU1 / interval, 1e-6 * numOperationsFastCPU2 / interval, opCount / interval);

    begin_time = high_resolution_clock::now();
    gradientMatrices = testCoilFast.computeAllBGradientMatrices(positionValues, CPU_ST);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G CPU_ST fast : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * numOperationsFastCPU1 / interval, 1e-6 * numOperationsFastCPU2 / interval, opCount / interval);

    printf("\n");

    // MT methods slow
    begin_time = high_resolution_clock::now();
    for (int i = 0; i < repeatCount; i++)
        potentialVectors = testCoilSlow.computeAllAPotentialVectors(positionValues, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_MT slow : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsSlowCPU * repeatCount) / interval,
           1e-6 * (numOperationsSlowCPU * repeatCount) / interval,
           (opCount * repeatCount) / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < repeatCount; i++)
        fieldVectors = testCoilSlow.computeAllBFieldVectors(positionValues, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Field     B CPU_MT slow : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsSlowCPU * repeatCount) / interval,
           1e-6 * (numOperationsSlowCPU * repeatCount) / interval,
           (opCount * repeatCount) / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < repeatCount; i++)
        gradientMatrices = testCoilSlow.computeAllBGradientMatrices(positionValues, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G CPU_MT slow : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsSlowCPU * repeatCount) / interval,
           1e-6 * (numOperationsSlowCPU * repeatCount) / interval,
           (opCount * repeatCount) / interval);

    printf("\n");
    // MT methods fast
    begin_time = high_resolution_clock::now();
    for (int i = 0; i < repeatCount; i++)
        potentialVectors = testCoilFast.computeAllAPotentialVectors(positionValues, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_MT fast : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsFastCPU1 * repeatCount) / interval,
           1e-6 * (numOperationsFastCPU2 * repeatCount) / interval,
           (opCount * repeatCount) / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < repeatCount; i++)
        fieldVectors = testCoilFast.computeAllBFieldVectors(positionValues, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Field     B CPU_MT fast : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsFastCPU1 * repeatCount) / interval,
           1e-6 * (numOperationsFastCPU2 * repeatCount) / interval,
           (opCount * repeatCount) / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < repeatCount; i++)
        gradientMatrices = testCoilFast.computeAllBGradientMatrices(positionValues, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G CPU_MT fast : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsFastCPU1 * repeatCount) / interval,
           1e-6 * (numOperationsFastCPU2 * repeatCount) / interval,
           (opCount * repeatCount) / interval);

    printf("\n");

    potentialVectors = testCoilSlow.computeAllAPotentialVectors(positionValues, GPU); // warmup for the GPU

    // GPU methods slow
    begin_time = high_resolution_clock::now();
    for (int i = 0; i < repeatCount; i++)
        potentialVectors = testCoilSlow.computeAllAPotentialVectors(positionValues, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A GPU    slow : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsSlowGPU * repeatCount) / interval,
           1e-6 * (numOperationsSlowGPU * repeatCount) / interval,
           (opCount * repeatCount) / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < repeatCount; i++)
        fieldVectors = testCoilSlow.computeAllBFieldVectors(positionValues, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Field     B GPU    slow : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsSlowGPU * repeatCount) / interval,
           1e-6 * (numOperationsSlowGPU * repeatCount) / interval,
           (opCount * repeatCount) / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < repeatCount; i++)
        gradientMatrices = testCoilSlow.computeAllBGradientMatrices(positionValues, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G GPU    slow : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsSlowGPU * repeatCount) / interval,
           1e-6 * (numOperationsSlowGPU * repeatCount) / interval,
           (opCount * repeatCount) / interval);

    printf("\n");
    // GPU methods fast

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < repeatCount; i++)
        potentialVectors = testCoilFast.computeAllAPotentialVectors(positionValues, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A GPU    fast : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsFastGPU1 * repeatCount) / interval,
           1e-6 * (numOperationsFastGPU2 * repeatCount) / interval,
           (opCount * repeatCount) / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < repeatCount; i++)
        fieldVectors = testCoilFast.computeAllBFieldVectors(positionValues, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Field     B GPU    fast : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsFastGPU1 * repeatCount) / interval,
           1e-6 * (numOperationsFastGPU2 * repeatCount) / interval,
           (opCount * repeatCount) / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < repeatCount; i++)
        gradientMatrices = testCoilFast.computeAllBGradientMatrices(positionValues, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G GPU    fast : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsFastGPU1 * repeatCount) / interval,
           1e-6 * (numOperationsFastGPU2 * repeatCount) / interval,
           (opCount * repeatCount) / interval);

    printf("\n");
}
#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"
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
    vec3::Vector3Array positionValues(opCount);

    for (int i = 0; i < opCount; ++i)
        positionValues[i] = vec3::Vector3::getFromSphericalCoords(1.0, M_PI * i / opCount, 0.0);

    vec3::Vector3Array potentialArr;
    vec3::Vector3Array fieldArr;
    vec3::Matrix3Array gradientArr;

    printf("This test is created for the purpose of generating performance charts\n\n");

    printf("Loop potential performance:\n");
    for(int i = 1; i <= 8; ++i)
    {
        auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(loop, PrecisionFactor(i));
        int numIterations = precision.angularBlockCount * precision.angularIncrementCount;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        potentialArr = loop.computeAllAPotentialVectors(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
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
            int numIterations = precision.angularBlockCount * precision.angularIncrementCount;
            long long totalIterations = (long long) opCount * numIterations;
            loop.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            potentialArr = loop.computeAllAPotentialVectors(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations / interval;
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
        auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(loop, PrecisionFactor(i));
        int numIterations = precision.angularBlockCount * precision.angularIncrementCount;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        fieldArr = loop.computeAllBFieldVectors(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
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
            int numIterations = precision.angularBlockCount * precision.angularIncrementCount;
            long long totalIterations = (long long) opCount * numIterations;
            loop.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            fieldArr = loop.computeAllBFieldVectors(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations / interval;
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
        auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(loop, PrecisionFactor(i));
        int numIterations = precision.angularBlockCount * precision.angularIncrementCount;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        gradientArr = loop.computeAllBGradientMatrices(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
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
            int numIterations = precision.angularBlockCount * precision.angularIncrementCount;
            long long totalIterations = (long long) opCount * numIterations;
            loop.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            gradientArr = loop.computeAllBGradientMatrices(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations / interval;
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
        auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(pancake, PrecisionFactor(i));
        int numIterations =
                precision.angularBlockCount * precision.angularIncrementCount
                * precision.thicknessBlockCount * precision.thicknessIncrementCount;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        potentialArr = pancake.computeAllAPotentialVectors(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
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
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(pancake, PrecisionFactor(i));
            int numIterations =
                    precision.angularBlockCount * precision.angularIncrementCount
                    * precision.thicknessBlockCount * precision.thicknessIncrementCount;
            long long totalIterations = (long long) opCount * numIterations;
            pancake.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            potentialArr = pancake.computeAllAPotentialVectors(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations/ interval;
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
        auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(pancake, PrecisionFactor(i));
        int numIterations =
                precision.angularBlockCount * precision.angularIncrementCount
                * precision.thicknessBlockCount * precision.thicknessIncrementCount;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        fieldArr = pancake.computeAllBFieldVectors(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
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
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(pancake, PrecisionFactor(i));
            int numIterations =
                    precision.angularBlockCount * precision.angularIncrementCount
                    * precision.thicknessBlockCount * precision.thicknessIncrementCount;
            long long totalIterations = (long long) opCount * numIterations;
            pancake.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            fieldArr = pancake.computeAllBFieldVectors(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations / interval;
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
        auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(pancake, PrecisionFactor(i));
        int numIterations =
                precision.angularBlockCount * precision.angularIncrementCount
                * precision.thicknessBlockCount * precision.thicknessIncrementCount;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        gradientArr = pancake.computeAllBGradientMatrices(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
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
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(pancake, PrecisionFactor(i));
            int numIterations =
                    precision.angularBlockCount * precision.angularIncrementCount
                    * precision.thicknessBlockCount * precision.thicknessIncrementCount;
            long long totalIterations = (long long) opCount * numIterations;
            pancake.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            gradientArr = pancake.computeAllBGradientMatrices(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations / interval;
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
        auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(thin, PrecisionFactor(i));
        int numIterations =
                precision.angularBlockCount * precision.angularIncrementCount
                * precision.thicknessBlockCount * precision.thicknessIncrementCount;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        potentialArr = thin.computeAllAPotentialVectors(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
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
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(thin, PrecisionFactor(i));
            int numIterations =
                    precision.angularBlockCount * precision.angularIncrementCount
                    * precision.thicknessBlockCount * precision.thicknessIncrementCount;
            long long totalIterations = (long long) opCount * numIterations;
            thin.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            potentialArr = thin.computeAllAPotentialVectors(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations/ interval;
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
        auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(thin, PrecisionFactor(i));
        int numIterations =
                precision.angularBlockCount * precision.angularIncrementCount
                * precision.thicknessBlockCount * precision.thicknessIncrementCount;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        fieldArr = thin.computeAllBFieldVectors(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
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
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(thin, PrecisionFactor(i));
            int numIterations =
                    precision.angularBlockCount * precision.angularIncrementCount
                    * precision.thicknessBlockCount * precision.thicknessIncrementCount;
            long long totalIterations = (long long) opCount * numIterations;
            thin.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            fieldArr = thin.computeAllBFieldVectors(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations / interval;
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
        auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(thin, PrecisionFactor(i));
        int numIterations =
                precision.angularBlockCount * precision.angularIncrementCount
                * precision.thicknessBlockCount * precision.thicknessIncrementCount;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        gradientArr = thin.computeAllBGradientMatrices(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
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
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(thin, PrecisionFactor(i));
            int numIterations =
                    precision.angularBlockCount * precision.angularIncrementCount
                    * precision.thicknessBlockCount * precision.thicknessIncrementCount;
            long long totalIterations = (long long) opCount * numIterations;
            thin.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            gradientArr = thin.computeAllBGradientMatrices(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations / interval;
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
        auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(thick, PrecisionFactor(i));
        int numIterations =
                precision.angularBlockCount * precision.angularIncrementCount
                * precision.thicknessBlockCount * precision.thicknessIncrementCount;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        potentialArr = thick.computeAllAPotentialVectors(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
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
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(thick, PrecisionFactor(i));
            int numIterations =
                    precision.angularBlockCount * precision.angularIncrementCount
                    * precision.thicknessBlockCount * precision.thicknessIncrementCount;
            long long totalIterations = (long long) opCount * numIterations;
            thick.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            potentialArr = thick.computeAllAPotentialVectors(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations/ interval;
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
        auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(thick, PrecisionFactor(i));
        int numIterations =
                precision.angularBlockCount * precision.angularIncrementCount
                * precision.thicknessBlockCount * precision.thicknessIncrementCount;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        fieldArr = thick.computeAllBFieldVectors(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
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
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(thick, PrecisionFactor(i));
            int numIterations =
                    precision.angularBlockCount * precision.angularIncrementCount
                    * precision.thicknessBlockCount * precision.thicknessIncrementCount;
            long long totalIterations = (long long) opCount * numIterations;
            thick.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            fieldArr = thick.computeAllBFieldVectors(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations / interval;
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
        auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(thick, PrecisionFactor(i));
        int numIterations =
                precision.angularBlockCount * precision.angularIncrementCount
                * precision.thicknessBlockCount * precision.thicknessIncrementCount;
        long long totalIterations = (long long) opCount * numIterations;

        beginTime = high_resolution_clock::now();
        gradientArr = thick.computeAllBGradientMatrices(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
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
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(thick, PrecisionFactor(i));
            int numIterations =
                    precision.angularBlockCount * precision.angularIncrementCount
                    * precision.thicknessBlockCount * precision.thicknessIncrementCount;
            long long totalIterations = (long long) opCount * numIterations;
            thick.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            gradientArr = thick.computeAllBGradientMatrices(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations / interval;
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
#pragma clang diagnostic pop

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
        int numPoints = int(std::pow(2, i));

        vec3::Vector3Array positions(numPoints);
        vec3::Vector3Array results;
        for (int j = 0; j < numPoints; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        results = flat.computeAllAPotentialVectors(positions, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = numPoints / interval;

        printf("%8d : %.1f\n", numPoints, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", numPoints, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic field performance slow for precision factor %.1f and %d threads\n",
           precisionFactor.relativePrecision, threadCount);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int numPoints = int(std::pow(2, i));

        vec3::Vector3Array positions(numPoints);
        vec3::Vector3Array results;
        for (int j = 0; j < numPoints; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        results = flat.computeAllBFieldVectors(positions, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = numPoints / interval;

        printf("%8d : %.1f\n", numPoints, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", numPoints, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic gradient performance slow for precision factor %.1f and %d threads\n",
           precisionFactor.relativePrecision, threadCount);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int numPoints = int(std::pow(2, i));

        vec3::Vector3Array positions(numPoints);
        vec3::Matrix3Array results;
        for (int j = 0; j < numPoints; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        results = flat.computeAllBGradientMatrices(positions, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = numPoints / interval;

        printf("%8d : %.1f\n", numPoints, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", numPoints, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Vector potential performance fast for precision factor %.1f and %d threads\n",
           precisionFactor.relativePrecision, threadCount);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int numPoints = int(std::pow(2, i));

        vec3::Vector3Array positions(numPoints);
        vec3::Vector3Array results;
        for (int j = 0; j < numPoints; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        results = coil.computeAllAPotentialVectors(positions, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = numPoints / interval;

        printf("%8d : %.1f\n", numPoints, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", numPoints, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic field performance fast for precision factor %.1f and %d threads\n",
           precisionFactor.relativePrecision, threadCount);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int numPoints = int(std::pow(2, i));

        vec3::Vector3Array positions(numPoints);
        vec3::Vector3Array results;
        for (int j = 0; j < numPoints; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        results = coil.computeAllBFieldVectors(positions, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = numPoints / interval;

        printf("%8d : %.1f\n", numPoints, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", numPoints, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic gradient performance fast for precision factor %.1f and %d threads\n",
           precisionFactor.relativePrecision, threadCount);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int numPoints = int(std::pow(2, i));

        vec3::Vector3Array positions(numPoints);
        vec3::Matrix3Array results;
        for (int j = 0; j < numPoints; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        results = coil.computeAllBGradientMatrices(positions, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = numPoints / interval;

        printf("%8d : %.1f\n", numPoints, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", numPoints, 0.001 * pointsPerSec);
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
        int numPoints = int(std::pow(2, i));

        vec3::Vector3Array positions(numPoints);
        vec3::Vector3Array potentialArr;
        for (int j = 0; j < numPoints; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        potentialArr = flat.computeAllAPotentialVectors(positions, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = numPoints / interval;

        printf("%8d : %.1f\n", numPoints, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", numPoints, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic field performance slow for precision factor %.1f\n",
           precisionFactor.relativePrecision);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int numPoints = int(std::pow(2, i));

        vec3::Vector3Array positions(numPoints);
        vec3::Vector3Array fieldArr;
        for (int j = 0; j < numPoints; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        fieldArr = flat.computeAllBFieldVectors(positions, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = numPoints / interval;

        printf("%8d : %.1f\n", numPoints, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", numPoints, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic gradient performance slow for precision factor %.1f\n",
           precisionFactor.relativePrecision);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int numPoints = int(std::pow(2, i));

        vec3::Vector3Array positions(numPoints);
        vec3::Matrix3Array gradientArr;
        for (int j = 0; j < numPoints; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        gradientArr = flat.computeAllBGradientMatrices(positions, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = numPoints / interval;

        printf("%8d : %.1f\n", numPoints, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", numPoints, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Vector potential performance fast for precision factor %.1f\n",
           precisionFactor.relativePrecision);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int numPoints = int(std::pow(2, i));

        vec3::Vector3Array positions(numPoints);
        vec3::Vector3Array potentialArr;
        for (int j = 0; j < numPoints; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        potentialArr = coil.computeAllAPotentialVectors(positions, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = numPoints / interval;

        printf("%8d : %.1f\n", numPoints, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", numPoints, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic field performance fast for precision factor %.1f\n",
           precisionFactor.relativePrecision);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int numPoints = int(std::pow(2, i));

        vec3::Vector3Array positions(numPoints);
        vec3::Vector3Array fieldArr;
        for (int j = 0; j < numPoints; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        fieldArr = coil.computeAllBFieldVectors(positions, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = numPoints / interval;

        printf("%8d : %.1f\n", numPoints, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", numPoints, 0.001 * pointsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic gradient performance fast for precision factor %.1f\n",
           precisionFactor.relativePrecision);

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int numPoints = int(std::pow(2, i));

        vec3::Vector3Array positions(numPoints);
        vec3::Matrix3Array gradientArr;
        for (int j = 0; j < numPoints; ++j)
            positions[j] = vec3::Vector3(0.1, 0.1, double(j));

        beginTime = high_resolution_clock::now();
        gradientArr = coil.computeAllBGradientMatrices(positions, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = numPoints / interval;

        printf("%8d : %.1f\n", numPoints, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", numPoints, 0.001 * pointsPerSec);
    }
    printf("\n");

    fclose(output);
}