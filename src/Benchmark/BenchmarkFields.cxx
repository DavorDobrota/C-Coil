#include "Benchmark.h"
#include "Coil.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <chrono>


void benchmarkCPU_ST(int nOps)
{
    using namespace std::chrono;

    double temp1;
    vec3::FieldVector3 temp2;
    vec3::Matrix3 temp3;

    high_resolution_clock::time_point begin_time;
    double interval;

    Coil testCoil = Coil(0.03, 0.03, 0.12, 3600);

    PrecisionArguments precision = testCoil.getPrecisionSettings();

    int numOperations = nOps *
                        precision.thicknessBlockCount * precision.thicknessIncrementCount *
                        precision.lengthBlockCount * precision.lengthIncrementCount *
                        precision.angularBlockCount * precision.angularIncrementCount;

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nOps; ++i){
        temp1 = testCoil.computeAPotentialAbs(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, i*0.000001));
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("potential A : %.1f MInc/s\n", 1e-6 * numOperations / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nOps; ++i){
        temp2 = testCoil.computeBFieldVector(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, i*0.000001));
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("combined B  : %.1f MInc/s\n", 1e-6 * numOperations / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nOps; ++i){
        temp3 = testCoil.computeBGradientTensor(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, i*0.000001));
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("gradient G  : %.1f MInc/s\n", 1e-6 * numOperations / interval);
}

void benchmarkComputeAll(PrecisionFactor precisionFactor, int nOps, int nRepeats, int nThreads)
{
    using namespace std::chrono;

    Coil testCoilFast = Coil(0.03, 0.03, 0.12, 3600, precisionFactor, nThreads);
    Coil testCoilSlow = Coil(0.03, 0.03, 0.0, 30, precisionFactor, nThreads);

    const double radius = 0.1;

    PrecisionArguments precisionSlow = testCoilSlow.getPrecisionSettings();
    PrecisionArguments precisionFast = testCoilFast.getPrecisionSettings();

    const long long numOperationsSlow = (long long) nOps *
                                        precisionSlow.thicknessBlockCount * precisionSlow.thicknessIncrementCount *
                                        precisionSlow.angularBlockCount * precisionSlow.angularIncrementCount;

    const long long numOperationsFast1 = (long long) nOps *
                                         precisionFast.thicknessBlockCount * precisionFast.thicknessIncrementCount *
                                         precisionFast.angularBlockCount * precisionFast.angularIncrementCount;

    const long long numOperationsFast2 = numOperationsFast1 * precisionFast.lengthBlockCount * precisionFast.lengthIncrementCount;

    const long long numOperationsGpu = nOps * 48 * 16 * 16;

    std::vector<vec3::CoordVector3> positionValues(nOps);

    for (int i = 0; i < nOps; i++)
        positionValues[i] = vec3::CoordVector3(vec3::CYLINDRICAL, 0.1, i * M_PI / nOps, 0.0);

    std::vector<double> cpuPotential;
    std::vector<double> gpuPotential;

    std::vector<vec3::FieldVector3> cpuFieldVectors;
    std::vector<vec3::FieldVector3> gpuFieldVectors;

    std::vector<vec3::Matrix3> cpuGradientMatrices;

    high_resolution_clock::time_point begin_time;
    double interval;

    // ST methods slow - repetitions removed for single core, it is just too slow compared to the rest
    begin_time = high_resolution_clock::now();
    cpuPotential = testCoilSlow.computeAllAPotentialAbs(positionValues, CPU_ST);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_ST slow : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * numOperationsSlow / interval, 1e-6 * numOperationsSlow / interval, nOps / interval);

    begin_time = high_resolution_clock::now();
    cpuFieldVectors = testCoilSlow.computeAllBFieldComponents(positionValues, CPU_ST);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("combined  B CPU_ST slow : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * numOperationsSlow / interval, 1e-6 * numOperationsSlow / interval, nOps / interval);

    begin_time = high_resolution_clock::now();
    cpuGradientMatrices = testCoilSlow.computeAllBGradientTensors(positionValues, CPU_ST);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G CPU_ST slow : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * numOperationsSlow / interval, 1e-6 * numOperationsSlow / interval, nOps / interval);


    // ST methods fast - repetitions removed for single core, it is just too slow compared to the rest
    begin_time = high_resolution_clock::now();
    cpuPotential = testCoilFast.computeAllAPotentialAbs(positionValues, CPU_ST);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_ST fast : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * numOperationsFast1 / interval, 1e-6 * numOperationsFast2 / interval, nOps / interval);

    begin_time = high_resolution_clock::now();
    cpuFieldVectors = testCoilFast.computeAllBFieldComponents(positionValues, CPU_ST);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("combined  B CPU_ST fast : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * numOperationsFast1 / interval, 1e-6 * numOperationsFast2 / interval, nOps / interval);

    begin_time = high_resolution_clock::now();
    cpuGradientMatrices = testCoilFast.computeAllBGradientTensors(positionValues, CPU_ST);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G CPU_ST fast : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * numOperationsFast1 / interval, 1e-6 * numOperationsFast2 / interval, nOps / interval);

    // MT methods slow
    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nRepeats; i++)
        cpuPotential = testCoilSlow.computeAllAPotentialAbs(positionValues, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_MT slow : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsSlow * nRepeats) / interval, 1e-6 * (numOperationsSlow * nRepeats) / interval, (nOps * nRepeats) / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nRepeats; i++)
        cpuFieldVectors = testCoilSlow.computeAllBFieldComponents(positionValues, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("combined  B CPU_MT slow : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsSlow * nRepeats) / interval, 1e-6 * (numOperationsSlow * nRepeats) / interval, (nOps * nRepeats) / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nRepeats; i++)
        cpuGradientMatrices = testCoilSlow.computeAllBGradientTensors(positionValues, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G CPU_MT slow : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsSlow * nRepeats) / interval, 1e-6 * (numOperationsSlow * nRepeats) / interval, (nOps * nRepeats) / interval);

    // MT methods fast
    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nRepeats; i++)
        cpuPotential = testCoilFast.computeAllAPotentialAbs(positionValues, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_MT fast : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsFast1 * nRepeats) / interval, 1e-6 * (numOperationsFast2 * nRepeats) / interval, (nOps * nRepeats) / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nRepeats; i++)
        cpuFieldVectors = testCoilFast.computeAllBFieldComponents(positionValues, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("combined  B CPU_MT fast : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsFast1 * nRepeats) / interval, 1e-6 * (numOperationsFast2 * nRepeats) / interval, (nOps * nRepeats) / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nRepeats; i++)
        cpuGradientMatrices = testCoilFast.computeAllBGradientTensors(positionValues, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G CPU_MT fast : %.1f MInc/s | eff: %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsFast1 * nRepeats) / interval, 1e-6 * (numOperationsFast2 * nRepeats) / interval, (nOps * nRepeats) / interval);

    // GPU methods
    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nRepeats; i++)
        gpuPotential = testCoilFast.computeAllAPotentialAbs(positionValues, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A GPU : %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsGpu * nRepeats) / interval, (nOps * nRepeats) / interval);


    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nRepeats; i++)
        gpuFieldVectors = testCoilFast.computeAllBFieldComponents(positionValues, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("combined  B GPU : %.1f MInc/s | %.0f points/s\n",
           1e-6 * (numOperationsGpu * nRepeats) / interval, (nOps * nRepeats) / interval);

    printf("\n");
}

void benchmarkVariousCoilTypes(int nOps, int nThreads)
{
    using namespace std::chrono;

    FILE *output = fopen("output.txt", "w");

    Coil loop = Coil(0.1, 0.0, 0.0, 1);
    Coil pancake = Coil(0.1, 0.4, 0.0, 40);
    Coil thin = Coil(0.1, 0.0, 0.1, 100);
    Coil thick = Coil(0.1, 0.4, 0.4, 1600);

    PrecisionArguments thinPrecision = pancake.getPrecisionSettings();
    PrecisionArguments thickPrecision = pancake.getPrecisionSettings();

    int thinNumIterations = thinPrecision.angularBlockCount * thinPrecision.angularIncrementCount;
    int thickNumIterations = thickPrecision.angularBlockCount * thickPrecision.angularIncrementCount
                             * thickPrecision.thicknessBlockCount * thickPrecision.thicknessIncrementCount;


    high_resolution_clock::time_point beginTime;
    double interval, incrementsPerSec, pointsPerSec;
    std::vector<vec3::CoordVector3> positionValues(nOps);

    for (int i = 0; i < nOps; ++i)
        positionValues[i] = vec3::CoordVector3(vec3::SPHERICAL, 1.0, M_PI * i / nOps, 0.0);

    std::vector<vec3::FieldVector3> potentialArr;
    std::vector<vec3::FieldVector3> fieldArr;
    std::vector<vec3::Matrix3> gradientArr;

    printf("This test is created for the purpose of generating performance charts\n\n");

    printf("Loop potential performance:\n");
    for(int i = 1; i <= 8; ++i)
    {
        auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(loop, PrecisionFactor(i));
        int numIterations = precision.angularBlockCount * precision.angularIncrementCount;
        long long totalIterations = (long long) nOps * numIterations;

        beginTime = high_resolution_clock::now();
        potentialArr = loop.computeAllAPotentialComponents(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
        pointsPerSec = nOps / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= nThreads; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(loop, PrecisionFactor(i));
            int numIterations = precision.angularBlockCount * precision.angularIncrementCount;
            long long totalIterations = (long long) nOps * numIterations;
            loop.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            potentialArr = loop.computeAllAPotentialComponents(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations / interval;
            pointsPerSec = nOps / interval;
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
        long long totalIterations = (long long) nOps * numIterations;

        beginTime = high_resolution_clock::now();
        fieldArr = loop.computeAllBFieldComponents(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
        pointsPerSec = nOps / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= nThreads; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(loop, PrecisionFactor(i));
            int numIterations = precision.angularBlockCount * precision.angularIncrementCount;
            long long totalIterations = (long long) nOps * numIterations;
            loop.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            fieldArr = loop.computeAllBFieldComponents(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations / interval;
            pointsPerSec = nOps / interval;
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
        long long totalIterations = (long long) nOps * numIterations;

        beginTime = high_resolution_clock::now();
        gradientArr = loop.computeAllBGradientTensors(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
        pointsPerSec = nOps / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= nThreads; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(loop, PrecisionFactor(i));
            int numIterations = precision.angularBlockCount * precision.angularIncrementCount;
            long long totalIterations = (long long) nOps * numIterations;
            loop.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            gradientArr = loop.computeAllBGradientTensors(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations / interval;
            pointsPerSec = nOps / interval;
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
        long long totalIterations = (long long) nOps * numIterations;

        beginTime = high_resolution_clock::now();
        potentialArr = pancake.computeAllAPotentialComponents(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
        pointsPerSec = nOps / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= nThreads; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(pancake, PrecisionFactor(i));
            int numIterations =
                    precision.angularBlockCount * precision.angularIncrementCount
                    * precision.thicknessBlockCount * precision.thicknessIncrementCount;
            long long totalIterations = (long long) nOps * numIterations;
            pancake.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            potentialArr = pancake.computeAllAPotentialComponents(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations/ interval;
            pointsPerSec = nOps / interval;
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
        long long totalIterations = (long long) nOps * numIterations;

        beginTime = high_resolution_clock::now();
        fieldArr = pancake.computeAllBFieldComponents(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
        pointsPerSec = nOps / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= nThreads; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(pancake, PrecisionFactor(i));
            int numIterations =
                    precision.angularBlockCount * precision.angularIncrementCount
                    * precision.thicknessBlockCount * precision.thicknessIncrementCount;
            long long totalIterations = (long long) nOps * numIterations;
            pancake.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            fieldArr = pancake.computeAllBFieldComponents(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations / interval;
            pointsPerSec = nOps / interval;
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
        long long totalIterations = (long long) nOps * numIterations;

        beginTime = high_resolution_clock::now();
        gradientArr = pancake.computeAllBGradientTensors(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
        pointsPerSec = nOps / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= nThreads; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(pancake, PrecisionFactor(i));
            int numIterations =
                    precision.angularBlockCount * precision.angularIncrementCount
                    * precision.thicknessBlockCount * precision.thicknessIncrementCount;
            long long totalIterations = (long long) nOps * numIterations;
            pancake.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            gradientArr = pancake.computeAllBGradientTensors(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations / interval;
            pointsPerSec = nOps / interval;
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
        long long totalIterations = (long long) nOps * numIterations;

        beginTime = high_resolution_clock::now();
        potentialArr = thin.computeAllAPotentialComponents(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
        pointsPerSec = nOps / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= nThreads; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(thin, PrecisionFactor(i));
            int numIterations =
                    precision.angularBlockCount * precision.angularIncrementCount
                    * precision.thicknessBlockCount * precision.thicknessIncrementCount;
            long long totalIterations = (long long) nOps * numIterations;
            thin.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            potentialArr = thin.computeAllAPotentialComponents(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations/ interval;
            pointsPerSec = nOps / interval;
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
        long long totalIterations = (long long) nOps * numIterations;

        beginTime = high_resolution_clock::now();
        fieldArr = thin.computeAllBFieldComponents(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
        pointsPerSec = nOps / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= nThreads; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(thin, PrecisionFactor(i));
            int numIterations =
                    precision.angularBlockCount * precision.angularIncrementCount
                    * precision.thicknessBlockCount * precision.thicknessIncrementCount;
            long long totalIterations = (long long) nOps * numIterations;
            thin.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            fieldArr = thin.computeAllBFieldComponents(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations / interval;
            pointsPerSec = nOps / interval;
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
        long long totalIterations = (long long) nOps * numIterations;

        beginTime = high_resolution_clock::now();
        gradientArr = thin.computeAllBGradientTensors(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
        pointsPerSec = nOps / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= nThreads; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(thin, PrecisionFactor(i));
            int numIterations =
                    precision.angularBlockCount * precision.angularIncrementCount
                    * precision.thicknessBlockCount * precision.thicknessIncrementCount;
            long long totalIterations = (long long) nOps * numIterations;
            thin.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            gradientArr = thin.computeAllBGradientTensors(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations / interval;
            pointsPerSec = nOps / interval;
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
        long long totalIterations = (long long) nOps * numIterations;

        beginTime = high_resolution_clock::now();
        potentialArr = thick.computeAllAPotentialComponents(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
        pointsPerSec = nOps / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= nThreads; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(thick, PrecisionFactor(i));
            int numIterations =
                    precision.angularBlockCount * precision.angularIncrementCount
                    * precision.thicknessBlockCount * precision.thicknessIncrementCount;
            long long totalIterations = (long long) nOps * numIterations;
            thick.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            potentialArr = thick.computeAllAPotentialComponents(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations/ interval;
            pointsPerSec = nOps / interval;
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
        long long totalIterations = (long long) nOps * numIterations;

        beginTime = high_resolution_clock::now();
        fieldArr = thick.computeAllBFieldComponents(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
        pointsPerSec = nOps / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= nThreads; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(thick, PrecisionFactor(i));
            int numIterations =
                    precision.angularBlockCount * precision.angularIncrementCount
                    * precision.thicknessBlockCount * precision.thicknessIncrementCount;
            long long totalIterations = (long long) nOps * numIterations;
            thick.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            fieldArr = thick.computeAllBFieldComponents(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations / interval;
            pointsPerSec = nOps / interval;
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
        long long totalIterations = (long long) nOps * numIterations;

        beginTime = high_resolution_clock::now();
        gradientArr = thick.computeAllBGradientTensors(positionValues, precision, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

        incrementsPerSec = totalIterations / interval;
        pointsPerSec = nOps / interval;
        printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
        fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
    }
    printf("\n");
    fprintf(output, "\n");

    for (int t = 2; t <= nThreads; ++t)
    {
        for (int i = 1; i <= 8; ++i)
        {
            auto precision = PrecisionArguments::getCoilPrecisionArgumentsCPU(thick, PrecisionFactor(i));
            int numIterations =
                    precision.angularBlockCount * precision.angularIncrementCount
                    * precision.thicknessBlockCount * precision.thicknessIncrementCount;
            long long totalIterations = (long long) nOps * numIterations;
            thick.setThreadCount(t);

            beginTime = high_resolution_clock::now();
            gradientArr = thick.computeAllBGradientTensors(positionValues, precision, CPU_MT);
            interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();

            incrementsPerSec = totalIterations / interval;
            pointsPerSec = nOps / interval;
            printf("%.1f MInc/s | %.0f kPoints/s\n", 1e-6 * incrementsPerSec, 1e-3 * pointsPerSec);
            fprintf(output, "%.3f\t", 1e-6 * incrementsPerSec);
        }
        printf("\n");
        fprintf(output, "\n");
    }
    printf("\n");
    fprintf(output, "\n");

    fclose(output);
}