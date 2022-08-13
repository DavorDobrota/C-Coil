#include "Benchmark.h"

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
                        precision.thicknessBlocks * precision.thicknessIncrements *
                        precision.lengthBlocks * precision.lengthIncrements *
                        precision.angularBlocks * precision.angularIncrements;

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
                                           precisionSlowCPU.thicknessBlocks * precisionSlowCPU.thicknessIncrements *
                                           precisionSlowCPU.angularBlocks * precisionSlowCPU.angularIncrements;

    const long long numOperationsFastCPU1 = (long long) opCount *
                                            precisionFastCPU.thicknessBlocks * precisionFastCPU.thicknessIncrements *
                                            precisionFastCPU.angularBlocks * precisionFastCPU.angularIncrements;
    const long long numOperationsFastCPU2 = numOperationsFastCPU1 *
                                            precisionFastCPU.lengthBlocks * precisionFastCPU.lengthIncrements;

    const long long numOperationsSlowGPU = (long long) opCount *
                                           precisionSlowGPU.thicknessIncrements * precisionSlowGPU.angularIncrements;

    const long long numOperationsFastGPU1 = (long long)opCount *
                                            precisionFastGPU.thicknessIncrements * precisionFastGPU.angularIncrements;
    const long long numOperationsFastGPU2 = numOperationsFastGPU1 * precisionFastGPU.lengthIncrements;

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
