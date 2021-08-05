#include "Test.h"
#include "Coil.h"
#include "ctpl.h"

#include <cmath>

void testPerformanceCPU_ST(int nOps)
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

void testPerformanceForComputeAll(int nOps, int nRepeats, int nThreads)
{
    using namespace std::chrono;

    Coil testCoil = Coil(0.03, 0.03, 0.12, 3600);
    testCoil.setThreadCount(nThreads);

    const double radius = 0.1;

    PrecisionArguments precision = testCoil.getPrecisionSettings();

    const long long numOperations1 = (long long) nOps *
                                     precision.thicknessBlockCount * precision.thicknessIncrementCount  *
                                     precision.angularBlockCount * precision.angularIncrementCount;

    const long long numOperations2 = numOperations1 * precision.lengthBlockCount * precision.lengthIncrementCount;

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

    // ST methods - repetitions removed for single core, it is just too slow compared to the rest
    begin_time = high_resolution_clock::now();
    cpuPotential = testCoil.computeAllAPotentialAbs(positionValues, CPU_ST);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_ST : %.1f MInc/s | eff: %.1f MInc/s\n",
           1e-6 * numOperations1 / interval, 1e-6 * numOperations2 / interval);

    begin_time = high_resolution_clock::now();
    cpuFieldVectors = testCoil.computeAllBFieldComponents(positionValues, CPU_ST);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("combined  B CPU_ST : %.1f MInc/s | eff: %.1f MInc/s\n",
           1e-6 * numOperations1 / interval, 1e-6 * numOperations2 / interval);

    begin_time = high_resolution_clock::now();
    cpuGradientMatrices = testCoil.computeAllBGradientTensors(positionValues, CPU_ST);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G CPU_ST : %.1f MInc/s | eff: %.1f MInc/s\n",
           1e-6 * numOperations1 / interval, 1e-6 * numOperations2 / interval);

    // MT methods
    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nRepeats; i++)
        cpuPotential = testCoil.computeAllAPotentialAbs(positionValues, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_MT : %.1f MInc/s | eff: %.1f MInc/s\n",
           1e-6 * (numOperations1 * nRepeats) / interval, 1e-6 * (numOperations2 * nRepeats) / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nRepeats; i++)
        cpuFieldVectors = testCoil.computeAllBFieldComponents(positionValues, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("combined  B CPU_MT : %.1f MInc/s | eff: %.1f MInc/s\n",
           1e-6 * (numOperations1 * nRepeats) / interval, 1e-6 * (numOperations2 * nRepeats) / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nRepeats; i++)
        cpuGradientMatrices = testCoil.computeAllBGradientTensors(positionValues, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G CPU_MT : %.1f MInc/s | eff: %.1f MInc/s\n",
           1e-6 * (numOperations1 * nRepeats) / interval, 1e-6 * (numOperations2 * nRepeats) / interval);

    // GPU methods
    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nRepeats; i++)
        gpuPotential = testCoil.computeAllAPotentialAbs(positionValues, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A GPU : %.1f MInc/s\n", 1e-6 * (numOperationsGpu * nRepeats) / interval);


    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nRepeats; i++)
        gpuFieldVectors = testCoil.computeAllBFieldComponents(positionValues, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("combined  B GPU : %.1f MInc/s\n", 1e-6 * (numOperationsGpu * nRepeats) / interval);

    printf("\n");
}