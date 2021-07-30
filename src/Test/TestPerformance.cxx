#include "Test.h"
#include "Coil.h"
#include "ctpl.h"

#include <cmath>

void testPerformanceCPU_ST(int nOps)
{
    using namespace std::chrono;

    double temp1;
    std::vector<double> temp2;
    std::vector<double> temp3;

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
        temp1 = testCoil.computeAPotentialAbs(i*0.000001, 0.0);
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("potential A : %.1f MInc/s\n", 1e-6 * numOperations / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nOps; ++i){
        temp2 = testCoil.computeBFieldVector(i*0.000001, 0.0, 0.0);
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("combined B  : %.1f MInc/s\n", 1e-6 * numOperations / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nOps; ++i){
        temp3 = testCoil.computeBGradientTensor(i*0.000001, 0.0, 0.0);
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

    const long long numOperations = (long long) nOps *
                                    precision.thicknessBlockCount * precision.thicknessIncrementCount *
                                    precision.lengthBlockCount * precision.lengthIncrementCount *
                                    precision.angularBlockCount * precision.angularIncrementCount;

    const long long numOperationsGpu = nOps * 48 * 16 * 16;

    std::vector<double> cylindricalZArr;
    std::vector<double> cylindricalRArr;
    std::vector<double> cylindricalPhiArr;

    for (int i = 0; i < nOps ; i++)
    {
        cylindricalZArr.push_back(radius * cos(i * 2*M_PI / nOps));
        cylindricalRArr.push_back(radius * sin(i * 2*M_PI / nOps));
        cylindricalPhiArr.push_back(0.0);
    }

    std::vector<double> singleResultsX, singleResultsY, singleResultsZ;

    std::vector<double> acceleratedResultsX, acceleratedResultsY, acceleratedResultsZ;

    std::vector<double> singlePotential;
    std::vector<double> acceleratedPotential;

    std::vector<double> tensorXX, tensorXY, tensorXZ, tensorYX, tensorYY, tensorYZ, tensorZX, tensorZY, tensorZZ;

    high_resolution_clock::time_point begin_time;
    double interval;

    // ST methods - repetitions removed for single core, it is just too slow compared to the rest
    begin_time = high_resolution_clock::now();
    testCoil.computeAllAPotentialAbs(cylindricalZArr, cylindricalRArr,
                                     singlePotential, CPU_ST);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_ST : %.1f MInc/s\n", 1e-6 * numOperations / interval);

    begin_time = high_resolution_clock::now();
    testCoil.computeAllBFieldComponents(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                                        singleResultsX, singleResultsY, singleResultsZ,
                                        CPU_ST);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("combined  B CPU_ST : %.1f MInc/s\n", 1e-6 * numOperations / interval);

    begin_time = high_resolution_clock::now();
    testCoil.computeAllBGradientTensors(cylindricalZArr, cylindricalRArr, singlePotential,
                                        tensorXX, tensorXY, tensorXZ,
                                        tensorYX, tensorYY, tensorYZ,
                                        tensorZX, tensorZY, tensorZZ,
                                        CPU_ST);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G CPU_ST : %.1f MInc/s\n", 1e-6 * numOperations / interval);

    // MT methods
    begin_time = high_resolution_clock::now();
    for(int i = 0; i < nRepeats; i++)
    {
        testCoil.computeAllAPotentialAbs(cylindricalZArr, cylindricalRArr,
                                         acceleratedPotential, CPU_MT);
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_MT : %.1f MInc/s\n", 1e-6 * (numOperations * nRepeats) / interval);

    begin_time = high_resolution_clock::now();
    for(int i = 0; i < nRepeats; i++)
    {
        testCoil.computeAllBFieldComponents(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                                            acceleratedResultsX, acceleratedResultsY, acceleratedResultsZ,
                                            CPU_MT);
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("combined  B CPU_MT : %.1f MInc/s\n", 1e-6 * (numOperations * nRepeats) / interval);

    begin_time = high_resolution_clock::now();
    testCoil.computeAllBGradientTensors(cylindricalZArr, cylindricalRArr, singlePotential,
                                        tensorXX, tensorXY, tensorXZ,
                                        tensorYX, tensorYY, tensorYZ,
                                        tensorZX, tensorZY, tensorZZ,
                                        CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G CPU_MT : %.1f MInc/s\n", 1e-6 * numOperations / interval);

    // GPU methods
    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nRepeats; i++)
    {
        testCoil.computeAllBFieldComponents(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                                            acceleratedResultsX, acceleratedResultsY, acceleratedResultsZ,
                                            GPU);
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("combined  B GPU : %.1f MInc/s\n", 1e-6 * (numOperationsGpu * nRepeats) / interval);

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < nRepeats; i++)
    {
        testCoil.computeAllAPotentialAbs(cylindricalZArr, cylindricalRArr,
                                         acceleratedPotential, GPU);
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A GPU : %.1f MInc/s\n", 1e-6 * (numOperationsGpu * nRepeats) / interval);

    printf("\n");
}