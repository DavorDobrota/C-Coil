#include "Test.h"
#include "Coil.h"
#include "ctpl.h"




void testPerformanceCPU_ST()
{
    int nOp = 100000;
    std::vector<double> temp1;
    Coil testCoil = Coil(0.03, 0.03, 0.12, 3600);

    PrecisionArguments precision = testCoil.getPrecisionSettings();

    int numOperations = nOp *
                        precision.thicknessBlockCount * precision.thicknessIncrementCount *
                        precision.lengthBlockCount * precision.lengthIncrementCount *
                        precision.angularBlockCount * precision.angularIncrementCount;

    clock_t begin_time1 = clock();
    for (int i = 0; i < nOp; ++i){
        temp1 = testCoil.computeBFieldVector(i*0.000001, 0.0, 0.0);
    }
    printf("combined B  : %.1f MInc/s\n", 1e-6/ (float(clock() - begin_time1) / CLOCKS_PER_SEC / numOperations));

    double temp2;
    clock_t begin_time2 = clock();
    for (int i = 0; i < nOp; ++i){
        temp2 = testCoil.computeAPotentialAbs(i*0.000001, 0.0);
    }
    printf("potential A : %.1f MInc/s\n", 1e-6 / (float(clock() - begin_time2) / CLOCKS_PER_SEC / numOperations));

}

void testPerformanceForComputeAll(int nOps, int nRepeats, int nThreads)
{
    using namespace std::chrono;

    Coil testCoil = Coil(0.03, 0.03, 0.12, 3600);
    testCoil.setThreadCount(nThreads);

    const double radius = 0.1;

    PrecisionArguments precision = testCoil.getPrecisionSettings();

    const int numOperations = (long long) nOps *
                              precision.thicknessBlockCount * precision.thicknessIncrementCount *
                              precision.lengthBlockCount * precision.lengthIncrementCount *
                              precision.angularBlockCount * precision.angularIncrementCount;

    const int numOperationsGpu = (long long) nOps * 48 * 16 * 16;

    std::vector<double> cylindricalZArr;
    std::vector<double> cylindricalRArr;
    std::vector<double> cylindricalPhiArr;

    for (int i = 0; i < nOps ; i++)
    {
        cylindricalZArr.push_back(radius * cos(i * 2*Pi / nOps));
        cylindricalRArr.push_back(radius * sin(i * 2*Pi / nOps));
        cylindricalPhiArr.push_back(0.0);
    }

    std::vector<double> singleResultsX;
    std::vector<double> singleResultsY;
    std::vector<double> singleResultsZ;

    std::vector<double> acceleratedResultsX;
    std::vector<double> acceleratedResultsY;
    std::vector<double> acceleratedResultsZ;

    std::vector<double> singlePotential;
    std::vector<double> acceleratedPotential;

    high_resolution_clock::time_point begin_time = high_resolution_clock::now();
    for(int i = 0; i < nRepeats; i++)
    {
        testCoil.computeAllBFieldComponents(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                                            singleResultsX, singleResultsY, singleResultsZ,
                                            CPU_ST);
    }
    double interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("combined  B CPU : %.1f MInc/s\n", 1e-6 * (numOperations * nRepeats) / interval);

    begin_time = high_resolution_clock::now();
    for(int i = 0; i < nRepeats; i++)
    {
        testCoil.computeAllAPotentialAbs(cylindricalZArr, cylindricalRArr,
                                         singlePotential, CPU_ST);
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU : %.1f MInc/s\n", 1e-6 * (numOperations * nRepeats) / interval);


    testCoil.computeAllBFieldComponents(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                                        acceleratedResultsX, acceleratedResultsY, acceleratedResultsZ,
                                        GPU);

    begin_time = high_resolution_clock::now();
    for(int i = 0; i < nRepeats; i++)
    {
        testCoil.computeAllBFieldComponents(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                                            acceleratedResultsX, acceleratedResultsY, acceleratedResultsZ,
                                            GPU);
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("combined  B GPU : %.1f MInc/s\n", 1e-6 * (numOperationsGpu * nRepeats) / interval);

    begin_time = high_resolution_clock::now();
    for(int i = 0; i < nRepeats; i++)
    {
        testCoil.computeAllAPotentialAbs(cylindricalZArr, cylindricalRArr,
                                         acceleratedPotential, GPU);
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A GPU : %.1f MInc/s\n", 1e-6 * (numOperationsGpu * nRepeats) / interval);

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
    for(int i = 0; i < nRepeats; i++)
    {
        testCoil.computeAllAPotentialAbs(cylindricalZArr, cylindricalRArr,
                                         acceleratedPotential, CPU_MT);
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_MT : %.1f MInc/s\n", 1e-6 * (numOperations * nRepeats) / interval);
}