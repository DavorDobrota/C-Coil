#ifndef GENERAL_COIL_PROGRAM_TEST_H
#define GENERAL_COIL_PROGRAM_TEST_H

void testLegendrePolynomials();

void testNewCoilParameters();

void testPerformanceCPU_ST();
void testPerformanceForComputeAll(int nOps = 80'000, int nRepeats = 1, int nThreads = 8);

void testMethodPrecisionCompareCPUvsGPU();

void testCoilMutualInductanceZAxis();
void testCoilMutualInductanceZAxisArgumentGeneration();
void testCoilMutualInductanceZAxisDifferentGeometries();
void testCoilMutualInductanceZAxisPerformance();
void testCoilSelfInductance();

void testOldCoilMutualInductanceZAxis();
void testOldCoilMutualInductanceZAxisPerformance();
void testOldCoilMutualInductanceGeneralPerformance();
void testOldCoilSelfInductance();

void testCoilMutualInductanceGeneralForZAxis();
void testCoilMutualInductanceGeneralArgumentGeneration();
void testCoilMutualInductanceGeneralDifferentGeometries();
void testCoilMutualInductanceForSpecialCase();

void testCoilMutualInductanceGeneralThinCoilAndFilament();
void testCoilMutualInductanceGeneralThinCoilAndThinCoil();
void testCoilMutualInductanceGeneralPancakeAndPancake();
void testCoilMutualInductanceGeneralRectangularAndFilament();


#endif //GENERAL_COIL_PROGRAM_TEST_H