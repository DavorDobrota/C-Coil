#ifndef GENERAL_COIL_PROGRAM_TEST_H
#define GENERAL_COIL_PROGRAM_TEST_H

#include "ComputeMethod.h"

void testLegendrePolynomials();
void testNewCoilParameters();
void testMethodPrecisionCompareCPUvsGPU();
void testCoilMutualInductanceForSpecialCase();

void testPerformanceCPU_ST();
void testPerformanceForComputeAll(int nOps = 80'000, int nRepeats = 1, int nThreads = 16);

void testCoilMutualInductanceZAxis();
void testCoilMutualInductanceZAxisArgumentGeneration();
void testCoilMutualInductanceZAxisDifferentGeometries();
void testCoilMutualInductanceZAxisPerformance(ComputeMethod method = CPU_ST);
void testCoilSelfInductance();

void testOldCoilMutualInductanceZAxis();
void testOldCoilMutualInductanceZAxisPerformance();
void testOldCoilMutualInductanceGeneralPerformance();
void testOldCoilSelfInductance();

void testCoilMutualInductanceGeneralForZAxis();
void testCoilMutualInductanceGeneralPerformance();
void testCoilMutualInductanceGeneralArgumentGeneration();
void testCoilMutualInductanceGeneralDifferentGeometries();

void testCoilMutualInductanceGeneralThinCoilAndFilament();
void testCoilMutualInductanceGeneralThinCoilAndThinCoil();
void testCoilMutualInductanceGeneralPancakeAndPancake();
void testCoilMutualInductanceGeneralRectangularAndFilament();

void testCoilAmpereForceZAxis();
void testCoilAmpereForceGeneralForZAxis();

#endif //GENERAL_COIL_PROGRAM_TEST_H
