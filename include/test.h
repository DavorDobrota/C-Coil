
#ifndef GENERAL_COIL_PROGRAM_TEST_H
#define GENERAL_COIL_PROGRAM_TEST_H

void testLegendrePolynomials();

void testNewCoilParameters();

void testPerformanceCPU_ST();
void testPerformanceForComputeAll();

void testMethodPrecisionCompareCPUvsGPU();

void testCoilMutualInductanceZAxis();
void testCoilMutualInductanceZAxisPerformance();

void testOldCoilMutualInductanceZAxis();
void testOldCoilMutualInductanceZAxisPerformance();
void testOldCoilMutualInductanceGeneralPerformance();
void testOldCoilSelfInductance();

void testCoilMutualInductanceGeneralForZAxis();
void testCoilMutualInductanceForSpecialCase();
void testCoilMutualInductanceGeneralThinCoilAndFilament();
void testCoilMutualInductanceGeneralThinCoilAndThinCoil();


#endif //GENERAL_COIL_PROGRAM_TEST_H
