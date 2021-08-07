#ifndef GENERAL_COIL_PROGRAM_TEST_H
#define GENERAL_COIL_PROGRAM_TEST_H

#include "ComputeMethod.h"
#include "Coil.h"

void testLegendrePolynomials();
void testNewCoilParameters();
void testMethodPrecisionCompareCPUvsGPU();
void testCoilMutualInductanceForSpecialCase();
void testVector3();

void testPerformanceCPU_ST(int nOps = 50'000);
void testPerformanceForComputeAll(PrecisionFactor precisionFactor = PrecisionFactor(),
                                  int nOps = 60'000, int nRepeats = 1, int nThreads = 12);

void testCoilMutualInductanceZAxis();
void testCoilMutualInductanceZAxisArgumentGeneration();
void testCoilMutualInductanceZAxisDifferentGeometries();
void testCoilMutualInductanceZAxisPerformance(ComputeMethod method = CPU_ST, int nThreads = 12);
void testCoilSelfInductance();

void testOldCoilMutualInductanceZAxis();
void testOldCoilMutualInductanceZAxisPerformance();
void testOldCoilMutualInductanceGeneralPerformance();
void testOldCoilSelfInductance();

void testCoilMutualInductanceGeneralForZAxis();
void testCoilMutualInductanceGeneralPerformance();
void testCoilMutualInductanceGeneralArgumentGeneration();
void testCoilMutualInductanceGeneralDifferentGeometries();
void testCoilMutualInductanceGeneralGraphs();

void testCoilMutualInductanceGeneralThinCoilAndFilament();
void testCoilMutualInductanceGeneralThinCoilAndThinCoil();
void testCoilMutualInductanceGeneralPancakeAndPancake();
void testCoilMutualInductanceGeneralRectangularAndFilament();

void testCoilAmpereForceZAxis();
void testCoilAmpereForceGeneralForZAxis();
void testCoilGradientTensor();

void testCoilAmpereForceForFilamentsZAxis();
void testCoilAmpereForceGeneralCase();
void testCoilAmpereForceThinCoils();


#endif //GENERAL_COIL_PROGRAM_TEST_H
