#ifndef GENERAL_COIL_PROGRAM_TEST_H
#define GENERAL_COIL_PROGRAM_TEST_H

#include "ComputeMethod.h"
#include "Coil.h"


const int g_default_nThreads = 8;

void testFunctionPerformance();
void testNewCoilParameters();
void testMethodPrecisionCompareCPUvsGPU();
void testCoilMutualInductanceForSpecialCase();
void testVector3();
void testCoilPositionAndRotation();
void testCoilGroupMTD(int numCoils = 100, int numPoints = 10'000, int threadCount = 8, bool print = true);
void testCoilGroupMTvsMTD(int threadCount = 8, int numPoints = 20'000);
void testCoilGroupMTDInductanceAndForce(int threadCount = 8);

void testPerformanceCPU_ST(int nOps = 50'000);
void testPerformanceForComputeAll(PrecisionFactor precisionFactor = PrecisionFactor(),
                                  int nOps = 20'000, int nRepeats = 1, int nThreads = g_default_nThreads);
void testPerformanceForVariousCoilTypes(int nOps = 100'000, int nThreads = g_default_nThreads);

void testMutualInductanceZAxis();
void testMutualInductanceZAxisArgumentGeneration();
void testMutualInductanceZAxisDifferentGeometries();
void testCoilMutualInductanceZAxisPerformance(ComputeMethod method = CPU_ST, int nThreads = g_default_nThreads);
void testMutualInductanceZAxisMTScaling(int maxThreads = g_default_nThreads);
void testSelfInductance();

void testMutualInductanceGeneralForZAxis(ComputeMethod method = CPU_ST, int nThreads = g_default_nThreads);
void testMutualInductanceGeneralPerformance(ComputeMethod method = CPU_ST, int nThreads = g_default_nThreads);
void testMutualInductanceGeneralMTScaling(int maxThreads = g_default_nThreads);
void testMutualInductanceGeneralArgumentGeneration();
void testMutualInductanceGeneralDifferentGeometries();
void testMutualInductanceGeneralGraphs();
void testMutualInductanceGeneralEdgeCases();

void testMutualInductanceGeneralMisalignedCoils();
void testMutualInductanceGeneralParallelAxes();
void testMutualInductanceGeneralConway();

void testAmpereForceZAxis();
void testAmpereForceZAxisPerformance(ComputeMethod method = CPU_ST, int nThreads = g_default_nThreads);
void testAmpereForceZAxisMTScaling(int maxThreads = g_default_nThreads);
void testAmpereForceGeneralForZAxis();
void testAmpereForceGeneralPerformance(ComputeMethod method = CPU_ST, int nThreads = g_default_nThreads);
void testAmpereForceGeneralMTScaling(int maxThreads = g_default_nThreads);
void testGradientTensor();
void testForceOnDipole();

void testAmpereForceForFilamentsZAxis();
void testAmpereForceGeneralCase();
void testAmpereForceThinCoils();
void testAmpereForceFilamentsGeneral();


#endif //GENERAL_COIL_PROGRAM_TEST_H
