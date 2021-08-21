#ifndef GENERAL_COIL_PROGRAM_TEST_H
#define GENERAL_COIL_PROGRAM_TEST_H

#include "ComputeMethod.h"
#include "Coil.h"


const int g_default_nThreads = 8;

void testNewCoilParameters();
void testMethodPrecisionCompareCPUvsGPU();
void testCoilPositionAndRotation();
void testCoilMutualInductanceForSpecialCase();
void testVector3();

void testPerformanceCPU_ST(int nOps = 50'000);
void testPerformanceForComputeAll(PrecisionFactor precisionFactor = PrecisionFactor(),
                                  int nOps = 60'000, int nRepeats = 1, int nThreads = g_default_nThreads);

void testCoilMutualInductanceZAxis();
void testCoilMutualInductanceZAxisArgumentGeneration();
void testCoilMutualInductanceZAxisDifferentGeometries();
void testCoilMutualInductanceZAxisPerformance(ComputeMethod method = CPU_ST, int nThreads = g_default_nThreads);
void testCoilMutualInductanceZAxisMTScaling(int maxThreads = g_default_nThreads);
void testCoilSelfInductance();

void testCoilMutualInductanceGeneralForZAxis(ComputeMethod method = CPU_ST, int nThreads = g_default_nThreads);
void testCoilMutualInductanceGeneralPerformance(ComputeMethod method = CPU_ST, int nThreads = g_default_nThreads);
void testCoilMutualInductanceGeneralMTScaling(int maxThreads = g_default_nThreads);
void testCoilMutualInductanceGeneralArgumentGeneration();
void testCoilMutualInductanceGeneralDifferentGeometries();
void testCoilMutualInductanceGeneralGraphs();

void testCoilMutualInductanceGeneralThinCoilAndFilament();
void testCoilMutualInductanceGeneralThinCoilAndThinCoil();
void testCoilMutualInductanceGeneralPancakeAndPancake();
void testCoilMutualInductanceGeneralRectangularAndFilament();

void testCoilAmpereForceZAxis();
void testCoilAmpereForceZAxisPerformance(ComputeMethod method = CPU_ST, int nThreads = g_default_nThreads);
void testCoilAmpereForceZAxisMTScaling(int maxThreads = g_default_nThreads);
void testCoilAmpereForceGeneralForZAxis();
void testCoilAmpereForceGeneralPerformance(ComputeMethod method = CPU_ST, int nThreads = g_default_nThreads);
void testCoilAmpereForceZGeneralMTScaling(int maxThreads = g_default_nThreads);
void testCoilGradientTensor();

void testCoilAmpereForceForFilamentsZAxis();
void testCoilAmpereForceGeneralCase();
void testCoilAmpereForceThinCoils();
void testCoilAmpereForceFilamentsGeneral();


#endif //GENERAL_COIL_PROGRAM_TEST_H
