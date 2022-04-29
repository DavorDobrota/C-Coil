#ifndef GENERAL_COIL_PROGRAM_BENCHMARK_H
#define GENERAL_COIL_PROGRAM_BENCHMARK_H

#include "Coil.h"


const int g_default_nThreads = 8;

void benchmarkCPU_ST(int nOps = 50'000);
void benchmarkComputeAll(PrecisionFactor precisionFactor = PrecisionFactor(),
                         int nOps = 20'000, int nRepeats = 1, int nThreads = g_default_nThreads);
void benchmarkVariousCoilTypes(int nOps = 100'000, int nThreads = g_default_nThreads);

void benchmarkMutualInductanceGeneral(ComputeMethod computeMethod = CPU_ST, int nThreads = g_default_nThreads);
void benchmarkMutualInductanceGeneralForZAxis(ComputeMethod computeMethod = CPU_ST, int nThreads = g_default_nThreads);
void benchmarkMutualInductanceGeneralMTScaling(int maxThreads = g_default_nThreads);

void benchmarkCoilMutualInductanceZAxis(ComputeMethod computeMethod = CPU_ST, int nThreads = g_default_nThreads);
void benchmarkMutualInductanceZAxisMTScaling(int maxThreads = g_default_nThreads);
void benchmarkSelfInductance();

void benchmarkAmpereForceZAxis(ComputeMethod computeMethod = CPU_ST, int nThreads = g_default_nThreads);
void benchmarkAmpereForceGeneral(ComputeMethod computeMethod = CPU_ST, int nThreads = g_default_nThreads);
void benchmarkAmpereForceZAxisMTScaling(int maxThreads = g_default_nThreads);
void benchmarkAmpereForceGeneralMTScaling(int maxThreads = g_default_nThreads);

void benchmarkMathFunctions();

void benchmarkCoilGroupMTvsMTD(int threadCount = g_default_nThreads, int numPoints = 20'000);
void benchmarkCoilGroupMTDFields(int threadCount = g_default_nThreads);
void benchmarkCoilGroupMTDInductanceAndForce(int threadCount = g_default_nThreads);
void benchmarkMInductanceAndForceComputeAll(PrecisionFactor precisionFactor = PrecisionFactor(), int threadCount = g_default_nThreads);

#endif //GENERAL_COIL_PROGRAM_BENCHMARK_H
