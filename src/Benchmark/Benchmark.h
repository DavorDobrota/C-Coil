#ifndef GENERAL_COIL_PROGRAM_BENCHMARK_H
#define GENERAL_COIL_PROGRAM_BENCHMARK_H

#include "Coil.h"


const int g_default_nThreads = 8;

void benchMathFunctions();

void benchComputeFieldsST(int nOps = 50'000);
void benchComputeAllFields(PrecisionFactor precisionFactor = PrecisionFactor(),
                           int nOps = 20'000, int nRepeats = 1, int nThreads = g_default_nThreads);
void benchComputeAllFieldsEveryCoilType(int nOps = 100'000, int nThreads = g_default_nThreads);


void benchMInductanceZAxis(ComputeMethod computeMethod = CPU_ST, int nThreads = g_default_nThreads);
void benchMInductanceZAxisMTScaling(int maxThreads = g_default_nThreads);
void benchSelfInductance();

void benchMInductanceGeneral(ComputeMethod computeMethod = CPU_ST, int nThreads = g_default_nThreads);
void benchMInductanceGeneralMTScaling(int maxThreads = g_default_nThreads);

void benchForceZAxis(ComputeMethod computeMethod = CPU_ST, int nThreads = g_default_nThreads);
void benchForceZAxisMTScaling(int maxThreads = g_default_nThreads);

void benchForceGeneral(ComputeMethod computeMethod = CPU_ST, int nThreads = g_default_nThreads);
void benchForceGeneralMTScaling(int maxThreads = g_default_nThreads);

void benchCoilGroupMTvsMTD(int threadCount = g_default_nThreads, int numPoints = 20'000);
void benchCoilGroupComputeAllFieldsMTD(int threadCount = g_default_nThreads);
void benchCoilGroupMInductanceAndForceMTD(int threadCount = g_default_nThreads);
void benchMInductanceAndForceComputeAll(PrecisionFactor precisionFactor = PrecisionFactor(), int threadCount = g_default_nThreads);

#endif //GENERAL_COIL_PROGRAM_BENCHMARK_H
