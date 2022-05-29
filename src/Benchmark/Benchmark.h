#ifndef GENERAL_COIL_PROGRAM_BENCHMARK_H
#define GENERAL_COIL_PROGRAM_BENCHMARK_H

#include "Coil.h"


const int g_defaultThreadCount = 8;
const int maxPot = 25;

void benchMathFunctions();

void benchComputeFieldsST(int opCount = 50'000);
void benchComputeAllFields(PrecisionFactor precisionFactor = PrecisionFactor(),
                           int opCount = 20'000, int repeatCount = 1, int threadCount = g_defaultThreadCount);
void benchComputeAllFieldsEveryCoilType(int opCount = 100'000, int threadCount = g_defaultThreadCount);
void benchComputeAllFieldsWorkloadScalingMT(PrecisionFactor precisionFactor = PrecisionFactor(),
                                            int threadCount = g_defaultThreadCount, int maxPointsLog2 = maxPot);

void benchMInductanceZAxis(ComputeMethod computeMethod = CPU_ST, int threadCount = g_defaultThreadCount);
void benchMInductanceZAxisMTScaling(int maxThreadCount = g_defaultThreadCount);
void benchSelfInductance();

void benchMInductanceGeneral(ComputeMethod computeMethod = CPU_ST, int threadCount = g_defaultThreadCount);
void benchMInductanceGeneralMTScaling(int maxThreadCount = g_defaultThreadCount);

void benchForceZAxis(ComputeMethod computeMethod = CPU_ST, int threadCount = g_defaultThreadCount);
void benchForceZAxisMTScaling(int maxThreadCount = g_defaultThreadCount);

void benchForceGeneral(ComputeMethod computeMethod = CPU_ST, int threadCount = g_defaultThreadCount);
void benchForceGeneralMTScaling(int maxThreadCount = g_defaultThreadCount);

void benchCoilGroupMTvsMTD(int threadCount = g_defaultThreadCount, int pointCount = 20'000);
void benchCoilGroupComputeAllFieldsMTD(int threadCount = g_defaultThreadCount);
void benchCoilGroupMInductanceAndForceMTD(int threadCount = g_defaultThreadCount);
void benchMInductanceAndForceComputeAll(PrecisionFactor precisionFactor = PrecisionFactor(), int threadCount = g_defaultThreadCount);

#endif //GENERAL_COIL_PROGRAM_BENCHMARK_H
