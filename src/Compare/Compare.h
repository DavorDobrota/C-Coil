#ifndef GENERAL_COIL_PROGRAM_COMPARE_H
#define GENERAL_COIL_PROGRAM_COMPARE_H

#include "Coil.h"
#include "Benchmark.h"


void compMethodPrecisionCPUvsGPU();
void compMInductanceForSpecialCase();

void compAmpereForceFilamentsZAxis();
void compAmpereForceThickCoilsGeneral();
void compAmpereForceThinCoilsZAxis();
void compAmpereForceFilamentsGeneral();
void compAmpereForceZAxis();
void compForceOnDipoleVsAmpereForce();

void compMInductanceGeneralMisalignedCoils();
void compMInductanceGeneralParallelAxes();
void compMInductanceGeneralCase();
void compMInductanceGeneralGraphs();
void compMInductanceGeneralEdgeCases();

void compMutualInductanceZAxis();
void compPrecisionCPUvsGPU(); // TODO - implement binding
void compSelfInductance();

void compCoilGroupMTD(int coilCount = 100, int pointCount = 10'000, int threadCount = g_defaultThreadCount, bool print = true);

#endif //GENERAL_COIL_PROGRAM_COMPARE_H
