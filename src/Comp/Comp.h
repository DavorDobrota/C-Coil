#ifndef GENERAL_COIL_PROGRAM_COMP_H
#define GENERAL_COIL_PROGRAM_COMP_H

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
void compSelfInductance();

void compCoilGroupMTD(int numCoils = 100, int numPoints = 10'000, int threadCount = g_default_nThreads, bool print = true);

#endif //GENERAL_COIL_PROGRAM_COMP_H
