#ifndef GENERAL_COIL_PROGRAM_TEST_H
#define GENERAL_COIL_PROGRAM_TEST_H

#include "ComputeMethod.h"
#include "Coil.h"


void testNewCoilParameters();
void compMethodPrecisionCPUvsGPU();
void compMInductanceForSpecialCase();
void testVector3();
void testCoilPositionAndRotation();

void compCoilGroupMTD(int numCoils = 100, int numPoints = 10'000, int threadCount = 8, bool print = true);

void compMutualInductanceZAxis();
void testMInductanceZAxisArgumentGeneration();
void testMInductanceZAxisDifferentGeometries();
void compSelfInductance();

void testMInductanceGeneralArgumentGeneration();
void testMInductanceGeneralDifferentGeometries();
void compMInductanceGeneralGraphs();
void compMInductanceGeneralEdgeCases();

void compMInductanceGeneralMisalignedCoils();
void compMInductanceGeneralParallelAxes();
void compMInductanceGeneralCase();

void compAmpereForceZAxis();

void testAmpereForceGeneralForZAxis();
void testGradientTensor();
void compForceOnDipoleVsAmpereForce();

void compAmpereForceFilamentsZAxis();
void compAmpereForceThinCoilsZAxis();
void compAmpereForceFilamentsGeneral();
void compAmpereForceThickCoilsGeneral();


#endif //GENERAL_COIL_PROGRAM_TEST_H
