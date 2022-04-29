#ifndef GENERAL_COIL_PROGRAM_TEST_H
#define GENERAL_COIL_PROGRAM_TEST_H

#include "ComputeMethod.h"
#include "Coil.h"


void testNewCoilParameters();
void testMethodPrecisionCompareCPUvsGPU();
void testCoilMutualInductanceForSpecialCase();
void testVector3();
void testCoilPositionAndRotation();

void testCoilGroupMTD(int numCoils = 100, int numPoints = 10'000, int threadCount = 8, bool print = true);

void testMutualInductanceZAxis();
void testMutualInductanceZAxisArgumentGeneration();
void testMutualInductanceZAxisDifferentGeometries();
void testSelfInductance();

void testMutualInductanceGeneralArgumentGeneration();
void testMutualInductanceGeneralDifferentGeometries();
void testMutualInductanceGeneralGraphs();
void testMutualInductanceGeneralEdgeCases();

void testMutualInductanceGeneralMisalignedCoils();
void testMutualInductanceGeneralParallelAxes();
void testMutualInductanceGeneralConway();

void testAmpereForceZAxis();

void testAmpereForceGeneralForZAxis();
void testGradientTensor();
void testForceOnDipole();

void testAmpereForceForFilamentsZAxis();
void testAmpereForceGeneralCase();
void testAmpereForceThinCoils();
void testAmpereForceFilamentsGeneral();


#endif //GENERAL_COIL_PROGRAM_TEST_H
