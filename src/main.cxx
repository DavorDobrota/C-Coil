#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>

#include "Benchmark.h"
#include "Coil.h"
#include "CoilGroup.h"
#include "Test.h"
#include "Compare.h"

//    benchMathFunctions();

#pragma clang diagnostic push
#pragma ide diagnostic ignored "Simplify"
int main()
{
//    benchComputeAllFieldsWorkloadScalingMT(PrecisionFactor(1.0), 16, 25);
//    benchComputeAllFieldsWorkloadScalingGPU(PrecisionFactor(3.0), 27);
//    benchComputeAllFields();

//    benchComputeAllFieldsEveryCoilType(1'000'000, 8);

//    benchMInductanceZAxisMTScaling(16);
//    benchMInductanceGeneralMTScaling(16);
//    benchSelfInductance();
//    benchForceZAxisMTScaling(16);
//    benchForceGeneralMTScaling(16);

//    compMInductanceGeneralMisalignedCoils();
//    compMInductanceGeneralCase();
//    compMInductanceGeneralParallelAxes();

//    benchComputeAllFieldsWorkloadScalingMT(PrecisionFactor(5.0), 16);

//    compPrecisionCPUvsGPU();

//    compMInductanceGeneralGraphs();

//    benchComputeAllFields();

//    compMutualInductanceZAxis();

//    testAmpereForceGeneralForZAxis();

//    compAmpereForceZAxis();
//    compAmpereForceThinCoilsZAxis();
//    compAmpereForceThickCoilsGeneral();
//    compAmpereForceFilamentsGeneral();
//    benchCoilGroupMInductanceAndForce(16);

//    benchCoilGroupComputeAllFieldsMTvsMTD(16, 100'000);
//    benchCoilMInductanceAndForceComputeAllMTvsMTD(PrecisionFactor(), 16);

//    benchCoilGroupComputeAllFields(PrecisionFactor(3.0), 1, 1'048'576, 16);
//    benchCoilGroupComputeAllFields(PrecisionFactor(1.0), 100, 131'072, 16);

//    benchCoilGroupComputeAllFieldsGPU(1000, 1000000);

//    benchCoilGroupComputeAllFieldsMTvsMTD();
    benchCoilGroupMInductanceAndForce(1, 12);

//    benchCoilMInductanceAndForceComputeAll();
//    benchCoilMInductanceAndForceComputeAllGPU();

//    benchCoilGroupMInductanceAndForceAll();
//    benchCoilGroupMInductanceAndForceAllGPU();

//    testCoilMInductanceArrangements();
//    testCoilForceArrangements();

//    testGroupMInductanceArrangements();
//    testGroupForceArrangements();

//    testMInductanceZAxisArgumentGeneration();
//    testMInductanceGeneralArgumentGeneration();

//    benchCoilGroupComputeAllFieldsMTScaling(PrecisionFactor(1.0), 16, 100, 20);
//    benchCoilGroupComputeAllFieldsGPUScaling(PrecisionFactor(3.0), 100, 20);

//    benchMathFunctions();

//    testCoilGroupFieldsMTD();


    return 0;
}
#pragma clang diagnostic pop
