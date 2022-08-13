#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>

#include "Benchmark.h"
#include "Coil.h"
#include "CoilGroup.h"
#include "Test.h"
#include "Compare.h"


#pragma clang diagnostic push
#pragma ide diagnostic ignored "Simplify"
int main()
{
//    benchComputeAllFieldsWorkloadScalingMT(PrecisionFactor(1.0), 16, 22);
//    benchComputeAllFieldsWorkloadScalingGPU(PrecisionFactor(3.0), 25);
//    benchComputeAllFields();

//    benchComputeAllFieldsEveryCoilType(200003, 8);

//    benchMInductanceZAxisMTScaling(16);
//    benchMInductanceGeneralMTScaling(16);
//    benchSelfInductance();
//    benchForceZAxisMTScaling(16);
//    benchForceGeneralMTScaling(16);

//    compMInductanceGeneralMisalignedCoils();
//    compMInductanceGeneralCase();

//    benchComputeAllFieldsWorkloadScalingMT(PrecisionFactor(5.0), 16);

//    compPrecisionCPUvsGPU();

//    compMInductanceGeneralGraphs();

//    benchComputeAllFields();

//    compMutualInductanceZAxis();

//    testAmpereForceGeneralForZAxis();

//    compAmpereForceZAxis();
//    compAmpereForceThinCoilsZAxis();

//    benchCoilGroupMInductanceAndForce(16);

//    benchCoilGroupComputeAllFieldsMTvsMTD(16, 100'000);
//    benchCoilMInductanceAndForceComputeAllMTvsMTD(PrecisionFactor(), 16);

//    benchCoilGroupComputeAllFields(PrecisionFactor(3.0), 1, 1'048'576, 16);
//    benchCoilGroupComputeAllFields(PrecisionFactor(1.0), 100, 131'072, 16);

//    benchCoilGroupComputeAllFieldsGPU(1000, 1000000);

//    benchCoilGroupComputeAllFieldsMTvsMTD();
//    benchCoilGroupMInductanceAndForce();

//    benchCoilMInductanceAndForceComputeAll();
//    benchCoilMInductanceAndForceComputeAllGPU();

    benchCoilGroupMInductanceAndForceAll();
    benchCoilGroupMInductanceAndForceAllGPU();

//    testCoilMInductanceArrangements();
//    testCoilForceArrangements();

//    testGroupMInductanceArrangements();
//    testGroupForceArrangements();

//    testMInductanceZAxisArgumentGeneration();
//    testMInductanceGeneralArgumentGeneration();


    return 0;
}
#pragma clang diagnostic pop
