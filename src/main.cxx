#include <cstdio>
#include <cmath>
#include <cstdint>
#include <cstring>

#include "Coil.h"
#include "Test.h"


int main()
{

//    testMutualInductanceGeneralThinCoilAndThinCoil();
//    testAmpereForceThinCoils();
//    testAmpereForceFilamentsGeneral();

//    testAmpereForceZAxisMTScaling(16);
//    testAmpereForceGeneralMTScaling(16);

//    testFunctionPerformance();

//    testPerformanceForComputeAll(PrecisionFactor(1.0), 50000, 5, 2);
//    testPerformanceForComputeAll(PrecisionFactor(1.0), 50000, 5, 3);
//    testPerformanceForComputeAll(PrecisionFactor(1.0), 50000, 5, 4);
//    testPerformanceForComputeAll(PrecisionFactor(1.0), 50000, 5, 5);
//    testPerformanceForComputeAll(PrecisionFactor(1.0), 50000, 5, 6);
//    testPerformanceForComputeAll(PrecisionFactor(1.0), 50000, 5, 7);
//    testPerformanceForComputeAll(PrecisionFactor(1.0), 50000, 5, 8);

//    testMutualInductanceZAxisMTScaling(16);
//    testMutualInductanceGeneralMTScaling(16);

//    testMutualInductanceGeneralEdgeCases();
//    testMutualInductanceZAxis();

//    testSelfInductance();
//    testMutualInductanceGeneralConway();

//testAmpereForceGeneralCase();
//testAmpereForceThinCoils();
//    testAmpereForceFilamentsGeneral();

testPerformanceForVariousCoilTypes(80'640);

    return 0;
}
