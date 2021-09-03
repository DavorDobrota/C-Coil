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

//    testFunctionPerformance();

//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 2);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 3);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 4);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 5);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 6);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 7);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 8);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 9);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 10);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 11);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 12);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 13);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 14);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 15);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 16);
    testMutualInductanceZAxisMTScaling(16);

//   testMutualInductanceGeneralMTScaling(24);

//    testMutualInductanceGeneralEdgeCases();
//    testMutualInductanceZAxis();

    testMutualInductanceGeneralConway();

    return 0;
}
