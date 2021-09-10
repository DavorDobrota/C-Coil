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

    testFunctionPerformance();

//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 2);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 3);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 4);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 5);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 6);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 7);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 8);

//    testMutualInductanceZAxisMTScaling(16);
//    testMutualInductanceGeneralMTScaling(16);
//   testAmpereForceGeneralMTScaling(16);

//    testMutualInductanceGeneralEdgeCases();
//    testMutualInductanceZAxis();

//    testPerformanceForVariousCoilTypes(56243);
//    testMutualInductanceGeneralConway();

    testMutualInductanceGeneralParallelAxes();

//

    return 0;
}
