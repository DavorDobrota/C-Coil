#include <iostream>

#include "OldCoil.h"
#include "Coil.h"
#include "Test.h"

extern thread_pool tp;

int main()
{

//    testCoilMutualInductanceZAxisDifferentGeometries();
//
//    testCoilMutualInductanceGeneralThinCoilAndFilament();
//    testCoilMutualInductanceGeneralThinCoilAndThinCoil();
//    testCoilMutualInductanceGeneralPancakeAndPancake();
//    testCoilMutualInductanceGeneralRectangularAndFilament();
//
//    testCoilSelfInductance();

    testCoilPositionAndRotation();

//    testPerformanceForComputeAll(PrecisionFactor(6.0), 50000, 10, 16);

    return 0;
}
