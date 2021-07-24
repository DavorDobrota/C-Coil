#include <iostream>
#include "OldCoil.h"
#include "Polynomial.h"
#include "Coil.h"
#include "Test.h"

extern thread_pool tp;

int main()
{
//    testCoilMutualInductanceZAxis();

//    testPerformanceForComputeAll(10'000, 1, 16);
//    testCoilMutualInductanceZAxisPerformance();

//    testCoilMutualInductanceZAxisArgumentGeneration();

//    testCoilMutualInductanceZAxisPerformance();
    testCoilMutualInductanceZAxisPerformance(CPU_MT);


//    testCoilMutualInductanceGeneralThinCoilAndFilament();
//    testCoilMutualInductanceGeneralThinCoilAndThinCoil();
//    testCoilMutualInductanceGeneralPancakeAndPancake();
//    testCoilMutualInductanceGeneralRectangularAndFilament();


    return 0;
}
