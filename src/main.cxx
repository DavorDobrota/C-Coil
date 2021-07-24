#include <iostream>
#include "OldCoil.h"
#include "Polynomial.h"
#include "Coil.h"
#include "Test.h"

extern thread_pool tp;

int main()
{
//    testCoilMutualInductanceZAxis();
//    testCoilMutualInductanceZAxisPerformance();

//    testCoilMutualInductanceZAxisArgumentGeneration();

    testCoilSelfInductance();

    testOldCoilSelfInductance();


//    testCoilMutualInductanceGeneralThinCoilAndFilament();
//    testCoilMutualInductanceGeneralThinCoilAndThinCoil();
//    testCoilMutualInductanceGeneralPancakeAndPancake();
//    testCoilMutualInductanceGeneralRectangularAndFilament();


    return 0;
}
