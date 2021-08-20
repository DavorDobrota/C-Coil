#include <iostream>

#include "OldCoil.h"
#include "Coil.h"
#include "Test.h"

extern thread_pool tp;

int main()
{
//    testCoilMutualInductanceGeneralThinCoilAndFilament();
//    testCoilMutualInductanceGeneralThinCoilAndThinCoil();
//    testCoilMutualInductanceGeneralPancakeAndPancake();
//    testCoilMutualInductanceGeneralRectangularAndFilament();

    testCoilSelfInductance();

    return 0;
}
