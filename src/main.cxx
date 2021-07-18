#include <iostream>
#include "../include/OldCoil.h"
#include "../include/Polynomial.h"
#include "../include/Coil.h"
#include "test.h"

extern thread_pool tp;

int main(){

    testCoilMutualInductanceGeneralThinCoilAndFilament();
    testCoilMutualInductanceGeneralThinCoilAndThinCoil();
    testCoilMutualInductanceGeneralPancakeAndPancake();
    testCoilMutualInductanceGeneralRectangularAndFilament();

    return 0;
}
