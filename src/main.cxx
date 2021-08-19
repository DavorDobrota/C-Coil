#include <iostream>

#include "OldCoil.h"
#include "Coil.h"
#include "Test.h"

extern thread_pool tp;

int main()
{
    testCoilMutualInductanceZAxis();
    testCoilMutualInductanceGeneralForZAxis(CPU_MT);

    return 0;
}
