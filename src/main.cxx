#include <cstdio>
#include <cmath>

#include "Coil.h"
#include "CoilGroup.h"
#include "Test.h"


int main()
{
//    testMutualInductanceZAxis();
    testSelfInductance();
    testSelfInductancePerformance();
    testMutualInductanceZAxisMTScaling();

    return 0;
}
