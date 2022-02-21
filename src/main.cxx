#include <cstdio>
#include <cmath>

#include "Coil.h"
#include "CoilGroup.h"
#include "Test.h"


int main()
{

//    testMutualInductanceZAxisMTScaling(16);
//    testMutualInductanceGeneralGraphs();
    testMutualInductanceGeneralComputeAll(PrecisionFactor(8.0), 24);
    return 0;
}
