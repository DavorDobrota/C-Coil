#include <iostream>

#include "OldCoil.h"
#include "Coil.h"
#include "Test.h"

extern thread_pool tp;

int main()
{

 //   testCoilAmpereForceForFilamentsZAxis();
//    testCoilAmpereForceGeneralCase();

    testPerformanceForComputeAll(PrecisionFactor(7.0), 80000, 10, 16);

//testMethodPrecisionCompareCPUvsGPU();

//    for (int i = 1; i <= 16; ++i)
//    {
//        printf("%d threads:\n", i);
//        testCoilMutualInductanceGeneralPerformance(CPU_MT, i);
//        printf("\n");
//    }


//    testCoilMutualInductanceZAxis();
//    testCoilMutualInductanceGeneralForZAxis(CPU_MT, 8);

//    testCoilGradientTensor();


    return 0;
}
