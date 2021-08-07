#include <iostream>

#include "OldCoil.h"
#include "Coil.h"
#include "Test.h"

extern thread_pool tp;

int main()
{

 //   testCoilAmpereForceForFilamentsZAxis();
//    testCoilAmpereForceGeneralCase();

    testPerformanceForComputeAll(PrecisionFactor(8.0), 60000, 10, 32);

//testMethodPrecisionCompareCPUvsGPU();

//    for (int i = 1; i <= 32; ++i)
//    {
//        printf("%d threads:\n", i);
//        testCoilMutualInductanceZAxisPerformance(CPU_MT, i);
//        printf("\n");
//    }

//    testCoilGradientTensor();


    return 0;
}
