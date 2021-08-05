#include <iostream>

#include "OldCoil.h"
#include "Coil.h"
#include "Test.h"

extern thread_pool tp;

int main()
{

 //   testCoilAmpereForceForFilamentsZAxis();
//    testCoilAmpereForceGeneralCase();

//    testPerformanceForComputeAll(60000, 20, 32);

//testMethodPrecisionCompareCPUvsGPU();

    testCoilMutualInductanceZAxisPerformance(CPU_MT);


    return 0;
}
