#include <iostream>
#include "OldCoil.h"
#include "Polynomial.h"
#include "Coil.h"
#include "Test.h"

extern thread_pool tp;

int main()
{

//    testPerformanceCPU_ST();
    testPerformanceForComputeAll(20000, 10, 32);

    return 0;
}
