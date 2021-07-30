#include <iostream>
#include "OldCoil.h"
#include "Polynomial.h"
#include "Coil.h"
#include "Test.h"

extern thread_pool tp;

int main()
{

    testPerformanceForComputeAll(60000, 5, 32);

    return 0;
}
