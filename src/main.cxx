#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>

#include "Benchmark.h"
#include "Coil.h"
#include "CoilGroup.h"
#include "Test.h"
#include "Compare.h"


#pragma clang diagnostic push
#pragma ide diagnostic ignored "Simplify"
int main()
{

//    Benchmark::coilGroupComputeAllFieldsGPU(100, 1'048'576);
    Test::testCoilGroupFieldsMTD();
    return 0;
}
#pragma clang diagnostic pop
