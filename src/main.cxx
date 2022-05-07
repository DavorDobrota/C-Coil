#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>

#include "Benchmark.h"
#include "Coil.h"
#include "CoilGroup.h"
#include "Test.h"


#pragma clang diagnostic push
#pragma ide diagnostic ignored "Simplify"
int main()
{

//    benchComputeAllFields(PrecisionFactor(9.0), 10'000, 2, 8);
    benchComputeAllFieldsEveryCoilType(415'701, 8);

//    benchMInductanceZAxisMTScaling(12);
//    benchMInductanceGeneralMTScaling(16);
//    benchForceZAxisMTScaling(12);
//    benchForceGeneralMTScaling(16);
//    benchSelfInductance();
//
//    Coil coil = Coil(0.1, 0.1, 0.1, 10000);
//
//    for (int i = 0; i <= 140; ++i)
//    {
//        printf("%.15g\n", coil.computeAndSetSelfInductance(PrecisionFactor(1.0 + i * 0.1 )));
//    }

    return 0;
}
#pragma clang diagnostic pop
