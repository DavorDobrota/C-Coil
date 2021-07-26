#include "Test.h"
#include "Coil.h"
#include "OldCoil.h"

#include <cstdio>
#include <cmath>


void testCoilAmpereForceZAxis()
{
    Coil prim1 = Coil(0.03, 0.03, 0.12, 3600, PrecisionFactor(6.0), 12);
    Coil sec1 = Coil(0.02, 0.025, 0.04, 1000, PrecisionFactor(6.0), 12);

    for (int i = 0; i < 500; ++i)
    {
        printf("%.15f\n",
               Coil::computeAmpereForce(prim1, sec1, 0.08 + i*0.001, PrecisionFactor(7.0), CPU_MT));
    }
    printf("\n");
}