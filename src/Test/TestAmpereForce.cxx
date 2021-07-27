#include "Test.h"
#include "Coil.h"
#include "OldCoil.h"

#include <cstdio>
#include <cmath>


void testCoilAmpereForceZAxis()
{
    Coil prim1 = Coil(0.03, 0.03, 0.12, 3600, PrecisionFactor(6.0), 16);
    Coil sec1 = Coil(0.02, 0.025, 0.04, 1000, PrecisionFactor(6.0), 16);

    for (int i = 0; i < 500; ++i)
    {
        printf("%.15f\n",
               Coil::computeAmpereForceZAxis(prim1, sec1, 0.08 + i*0.001, PrecisionFactor(8.0), CPU_MT));
    }
    printf("\n");
}

void testCoilAmpereForceGeneralForZAxis()
{
    Coil prim1 = Coil(0.03, 0.03, 0.12, 3600, PrecisionFactor(6.0), 16);
    Coil sec1 = Coil(0.02, 0.025, 0.04, 1000, PrecisionFactor(6.0), 16);

    std::vector<double> tempVector;

    for (int i = 0; i < 100; ++i)
    {
        tempVector = Coil::computeAmpereForceGeneral(prim1, sec1, 0.08 + i*0.001, 1e-18,
                                                     PrecisionFactor(8.0), CPU_MT);
        printf("%.15f %.15f %.15f %.15f %.15f %.15f\n",
               tempVector[0], tempVector[1], tempVector[2], tempVector[3], tempVector[4], tempVector[5]);
    }
    printf("\n");
}