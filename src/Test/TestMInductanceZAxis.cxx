#include "Test.h"
#include "Coil.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>


void testMInductanceZAxisArgumentGeneration()
{
    Coil coil1 = Coil(0.05, 0.1, 0.1, 100);
    Coil coil2 = Coil(0.05, 0.1, 0.0, 10);
    Coil coil3 = Coil(0.05, 0.0, 0.1, 10);
    Coil coil4 = Coil(0.05, 0.0, 0.0, 1);

    CoilPairArguments args;

    printf("Testing mutual inductance z-axis argument generation\n\n");

    printf("CPU z axis\n");
    for (double i = 1.0; i <= 9.0; i += 1.0)
    {
        printf("precision = %.1f\n", i);
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil1, PrecisionFactor(i), CPU_ST, true);
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil2, PrecisionFactor(i), CPU_ST, true);
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil3, PrecisionFactor(i), CPU_ST, true);
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil4, PrecisionFactor(i), CPU_ST, true);
        printf("\n");
    }

    printf("GPU z axis\n");
    for (double i = 1.0; i <= 9.0; i += 1.0)
    {
        printf("precision = %.1f\n", i);
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil1, PrecisionFactor(i), GPU, true);
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil2, PrecisionFactor(i), GPU, true);
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil3, PrecisionFactor(i), GPU, true);
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil4, PrecisionFactor(i), GPU, true);
        printf("\n");
    }

}
