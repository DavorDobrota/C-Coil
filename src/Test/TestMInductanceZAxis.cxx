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
    for (int i = 1; i <= 15; ++i)
    {
        auto precisionFactor = PrecisionFactor(double(i));
        printf("precision = %.1f\n", double(i));

        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil1, precisionFactor, CPU_ST, true, false
        );
        printf("%d %d %d | %d %d %d\n",
               args.primaryPrecision.lengthBlockCount * args.primaryPrecision.lengthIncrementCount,
               args.primaryPrecision.thicknessBlockCount * args.primaryPrecision.thicknessIncrementCount,
               args.primaryPrecision.angularBlockCount * args.primaryPrecision.angularIncrementCount,
               args.secondaryPrecision.lengthBlockCount * args.secondaryPrecision.lengthIncrementCount,
               args.secondaryPrecision.thicknessBlockCount * args.secondaryPrecision.thicknessIncrementCount,
               args.secondaryPrecision.angularBlockCount * args.secondaryPrecision.angularIncrementCount
        );
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil2, precisionFactor, CPU_ST, true, false
        );
        printf("%d %d %d | %d %d %d\n",
               args.primaryPrecision.lengthBlockCount * args.primaryPrecision.lengthIncrementCount,
               args.primaryPrecision.thicknessBlockCount * args.primaryPrecision.thicknessIncrementCount,
               args.primaryPrecision.angularBlockCount * args.primaryPrecision.angularIncrementCount,
               args.secondaryPrecision.lengthBlockCount * args.secondaryPrecision.lengthIncrementCount,
               args.secondaryPrecision.thicknessBlockCount * args.secondaryPrecision.thicknessIncrementCount,
               args.secondaryPrecision.angularBlockCount * args.secondaryPrecision.angularIncrementCount
        );
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil3, precisionFactor, CPU_ST, true, false
        );
        printf("%d %d %d | %d %d %d\n",
               args.primaryPrecision.lengthBlockCount * args.primaryPrecision.lengthIncrementCount,
               args.primaryPrecision.thicknessBlockCount * args.primaryPrecision.thicknessIncrementCount,
               args.primaryPrecision.angularBlockCount * args.primaryPrecision.angularIncrementCount,
               args.secondaryPrecision.lengthBlockCount * args.secondaryPrecision.lengthIncrementCount,
               args.secondaryPrecision.thicknessBlockCount * args.secondaryPrecision.thicknessIncrementCount,
               args.secondaryPrecision.angularBlockCount * args.secondaryPrecision.angularIncrementCount
        );
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil4, precisionFactor, CPU_ST, true, false
        );
        printf("%d %d %d | %d %d %d\n",
               args.primaryPrecision.lengthBlockCount * args.primaryPrecision.lengthIncrementCount,
               args.primaryPrecision.thicknessBlockCount * args.primaryPrecision.thicknessIncrementCount,
               args.primaryPrecision.angularBlockCount * args.primaryPrecision.angularIncrementCount,
               args.secondaryPrecision.lengthBlockCount * args.secondaryPrecision.lengthIncrementCount,
               args.secondaryPrecision.thicknessBlockCount * args.secondaryPrecision.thicknessIncrementCount,
               args.secondaryPrecision.angularBlockCount * args.secondaryPrecision.angularIncrementCount
        );
        printf("\n");
    }

    printf("GPU z axis\n");
    for (int i = 1; i <= 15; ++i)
    {
        auto precisionFactor = PrecisionFactor(double(i));
        printf("precision = %.1f\n", double(i));

        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil1, precisionFactor, GPU, true, false
        );
        printf("%d %d %d | %d %d %d\n",
               args.primaryPrecision.lengthBlockCount * args.primaryPrecision.lengthIncrementCount,
               args.primaryPrecision.thicknessBlockCount * args.primaryPrecision.thicknessIncrementCount,
               args.primaryPrecision.angularBlockCount * args.primaryPrecision.angularIncrementCount,
               args.secondaryPrecision.lengthBlockCount * args.secondaryPrecision.lengthIncrementCount,
               args.secondaryPrecision.thicknessBlockCount * args.secondaryPrecision.thicknessIncrementCount,
               args.secondaryPrecision.angularBlockCount * args.secondaryPrecision.angularIncrementCount
        );
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil2, precisionFactor, GPU, true, false
        );
        printf("%d %d %d | %d %d %d\n",
               args.primaryPrecision.lengthBlockCount * args.primaryPrecision.lengthIncrementCount,
               args.primaryPrecision.thicknessBlockCount * args.primaryPrecision.thicknessIncrementCount,
               args.primaryPrecision.angularBlockCount * args.primaryPrecision.angularIncrementCount,
               args.secondaryPrecision.lengthBlockCount * args.secondaryPrecision.lengthIncrementCount,
               args.secondaryPrecision.thicknessBlockCount * args.secondaryPrecision.thicknessIncrementCount,
               args.secondaryPrecision.angularBlockCount * args.secondaryPrecision.angularIncrementCount
        );
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil3, precisionFactor, GPU, true, false
        );
        printf("%d %d %d | %d %d %d\n",
               args.primaryPrecision.lengthBlockCount * args.primaryPrecision.lengthIncrementCount,
               args.primaryPrecision.thicknessBlockCount * args.primaryPrecision.thicknessIncrementCount,
               args.primaryPrecision.angularBlockCount * args.primaryPrecision.angularIncrementCount,
               args.secondaryPrecision.lengthBlockCount * args.secondaryPrecision.lengthIncrementCount,
               args.secondaryPrecision.thicknessBlockCount * args.secondaryPrecision.thicknessIncrementCount,
               args.secondaryPrecision.angularBlockCount * args.secondaryPrecision.angularIncrementCount
        );
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil4, precisionFactor, GPU, true, false
        );
        printf("%d %d %d | %d %d %d\n",
               args.primaryPrecision.lengthBlockCount * args.primaryPrecision.lengthIncrementCount,
               args.primaryPrecision.thicknessBlockCount * args.primaryPrecision.thicknessIncrementCount,
               args.primaryPrecision.angularBlockCount * args.primaryPrecision.angularIncrementCount,
               args.secondaryPrecision.lengthBlockCount * args.secondaryPrecision.lengthIncrementCount,
               args.secondaryPrecision.thicknessBlockCount * args.secondaryPrecision.thicknessIncrementCount,
               args.secondaryPrecision.angularBlockCount * args.secondaryPrecision.angularIncrementCount
        );
        printf("\n");
    }

    printf("Pure GPU z axis\n");
    for (int i = 1; i <= 15; ++i)
    {
        auto precisionFactor = PrecisionFactor(double(i));
        printf("precision = %.1f\n", double(i));

        args = CoilPairArguments::getAppropriateCoilPairArguments(
                coil1, coil1, precisionFactor, GPU, true, true
        );
        printf("%d %d %d | %d %d %d\n",
               args.primaryPrecision.lengthBlockCount * args.primaryPrecision.lengthIncrementCount,
               args.primaryPrecision.thicknessBlockCount * args.primaryPrecision.thicknessIncrementCount,
               args.primaryPrecision.angularBlockCount * args.primaryPrecision.angularIncrementCount,
               args.secondaryPrecision.lengthBlockCount * args.secondaryPrecision.lengthIncrementCount,
               args.secondaryPrecision.thicknessBlockCount * args.secondaryPrecision.thicknessIncrementCount,
               args.secondaryPrecision.angularBlockCount * args.secondaryPrecision.angularIncrementCount
        );
        args = CoilPairArguments::getAppropriateCoilPairArguments(
                coil1, coil2, precisionFactor, GPU, true, true
        );
        printf("%d %d %d | %d %d %d\n",
               args.primaryPrecision.lengthBlockCount * args.primaryPrecision.lengthIncrementCount,
               args.primaryPrecision.thicknessBlockCount * args.primaryPrecision.thicknessIncrementCount,
               args.primaryPrecision.angularBlockCount * args.primaryPrecision.angularIncrementCount,
               args.secondaryPrecision.lengthBlockCount * args.secondaryPrecision.lengthIncrementCount,
               args.secondaryPrecision.thicknessBlockCount * args.secondaryPrecision.thicknessIncrementCount,
               args.secondaryPrecision.angularBlockCount * args.secondaryPrecision.angularIncrementCount
        );
        args = CoilPairArguments::getAppropriateCoilPairArguments(
                coil1, coil3, precisionFactor, GPU, true, true
        );
        printf("%d %d %d | %d %d %d\n",
               args.primaryPrecision.lengthBlockCount * args.primaryPrecision.lengthIncrementCount,
               args.primaryPrecision.thicknessBlockCount * args.primaryPrecision.thicknessIncrementCount,
               args.primaryPrecision.angularBlockCount * args.primaryPrecision.angularIncrementCount,
               args.secondaryPrecision.lengthBlockCount * args.secondaryPrecision.lengthIncrementCount,
               args.secondaryPrecision.thicknessBlockCount * args.secondaryPrecision.thicknessIncrementCount,
               args.secondaryPrecision.angularBlockCount * args.secondaryPrecision.angularIncrementCount
        );
        args = CoilPairArguments::getAppropriateCoilPairArguments(
                coil1, coil4, precisionFactor, GPU, true, true
        );
        printf("%d %d %d | %d %d %d\n",
               args.primaryPrecision.lengthBlockCount * args.primaryPrecision.lengthIncrementCount,
               args.primaryPrecision.thicknessBlockCount * args.primaryPrecision.thicknessIncrementCount,
               args.primaryPrecision.angularBlockCount * args.primaryPrecision.angularIncrementCount,
               args.secondaryPrecision.lengthBlockCount * args.secondaryPrecision.lengthIncrementCount,
               args.secondaryPrecision.thicknessBlockCount * args.secondaryPrecision.thicknessIncrementCount,
               args.secondaryPrecision.angularBlockCount * args.secondaryPrecision.angularIncrementCount
        );
        printf("\n");
    }

}
