#include "Test.h"

#include <cstdio>
#include <chrono>


void Test::testMInductanceZAxisArgumentGeneration()
{
    using namespace std::chrono;

    high_resolution_clock::time_point begin_time;
    double interval;

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

        begin_time = high_resolution_clock::now();
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil1, precisionFactor, CPU_ST, true, false
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("%d %d %d | %d %d %d Time: %.0f ns\n",
               args.primaryPrecision.angularBlocks * args.primaryPrecision.angularIncrements,
               args.primaryPrecision.thicknessBlocks * args.primaryPrecision.thicknessIncrements,
               args.primaryPrecision.lengthBlocks * args.primaryPrecision.lengthIncrements,
               args.secondaryPrecision.angularBlocks * args.secondaryPrecision.angularIncrements,
               args.secondaryPrecision.thicknessBlocks * args.secondaryPrecision.thicknessIncrements,
               args.secondaryPrecision.lengthBlocks * args.secondaryPrecision.lengthIncrements,
               1e9 * interval
        );

        begin_time = high_resolution_clock::now();
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil2, precisionFactor, CPU_ST, true, false
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("%d %d %d | %d %d %d Time: %.0f ns\n",
               args.primaryPrecision.angularBlocks * args.primaryPrecision.angularIncrements,
               args.primaryPrecision.thicknessBlocks * args.primaryPrecision.thicknessIncrements,
               args.primaryPrecision.lengthBlocks * args.primaryPrecision.lengthIncrements,
               args.secondaryPrecision.angularBlocks * args.secondaryPrecision.angularIncrements,
               args.secondaryPrecision.thicknessBlocks * args.secondaryPrecision.thicknessIncrements,
               args.secondaryPrecision.lengthBlocks * args.secondaryPrecision.lengthIncrements,
               1e9 * interval
        );

        begin_time = high_resolution_clock::now();
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil3, precisionFactor, CPU_ST, true, false
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("%d %d %d | %d %d %d Time: %.0f ns\n",
               args.primaryPrecision.angularBlocks * args.primaryPrecision.angularIncrements,
               args.primaryPrecision.thicknessBlocks * args.primaryPrecision.thicknessIncrements,
               args.primaryPrecision.lengthBlocks * args.primaryPrecision.lengthIncrements,
               args.secondaryPrecision.angularBlocks * args.secondaryPrecision.angularIncrements,
               args.secondaryPrecision.thicknessBlocks * args.secondaryPrecision.thicknessIncrements,
               args.secondaryPrecision.lengthBlocks * args.secondaryPrecision.lengthIncrements,
               1e9 * interval
        );

        begin_time = high_resolution_clock::now();
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil4, precisionFactor, CPU_ST, true, false
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("%d %d %d | %d %d %d Time: %.0f ns\n",
               args.primaryPrecision.angularBlocks * args.primaryPrecision.angularIncrements,
               args.primaryPrecision.thicknessBlocks * args.primaryPrecision.thicknessIncrements,
               args.primaryPrecision.lengthBlocks * args.primaryPrecision.lengthIncrements,
               args.secondaryPrecision.angularBlocks * args.secondaryPrecision.angularIncrements,
               args.secondaryPrecision.thicknessBlocks * args.secondaryPrecision.thicknessIncrements,
               args.secondaryPrecision.lengthBlocks * args.secondaryPrecision.lengthIncrements,
               1e9 * interval
        );

        printf("\n");
    }

    printf("GPU z axis\n");
    for (int i = 1; i <= 15; ++i)
    {
        auto precisionFactor = PrecisionFactor(double(i));
        printf("precision = %.1f\n", double(i));

        begin_time = high_resolution_clock::now();
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil1, precisionFactor, GPU, true, false
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("%d %d %d | %d %d %d Time: %.0f ns\n",
               args.primaryPrecision.angularBlocks * args.primaryPrecision.angularIncrements,
               args.primaryPrecision.thicknessBlocks * args.primaryPrecision.thicknessIncrements,
               args.primaryPrecision.lengthBlocks * args.primaryPrecision.lengthIncrements,
               args.secondaryPrecision.angularBlocks * args.secondaryPrecision.angularIncrements,
               args.secondaryPrecision.thicknessBlocks * args.secondaryPrecision.thicknessIncrements,
               args.secondaryPrecision.lengthBlocks * args.secondaryPrecision.lengthIncrements,
               1e9 * interval
        );

        begin_time = high_resolution_clock::now();
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil2, precisionFactor, GPU, true, false
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("%d %d %d | %d %d %d Time: %.0f ns\n",
               args.primaryPrecision.angularBlocks * args.primaryPrecision.angularIncrements,
               args.primaryPrecision.thicknessBlocks * args.primaryPrecision.thicknessIncrements,
               args.primaryPrecision.lengthBlocks * args.primaryPrecision.lengthIncrements,
               args.secondaryPrecision.angularBlocks * args.secondaryPrecision.angularIncrements,
               args.secondaryPrecision.thicknessBlocks * args.secondaryPrecision.thicknessIncrements,
               args.secondaryPrecision.lengthBlocks * args.secondaryPrecision.lengthIncrements,
               1e9 * interval
        );

        begin_time = high_resolution_clock::now();
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil3, precisionFactor, GPU, true, false
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("%d %d %d | %d %d %d Time: %.0f ns\n",
               args.primaryPrecision.angularBlocks * args.primaryPrecision.angularIncrements,
               args.primaryPrecision.thicknessBlocks * args.primaryPrecision.thicknessIncrements,
               args.primaryPrecision.lengthBlocks * args.primaryPrecision.lengthIncrements,
               args.secondaryPrecision.angularBlocks * args.secondaryPrecision.angularIncrements,
               args.secondaryPrecision.thicknessBlocks * args.secondaryPrecision.thicknessIncrements,
               args.secondaryPrecision.lengthBlocks * args.secondaryPrecision.lengthIncrements,
               1e9 * interval
        );

        begin_time = high_resolution_clock::now();
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil4, precisionFactor, GPU, true, false
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("%d %d %d | %d %d %d Time: %.0f ns\n",
               args.primaryPrecision.angularBlocks * args.primaryPrecision.angularIncrements,
               args.primaryPrecision.thicknessBlocks * args.primaryPrecision.thicknessIncrements,
               args.primaryPrecision.lengthBlocks * args.primaryPrecision.lengthIncrements,
               args.secondaryPrecision.angularBlocks * args.secondaryPrecision.angularIncrements,
               args.secondaryPrecision.thicknessBlocks * args.secondaryPrecision.thicknessIncrements,
               args.secondaryPrecision.lengthBlocks * args.secondaryPrecision.lengthIncrements,
               1e9 * interval
        );

        printf("\n");
    }

    printf("Pure GPU z axis\n");
    for (int i = 1; i <= 15; ++i)
    {
        auto precisionFactor = PrecisionFactor(double(i));
        printf("precision = %.1f\n", double(i));

        begin_time = high_resolution_clock::now();
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil1, precisionFactor, GPU, true, true
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("%d %d %d | %d %d %d Time: %.0f ns\n",
               args.primaryPrecision.angularBlocks * args.primaryPrecision.angularIncrements,
               args.primaryPrecision.thicknessBlocks * args.primaryPrecision.thicknessIncrements,
               args.primaryPrecision.lengthBlocks * args.primaryPrecision.lengthIncrements,
               args.secondaryPrecision.angularBlocks * args.secondaryPrecision.angularIncrements,
               args.secondaryPrecision.thicknessBlocks * args.secondaryPrecision.thicknessIncrements,
               args.secondaryPrecision.lengthBlocks * args.secondaryPrecision.lengthIncrements,
               1e9 * interval
        );

        begin_time = high_resolution_clock::now();
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil2, precisionFactor, GPU, true, true
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("%d %d %d | %d %d %d Time: %.0f ns\n",
               args.primaryPrecision.angularBlocks * args.primaryPrecision.angularIncrements,
               args.primaryPrecision.thicknessBlocks * args.primaryPrecision.thicknessIncrements,
               args.primaryPrecision.lengthBlocks * args.primaryPrecision.lengthIncrements,
               args.secondaryPrecision.angularBlocks * args.secondaryPrecision.angularIncrements,
               args.secondaryPrecision.thicknessBlocks * args.secondaryPrecision.thicknessIncrements,
               args.secondaryPrecision.lengthBlocks * args.secondaryPrecision.lengthIncrements,
               1e9 * interval
        );

        begin_time = high_resolution_clock::now();
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil3, precisionFactor, GPU, true, true
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("%d %d %d | %d %d %d Time: %.0f ns\n",
               args.primaryPrecision.angularBlocks * args.primaryPrecision.angularIncrements,
               args.primaryPrecision.thicknessBlocks * args.primaryPrecision.thicknessIncrements,
               args.primaryPrecision.lengthBlocks * args.primaryPrecision.lengthIncrements,
               args.secondaryPrecision.angularBlocks * args.secondaryPrecision.angularIncrements,
               args.secondaryPrecision.thicknessBlocks * args.secondaryPrecision.thicknessIncrements,
               args.secondaryPrecision.lengthBlocks * args.secondaryPrecision.lengthIncrements,
               1e9 * interval
        );

        begin_time = high_resolution_clock::now();
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil4, precisionFactor, GPU, true, true
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("%d %d %d | %d %d %d Time: %.0f ns\n",
               args.primaryPrecision.angularBlocks * args.primaryPrecision.angularIncrements,
               args.primaryPrecision.thicknessBlocks * args.primaryPrecision.thicknessIncrements,
               args.primaryPrecision.lengthBlocks * args.primaryPrecision.lengthIncrements,
               args.secondaryPrecision.angularBlocks * args.secondaryPrecision.angularIncrements,
               args.secondaryPrecision.thicknessBlocks * args.secondaryPrecision.thicknessIncrements,
               args.secondaryPrecision.lengthBlocks * args.secondaryPrecision.lengthIncrements,
               1e9 * interval
        );

        printf("\n");
    }

}
