#include "Test.h"

#include <cstdio>
#include <chrono>


void Test::testMutualInductanceGeneralForZAxis(ComputeMethod computeMethod)
{
    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);

    primary.setThreadCount(g_defaultThreadCount);
    secondary.setThreadCount(g_defaultThreadCount);
    primary.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.0), 0.0, 0.0);
    secondary.setPositionAndOrientation(vec3::Vector3(1e-10, 0.0, 0.2), 0.0, 0.0);

    printf("%.20f\n\n", Coil::computeMutualInductance(primary, secondary));

    FILE *input = fopen("../data/values_MInductance_zAxis.txt", "r");
    FILE *output = fopen("output.txt", "w");

    double Rt1, at1, bt1; int Nt1;
    double Rt2, at2, bt2; int Nt2;
    double distance;
    double temp;

    while (fscanf(input, "%lf %lf %lf %d %lf %lf %lf %d %lf", &Rt1, &at1, &bt1, &Nt1, &Rt2, &at2, &bt2, &Nt2, &distance) == 9)
    {
        printf("%f %f %f %d %f %f %f %d %f\n", Rt1, at1, bt1, Nt1, Rt2, at2, bt2, Nt2, distance);

        Coil prim = Coil(Rt1, at1, bt1, Nt1);
        Coil sec = Coil(Rt2, at2, bt2, Nt2);
        prim.setPositionAndOrientation(vec3::Vector3(0.1, 0.0, 0.0), 0.0, 0.0);
        sec.setPositionAndOrientation(vec3::Vector3(0.1, 0.0, distance), 0.0, 0.0);


        for (int i = 1.0; i <= 8; ++i)
        {
            temp = Coil::computeMutualInductance(prim, sec, PrecisionFactor(i), computeMethod);
            printf("%.18f\n", temp);
            fprintf(output, "%.20f\t", temp);
        }

        printf("====================================================================================\n");
        fprintf(output, "\n");
    }
    fclose(input);
    fclose(output);
}

void Test::testMInductanceGeneralArgumentGeneration()
{
    using namespace std::chrono;

    high_resolution_clock::time_point begin_time;
    double interval;

    Coil coil1 = Coil(0.05, 0.1, 0.1, 100);
    Coil coil2 = Coil(0.05, 0.1, 0.1, 100);
    Coil coil3 = Coil(0.05, 0.1, 0.0, 10);
    Coil coil4 = Coil(0.05, 0.0, 0.1, 10);
    Coil coil5 = Coil(0.05, 0.0, 0.0, 1);

    vec3::Vector3 defaultVec = vec3::Vector3(0.1, 0.0, 0.2);

    coil2.setPositionAndOrientation(defaultVec);
    coil3.setPositionAndOrientation(defaultVec);
    coil4.setPositionAndOrientation(defaultVec);
    coil5.setPositionAndOrientation(defaultVec);

    CoilPairArguments args;

    printf("Testing mutual inductance general argument generation\n\n");

    printf("CPU general\n");
    for (int i = 1; i <= 15; ++i)
    {
        auto precisionFactor = PrecisionFactor(double(i));
        printf("precision = %.1f\n", double(i));

        begin_time = high_resolution_clock::now();
        args = CoilPairArguments::getAppropriateCoilPairArguments(
            coil1, coil2, precisionFactor, CPU_ST, false, false
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
                coil1, coil3, precisionFactor, CPU_ST, false, false
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
                coil1, coil4, precisionFactor, CPU_ST, false, false
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
                coil1, coil5, precisionFactor, CPU_ST, false, false
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

    printf("GPU general\n");
    for (int i = 1; i <= 15; ++i)
    {
        auto precisionFactor = PrecisionFactor(double(i));
        printf("precision = %.1f\n", double(i));

        begin_time = high_resolution_clock::now();
        args = CoilPairArguments::getAppropriateCoilPairArguments(
                coil1, coil2, precisionFactor, GPU, false, false
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
                coil1, coil3, precisionFactor, GPU, false, false
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
                coil1, coil4, precisionFactor, GPU, false, false
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
                coil1, coil5, precisionFactor, GPU, false, false
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

    printf("Pure GPU general\n");
    for (int i = 1; i <= 15; ++i)
    {
        auto precisionFactor = PrecisionFactor(double(i));
        printf("precision = %.1f\n", double(i));

        begin_time = high_resolution_clock::now();
        args = CoilPairArguments::getAppropriateCoilPairArguments(
                coil1, coil2, precisionFactor, GPU, false, true
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
                coil1, coil3, precisionFactor, GPU, false, true
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
                coil1, coil4, precisionFactor, GPU, false, true
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
                coil1, coil5, precisionFactor, GPU, false, true
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
