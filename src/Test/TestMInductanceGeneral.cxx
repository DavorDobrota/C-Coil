#include "Test.h"
#include "Coil.h"
#include "ctpl.h"

#include <cstdio>
#include <cmath>

void testMutualInductanceGeneralForZAxis(ComputeMethod method, int nThreads)
{
    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);

    primary.setThreadCount(nThreads);
    secondary.setThreadCount(nThreads);
    primary.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.0), 0.0, 0.0);
    secondary.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 1e-10, 0.0, 0.2), 0.0, 0.0);


    printf("%.20f\n\n", Coil::computeMutualInductance(primary, secondary));

    FILE *input = fopen("values_MInductance_zAxis.txt", "r");
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
        prim.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.1, 0.0, 0.0), 0.0, 0.0);
        sec.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.1, 0.0, distance), 0.0, 0.0);


        for (int i = 1.0; i <= 8; ++i)
        {
            temp = Coil::computeMutualInductance(prim, sec, PrecisionFactor(i), method);
            printf("%.18f\n", temp);
            fprintf(output, "%.20f\t", temp);
        }

        printf("====================================================================================\n");
        fprintf(output, "\n");
    }
    fclose(input);
    fclose(output);
}

void testMutualInductanceGeneralPerformance(ComputeMethod method, int nThreads)
{
    using namespace std::chrono;

    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);

    primary.setThreadCount(nThreads);
    primary.setPositionAndOrientation();
    secondary.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.1, 0.0, 0.2));

    int nOps = 1024;
    double temp;

    printf("Expected execution time for one MInductance general case calculation of specified precision\n");

    for (int i = 1; i <= 9; ++i)
    {
        int currentOperations = nOps / (int) pow(2, i);

        high_resolution_clock::time_point begin_time = high_resolution_clock::now();
        for (int j = 0; j < currentOperations; ++j)
            temp = Coil::computeMutualInductance(primary, secondary, PrecisionFactor(i), method);
        double interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("precisionFactor(%.1f) : %6.2f ms/op\n", (double) i, 1'000.0 * interval / currentOperations);

    }
}

void testMutualInductanceGeneralMTScaling(int maxThreads)
{
    printf("Performance comparison between different numbers of threads:\n");

    printf(" -> single thread:\n");
    testMutualInductanceGeneralPerformance(CPU_ST);
    printf("\n");

    for (int i = 2; i <= maxThreads; ++i)
    {
        printf(" -> %2d threads:\n", i);
        testMutualInductanceGeneralPerformance(CPU_MT, i);
        printf("\n");
    }
}


void testMutualInductanceGeneralArgumentGeneration()
{
    Coil coil1 = Coil(0.05, 0.1, 0.1, 100);
    Coil coil2 = Coil(0.05, 0.1, 0.0, 10);
    Coil coil3 = Coil(0.05, 0.0, 0.1, 10);
    Coil coil4 = Coil(0.05, 0.0, 0.0, 1);

    CoilPairArguments args;

    printf("CPU general\n");
    for (double i = 1.0; i <= 8.0; i += 1.0)
    {
        printf("precision = %.1f\n", i);
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil2, PrecisionFactor(i));
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil3, PrecisionFactor(i));
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil4, PrecisionFactor(i));
        printf("\n");
    }

    printf("GPU general\n");
    for (double i = 1.0; i <= 8.0; i += 1.0)
    {
        printf("precision = %.1f\n", i);
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil2, PrecisionFactor(i), GPU);
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil3, PrecisionFactor(i), GPU);
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil4, PrecisionFactor(i), GPU);
        printf("\n");
    }
}

void testMutualInductanceGeneralGraphs()
{
    FILE *output = fopen("output.txt", "w");

    Coil prim = Coil(0.01022, 0.011, 0.0022, 20, PrecisionFactor(6.0), 12);
    Coil sec = Coil(0.01022, 0.011, 0.0022, 20, PrecisionFactor(6.0), 12);

    auto precision = PrecisionFactor(6.0);

    for (int i = 0; i <= 200; ++i)
    {
        sec.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, i * 0.0005));
        fprintf(output, "%f\t%.10g\n", i * 0.0005, Coil::computeMutualInductance(prim, sec, precision, CPU_MT));
    }
    fprintf(output,"\n");
    printf("1/6 tests Done\n");

    for (int i = 0; i <= 200; ++i)
    {
        sec.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, i * 0.001, 0.0, 0.02));
        fprintf(output, "%f\t%.10g\n", i * 0.001, Coil::computeMutualInductance(prim, sec,precision, CPU_MT));
    }
    fprintf(output,"\n");
    printf("2/6 tests done\n");

    for (int i = 0; i <= 200; ++i)
    {
        sec.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.04), M_PI/100 * i);
        fprintf(output, "%f\t%.10g\n", M_PI/100 * i, Coil::computeMutualInductance(prim, sec, precision, CPU_MT));
    }
    fprintf(output,"\n");
    printf("3/6 tests done\n");

    for (int i = 0; i <= 100; ++i)
    {
        for (int j = 0; j <= 100; ++j)
        {
            sec.setPositionAndOrientation(
                    vec3::CoordVector3(vec3::CARTESIAN, j * 0.001, 0.0, i * 0.001));
            fprintf(output, "%.10g\t", Coil::computeMutualInductance(prim, sec, precision, CPU_MT));
        }
        fprintf(output,"\n");
        printf("test progress: %d / 100\n", i);
    }
    fprintf(output,"\n");
    printf("4/6 tests done\n");

    for (int i = 0; i <= 100; ++i)
    {
        for (int j = 0; j <= 100; ++j)
        {
            sec.setPositionAndOrientation(
                    vec3::CoordVector3(vec3::CARTESIAN, i * 0.002, 0.0, 0.04), -M_PI/50 * j);
            fprintf(output, "%.10g\t", Coil::computeMutualInductance(prim, sec, precision, CPU_MT));
        }

        fprintf(output,"\n");
        printf("test progress: %d / 100\n", i);
    }
    fprintf(output,"\n");
    printf("5/6 tests done\n");

    for (int i = 0; i <= 100; ++i)
    {
        for (int j = 0; j <= 100; ++j)
        {
            sec.setPositionAndOrientation(
                    vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.025 + 0.001 * i), -M_PI/50 * j, 0.0);
            fprintf(output, "%.10g\t", Coil::computeMutualInductance(prim, sec, precision, CPU_MT));
        }
        fprintf(output,"\n");
        printf("test progress: %d / 100\n", i);
    }
    fprintf(output,"\n");
    printf("6/6 tests done\n\n");

    fclose(output);
}

void testMutualInductanceGeneralEdgeCases()
{
    Coil coil1 = Coil(0.03, 0.12, 0.12, 3600);
    Coil coil2 = Coil(0.03, 0.12, 0.12, 3600);
    coil2.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.12));
    for (int i = 1; i <= 15; ++i)
        printf("%.15g\n", Coil::computeMutualInductance(coil1, coil2, PrecisionFactor(i), CPU_MT));
    printf("\n");

    Coil coil3 = Coil(0.15, 0.12, 0.12, 3600);
    Coil coil4 = Coil(0.03, 0.12, 0.12, 3600);
    for (int i = 1; i <= 15; ++i)
        printf("%.15g\n", Coil::computeMutualInductance(coil3, coil4, PrecisionFactor(i), CPU_MT));
    printf("\n");

    Coil coil5 = Coil(0.03, 0.12, 0.12, 3600);
    Coil coil6 = Coil(0.03, 0.12, 0.12, 3600);
    coil6.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.00001, 0.0, 0.12));
    for (int i = 1; i <= 15; ++i)
        printf("%.15g\n", Coil::computeMutualInductance(coil5, coil6, PrecisionFactor(i), CPU_MT));
    printf("\n");

    Coil coil7 = Coil(3.0, 0.2, 0.1, 1250);
    Coil coil8 = Coil(2.8, 0.2, 0.1, 1250);
    for (int i = 1; i <= 15; ++i)
        printf("%.15g\n", Coil::computeMutualInductance(coil7, coil8, PrecisionFactor(i), CPU_MT));
    printf("\n");

    Coil coil9 = Coil(3.0, 0.2, 0.1, 1250);
    Coil coil10 = Coil(2.7, 0.2, 0.1, 1250);
    coil10.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.1, 0.0, 0.0));
    for (int i = 1; i <= 15; ++i)
        printf("%.15g\n", Coil::computeMutualInductance(coil9, coil10, PrecisionFactor(i), CPU_MT));
    printf("\n");
}