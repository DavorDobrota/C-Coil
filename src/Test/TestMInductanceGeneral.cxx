#include "Test.h"
#include "Coil.h"
#include "ctpl.h"

#include <cstdio>
#include <cmath>

void testCoilMutualInductanceGeneralForZAxis(ComputeMethod method, int nThreads)
{
    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);

    primary.setThreadCount(nThreads);
    secondary.setThreadCount(nThreads);

    printf("%.20f\n\n", Coil::computeMutualInductance(primary, secondary, 0.2, 1e-15));

    FILE *input = fopen("values.txt", "r");
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

        for (double i = 1.0; i <= 8.0; i += 1.0)
        {
            temp = Coil::computeMutualInductance(prim, sec, distance, 1e-15, PrecisionFactor(i), method);
            printf("%.18f\n", temp);
            fprintf(output, "%.20f\t", temp);
        }

        printf("====================================================================================\n");
        fprintf(output, "\n");
    }
    fclose(input);
    fclose(output);
}

void testCoilMutualInductanceGeneralPerformance(ComputeMethod method, int nThreads)
{
    using namespace std::chrono;

    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);

    primary.setThreadCount(nThreads);

    int nOps = 1024;
    double temp;

    printf("Expected execution time for one MInductance general case calculation of specified precision\n");

    for (int i = 1; i <= 9; ++i)
    {
        int currentOperations = nOps / (int) pow(2, i);

        high_resolution_clock::time_point begin_time = high_resolution_clock::now();
        for (int j = 0; j < currentOperations; ++j)
            temp = Coil::computeMutualInductance(primary, secondary, 0.2, 1e-10, PrecisionFactor(i), method);
        double interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("precisionFactor(%.1f) : %6.2f ms/op\n", (double) i, 1'000.0 * interval / currentOperations);

    }
}

void testCoilMutualInductanceGeneralMTScaling(int maxThreads)
{
    printf("Performance comparison between different numbers of threads:\n");

    printf(" -> single thread:\n");
    testCoilMutualInductanceGeneralPerformance(CPU_ST);
    printf("\n");

    for (int i = 2; i <= maxThreads; ++i)
    {
        printf(" -> %2d threads:\n", i);
        testCoilMutualInductanceGeneralPerformance(CPU_MT, i);
        printf("\n");
    }
}


void testCoilMutualInductanceGeneralArgumentGeneration()
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

void testCoilMutualInductanceGeneralGraphs()
{
    FILE *output = fopen("output.txt", "w");

    Coil prim = Coil(0.01022, 0.011, 0.0022, 20, PrecisionFactor(6.0), 16);
    Coil sec = Coil(0.01022, 0.011, 0.0022, 20, PrecisionFactor(6.0), 16);

    double temp;

    for (int i = 0; i <= 200; ++i)
        fprintf(output, "%f\t%.10g\n", i * 0.0005,
                Coil::computeMutualInductance(prim, sec, i * 0.0005,
                                              PrecisionFactor(6.0), CPU_MT));
    fprintf(output,"\n");
    printf("1/6 tests Done\n");

    for (int i = 0; i <= 200; ++i)
        fprintf(output, "%f\t%.10g\n", i * 0.001,
               Coil::computeMutualInductance(prim, sec, 0.02, i * 0.001,
                                             PrecisionFactor(6.0), CPU_MT));
    fprintf(output,"\n");
    printf("2/6 tests done\n");

    for (int i = 0; i <= 200; ++i)
        fprintf(output, "%f\t%.10g\n", M_PI/100 * i,
               Coil::computeMutualInductance(prim, sec, 0.04, 0.0, M_PI/50 * i,
                                             PrecisionFactor(6.0), CPU_MT));
    fprintf(output,"\n");
    printf("3/6 tests done\n");

    for (int i = 0; i <= 100; ++i)
    {
        for (int j = 0; j <= 100; ++j)
            fprintf(output, "%.10g\t",
                    Coil::computeMutualInductance(prim, sec, i * 0.001, j * 0.001,
                                                  PrecisionFactor(6.0), CPU_MT));
        fprintf(output,"\n");
        printf("test progress: %d / 100\n", i);
    }
    fprintf(output,"\n");
    printf("4/6 tests done\n");

    for (int i = 0; i <= 100; ++i)
    {
        for (int j = 0; j <= 100; ++j)
            fprintf(output, "%.10g\t",
                    Coil::computeMutualInductance(prim, sec, 0.04, i * 0.002, -M_PI/50 * j,
                                                  PrecisionFactor(6.0), CPU_MT));
        fprintf(output,"\n");
        printf("test progress: %d / 100\n", i);
    }
    fprintf(output,"\n");
    printf("5/6 tests done\n");

    for (int i = 0; i <= 100; ++i)
    {
        for (int j = 0; j <= 100; ++j)
            fprintf(output, "%.10g\t",
                    Coil::computeMutualInductance(prim, sec, 0.025 + 0.001 * i, 0.0, -M_PI/50 * j,
                                                  PrecisionFactor(6.0), CPU_MT));
        fprintf(output,"\n");
        printf("test progress: %d / 100\n", i);
    }
    fprintf(output,"\n");
    printf("6/6 tests done\n\n");

    fclose(output);
}
