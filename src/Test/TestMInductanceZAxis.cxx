#include "Test.h"
#include "Coil.h"
#include "ctpl.h"

#include <cstdio>
#include <cmath>


void testMutualInductanceZAxis()
{
    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);

    primary.setPositionAndOrientation();
    secondary.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.2));

    printf("%.20f\n\n", Coil::computeMutualInductance(primary, secondary));

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

        prim.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.0));
        sec.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, distance));

        for (int i = 1; i <= 8; i++)
        {
            temp = Coil::computeMutualInductance(prim, sec, PrecisionFactor(i));
            printf("%.16g\n", temp);
            fprintf(output, "%.16g\t", temp);
        }

        printf("====================================================================================\n");
        fprintf(output, "\n");
    }

    fclose(input);
    fclose(output);
}

void testMutualInductanceZAxisArgumentGeneration()
{
    Coil coil1 = Coil(0.05, 0.1, 0.1, 100);
    Coil coil2 = Coil(0.05, 0.1, 0.0, 10);
    Coil coil3 = Coil(0.05, 0.0, 0.1, 10);
    Coil coil4 = Coil(0.05, 0.0, 0.0, 1);

    CoilPairArguments args;

    printf("CPU z axis\n");
    for (double i = 1.0; i <= 8.0; i += 1.0)
    {
        printf("precision = %.1f\n", i);
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil2, PrecisionFactor(i), CPU_ST, false);
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil3, PrecisionFactor(i), CPU_ST, false);
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil4, PrecisionFactor(i), CPU_ST, false);
        printf("\n");
    }

    printf("GPU z axis\n");
    for (double i = 1.0; i <= 8.0; i += 1.0)
    {
        printf("precision = %.1f\n", i);
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil2, PrecisionFactor(i), GPU, false);
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil3, PrecisionFactor(i), GPU, false);
        args = CoilPairArguments::getAppropriateCoilPairArguments(coil1, coil4, PrecisionFactor(i), GPU, false);
        printf("\n");
    }
}

void testCoilMutualInductanceZAxisPerformance(ComputeMethod method, int nThreads)
{
    using namespace std::chrono;

    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);

    primary.setThreadCount(nThreads);
    primary.setPositionAndOrientation();
    secondary.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.2));

    int nOps = 8192;
    double temp;

    printf("Expected execution time for one MInductance z-axis calculation of specified precision\n");

    for (int i = 1; i <= 9; ++i)
    {
        int currentOperations = nOps / (int) pow(2, i);

        high_resolution_clock::time_point begin_time = high_resolution_clock::now();
        for (int j = 0; j < currentOperations; ++j)
            temp = Coil::computeMutualInductance(primary, secondary, PrecisionFactor(i), method);
        double interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("precisionFactor(%.1f) : %6.3f ms/op\n", (double) i, 1'000.0 * interval / currentOperations);
    }
}

void testMutualInductanceZAxisMTScaling(int maxThreads)
{
    printf("Performance comparison between different numbers of threads:\n");

    printf(" -> single thread:\n");
    testCoilMutualInductanceZAxisPerformance(CPU_ST);
    printf("\n");

    for (int i = 2; i <= maxThreads; ++i)
    {
        printf(" -> %2d threads:\n", i);
        testCoilMutualInductanceZAxisPerformance(CPU_MT, i);
        printf("\n");
    }
}

void testSelfInductance()
{
    Coil coil1 = Coil(0.03, 0.03, 0.12, 3600, PrecisionFactor(), 12);
    Coil coil2 = Coil(0.03, 0.0, 0.12, 120);
    Coil coil3 = Coil(0.03, 0.03, 0.0, 30);
    Coil coil4 = Coil(0.03, 0.0, 0.0, 1);

    coil1.setThreadCount(8);

    for (int i = 1; i <= 12; ++i)
        printf("%.15g\n", coil1.computeAndSetSelfInductance(PrecisionFactor(i), CPU_MT));
    printf("\n");

    for (int i = 1; i <= 12; ++i)
        printf("%.15g\n", coil2.computeAndSetSelfInductance(PrecisionFactor(i)));
    printf("\n");

    for (int i = 1; i <= 12; ++i)
        printf("%.15g\n", coil3.computeAndSetSelfInductance(PrecisionFactor(i)));
    printf("\n");

    for (int i = 1; i <= 12; ++i)
        printf("%.15g\n", coil4.computeAndSetSelfInductance(PrecisionFactor(i)));
    printf("\n");
}