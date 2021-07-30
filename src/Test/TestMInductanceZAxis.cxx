#include "Test.h"
#include "Coil.h"
#include "ctpl.h"

#include <cstdio>
#include <cmath>


void testCoilMutualInductanceZAxis()
{
    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);

    printf("%.20f\n\n", Coil::computeMutualInductance(primary, secondary, 0.2));

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
            temp = Coil::computeMutualInductance(prim, sec, distance, PrecisionFactor(i));
            printf("%.18f\n", temp);
            fprintf(output, "%.15g\t", temp);
        }

        printf("====================================================================================\n");
        fprintf(output, "\n");
    }

    fclose(input);
    fclose(output);
}

void testCoilMutualInductanceZAxisArgumentGeneration()
{
    Coil coil1 = Coil(0.05, 0.1, 0.1, 100);
    Coil coil2 = Coil(0.05, 0.1, 1e-15, 10);
    Coil coil3 = Coil(0.05, 1e-15, 0.1, 10);
    Coil coil4 = Coil(0.05, 1e-15, 1e-15, 1);

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

void testCoilMutualInductanceZAxisDifferentGeometries()
{
    Coil prim1 = Coil(0.03, 1e-15, 1e-15, 1);
    Coil sec1 = Coil(0.03, 1e-15, 1e-15, 1);
    printf("%.12g\n", Coil::computeMutualInductance(prim1, sec1, 0.2));

    Coil prim2 = Coil(0.03, 0.03, 1e-15, 30);
    Coil sec2 = Coil(0.03, 1e-15, 1e-15, 1);
    printf("%.12g\n", Coil::computeMutualInductance(prim2, sec2, 0.2));

    Coil prim3 = Coil(0.03, 1e-15, 1e-15, 1);
    Coil sec3 = Coil(0.03, 0.03, 1e-15, 30);
    printf("%.12g\n", Coil::computeMutualInductance(prim3, sec3, 0.2));

    Coil prim4 = Coil(0.03, 0.03, 1e-15, 30);
    Coil sec4 = Coil(0.03, 0.03, 1e-15, 30);
    printf("%.12g\n", Coil::computeMutualInductance(prim4, sec4, 0.2));

    Coil prim5 = Coil(0.03, 1e-15, 0.12, 120);
    Coil sec5 = Coil(0.03,  1e-15, 1e-15, 1);
    printf("%.12g\n", Coil::computeMutualInductance(prim5, sec5, 0.2));

    Coil prim6 = Coil(0.03, 1e-15, 1e-15, 1);
    Coil sec6 = Coil(0.03,  1e-15, 0.12, 120);
    printf("%.12g\n", Coil::computeMutualInductance(prim6, sec6, 0.2));

    Coil prim7 = Coil(0.03, 1e-15, 0.12, 120);
    Coil sec7 = Coil(0.03,  1e-15, 0.12, 120);
    printf("%.12g\n", Coil::computeMutualInductance(prim7, sec7, 0.2));

    Coil prim8 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec8 = Coil(0.03,  1e-15, 1e-15, 1);
    printf("%.12g\n", Coil::computeMutualInductance(prim8, sec8, 0.2));

    Coil prim9 = Coil(0.03,  1e-15, 1e-15, 1);
    Coil sec9 = Coil(0.03, 0.03, 0.12, 3600);
    printf("%.12g\n", Coil::computeMutualInductance(prim9, sec9, 0.2));

    Coil prim10 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec10 = Coil(0.03, 0.03, 0.12, 3600);
    printf("%.12g\n", Coil::computeMutualInductance(prim10, sec10, 0.2));

    Coil prim11 = Coil(0.03, 0.03, 1e-15, 30);
    Coil sec11 = Coil(0.03, 1e-15, 0.12, 120);
    printf("%.12g\n", Coil::computeMutualInductance(prim11, sec11, 0.2));

    Coil prim12 = Coil(0.03, 1e-15, 0.12, 120);
    Coil sec12 = Coil(0.03, 0.03, 1e-15, 30);
    printf("%.12g\n", Coil::computeMutualInductance(prim12, sec12, 0.2));

    Coil prim13 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec13 = Coil(0.03, 1e-15, 0.12, 120);
    printf("%.12g\n", Coil::computeMutualInductance(prim13, sec13, 0.2));

    Coil prim14 = Coil(0.03, 1e-15, 0.12, 120);
    Coil sec14 = Coil(0.03, 0.03, 0.12, 3600);
    printf("%.12g\n", Coil::computeMutualInductance(prim14, sec14, 0.2));

    Coil prim15 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec15 = Coil(0.03, 0.03, 1e-15, 30);
    printf("%.12g\n", Coil::computeMutualInductance(prim15, sec15, 0.2));

    Coil prim16 = Coil(0.03, 0.03, 1e-15, 30);
    Coil sec16 = Coil(0.03, 0.03, 0.12, 3600);
    printf("%.12g\n", Coil::computeMutualInductance(prim16, sec16, 0.2));

    printf("\n");
}

void testCoilMutualInductanceZAxisPerformance(ComputeMethod method)
{
    using namespace std::chrono;

    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);

    primary.setThreadCount(12);

    int nOps = 5120;
    int numIncrements[] = {78732, 147000, 263296, 547560, 1057500, 2247264, 4528384, 9168896};
    double temp;

    for (double i = 1.0; i <= 8.0; i += 1.0)
    {
        int currentOperations = nOps / pow(2, i);
        double relativeOperations = currentOperations * numIncrements[int(round(i - 1))] / pow(2, 15 + i);

        high_resolution_clock::time_point begin_time = high_resolution_clock::now();
        for (int j = 0; j < currentOperations; ++j)
            temp = Coil::computeMutualInductance(primary, secondary, 0.2, PrecisionFactor(i), method);
        double interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("inductance calc time for %.0f : %.2f ms/op\n", i, 1'000.0 * interval / currentOperations);

    }
}

void testCoilSelfInductance()
{
    Coil coil1 = Coil(0.03, 0.03, 0.12, 3600);

    for (double i = 1.0; i <= 8.0; i += 0.5)
    {
        printf("%.15g %.15g\n",
               coil1.computeAndSetSelfInductance(PrecisionFactor(i)),
               coil1.computeAndSetApproximateSelfInductance(PrecisionFactor(i)));
    }

}