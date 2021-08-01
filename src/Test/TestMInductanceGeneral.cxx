#include "Test.h"
#include "Coil.h"

#include <cstdio>
#include <cmath>

void testCoilMutualInductanceGeneralForZAxis()
{
    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);

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
            temp = Coil::computeMutualInductance(prim, sec, distance, 1e-15, PrecisionFactor(i));
            printf("%.18f\n", temp);
            fprintf(output, "%.20f\t", temp);
        }

        printf("====================================================================================\n");
        fprintf(output, "\n");
    }

    fclose(input);
    fclose(output);
}


void testCoilMutualInductanceGeneralPerformance()
{

}


void testCoilMutualInductanceGeneralArgumentGeneration()
{
    Coil coil1 = Coil(0.05, 0.1, 0.1, 100);
    Coil coil2 = Coil(0.05, 0.1, 1e-15, 10);
    Coil coil3 = Coil(0.05, 1e-15, 0.1, 10);
    Coil coil4 = Coil(0.05, 1e-15, 1e-15, 1);

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

void testCoilMutualInductanceGeneralDifferentGeometries()
{
    Coil prim1 = Coil(0.03, 1e-15, 1e-15, 1);
    Coil sec1 = Coil(0.03, 1e-15, 1e-15, 1);
    printf("%.12g\n", Coil::computeMutualInductance(prim1, sec1, 0.2, 0.0, 1e-15));

    Coil prim2 = Coil(0.03, 0.03, 1e-15, 30);
    Coil sec2 = Coil(0.03, 1e-15, 1e-15, 1);
    printf("%.12g\n", Coil::computeMutualInductance(prim2, sec2, 0.2, 0.0, 1e-15));

    Coil prim3 = Coil(0.03, 1e-15, 1e-15, 1);
    Coil sec3 = Coil(0.03, 0.03, 1e-15, 30);
    printf("%.12g\n", Coil::computeMutualInductance(prim3, sec3, 0.2, 0.0, 1e-15));

    Coil prim4 = Coil(0.03, 0.03, 1e-15, 30);
    Coil sec4 = Coil(0.03, 0.03, 1e-15, 30);
    printf("%.12g\n", Coil::computeMutualInductance(prim4, sec4, 0.2, 0.0, 1e-15));

    Coil prim5 = Coil(0.03, 1e-15, 0.12, 120);
    Coil sec5 = Coil(0.03,  1e-15, 1e-15, 1);
    printf("%.12g\n", Coil::computeMutualInductance(prim5, sec5, 0.2, 0.0, 1e-15));

    Coil prim6 = Coil(0.03, 1e-15, 1e-15, 1);
    Coil sec6 = Coil(0.03,  1e-15, 0.12, 120);
    printf("%.12g\n", Coil::computeMutualInductance(prim6, sec6, 0.2, 0.0, 1e-15));

    Coil prim7 = Coil(0.03, 1e-15, 0.12, 120);
    Coil sec7 = Coil(0.03,  1e-15, 0.12, 120);
    printf("%.12g\n", Coil::computeMutualInductance(prim7, sec7, 0.2, 0.0, 1e-15));

    Coil prim8 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec8 = Coil(0.03,  1e-15, 1e-15, 1);
    printf("%.12g\n", Coil::computeMutualInductance(prim8, sec8, 0.2, 0.0, 1e-15));

    Coil prim9 = Coil(0.03,  1e-15, 1e-15, 1);
    Coil sec9 = Coil(0.03, 0.03, 0.12, 3600);
    printf("%.12g\n", Coil::computeMutualInductance(prim9, sec9, 0.2, 0.0, 1e-15));

    Coil prim10 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec10 = Coil(0.03, 0.03, 0.12, 3600);
    printf("%.12g\n", Coil::computeMutualInductance(prim10, sec10, 0.2, 0.0, 1e-15));

    Coil prim11 = Coil(0.03, 0.03, 1e-15, 30);
    Coil sec11 = Coil(0.03, 1e-15, 0.12, 120);
    printf("%.12g\n", Coil::computeMutualInductance(prim11, sec11, 0.2, 0.0, 1e-15));

    Coil prim12 = Coil(0.03, 1e-15, 0.12, 120);
    Coil sec12 = Coil(0.03, 0.03, 1e-15, 30);
    printf("%.12g\n", Coil::computeMutualInductance(prim12, sec12, 0.2, 0.0, 1e-15));

    Coil prim13 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec13 = Coil(0.03, 1e-15, 0.12, 120);
    printf("%.12g\n", Coil::computeMutualInductance(prim13, sec13, 0.2, 0.0, 1e-15));

    Coil prim14 = Coil(0.03, 1e-15, 0.12, 120);
    Coil sec14 = Coil(0.03, 0.03, 0.12, 3600);
    printf("%.12g\n", Coil::computeMutualInductance(prim14, sec14, 0.2, 0.0, 1e-15));

    Coil prim15 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec15 = Coil(0.03, 0.03, 1e-15, 30);
    printf("%.12g\n", Coil::computeMutualInductance(prim15, sec15, 0.2, 0.0, 1e-15));

    Coil prim16 = Coil(0.03, 0.03, 1e-15, 30);
    Coil sec16 = Coil(0.03, 0.03, 0.12, 3600);
    printf("%.12g\n", Coil::computeMutualInductance(prim16, sec16, 0.2, 0.0, 1e-15));

    printf("\n");
}

void testCoilMutualInductanceGeneralGraphs()
{
    FILE *output = fopen("output.txt", "w");

    Coil prim = Coil(0.01022, 0.011, 0.0022, 20, PrecisionFactor(6.0), 16);
    Coil sec = Coil(0.01022, 0.011, 0.0022, 20, PrecisionFactor(6.0), 16);

    double temp;

//    for (int i = 0; i <= 200; ++i)
//        printf("%.10g\n",
//                Coil::computeMutualInductance(prim, sec, 0.04, i * 0.001, M_PI/2, PrecisionFactor(6.0), CPU_MT));
//    printf("\n");
//    printf("0. Done\n");

    for (int i = 0; i <= 100; ++i)
        fprintf(output, "%f\t%.10g\n", i * 0.001,
                Coil::computeMutualInductance(prim, sec, i * 0.001, PrecisionFactor(6.0), CPU_MT));
    fprintf(output,"\n");
    printf("1/6 tests Done\n");

    for (int i = 0; i <= 100; ++i)
        fprintf(output, "%f\t%.10g\n", i * 0.002,
               Coil::computeMutualInductance(prim, sec, 0.02, i * 0.002, PrecisionFactor(6.0), CPU_MT));
    fprintf(output,"\n");
    printf("2/6 tests done\n");

    for (int i = 0; i <= 100; ++i)
        fprintf(output, "%f\t%.10g\n", M_PI/50 * i,
               Coil::computeMutualInductance(prim, sec, 0.04, 0.0, M_PI/50 * i, PrecisionFactor(6.0), CPU_MT));
    fprintf(output,"\n");
    printf("3/6 tests done\n");

    for (int i = 0; i <= 50; ++i)
    {
        for (int j = 0; j <= 50; ++j)
            fprintf(output, "%.10g\t",
                    Coil::computeMutualInductance(prim, sec, i * 0.002, j * 0.002, PrecisionFactor(6.0), CPU_MT));
        fprintf(output,"\n");
        printf("test progress: %d / 50\n", i);
    }
    fprintf(output,"\n");
    printf("4/6 tests done\n");

    for (int i = 0; i <= 50; ++i)
    {
        for (int j = 0; j <= 50; ++j)
            fprintf(output, "%.10g\t",
                    Coil::computeMutualInductance(prim, sec, 0.04, i * 0.004, -M_PI/50 * j,
                                                  PrecisionFactor(6.0), CPU_MT));
        fprintf(output,"\n");
        printf("test progress: %d / 50\n", i);
    }
    fprintf(output,"\n");
    printf("5/6 tests done\n");

    for (int i = 0; i <= 50; ++i)
    {
        for (int j = 0; j <= 50; ++j)
            fprintf(output, "%.10g\t",
                    Coil::computeMutualInductance(prim, sec, 0.025 + 0.002 * i, 0.0, -M_PI/50 * j,
                                                  PrecisionFactor(6.0), CPU_MT));
        fprintf(output,"\n");
        printf("test progress: %d / 50\n", i);
    }
    fprintf(output,"\n");
    printf("6/6 tests done\n\n");

    fclose(output);
}
