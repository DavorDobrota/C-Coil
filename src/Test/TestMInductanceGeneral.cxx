#include "Test.h"
#include "Coil.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>


void testMutualInductanceGeneralForZAxis(ComputeMethod computeMethod, int nThreads)
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

void testMInductanceGeneralArgumentGeneration()
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
