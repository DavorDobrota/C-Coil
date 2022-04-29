#include "Test.h"
#include "Coil.h"

#include <cstdio>
#include <math.h>


void testMutualInductanceGeneralMisalignedCoils()
{
    Coil primary = Coil(0.06, 0.0, 0.12, 120);
    Coil secondary = Coil(0.05, 0.0, 0.0, 1);

    auto precision = PrecisionFactor(8.0);

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(vec3::CoordVector3(), std::acos(i * 0.1));
        printf("%.15g\n", 1e6 * Coil::computeMutualInductance(primary, secondary, precision));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.03), std::acos(i * 0.1));
        printf("%.15g\n", 1e6 * Coil::computeMutualInductance(primary, secondary, precision));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.06), std::acos(i * 0.1));
        printf("%.15g\n", 1e6 * Coil::computeMutualInductance(primary, secondary, precision));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.12), std::acos(i * 0.1));
        printf("%.15g\n", 1e6 * Coil::computeMutualInductance(primary, secondary, precision));
    }
    printf("\n");


    primary = Coil(0.06, 0.0, 0.12, 120);
    secondary = Coil(0.05, 0.0, 0.04, 60);

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(vec3::CoordVector3(), std::acos(i * 0.1));
        printf("%.15g\n", 1e6 * Coil::computeMutualInductance(primary, secondary, precision));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
        //   printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.03), std::acos(i * 0.1));
        printf("%.15g\n", 1e6 * Coil::computeMutualInductance(primary, secondary, precision));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.06), std::acos(i * 0.1));
        printf("%.15g\n", 1e6 * Coil::computeMutualInductance(secondary, primary, precision));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.12), std::acos(i * 0.1));
        printf("%.15g\n", 1e6 * Coil::computeMutualInductance(secondary, primary, precision));
    }
    printf("\n");


    primary = Coil(0.04, 0.02, 0.0, 200);
    secondary = Coil(0.015, 0.01, 0.0, 100);

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.05), std::acos(i * 0.1));
        printf("%.15g\n", 1e6 * Coil::computeMutualInductance(primary, secondary, precision));
    }
    printf("\n");


    primary = Coil(0.04, 0.02, 0.01, 100);
    secondary = Coil(0.02, 0.0, 0.0, 1);

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(vec3::CoordVector3(), std::acos(i * 0.1));
        printf("%.15g\n", 1e6 * Coil::computeMutualInductance(primary, secondary, precision));
    }
    printf("\n");
}

void testMutualInductanceGeneralParallelAxes()
{
    double tempInductance;
    auto precision = PrecisionFactor(6.0);

    Coil coil1 = Coil(0.071247, 0.01397, 0.142748, 1142);
    Coil coil2 = Coil(0.0969645, 0.041529, 0.02413, 516);

    double rArr1[] = {0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.0117475,
                      0.2237105, 0.224, 0.225, 0.23, 0.24, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 5.0, 10.0};

    for (double i : rArr1)
    {
        coil2.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, i, 0.0, 0.0));
        tempInductance = Coil::computeMutualInductance(coil1, coil2, precision, CPU_MT);
        printf("%9g : %.14g mH\n", i, 1e3 * tempInductance);
    }
    printf("\n");

    double rArr2[] = {0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.02, 0.02,
                      0.02, 0.02, 0.02, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
    double zArr2[] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.059309, 0.07, 0.083439, 0.09, 0.1, 0.6, 1.0, 0.083439, 0.09, 0.1,
                      0.6, 1.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.059309, 0.07, 0.083439, 0.09, 0.1, 0.6, 1.0};
    for (int i = 0; i < 29; ++i)
    {
        coil2.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, rArr2[i], 0.0, zArr2[i]));
        tempInductance = Coil::computeMutualInductance(coil1, coil2, precision, CPU_MT);
        printf("%8g %5g : %.14g mH\n", zArr2[i], rArr2[i], 1e3 * tempInductance);
    }
    printf("\n");

    FILE *input = fopen("values_MInductance_general.txt", "r");
    FILE *output = fopen("output.txt", "w");

    double Rt1, at1, bt1; int Nt1;
    double Rt2, at2, bt2; int Nt2;
    double distanceZ, distanceR;
    double temp;

    while (fscanf(input, "%lf %lf %lf %d %lf %lf %lf %d %lf %lf", &Rt1, &at1, &bt1, &Nt1, &Rt2, &at2, &bt2, &Nt2, &distanceZ, &distanceR) == 10)
    {
        printf("%g %g %g %d : %g %g %g %d | %g %g\n", Rt1, at1, bt1, Nt1, Rt2, at2, bt2, Nt2, distanceZ, distanceR);

        Coil prim = Coil(Rt1, at1, bt1, Nt1);
        Coil sec = Coil(Rt2, at2, bt2, Nt2);
        sec.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, distanceR, 0.0, distanceZ));

        for (int i = 1; i <= 9; i++)
        {
            temp = Coil::computeMutualInductance(prim, sec, PrecisionFactor(i), CPU_MT);
            printf("%.12f\n", 1e9 *  temp);
            fprintf(output, "%.12f\t", 1e9 * temp);
        }
        printf("===========================================================================\n");
        fprintf(output, "\n");
    }
    fclose(input);
    fclose(output);
}

void testMutualInductanceGeneralConway()
{
    FILE *output = fopen("output.txt", "w");

    Coil coil1 = Coil(0.1, 0.02, 0.24, 1200);
    Coil coil2 = Coil(0.04, 0.01, 0.06, 1200);
    double temp;

    for (int i = 10; i > 0; --i)
    {
        coil2.setPositionAndOrientation(vec3::CoordVector3(), std::acos(0.1 * i));
        for (int j = 1; j <= 9; ++j)
        {
            temp = Coil::computeMutualInductance(coil1, coil2, PrecisionFactor(j), CPU_MT);
            printf("%.16g\n", 1e3 * temp);
            fprintf(output, "%.16g\t", 1e3 * temp);
        }
        printf("\n");
        fprintf(output, "\n");
    }
    fprintf(output, "\n");
    printf("\n");

    for (int i = 10; i > 0; --i)
    {
        coil2.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.1, 0.0, 0.18),
                                        std::acos(0.1 * i), 0.0);
        for (int j = 1; j <= 9; ++j)
        {
            temp = Coil::computeMutualInductance(coil1, coil2, PrecisionFactor(j), CPU_MT);
            printf("%.16g\n", 1e3 * temp);
            fprintf(output, "%.16g\t", 1e3 * temp);
        }
        printf("\n");
        fprintf(output, "\n");
    }
    printf("\n");

    fclose(output);
}

