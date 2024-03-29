#include "Compare.h"

#include <cstdio>


void Compare::mutualInductanceZAxis()
{
    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);

    primary.setPositionAndOrientation();
    secondary.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.2));

    printf("%.20f\n\n", Coil::computeMutualInductance(primary, secondary));

    FILE *input = fopen("../data/values_MInductance_zAxis.txt", "r");
    FILE *output = fopen("output.txt", "w");

    if(!input || !output)
    {
        printf("Couldn't load file!\n");
        return;
    }

    double Rt1, at1, bt1; int Nt1;
    double Rt2, at2, bt2; int Nt2;
    double distance;
    double temp;

    while (fscanf(input, "%lf %lf %lf %d %lf %lf %lf %d %lf", &Rt1, &at1, &bt1, &Nt1, &Rt2, &at2, &bt2, &Nt2, &distance) == 9)
    {
        printf("%f %f %f %d %f %f %f %d %f\n", Rt1, at1, bt1, Nt1, Rt2, at2, bt2, Nt2, distance);

        Coil prim = Coil(Rt1, at1, bt1, Nt1);
        Coil sec = Coil(Rt2, at2, bt2, Nt2);

        prim.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.0));
        sec.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, distance));

        for (int i = 1; i <= 8; i++)
        {
            temp = Coil::computeMutualInductance(prim, sec, PrecisionFactor(i));
            printf("%.16g\n", temp);
            fprintf(output, "%.15g\t", temp);
        }

        printf("====================================================================================\n");
        fprintf(output, "\n");
    }

    fclose(input);
    fclose(output);
}

void Compare::selfInductance()
{
    Coil coil1 = Coil(0.03, 0.03, 0.12, 3600, PrecisionFactor(), 12);
    Coil coil2 = Coil(0.03, 0.0, 0.12, 120);
    Coil coil3 = Coil(0.03, 0.03, 0.0, 30);

    coil1.setThreadCount(8);

    for (int i = 1; i <= 12; ++i)
        printf("%.15g\n", coil1.computeAndSetSelfInductance(PrecisionFactor(i)));
    printf("\n");

    for (int i = 1; i <= 12; ++i)
        printf("%.15g\n", coil2.computeAndSetSelfInductance(PrecisionFactor(i)));
    printf("\n");

    for (int i = 1; i <= 12; ++i)
        printf("%.15g\n", coil3.computeAndSetSelfInductance(PrecisionFactor(i)));
    printf("\n");

    FILE *output = fopen("output.txt", "w");

    double r1 = 0.1;
    double alpha, beta;
    double temp;

    alpha = 1.5; beta = 0.25;
    Coil coil4 = Coil(r1, r1 * (alpha - 1), 2 * r1 * beta, 1);
    for (int i = 1; i <= 15; ++i)
    {
        temp = coil4.computeAndSetSelfInductance(PrecisionFactor(i)) / coil4.getInnerRadius();
        printf("%.11g\n", 1e6 * temp);
        fprintf(output, "%.11g\t", 1e6 * temp);
    }
    printf("\n");
    fprintf(output, "\n");

    alpha = 3.0; beta = 1.0;
    Coil coil5 = Coil(r1, r1 * (alpha - 1), 2 * r1 * beta, 1);
    for (int i = 1; i <= 15; ++i)
    {
        temp = coil5.computeAndSetSelfInductance(PrecisionFactor(i)) / coil4.getInnerRadius();
        printf("%.11g\n", 1e6 * temp);
        fprintf(output, "%.11g\t", 1e6 * temp);
    }
    printf("\n");
    fprintf(output, "\n");

    alpha = 4.0; beta = 3.0;
    Coil coil6 = Coil(r1, r1 * (alpha - 1), 2 * r1 * beta, 1);
    for (int i = 1; i <= 15; ++i)
    {
        temp = coil6.computeAndSetSelfInductance(PrecisionFactor(i)) / coil4.getInnerRadius();
        printf("%.11g\n", 1e6 * temp);
        fprintf(output, "%.11g\t", 1e6 * temp);
    }
    printf("\n");
    fprintf(output, "\n");

    alpha = 7.0; beta = 6.0;
    Coil coil7 = Coil(r1, r1 * (alpha - 1), 2 * r1 * beta, 1);
    for (int i = 1; i <= 15; ++i)
    {
        temp = coil7.computeAndSetSelfInductance(PrecisionFactor(i)) / coil4.getInnerRadius();
        printf("%.11g\n", 1e6 * temp);
        fprintf(output, "%.11g\t", 1e6 * temp);
    }
    printf("\n");
    fprintf(output, "\n");

    alpha = 9.0; beta = 4.0;
    Coil coil8 = Coil(r1, r1 * (alpha - 1), 2 * r1 * beta, 1);
    for (int i = 1; i <= 15; ++i)
    {
        temp = coil8.computeAndSetSelfInductance(PrecisionFactor(i)) / coil4.getInnerRadius();
        printf("%.11g\n", 1e6 * temp);
        fprintf(output, "%.11g\t", 1e6 * temp);
    }
    printf("\n");
    fprintf(output, "\n");

    fclose(output);
}
