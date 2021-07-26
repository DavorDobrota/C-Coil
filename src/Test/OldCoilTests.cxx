#include "Test.h"
#include "OldCoil.h"

#include <cstdio>

void testOldCoilMutualInductanceZAxis()
{
    FILE *input = fopen("values.txt", "r");
    FILE *output = fopen("output.txt", "w");

    Type Rt1, at1, bt1; int Nt1;
    Type Rt2, at2, bt2; int Nt2;
    Type distance;
    OldCoil prim1, sec1;
    Type Temp;

    while (fscanf(input, "%lf %lf %lf %d %lf %lf %lf %d %lf", &Rt1, &at1, &bt1, &Nt1, &Rt2, &at2, &bt2, &Nt2, &distance) == 9){
        printf("%f %f %f %d %f %f %f %d %f\n", Rt1, at1, bt1, Nt1, Rt2, at2, bt2, Nt2, distance);
        for (Type i = 1.0; i <= 9.0; i += 1.0){
            prim1 = OldCoil(1, Rt1, at1, bt1, 0.001, 12, 12, 80, Nt1, true, 100000, 1.63e-8);
            sec1 = OldCoil(1, Rt2, at2, bt2, 0.001, 8, 8, 12, Nt2, true, 100000, 1.63e-8);
            Temp = sec1.MutualInductanceCalc(distance, prim1, true, i);
            printf(" %.18f\n", Temp);
            fprintf(output, "%.20f\t", Temp);
        }
        printf("====================================================================================\n");
        fprintf(output, "\n");
    }

    fclose(input);
    fclose(output);
}

void testOldCoilMutualInductanceZAxisPerformance()
{
    OldCoil prim = OldCoil(1, 1.0, 0.1, 0.1, 0.001, 32, 32, 48, 100, true, 100000, 1.63e-8);
    OldCoil sec = OldCoil(1, 0.1, 0.1, 0.1, 0.001, 32, 32, 20, 100, true, 100000, 1.63e-8);
    int nOp = 1000;
    Type temp;

    for (Type p = 1.0; p <= 9.0; p += 1.0){
        clock_t begin_time = clock();
        for (int i = 0; i < nOp; ++i){
            temp = prim.MutualInductanceCalc(0.2+0.00001*i, sec, true, p);
        }
        printf("inductance calc time for %.0f : %.2f ms\n", p, 1000 * float(clock() - begin_time) / CLOCKS_PER_SEC / nOp);
    }
}

void testOldCoilMutualInductanceGeneralPerformance()
{

    OldCoil prim = OldCoil(1, 0.071335, 0.01397, 0.142748, 0.001, 20, 20, 60, 1142, true, 1, 1.63e-8);
    OldCoil sec = OldCoil(1, 0.096945, 0.041529, 0.02413, 0.001, 20, 20, 80, 516, true, 1, 1.63e-8);

    int nOp = 10;
    double temp;

    clock_t begin_time = clock();
    for (int i = 0; i < nOp; ++i){
        temp = prim.MutualInductanceGeneralCalc(sec, 0.5+0.0001*i, 0.0, 0.0, 0.0, true, 2.0);
    }
    printf("%f\n", float(clock() - begin_time) / CLOCKS_PER_SEC / nOp * 1000);
}

void testOldCoilSelfInductance()
{

    float temp, theta = Pi/2, dist = 0.1;
    OldCoil test = OldCoil(1, 1.0, 1.0, 2.0, 0.1, 24, 24, 64, 100, false, 0.0, 1.63e-8);
    test.SelfInductanceCalc();
    printf("%.15e", test.L);

}