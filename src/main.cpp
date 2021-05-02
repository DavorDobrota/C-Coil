#include <iostream>
#include "../include/Coil.h"
#include "../include/Polynomial.h"

Type ro1 = 1.63e-8;
extern thread_pool tp;

int main(){
    //setting the number of threads
    tp.resize(16);

    // polynomial testing
    int numPol = 26;

    std::vector<Polynomial> legendrePolynomials = Polynomial::getLegendreSequenceUpToN(numPol);

    for (int i = 0; i < numPol; i++)
    {
        legendrePolynomials[i].printForGrapher();
    }

    std::vector<double> zeros;
    std::vector<double> weights;

    for (int i = 1; i < numPol; ++i)
    {
        Polynomial::getLegendreParametersForN(i, zeros, weights);

        printf("%2d:\n", i);
        for (double j : zeros)
        {
            printf("%.12f, ", j);
        }
        printf("\n");

        for (double j : weights){
            printf("%.12f, ", j);
        }
        printf("\n");
    }

//	FILE *input = fopen("values.txt", "r");
//	FILE *output = fopen("output.txt", "w");
//
//	Type Rt1, at1, bt1; int Nt1;
//	Type Rt2, at2, bt2; int Nt2;
//	Type distance;
//	Coil prim1, sec1;
//	Type Temp;
//
//	while (fscanf(input, "%lf %lf %lf %d %lf %lf %lf %d %lf", &Rt1, &at1, &bt1, &Nt1, &Rt2, &at2, &bt2, &Nt2, &distance) == 9){
//		printf("%f %f %f %d %f %f %f %d %f\n", Rt1, at1, bt1, Nt1, Rt2, at2, bt2, Nt2, distance);
//		for (Type i = 1; i <= 9.0; i += 1.0){
//			prim1 = Coil(1, Rt1, at1, bt1, 0.001, 12, 12, 80, Nt1, true, 100000, ro1);
//			sec1 = Coil(1, Rt2, at2, bt2, 0.001, 8, 8, 12, Nt2, true, 100000, ro1);
//			Temp = sec1.MutualInductanceCalc(distance, prim1, true, i);
//			printf(" %.18f\n", Temp);
//			fprintf(output, "%.20f\t", Temp);
//		}
//		printf("====================================================================================\n");
//		fprintf(output, "\n");
//	}
//
//	fclose(input);
//	fclose(output);

    /*
    for (int i = 4; i <= 120; i += 4){
        for (int j = 4; j <= 120; j += 4){
            prim1 = Coil(1, Rt1, at1, bt1, 0.001, i, j, 128, Nt1, true, 100000, ro1);
            sec1 = Coil(1, Rt2, at2, bt2, 0.001, 20, 20, 128, Nt2, true, 100000, ro1);
            printf("%.20f\n", prim1.MutualInductanceCalc(distance, sec1, false, 4.0));
        }
    }
    */


    /*
    Coil prim = Coil(1, 1.0, 0.1, 0.1, 0.001, 32, 32, 48, 100, true, 100000, ro1);
    Coil sec = Coil(1, 0.1, 0.1, 0.1, 0.001, 32, 32, 20, 100, true, 100000, ro1);
    int nOp = 2000;
    Type temp;

    for (Type p = 1.0; p <= 9.0; p += 1.0){
        clock_t begin_time = clock();
        for (int i = 0; i < nOp; ++i){
            temp = prim.MutualInductanceCalc(0.2+0.00001*i, sec, true, p);
        }
        printf("%f\n", float(clock() - begin_time) / CLOCKS_PER_SEC / nOp * 1000);
    }
    */
/*
    float temp, theta = Pi/2, dist = 0.1;
    Coil test = Coil(1, 0.03, 0.03, 0.12, 0.001, 24, 24, 64, 3600, false, 0.0, 1.63e-8);
    test.newSelfInductanceCalc();
    printf("%.15e", test.L);
*/

    /*
    for (int i = 4; i <= 128; i += 4){

        prim = Coil(1, 0.2, 0.2, 0.2, 0.001, 8, 8, 48, 500, true, 100000, ro1);
        sec = Coil(1, 0.21, 0.2, 0.2, 0.001, i ,i, 24, 500, true, 100000, ro1);
        printf("%.15f\n", prim.MutualInductanceCalc(0.0, sec));
    }
    */

/*	Coil prim = Coil(1, 0.071335, 0.01397, 0.142748, 0.001, 20, 20, 60, 1142, true, 1, ro1);
	Coil sec = Coil(1, 0.096945, 0.041529, 0.02413, 0.001, 20, 20, 80, 516, true, 1, ro1);
	int nOp = 100;
	Type temp;*/
    /*
    clock_t begin_time = clock();
    for (int i = 0; i < nOp; ++i){
        temp = prim.MutualInductanceGeneralCalc(sec, 0.5+0.0001*i, 0.0, 0.0, 0.0, false, 1.0);
    }
    printf("%f\n", float(clock() - begin_time) / CLOCKS_PER_SEC / nOp * 1000);*/

//	printf("%.15f\n", prim.MutualInductanceGeneralCalc(sec, 0.07366, 0.30988, 0.0, 0.0, true, 4.5));
//	printf("%.15f\n", prim.MutualInductanceCalc(0.7, sec, true, 5.0));

//	Coil prim = Coil(1, 0.0375, 0.01, 0.01, 0.01, 16, 16, 60, 150, true, 1, ro1);
//	Coil sec = Coil(1, 0.018, 0.004, 0.004, 0.01, 16, 16, 60, 50, true, 1, ro1);

    /*
    for (int i = 0; i <= 300; ++i){
        printf("%.15f\n", prim.MutualInductanceGeneralCalc(sec, 0.02, i*0.001, 0.0, 0.0, true, 4.0));
    }
*/
//	Calculate_hardware_accelerated_e (1, &theta, &dist, 1000000.f, 0.03f, 0.12f, 0.03f, 0.001f, 0.001f, 0.2617993878f, nullptr, nullptr, &temp);
//	printf("%.30f\n", temp);

    return 0;
}
