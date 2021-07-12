#include <iostream>
#include "../include/OldCoil.h"
#include "../include/Polynomial.h"
#include "../include/Coil.h"
#include "test.h"

Type ro1 = 1.63e-8;
extern thread_pool tp;

int main(){

    testCoilMutualInductanceZAxis();

//    for (int i = 1; i <= 5000; ++i)
//        if (i % (i / 50 + 1) == 0 || i % 50 == 0 || i < 50)
//        {
//            if (i % 50 == 0)
//                printf("%d, ", i / (i / 50));
//            else
//                printf("%d, ", i / (i / 50 + 1));
//
//        }


    /*
    for (int i = 4; i <= 120; i += 4){
        for (int j = 4; j <= 120; j += 4){
            prim1 = OldCoil(1, Rt1, at1, bt1, 0.001, i, j, 128, Nt1, true, 100000, ro1);
            sec1 = OldCoil(1, Rt2, at2, bt2, 0.001, 20, 20, 128, Nt2, true, 100000, ro1);
            printf("%.20f\n", prim1.MutualInductanceCalc(distance, sec1, false, 4.0));
        }
    }
    */

//    OldCoil prim = OldCoil(1, 1.0, 0.1, 0.1, 0.001, 32, 32, 48, 100, true, 100000, ro1);
//    OldCoil sec = OldCoil(1, 0.1, 0.1, 0.1, 0.001, 32, 32, 20, 100, true, 100000, ro1);
//    int nOp = 2000;
//    Type temp;
//
//    for (Type p = 1.0; p <= 9.0; p += 1.0){
//        clock_t begin_time = clock();
//        for (int i = 0; i < nOp; ++i){
//            temp = prim.MutualInductanceCalc(0.2+0.00001*i, sec, true, p);
//        }
//        printf("%f\n", float(clock() - begin_time) / CLOCKS_PER_SEC / nOp * 1000);
//    }

/*
    float temp, theta = Pi/2, dist = 0.1;
    OldCoil test = OldCoil(1, 0.03, 0.03, 0.12, 0.001, 24, 24, 64, 3600, false, 0.0, 1.63e-8);
    test.newSelfInductanceCalc();
    printf("%.15e", test.L);
*/

    /*
    for (int i = 4; i <= 128; i += 4){

        prim = OldCoil(1, 0.2, 0.2, 0.2, 0.001, 8, 8, 48, 500, true, 100000, ro1);
        sec = OldCoil(1, 0.21, 0.2, 0.2, 0.001, i ,i, 24, 500, true, 100000, ro1);
        printf("%.15f\n", prim.MutualInductanceCalc(0.0, sec));
    }
    */

/*	OldCoil prim = OldCoil(1, 0.071335, 0.01397, 0.142748, 0.001, 20, 20, 60, 1142, true, 1, ro1);
	OldCoil sec = OldCoil(1, 0.096945, 0.041529, 0.02413, 0.001, 20, 20, 80, 516, true, 1, ro1);
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

//	OldCoil prim = OldCoil(1, 0.0375, 0.01, 0.01, 0.01, 16, 16, 60, 150, true, 1, ro1);
//	OldCoil sec = OldCoil(1, 0.018, 0.004, 0.004, 0.01, 16, 16, 60, 50, true, 1, ro1);

    /*
    for (int i = 0; i <= 300; ++i){
        printf("%.15f\n", prim.MutualInductanceGeneralCalc(sec, 0.02, i*0.001, 0.0, 0.0, true, 4.0));
    }
*/
//	Calculate_hardware_accelerated_e (1, &theta, &dist, 1000000.f, 0.03f, 0.12f, 0.03f, 0.001f, 0.001f, 0.2617993878f, nullptr, nullptr, &temp);
//	printf("%.30f\n", temp);

    return 0;
}
