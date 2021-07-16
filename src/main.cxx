#include <iostream>
#include "../include/OldCoil.h"
#include "../include/Polynomial.h"
#include "../include/Coil.h"
#include "test.h"

extern thread_pool tp;

int main(){

//    testCoilMutualInductanceZAxis();

    testPerformanceForComputeAll();

//    testCoilMutualInductanceZAxisPerformance();
//    testOldCoilMutualInductanceZAxisPerformance();

//    testCoilMutualInductanceGeneralForZAxis();


    /*
    for (int i = 4; i <= 120; i += 4){
        for (int j = 4; j <= 120; j += 4){
            prim1 = OldCoil(1, Rt1, at1, bt1, 0.001, i, j, 128, Nt1, true, 100000, ro1);
            sec1 = OldCoil(1, Rt2, at2, bt2, 0.001, 20, 20, 128, Nt2, true, 100000, ro1);
            printf("%.20f\n", prim1.MutualInductanceCalc(distance, sec1, false, 4.0));
        }
    }
    */

/*
    float temp, theta = Pi/2, dist = 0.1;
    OldCoil test = OldCoil(1, 0.03, 0.03, 0.12, 0.001, 24, 24, 64, 3600, false, 0.0, 1.63e-8);
    test.newSelfInductanceCalc();
    printf("%.15e", test.L);
*/


//	OldCoil prim = OldCoil(1, 0.071335, 0.01397, 0.142748, 0.001, 20, 20, 60, 1142, true, 1, 1.63e-8);
//	OldCoil sec = OldCoil(1, 0.096945, 0.041529, 0.02413, 0.001, 20, 20, 80, 516, true, 1, 1.63e-8);
//
//	printf("%.15f\n", prim.MutualInductanceGeneralCalc(sec, 0.07366, 0.30988, 0.0, 0.0, true, 4.5));
//
//	Coil primary = Coil(0.071335, 0.01397, 0.142748, 1142);
//	Coil secondary = Coil(0.096945, 0.041529, 0.02413, 516);
//
//	printf("%.15f\n", Coil::computeMutualInductance(primary, secondary,
//                                                    0.07366, 0.30988, MInductanceArguments(4.0)));
//
//	Coil primaryGeneral = Coil(0.06, 1e-18, 0.12, 120);
//	Coil secondaryGeneral = Coil(0.05, 1e-18, 1e-18, 1);
//
//	for (int i = 0; i <= 10; i++){
//        printf("%.15f\n", Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
//                                                        0.0, 0.0,
//                                                        acos(i * 0.1), MInductanceArguments(6.0)));
//	}
//
//    printf("%.15f\n", Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
//                                                    0.0, 0.0,  MInductanceArguments(7.0)));
//
//    printf("%.15f\n", Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
//                                                    0.03, 0.0, MInductanceArguments(7.0)));
//
//    printf("%.15f\n", Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
//                                                    0.06, 0.0, MInductanceArguments(7.0)));
//
//    printf("%.15f\n", Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
//                                                    0.12, 0.0, MInductanceArguments(7.0)));



//	int nOp = 100;
//	Type temp;

//    clock_t begin_time = clock();
//    for (int i = 0; i < nOp; ++i){
//        temp = prim.MutualInductanceGeneralCalc(sec, 0.5+0.0001*i, 0.0, 0.0, 0.0, false, 1.0);
//    }
//    printf("%f\n", float(clock() - begin_time) / CLOCKS_PER_SEC / nOp * 1000);


//	printf("%.15f\n", prim.MutualInductanceCalc(0.7, sec, true, 5.0));

//	OldCoil prim = OldCoil(1, 0.0375, 0.01, 0.01, 0.01, 16, 16, 60, 150, true, 1, ro1);
//	OldCoil sec = OldCoil(1, 0.018, 0.004, 0.004, 0.01, 16, 16, 60, 50, true, 1, ro1);

    /*
    for (int i = 0; i <= 300; ++i){
        printf("%.15f\n", prim.MutualInductanceGeneralCalc(sec, 0.02, i*0.001, 0.0, 0.0, true, 4.0));
    }
*/

    return 0;
}
