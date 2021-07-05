#include <iostream>
#include "../include/OldCoil.h"
#include "../include/Polynomial.h"
#include "../include/Coil.h"
#include "../include/Precision.h"

Type ro1 = 1.63e-8;
extern thread_pool tp;

int main(){
    //setting the number of threads
    tp.resize(4);

    // polynomial testing
    int numPol = 20;

//    std::vector<Polynomial> legendrePolynomials = Polynomial::getLegendreSequenceUpToN(numPol);
//
//    for (int i = 0; i < numPol; i++)
//    {
//        legendrePolynomials[i].printForGrapher();
//    }

//    std::vector<double> zeros;
//    std::vector<double> weights;
//
//    for (int i = 1; i < numPol; ++i)
//    {
//        Polynomial::getLegendreParametersForN(i,  zeros,  weights);
//
//        printf("%2d:\n", i);
//        for (double j : zeros)
//        {
//            printf("%.12f, ", j);
//        }
//        printf("\n");
//
//        for (double j : weights){
//            printf("%.12f, ", j);
//        }
//        printf("\n");
//    }

    Coil testCoil1 = Coil(0.03, 0.03, 0.12, 3600);
    OldCoil oldCoil = OldCoil(1, 0.03, 0.03, 0.12, 0.001, 16, 16, 32, 3600, true, 50, 1.63e-8);

    printf("%.15f, %.15f, %.15f\n\n", testCoil1.getCurrentDensity(), testCoil1.getWireResistivity(), testCoil1.getSineFrequency());

    printf("%.15f, %.15f, %.15f\n\n", testCoil1.getMagneticMoment(), testCoil1.getAverageWireThickness(), testCoil1.getResistance());
    printf("%.15f, %.15f, %.15f\n\n", oldCoil.mM, oldCoil.d, oldCoil.Res);

    testCoil1.setSineFrequency(100000);
    testCoil1.setCurrentDensity(500000);
    printf("%.15f, %.15f, %.15f\n\n", testCoil1.getMagneticMoment(), testCoil1.getAverageWireThickness(), testCoil1.getResistance());

    PrecisionArguments arguments = testCoil1.getPrecisionSettings();
//    printf("%.15f", arguments.angularIncrementPositions[0]);

    for (double value : arguments.angularIncrementPositions)
        printf("%.15f ", value);
    printf("\n");
    for (double value : arguments.angularIncrementWeights)
        printf("%.15f ", value);
    printf("\n");

    for (double value : arguments.lengthIncrementPositions)
        printf("%.15f ", value);
    printf("\n");
    for (double value : arguments.lengthIncrementWeights)
        printf("%.15f ", value);
    printf("\n\n");


    testCoil1.setSineFrequency(0);
    testCoil1.setCurrent(1);
    printf("%.20f\n", testCoil1.computeBFieldZ(0.0, 0.0));
    printf("%.20f\n", testCoil1.computeBFieldH(0.0, 0.0));

    int nOp = 20000;
    std::vector<double> temp1;

    clock_t begin_time1 = clock();
    for (int i = 0; i < nOp; ++i){
        temp1 = testCoil1.computeBFieldVector(i*0.000001, 0.0, 0.0);
    }
    printf("combined B : %.0f dots/s\n", 1.0 / (float(clock() - begin_time1) / CLOCKS_PER_SEC / nOp));

    double temp2;
    clock_t begin_time2 = clock();
    for (int i = 0; i < nOp; ++i){
        temp2 = testCoil1.computeBFieldH(0.0, 0.0);
    }
    printf("field Bh : %.0f dots/s\n", 1.0 / (float(clock() - begin_time2) / CLOCKS_PER_SEC / nOp));

    double temp3;
    clock_t begin_time3 = clock();
    for (int i = 0; i < nOp; ++i){
        temp3 = testCoil1.computeBFieldH(0.0, 0.0);
    }
    printf("field Bz : %.0f dots/s\n", 1.0 / (float(clock() - begin_time3) / CLOCKS_PER_SEC / nOp));

    double temp4;
    clock_t begin_time4 = clock();
    for (int i = 0; i < nOp; ++i){
        temp4 = testCoil1.computeAPotentialAbs(i*0.000001, 0.0);
    }
    printf("potential A : %.0f dots/s\n", 1.0 / (float(clock() - begin_time4) / CLOCKS_PER_SEC / nOp));


//    PrecisionArguments precisionArguments = PrecisionArguments(2, 1, 1, 16, 12, 12);


//	FILE *input = fopen("values.txt", "r");
//	FILE *output = fopen("output.txt", "w");
//
//	Type Rt1, at1, bt1; int Nt1;
//	Type Rt2, at2, bt2; int Nt2;
//	Type distance;
//	OldCoil prim1, sec1;
//	Type Temp;
//
//	while (fscanf(input, "%lf %lf %lf %d %lf %lf %lf %d %lf", &Rt1, &at1, &bt1, &Nt1, &Rt2, &at2, &bt2, &Nt2, &distance) == 9){
//		printf("%f %f %f %d %f %f %f %d %f\n", Rt1, at1, bt1, Nt1, Rt2, at2, bt2, Nt2, distance);
//		for (Type i = 1; i <= 9.0; i += 1.0){
//			prim1 = OldCoil(1, Rt1, at1, bt1, 0.001, 12, 12, 80, Nt1, true, 100000, ro1);
//			sec1 = OldCoil(1, Rt2, at2, bt2, 0.001, 8, 8, 12, Nt2, true, 100000, ro1);
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