
#include <cstdio>
#include <vector>

#include "test.h"
#include "Polynomial.h"
#include "Coil.h"
#include "OldCoil.h"
#include "Precision.h"
#include "ComputeMethod.h"

void testLegendrePolynomials()
{
    int numPol = 20;

    std::vector<Polynomial> legendrePolynomials = Polynomial::getLegendreSequenceUpToN(numPol);

    for (int i = 0; i < numPol; i++)
    {
        legendrePolynomials[i].printForGrapher();
    }

    std::vector<double> zeros;
    std::vector<double> weights;

    for (int i = 1; i < numPol; ++i)
    {
        Polynomial::getLegendreParametersForN(i,  zeros,  weights);

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
}


void testNewCoilParameters()
{
    Coil testCoil1 = Coil(0.03, 0.03, 0.12, 3600);
    OldCoil oldCoil = OldCoil(1, 0.03, 0.03, 0.12, 0.001, 16, 16, 32, 3600, true, 50, 1.63e-8);

    printf("%.15f, %.15f, %.15f\n\n", testCoil1.getCurrentDensity(), testCoil1.getWireResistivity(), testCoil1.getSineFrequency());

    printf("%.15f, %.15f, %.15f\n\n", testCoil1.getMagneticMoment(), testCoil1.getAverageWireThickness(), testCoil1.getResistance());
    printf("%.15f, %.15f, %.15f\n\n", oldCoil.mM, oldCoil.d, oldCoil.Res);

    testCoil1.setSineFrequency(100000);
    testCoil1.setCurrentDensity(500000);
    printf("%.15f, %.15f, %.15f\n\n", testCoil1.getMagneticMoment(), testCoil1.getAverageWireThickness(), testCoil1.getResistance());

    testCoil1.setSineFrequency(0);
    testCoil1.setCurrent(1);
    std::vector<double> fieldVector = testCoil1.computeBFieldVector(0.0, 0.0, 0.0);
    printf("%.25f %.25f\n", fieldVector[2], testCoil1.computeBFieldZ(0.0, 0.0));
    printf("%.25f %.25f\n", fieldVector[0], testCoil1.computeBFieldH(0.0, 0.0));
}

void testPerformanceCPU_ST()
{
    int nOp = 100000;
    std::vector<double> temp1;
    Coil testCoil = Coil(0.03, 0.03, 0.12, 3600);

    PrecisionArguments precision = testCoil.getPrecisionSettings();

    int numOperations = nOp *
            precision.thicknessBlockCount * precision.thicknessIncrementCount *
            precision.lengthBlockCount * precision.lengthIncrementCount *
            precision.angularBlockCount * precision.angularIncrementCount;

    clock_t begin_time1 = clock();
    for (int i = 0; i < nOp; ++i){
        temp1 = testCoil.computeBFieldVector(i*0.000001, 0.0, 0.0);
    }
    printf("combined B  : %.1f MInc/s\n", 1e-6/ (float(clock() - begin_time1) / CLOCKS_PER_SEC / numOperations));

    double temp2;
    clock_t begin_time2 = clock();
    for (int i = 0; i < nOp; ++i){
        temp2 = testCoil.computeAPotentialAbs(i*0.000001, 0.0);
    }
    printf("potential A : %.1f MInc/s\n", 1e-6 / (float(clock() - begin_time2) / CLOCKS_PER_SEC / numOperations));

}

void testPerformanceForComputeAll()
{
    Coil testCoil = Coil(0.03, 0.03, 0.12, 3600);

    int nOps = 80000;
    double radius = 0.1;

    PrecisionArguments precision = testCoil.getPrecisionSettings();

    int numOperations = nOps *
                        precision.thicknessBlockCount * precision.thicknessIncrementCount *
                        precision.lengthBlockCount * precision.lengthIncrementCount *
                        precision.angularBlockCount * precision.angularIncrementCount;

    int numOperationsGpu = nOps * 48 * 16 * 16;

    std::vector<double> cylindricalZArr;
    std::vector<double> cylindricalRArr;
    std::vector<double> cylindricalPhiArr;

    for (int i = 0; i < nOps ; i++)
    {
        cylindricalZArr.push_back(radius * cos(i * 2*Pi / nOps));
        cylindricalRArr.push_back(radius * sin(i * 2*Pi / nOps));
        cylindricalPhiArr.push_back(0.0);
    }

    std::vector<double> singleResultsX;
    std::vector<double> singleResultsY;
    std::vector<double> singleResultsZ;

    std::vector<double> acceleratedResultsX;
    std::vector<double> acceleratedResultsY;
    std::vector<double> acceleratedResultsZ;

    std::vector<double> singlePotential;
    std::vector<double> acceleratedPotential;

    clock_t begin_time11 = clock();
    testCoil.computeAllBFieldComponents(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                                        singleResultsX, singleResultsY, singleResultsZ,
                                        CPU_ST);
    printf("combined  B CPU : %.1f MInc/s\n", 1e-6 / (float(clock() - begin_time11) / CLOCKS_PER_SEC / numOperations));

    clock_t begin_time12 = clock();
    testCoil.computeAllAPotentialAbs(cylindricalZArr, cylindricalRArr,
                                     singlePotential, CPU_ST);
    printf("Potential A CPU : %.1f MInc/s\n", 1e-6 / (float(clock() - begin_time12) / CLOCKS_PER_SEC / numOperations));

    testCoil.computeAllBFieldComponents(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                                        acceleratedResultsX, acceleratedResultsY, acceleratedResultsZ,
                                        GPU);

    clock_t begin_time13 = clock();
    testCoil.computeAllBFieldComponents(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                                        acceleratedResultsX, acceleratedResultsY, acceleratedResultsZ,
                                        GPU);
    printf("combined  B GPU : %.1f MInc/s\n", 1e-6 / (float(clock() - begin_time13) / CLOCKS_PER_SEC / numOperationsGpu));

    clock_t begin_time14 = clock();
    testCoil.computeAllAPotentialAbs(cylindricalZArr, cylindricalRArr,
                                     acceleratedPotential, GPU);
    printf("Potential A GPU : %.1f MInc/s\n", 1e-6 / (float(clock() - begin_time14) / CLOCKS_PER_SEC / numOperationsGpu));
}

void testMethodPrecisionCompareCPUvsGPU()
{
    Coil testCoil = Coil(0.03, 0.03, 0.12, 3600);

    int pointCount = 2000;
    double radius = 0.075;

    std::vector<double> cylindricalZArr;
    std::vector<double> cylindricalRArr;
    std::vector<double> cylindricalPhiArr;

    for (int i = 0; i < pointCount ; i++)
    {
        cylindricalZArr.push_back(radius * cos(i * 2 * Pi / pointCount));
        cylindricalRArr.push_back(radius * sin(i * 2 * Pi / pointCount));
        cylindricalPhiArr.push_back(0.0);
    }

    std::vector<double> singleResultsX;
    std::vector<double> singleResultsY;
    std::vector<double> singleResultsZ;

    std::vector<double> acceleratedResultsX;
    std::vector<double> acceleratedResultsY;
    std::vector<double> acceleratedResultsZ;

    std::vector<double> singlePotential;
    std::vector<double> acceleratedPotential;

    testCoil.computeAllBFieldComponents(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                                        singleResultsX, singleResultsY, singleResultsZ,
                                        CPU_ST);
    testCoil.computeAllAPotentialAbs(cylindricalZArr, cylindricalRArr,
                                     singlePotential, CPU_ST);

    testCoil.computeAllBFieldComponents(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                                        acceleratedResultsX, acceleratedResultsY, acceleratedResultsZ,
                                        GPU);
    testCoil.computeAllAPotentialAbs(cylindricalZArr, cylindricalRArr,
                                     acceleratedPotential, GPU);

    FILE *output = fopen("output.txt", "w");

    for (int i = 0; i < pointCount; ++i)
    {
        fprintf(output, "%.20f\t%.20f\t%.20f\t%.20f\t%.20f\t%.20f\n",
               singleResultsX[i], acceleratedResultsX[i],
               singleResultsZ[i], acceleratedResultsZ[i],
               singlePotential[i], acceleratedPotential[i]);
    }
}

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
            fprintf(output, "%.7g\t", temp);
        }

        printf("====================================================================================\n");
        fprintf(output, "\n");
    }

    fclose(input);
    fclose(output);
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

void testCoilMutualInductanceZAxisPerformance()
{
    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);

    int nOps = 2560;
    int numIncrements[] = {78732, 147000, 263296, 547560, 1057500, 2247264, 4528384, 9168896};
    double temp;

    for (double i = 1.0; i <= 8.0; i += 1.0)
    {
        int currentOperations = nOps / pow(2, i);
        double relativeOperations = currentOperations * numIncrements[int(round(i - 1))] / pow(2, 15 + i);

        clock_t begin_time = clock();
        for (int j = 0; j < currentOperations; ++j)
            temp = Coil::computeMutualInductance(primary, secondary, 0.2, PrecisionFactor(i));
        printf("inductance calc time for %.0f : %.2f ms\n", i, 1000 * (float(clock() - begin_time) / CLOCKS_PER_SEC / relativeOperations));

    }
}

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
    OldCoil test = OldCoil(1, 0.03, 0.03, 0.12, 0.001, 24, 24, 64, 3600, false, 0.0, 1.63e-8);
    test.newSelfInductanceCalc();
    printf("%.15e", test.L);

}

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

void testCoilMutualInductanceForSpecialCase()
{
    OldCoil prim = OldCoil(1, 0.071335, 0.01397, 0.142748, 0.001, 20, 20, 60, 1142, true, 1, 1.63e-8);
	OldCoil sec = OldCoil(1, 0.096945, 0.041529, 0.02413, 0.001, 20, 20, 80, 516, true, 1, 1.63e-8);

	printf("%.15f\n", prim.MutualInductanceGeneralCalc(sec, 0.07366, 0.30988, 0.0, 0.0, true, 4.5));

	Coil primary = Coil(0.071335, 0.01397, 0.142748, 1142);
	Coil secondary = Coil(0.096945, 0.041529, 0.02413, 516);

	printf("%.15f\n", Coil::computeMutualInductance(primary, secondary,
                                                    0.07366, 0.30988,
                                                    PrecisionFactor(4.0)));
}

void testCoilMutualInductanceGeneralThinCoilAndFilament()
{
    Coil primaryGeneral = Coil(0.06, 1e-15, 0.12, 120);
    Coil secondaryGeneral = Coil(0.05, 1e-15, 1e-15, 1);

    MInductanceArguments inductanceArguments = MInductanceArguments(PrecisionArguments(4, 1, 1, 50, 1, 32),
                                                                    PrecisionArguments(4, 1, 1, 50, 1, 1));

    for (int i = 10; i >= 0; --i){
    //    printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                        0.00, 0.0,
                                                        acos(i * 0.1), 0.0, inductanceArguments));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
    //    printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                        0.03, 0.0,
                                                        acos(i * 0.1), 0.0, inductanceArguments));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
    //    printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                        0.06, 0.0,
                                                        acos(i * 0.1), 0.0, inductanceArguments));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
    //    printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                        0.12, 0.0,
                                                        acos(i * 0.1), 0.0, inductanceArguments));
    }
    printf("\n");
}

void testCoilMutualInductanceGeneralThinCoilAndThinCoil()
{
    Coil primaryGeneral = Coil(0.06, 1e-18, 0.12, 120);
    Coil secondaryGeneral = Coil(0.05, 1e-18, 0.04, 60);

    MInductanceArguments inductanceArguments = MInductanceArguments(PrecisionArguments(2, 1, 1, 50, 1, 32),
                                                                    PrecisionArguments(2, 1, 1, 50, 1, 20));

    for (int i = 10; i >= 0; --i){
    //    printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                              0.0, 0.0,
                                                              acos(i * 0.1), inductanceArguments));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
     //   printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                              0.03, 0.0,
                                                              acos(i * 0.1),  inductanceArguments));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
    //    printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                              0.06, 0.0,
                                                              acos(i * 0.1), inductanceArguments));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
    //    printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                              0.12, 0.0,
                                                              acos(i * 0.1), inductanceArguments));
    }
    printf("\n");
}

void testCoilMutualInductanceGeneralPancakeAndPancake()
{
    Coil primaryGeneral = Coil(0.04, 0.02, 1e-15, 200);
    Coil secondaryGeneral = Coil(0.015, 0.01, 1e-15, 100);

    MInductanceArguments inductanceArguments = MInductanceArguments(PrecisionArguments(2, 1, 1, 50, 32, 1),
                                                                    PrecisionArguments(2, 1, 1, 50, 20, 1));

    for (int i = 10; i >= 0; --i){
    //    printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                              0.05, 0.0,
                                                              acos(i * 0.1)));
    }
    printf("\n");
}

void testCoilMutualInductanceGeneralRectangularAndFilament()
{
    Coil primaryGeneral = Coil(0.04, 0.02, 0.01, 100);
    Coil secondaryGeneral = Coil(0.02, 1e-15, 1e-15, 1);

    MInductanceArguments inductanceArguments = MInductanceArguments(PrecisionArguments(2, 1, 1, 50, 24, 24),
                                                                    PrecisionArguments(2, 1, 1, 50, 1, 1));

    for (int i = 10; i >= 0; --i){
    //    printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                              0.0, 0.0,
                                                              acos(i * 0.1), inductanceArguments));
    }
    printf("\n");
}

