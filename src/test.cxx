
#include <cstdio>
#include <vector>

#include "test.h"
#include "Polynomial.h"
#include "Coil.h"
#include "OldCoil.h"
#include "Precision.h"
#include "LegendreMatrix.h"
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
    printf("combined B  : %.0f kInc/s\n", 0.001 / (float(clock() - begin_time1) / CLOCKS_PER_SEC / numOperations));

    double temp2;
    clock_t begin_time2 = clock();
    for (int i = 0; i < nOp; ++i){
        temp2 = testCoil.computeBFieldH(0.0, 0.0);
    }
    printf("field Bh    : %.0f kInc/s\n", 0.001 / (float(clock() - begin_time2) / CLOCKS_PER_SEC / numOperations));

    double temp3;
    clock_t begin_time3 = clock();
    for (int i = 0; i < nOp; ++i){
        temp3 = testCoil.computeBFieldH(0.0, 0.0);
    }
    printf("field Bz    : %.0f kInc/s\n", 0.001 / (float(clock() - begin_time3) / CLOCKS_PER_SEC / numOperations));

    double temp4;
    clock_t begin_time4 = clock();
    for (int i = 0; i < nOp; ++i){
        temp4 = testCoil.computeAPotentialAbs(i*0.000001, 0.0);
    }
    printf("potential A : %.0f kInc/s\n", 0.001 / (float(clock() - begin_time4) / CLOCKS_PER_SEC / numOperations));

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
    printf("combined  B CPU : %.0f kInc/s\n", 0.001 / (float(clock() - begin_time11) / CLOCKS_PER_SEC / numOperations));

    clock_t begin_time12 = clock();
    testCoil.computeAllAPotentialAbs(cylindricalZArr, cylindricalRArr,
                                     singlePotential, CPU_ST);
    printf("Potential A CPU : %.0f kInc/s\n", 0.001 / (float(clock() - begin_time12) / CLOCKS_PER_SEC / numOperations));

    testCoil.computeAllBFieldComponents(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                                        acceleratedResultsX, acceleratedResultsY, acceleratedResultsZ,
                                        GPU);

    clock_t begin_time13 = clock();
    testCoil.computeAllBFieldComponents(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                                        acceleratedResultsX, acceleratedResultsY, acceleratedResultsZ,
                                        GPU);
    printf("combined  B GPU : %.0f kInc/s\n", 0.001 / (float(clock() - begin_time13) / CLOCKS_PER_SEC / numOperationsGpu));

    clock_t begin_time14 = clock();
    testCoil.computeAllAPotentialAbs(cylindricalZArr, cylindricalRArr,
                                     acceleratedPotential, GPU);
    printf("Potential A GPU : %.0f kInc/s\n", 0.001 / (float(clock() - begin_time14) / CLOCKS_PER_SEC / numOperationsGpu));
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

    printf("%.20f\n\n", primary.computeMutualInductance(0.2, secondary));
    printf("%.20f\n\n", primary.computeMutualInductance(0.2, secondary, GPU));

    FILE *input = fopen("values.txt", "r");
    FILE *output = fopen("output.txt", "w");

    double Rt1, at1, bt1; int Nt1;
    double Rt2, at2, bt2; int Nt2;
    double distance;
    double Temp;

    while (fscanf(input, "%lf %lf %lf %d %lf %lf %lf %d %lf", &Rt1, &at1, &bt1, &Nt1, &Rt2, &at2, &bt2, &Nt2, &distance) == 9)
    {
        printf("%f %f %f %d %f %f %f %d %f\n", Rt1, at1, bt1, Nt1, Rt2, at2, bt2, Nt2, distance);

        for (int i = 1; i <= 9; ++i)
        {
            Coil prim1 = Coil(Rt1, at1, bt1, Nt1, PrecisionArguments(i, 1, 1, 12, 12, 12));
            Coil sec1 = Coil(Rt2, at2, bt2, Nt2);
            Temp = prim1.computeMutualInductance(distance, sec1);
            printf("%.18f\n", Temp);
            fprintf(output, "%.20f\t", Temp);
        }

        printf("====================================================================================\n");
        fprintf(output, "\n");
    }

    fclose(input);
    fclose(output);
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
		for (Type i = 1; i <= 9.0; i += 1.0){
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
