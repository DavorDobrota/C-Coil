#include <cstdio>
#include <vector>

#include "Test.h"
#include "Polynomial.h"
#include "Coil.h"
#include "OldCoil.h"
#include "ComputeMethod.h"
#include "Vector3.h"

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

    vec3::CoordVector3 positionVector = vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.0);

    vec3::FieldVector3 vector = testCoil1.computeBFieldVector(positionVector);
    printf("%.25f %.25f\n", vector.zComponent, testCoil1.computeBFieldZ(positionVector));
    printf("%.25f %.25f\n", vector.xComponent, testCoil1.computeBFieldH(positionVector));
}

void testMethodPrecisionCompareCPUvsGPU()
{
    Coil testCoil = Coil(0.03, 0.03, 0.12, 3600, PrecisionFactor(6.0), 12);

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
