#include <cstdio>
#include <vector>

#include "Test.h"
#include "Polynomial.h"
#include "Coil.h"
#include "OldCoil.h"
#include "ComputeMethod.h"
#include "Tensor.h"

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
    printf("%.25f %.25f\n", vector.xComponent, testCoil1.computeBFieldZ(positionVector));
    printf("%.25f %.25f\n", vector.zComponent, testCoil1.computeBFieldH(positionVector));
}

void testMethodPrecisionCompareCPUvsGPU()
{
    Coil testCoil = Coil(0.03, 0.03, 0.12, 3600, PrecisionFactor(6.0), 12);

    int pointCount = 2000;
    double radius = 0.075;

    std::vector<vec3::CoordVector3> positionValues(pointCount);

    for (int i = 0; i < pointCount; i++)
        positionValues[i] = vec3::CoordVector3(vec3::CYLINDRICAL, radius, i * M_PI / pointCount, 0.0);

    std::vector<double> cpuPotential;
    std::vector<double> gpuPotential;

    std::vector<vec3::FieldVector3> cpuFieldVectors;
    std::vector<vec3::FieldVector3> gpuFieldVectors;

    cpuPotential = testCoil.computeAllAPotentialAbs(positionValues, CPU_ST);
    cpuFieldVectors = testCoil.computeAllBFieldComponents(positionValues, CPU_ST);

    gpuPotential = testCoil.computeAllAPotentialAbs(positionValues, GPU);
    gpuFieldVectors = testCoil.computeAllBFieldComponents(positionValues, GPU);

    FILE *output = fopen("output.txt", "w");

    for (int i = 0; i < pointCount; ++i)
    {
        fprintf(output, "%.20f\t%.20f\t%.20f\t%.20f\t%.20f\t%.20f\n",
               cpuFieldVectors[i].xComponent, gpuFieldVectors[i].xComponent,
               cpuFieldVectors[i].zComponent, gpuFieldVectors[i].zComponent,
               cpuPotential[i], gpuPotential[i]);
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

void testVector3()
{
    for (int i = 0; i < 100; ++i)
    {
        vec3::CoordVector3 vector = vec3::CoordVector3(vec3::CYLINDRICAL, 0.0, 1, i * M_PI / 25);
        printf("%.6f %.6f %.6f | ", vector.component1, vector.component2, vector.component3);
        vector.convertToCartesian();
        printf("%.6f %.6f %.6f | ", vector.component1, vector.component2, vector.component3);
        vector.convertToSpherical();
        printf("%.6f %.6f %.6f\n", vector.component1, vector.component2, vector.component3);
    }
    printf("\n");

    for (int i = 0; i < 100; ++i)
    {
        vec3::CoordVector3 vector = vec3::CoordVector3(vec3::CYLINDRICAL, 1, 1, i * M_PI / 25);
        printf("%.6f %.6f %.6f | ", vector.component1, vector.component2, vector.component3);
        vector.convertToCartesian();
        printf("%.6f %.6f %.6f | ", vector.component1, vector.component2, vector.component3);
        vector.convertToSpherical();
        printf("%.6f %.6f %.6f\n", vector.component1, vector.component2, vector.component3);
    }
    printf("\n");

    for (int i = 0; i < 100; ++i)
    {
        vec3::CoordVector3 vector = vec3::CoordVector3(vec3::SPHERICAL, 1, M_PI/2, -M_PI + i * M_PI / 25);
        printf("%.6f %.6f %.6f | ", vector.component1, vector.component2, vector.component3);
        vector.convertToCartesian();
        printf("%.6f %.6f %.6f | ", vector.component1, vector.component2, vector.component3);
        vector.convertToCylindrical();
        printf("%.6f %.6f %.6f\n", vector.component1, vector.component2, vector.component3);
    }
    printf("\n");

    for (int i = 0; i < 100; ++i)
    {
        vec3::CoordVector3 vector = vec3::CoordVector3(vec3::SPHERICAL, 1, i * M_PI / 25, 0.0);
        printf("%.6f %.6f %.6f | ", vector.component1, vector.component2, vector.component3);
        vector.convertToCartesian();
        printf("%.6f %.6f %.6f | ", vector.component1, vector.component2, vector.component3);
        vector.convertToCylindrical();
        printf("%.6f %.6f %.6f\n", vector.component1, vector.component2, vector.component3);
    }
    printf("\n");
}

void testCoilPositionAndRotation()
{
    Coil coil1 = Coil(0.03, 0.03, 0.12, 3600);
    coil1.setPositionAndOrientation(vec3::CoordVector3(), 0.0, 0.0);

    Coil coil2 = Coil(0.03, 0.03, 0.12, 3600);
    coil2.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, -0.02), 0, 0);

    Coil coil3 = Coil(0.03, 0.03, 0.12, 3600);
    coil3.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.0), 0, M_PI_2);

    int numPoints = 100;
    std::vector<vec3::CoordVector3> pointPositions1(numPoints), pointPositions2(numPoints);

    for (int i = 0; i < numPoints; ++i)
        pointPositions1[i] = vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.1 + i * 0.001);

    std::vector<vec3::FieldVector3> fieldVectors1, fieldVectors2;

    fieldVectors1 = coil1.computeAllBFieldComponents(pointPositions1);
    fieldVectors2 = coil2.computeAllBFieldComponents(pointPositions1);

    for (int i = 0; i < numPoints; ++i)
        printf("%f : %.15g %.15g\n", 0.1 + i * 0.001, fieldVectors1[i].zComponent, fieldVectors2[i].zComponent);
    printf("\n");

    for (int i = 0; i < numPoints; ++i)
    {
        pointPositions1[i] = vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.001 * i);
        pointPositions2[i] = vec3::CoordVector3(vec3::CARTESIAN, 0.00 * i, 0.00 * i, 0.001 * i);
    }

    fieldVectors1 = coil1.computeAllBFieldComponents(pointPositions1);
    fieldVectors2 = coil3.computeAllBFieldComponents(pointPositions2);

    for (int i = 0; i < numPoints; ++i)
        printf("%15.10g %15.10g %.10g | %15.10g %15.10g %15.10g\n",
               fieldVectors1[i].xComponent, fieldVectors1[i].yComponent, fieldVectors1[i].zComponent,
               fieldVectors2[i].xComponent, fieldVectors2[i].yComponent, fieldVectors2[i].zComponent);
    printf("\n");
}
