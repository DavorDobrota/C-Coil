#include "Test.h"
#include "Tensor.h"
#include "CoilGroup.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>


void Test::testNewCoilParameters()
{
    Coil testCoil1 = Coil(0.03, 0.03, 0.12, 3600);

    printf("%.15f, %.15f, %.15f\n\n", testCoil1.getCurrentDensity(), testCoil1.getWireResistivity(), testCoil1.getSineFrequency());

    printf("%.15f, %.15f, %.15f\n\n", testCoil1.getMagneticMoment().abs(), testCoil1.getAverageWireThickness(), testCoil1.getResistance());

    testCoil1.setSineFrequency(100000);
    testCoil1.setCurrentDensity(500000);
    printf("%.15f, %.15f, %.15f\n\n", testCoil1.getMagneticMoment().abs(), testCoil1.getAverageWireThickness(), testCoil1.getResistance());

    testCoil1.setSineFrequency(0);
    testCoil1.setCurrent(1);

    vec3::Vector3 positionVector = vec3::Vector3(0.0, 0.0, 0.0);

    vec3::Vector3 vector = testCoil1.computeBFieldVector(positionVector);
    printf("%.25f %.25f\n", vector.x, testCoil1.computeBFieldVector(positionVector).x);
    printf("%.25f %.25f\n", vector.z, testCoil1.computeBFieldVector(positionVector).z);
}

void Test::testCoilPositionAndRotation()
{
    Coil coil1 = Coil(0.03, 0.03, 0.12, 3600);
    coil1.setPositionAndOrientation(vec3::Vector3(), 0.0, 0.0);

    Coil coil2 = Coil(0.03, 0.03, 0.12, 3600);
    coil2.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, -0.02), 0, 0);

    Coil coil3 = Coil(0.03, 0.03, 0.12, 3600);
    coil3.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.0), M_PI_2, M_PI_2);

    int pointCount = 100;
    vec3::Vector3Array pointPositions1(pointCount), pointPositions2(pointCount);

    for (int i = 0; i < pointCount; ++i)
        pointPositions1[i] = vec3::Vector3(0.0, 0.0, 0.1 + i * 0.001);

    vec3::Vector3Array fieldVectors1, fieldVectors2;

    fieldVectors1 = coil1.computeAllBFieldVectors(pointPositions1);
    fieldVectors2 = coil2.computeAllBFieldVectors(pointPositions1);

    for (int i = 0; i < pointCount; ++i)
        printf("%f : %.15g %.15g\n", 0.1 + i * 0.001, fieldVectors1[i].z, fieldVectors2[i].z);
    printf("\n");

    for (int i = 0; i < pointCount; ++i)
    {
        pointPositions1[i] = vec3::Vector3(0.0, 0.0, 0.001 * i);
        pointPositions2[i] = vec3::Vector3(0.0, 0.001 * i, 0.0);
    }

    fieldVectors1 = coil1.computeAllBFieldVectors(pointPositions1);
    fieldVectors2 = coil3.computeAllBFieldVectors(pointPositions2);

    for (int i = 0; i < pointCount; ++i)
        printf("%16.10g %16.10g %16.10g | %16.10g %16.10g %16.10g\n",
               fieldVectors1[i].x, fieldVectors1[i].y, fieldVectors1[i].z,
               fieldVectors2[i].x, fieldVectors2[i].y, fieldVectors2[i].z);
    printf("\n");
}

void Test::testCoilGroupComputeAllMTD()
{
    int coilCount = 8 * 8;
    const int pointCount = 1'000;

    CoilGroup coilGroup = CoilGroup();

    for (int i = 1; i <= coilCount; ++i)
    {
        coilGroup.addCoil(
            0.1, 0.1, 0.1, 10000, 10,
            PrecisionFactor(), 8,
            vec3::Vector3(0.0, 0.0, 0.15*i), 0.0, 0.0
        );
    }

    printf("Testing if ST and MTD methods for field computation return the same values\n\n");

    vec3::Vector3Array referencePoints(pointCount);
    vec3::Vector3Array computedAPotential;
    vec3::Vector3Array computedBField;
    vec3::Matrix3Array computedBGradient;

    for (int i = 0; i < pointCount; ++i)
        referencePoints[i] = vec3::Vector3(0.1, 1.0 * i / pointCount, -0.1);

    computedAPotential = coilGroup.computeAllAPotentialVectors(referencePoints, CPU_ST);
    printf("ST  potential: %.15g\n", computedAPotential[pointCount / 2].x);
    computedAPotential = coilGroup.computeAllAPotentialVectors(referencePoints, CPU_MT);
    printf("MTD potential: %.15g\n", computedAPotential[pointCount / 2].x);
    printf("\n");

    computedBField = coilGroup.computeAllBFieldVectors(referencePoints, CPU_ST);
    printf("ST  field      %.15g\n", computedBField[pointCount / 2].x);
    computedBField = coilGroup.computeAllBFieldVectors(referencePoints, CPU_MT);
    printf("MTD field      %.15g\n", computedBField[pointCount / 2].x);
    printf("\n");

    computedBGradient = coilGroup.computeAllBGradientMatrices(referencePoints, CPU_ST);
    printf("ST  gradient:  %.15g\n", computedBGradient[pointCount / 2].xy);
    computedBGradient = coilGroup.computeAllBGradientMatrices(referencePoints, CPU_MT);
    printf("MTD gradient:  %.15g\n", computedBGradient[pointCount / 2].xy);
    printf("\n");
}
