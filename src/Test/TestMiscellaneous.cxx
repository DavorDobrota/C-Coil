#include "Test.h"
#include "Coil.h"
#include "ComputeMethod.h"
#include "Tensor.h"
#include "hardware_acceleration.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>
#include <vector>
#include <chrono>


void testNewCoilParameters()
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

void testVector3()
{
    for (int i = 0; i < 100; ++i)
    {
        vec3::CoordVector3 vector = vec3::CoordVector3(vec3::CYLINDRICAL, 0.0, 1, i * M_PI / 25);
        printf("%.6f %.6f %.6f | ", vector.comp1, vector.comp2, vector.comp3);
        vector.convertToCartesian();
        printf("%.6f %.6f %.6f | ", vector.comp1, vector.comp2, vector.comp3);
        vector.convertToSpherical();
        printf("%.6f %.6f %.6f\n", vector.comp1, vector.comp2, vector.comp3);
    }
    printf("\n");

    for (int i = 0; i < 100; ++i)
    {
        vec3::CoordVector3 vector = vec3::CoordVector3(vec3::CYLINDRICAL, 1, 1, i * M_PI / 25);
        printf("%.6f %.6f %.6f | ", vector.comp1, vector.comp2, vector.comp3);
        vector.convertToCartesian();
        printf("%.6f %.6f %.6f | ", vector.comp1, vector.comp2, vector.comp3);
        vector.convertToSpherical();
        printf("%.6f %.6f %.6f\n", vector.comp1, vector.comp2, vector.comp3);
    }
    printf("\n");

    for (int i = 0; i < 100; ++i)
    {
        vec3::CoordVector3 vector = vec3::CoordVector3(vec3::SPHERICAL, 1, M_PI/2, -M_PI + i * M_PI / 25);
        printf("%.6f %.6f %.6f | ", vector.comp1, vector.comp2, vector.comp3);
        vector.convertToCartesian();
        printf("%.6f %.6f %.6f | ", vector.comp1, vector.comp2, vector.comp3);
        vector.convertToCylindrical();
        printf("%.6f %.6f %.6f\n", vector.comp1, vector.comp2, vector.comp3);
    }
    printf("\n");

    for (int i = 0; i < 100; ++i)
    {
        vec3::CoordVector3 vector = vec3::CoordVector3(vec3::SPHERICAL, 1, i * M_PI / 25, 0.0);
        printf("%.6f %.6f %.6f | ", vector.comp1, vector.comp2, vector.comp3);
        vector.convertToCartesian();
        printf("%.6f %.6f %.6f | ", vector.comp1, vector.comp2, vector.comp3);
        vector.convertToCylindrical();
        printf("%.6f %.6f %.6f\n", vector.comp1, vector.comp2, vector.comp3);
    }
    printf("\n");
}

void testCoilPositionAndRotation()
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

