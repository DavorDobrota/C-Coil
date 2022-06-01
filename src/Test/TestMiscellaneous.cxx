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

    printf("%.15f, %.15f, %.15f\n\n", testCoil1.getMagneticMoment().magnitude(), testCoil1.getAverageWireThickness(), testCoil1.getResistance());

    testCoil1.setSineFrequency(100000);
    testCoil1.setCurrentDensity(500000);
    printf("%.15f, %.15f, %.15f\n\n", testCoil1.getMagneticMoment().magnitude(), testCoil1.getAverageWireThickness(), testCoil1.getResistance());

    testCoil1.setSineFrequency(0);
    testCoil1.setCurrent(1);

    vec3::CoordVector3 positionVector = vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.0);

    vec3::FieldVector3 vector = testCoil1.computeBFieldVector(positionVector);
    printf("%.25f %.25f\n", vector.x, testCoil1.computeBFieldZ(positionVector));
    printf("%.25f %.25f\n", vector.z, testCoil1.computeBFieldX(positionVector));
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
    coil1.setPositionAndOrientation(vec3::CoordVector3(), 0.0, 0.0);

    Coil coil2 = Coil(0.03, 0.03, 0.12, 3600);
    coil2.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, -0.02), 0, 0);

    Coil coil3 = Coil(0.03, 0.03, 0.12, 3600);
    coil3.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.0), M_PI_2, M_PI_2);

    int pointCount = 100;
    std::vector<vec3::CoordVector3> pointPositions1(pointCount), pointPositions2(pointCount);

    for (int i = 0; i < pointCount; ++i)
        pointPositions1[i] = vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.1 + i * 0.001);

    std::vector<vec3::FieldVector3> fieldVectors1, fieldVectors2;

    fieldVectors1 = coil1.computeAllBFieldComponents(pointPositions1);
    fieldVectors2 = coil2.computeAllBFieldComponents(pointPositions1);

    for (int i = 0; i < pointCount; ++i)
        printf("%f : %.15g %.15g\n", 0.1 + i * 0.001, fieldVectors1[i].z, fieldVectors2[i].z);
    printf("\n");

    for (int i = 0; i < pointCount; ++i)
    {
        pointPositions1[i] = vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.001 * i);
        pointPositions2[i] = vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.001 * i, 0.0);
    }

    fieldVectors1 = coil1.computeAllBFieldComponents(pointPositions1);
    fieldVectors2 = coil3.computeAllBFieldComponents(pointPositions2);

    for (int i = 0; i < pointCount; ++i)
        printf("%16.10g %16.10g %16.10g | %16.10g %16.10g %16.10g\n",
               fieldVectors1[i].x, fieldVectors1[i].y, fieldVectors1[i].z,
               fieldVectors2[i].x, fieldVectors2[i].y, fieldVectors2[i].z);
    printf("\n");
}

void testRawGPUPerformance()
{
    using namespace std::chrono;

    FILE *output = fopen("output.txt", "w");

    printf("Testing raw GPU performance for a given number of points\n\n");

    high_resolution_clock::time_point beginTime;
    double interval;
    double pointsPerSec;

    int maxPointsLog2 = 28;

    std::vector<float> zCoords;
    std::vector<float> rCoords;

    std::vector<float> potentialArr;

    std::vector<float> fieldHArr;
    std::vector<float> fieldZArr;

    std::vector<float> gradientRPArr;
    std::vector<float> gradientRRArr;
    std::vector<float> gradientRZArr;
    std::vector<float> gradientZZArr;

    printf("Vector potential performance\n");

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int numPoints = int(std::pow(2, i));

        zCoords.resize(numPoints);
        rCoords.resize(numPoints);
        potentialArr.resize(numPoints);

        for (int j = 0; j < numPoints; ++j)
        {
            zCoords[j] = j / 1000.0;
            rCoords[j] = 0.1 + j / 1000000.0;
        }

        beginTime = high_resolution_clock::now();
        Calculate_hardware_accelerated_a(rCoords.size(), &zCoords[0], &rCoords[0],
                                         1000000, 0.03, 0.12, 0.03,
                                         16, 16, 16,
                                         &potentialArr[0]);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = numPoints / interval;

        printf("%8d : %.1f\n", numPoints, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", numPoints, 0.001 * pointsPerSec);

        rCoords.clear();
        zCoords.clear();
        potentialArr.clear();
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic field performance\n");

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int numPoints = int(std::pow(2, i));

        zCoords.resize(numPoints);
        rCoords.resize(numPoints);
        fieldHArr.resize(numPoints);
        fieldZArr.resize(numPoints);

        for (int j = 0; j < numPoints; ++j)
        {
            zCoords[j] = j / 1000.0;
            rCoords[j] = 0.1 + j / 1000000.0;
        }

        beginTime = high_resolution_clock::now();
        Calculate_hardware_accelerated_b(rCoords.size(), &zCoords[0], &rCoords[0],
                                         1000000, 0.03, 0.12, 0.03,
                                         16, 16, 16,
                                         &fieldHArr[0], &fieldZArr[0]);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = numPoints / interval;

        printf("%8d : %.1f\n", numPoints, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", numPoints, 0.001 * pointsPerSec);

        rCoords.clear();
        zCoords.clear();
        fieldHArr.clear();
        fieldZArr.clear();
    }
    printf("\n");
    fprintf(output, "\n");

    printf("Magnetic gradient performance\n");

    for (int i = 0; i <= maxPointsLog2; ++i)
    {
        int numPoints = int(std::pow(2, i));

        zCoords.resize(numPoints);
        rCoords.resize(numPoints);
        gradientRPArr.resize(numPoints);
        gradientRRArr.resize(numPoints);
        gradientRZArr.resize(numPoints);
        gradientZZArr.resize(numPoints);

        for (int j = 0; j < numPoints; ++j)
        {
            zCoords[j] = j / 1000.0;
            rCoords[j] = 0.1 + j / 1000000.0;
        }

        beginTime = high_resolution_clock::now();
        Calculate_hardware_accelerated_g(rCoords.size(), &zCoords[0], &rCoords[0],
                                         1000000, 0.03, 0.12, 0.03,
                                         16, 16, 16,
                                         &gradientRPArr[0], &gradientRRArr[0],
                                         &gradientRZArr[0], &gradientZZArr[0]);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - beginTime).count();
        pointsPerSec = numPoints / interval;

        printf("%8d : %.1f\n", numPoints, 0.001 * pointsPerSec);
        fprintf(output, "%8d\t%.7g\n", numPoints, 0.001 * pointsPerSec);

        rCoords.clear();
        zCoords.clear();
        gradientRPArr.clear();
        gradientRRArr.clear();
        gradientRZArr.clear();
        gradientZZArr.clear();
    }
    printf("\n");

    fclose(output);
}
