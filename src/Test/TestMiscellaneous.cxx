#include "Test.h"
#include "Coil.h"
#include "ComputeMethod.h"
#include "Tensor.h"
#include "CustomMath.h"
#include "CoilGroup.h"

#include <cmath>
#include <cstdio>
#include <vector>
#include <chrono>


void testFunctionPerformance()
{
    using namespace std::chrono;

    int nOps = 200'000'000;
    double temp, interval;
    high_resolution_clock::time_point begin_time;

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= nOps; ++i)
    {
        temp += std::sqrt(i);
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("sqrt  : %.1f MOps/s\n", 1e-6 * nOps / interval);

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= nOps; ++i)
    {
        temp += customMath::ln(i);
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("ln    : %.1f MOps/s\n", 1e-6 * nOps / interval);

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= nOps; ++i)
    {
        temp += customMath::lnf(i);
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("lnf   : %.1f MOps/s\n", 1e-6 * nOps / interval);

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= nOps; ++i)
    {
        temp += std::log10(400000.0 / i);
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("log10 : %.1f MOps/s\n", 1e-6 * nOps / interval);

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= nOps; ++i)
    {
        temp += std::log10((float) (i));
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("log10f: %.1f MOps/s\n", 1e-6 * nOps / interval);

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= nOps; ++i)
    {
        temp += std::log(i);
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("log    : %.1f MOps/s\n", 1e-6 * nOps / interval);

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= nOps; ++i)
    {
        temp += customMath::cos( 1.0 / i);
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("custom cos: %.1f MOps/s\n", 1e-6 * nOps / interval);

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= nOps; ++i)
    {
        temp += std::cos( 1.0 / i);
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("std cos   : %.1f MOps/s\n", 1e-6 * nOps / interval);
}

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
    printf("%.25f %.25f\n", vector.xComponent, testCoil1.computeBFieldZ(positionVector));
    printf("%.25f %.25f\n", vector.zComponent, testCoil1.computeBFieldX(positionVector));
}

void testMethodPrecisionCompareCPUvsGPU()
{
    Coil testCoil = Coil(0.03, 0.03, 0.12, 3600, PrecisionFactor(6.0), 12);

    int pointCount = 2000;
    double radius = 0.075;

    std::vector<vec3::CoordVector3> positionValues(pointCount);

    for (int i = 0; i < pointCount; i++)
        positionValues[i] = vec3::CoordVector3(vec3::SPHERICAL, radius, i * M_PI / pointCount, 0.0);

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
	Coil primary = Coil(0.071335, 0.01397, 0.142748, 1142);
	Coil secondary = Coil(0.096945, 0.041529, 0.02413, 516);

	secondary.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.30988, 0.0, 0.07366));

	printf("%.15f\n", Coil::computeMutualInductance(primary, secondary));
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
    coil3.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.0), M_PI_2, M_PI_2);

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
        pointPositions2[i] = vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.001 * i, 0.0);
    }

    fieldVectors1 = coil1.computeAllBFieldComponents(pointPositions1);
    fieldVectors2 = coil3.computeAllBFieldComponents(pointPositions2);

    for (int i = 0; i < numPoints; ++i)
        printf("%16.10g %16.10g %16.10g | %16.10g %16.10g %16.10g\n",
               fieldVectors1[i].xComponent, fieldVectors1[i].yComponent, fieldVectors1[i].zComponent,
               fieldVectors2[i].xComponent, fieldVectors2[i].yComponent, fieldVectors2[i].zComponent);
    printf("\n");
}

void testCoilGroupMTD(int numCoils, int numPoints, int threadCount, bool print)
{
    double torusRadius = 1.0;

    CoilGroup torusGroup = CoilGroup();
    std::vector<vec3::CoordVector3> fieldPoints(numPoints);
    std::vector<vec3::FieldVector3> computedBField;

    for (int i = 0; i < numPoints; ++i)
        fieldPoints[i] = vec3::CoordVector3(vec3::CYLINDRICAL, 0.0, torusRadius, 2*M_PI * i / numPoints);

    for (int i = 0; i < numCoils; ++i)
    {
        Coil tempCoil = Coil(torusRadius / 10.0, torusRadius / 100.0, torusRadius / 100.0, 10000, 10);
        tempCoil.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CYLINDRICAL, 0.0, torusRadius, 2*M_PI * i / numCoils),
                M_PI_2, 2*M_PI * i / numCoils + M_PI_2);
        torusGroup.addCoil(tempCoil);
    }

    torusGroup.setThreadCount(threadCount);

    computedBField = torusGroup.computeAllBFieldComponents(fieldPoints, CPU_MT);

    if (print)
        for (int i = 0; i < numPoints; ++i)
            printf("%.15g\n",
                   std::sqrt(computedBField[i].xComponent * computedBField[i].xComponent +
                   computedBField[i].yComponent * computedBField[i].yComponent));
}

void testCoilGroupMTvsMTD(int threadCount, int numPoints)
{
    using namespace std::chrono;

    high_resolution_clock::time_point begin_time;
    double interval;

    int coilCount1 = 2 * threadCount;
    int coilCount2 = 32 * threadCount;

    begin_time = high_resolution_clock::now();
    testCoilGroupMTD(coilCount1, numPoints, threadCount, false);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MT perf  : %.1f coils/s\n", coilCount1 / interval);

    begin_time = high_resolution_clock::now();
    testCoilGroupMTD(coilCount2, numPoints, threadCount, false);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MTD perf : %.1f coils/s", coilCount2 / interval);
}
