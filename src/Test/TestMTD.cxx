#include "Test.h"
#include "Coil.h"
#include "ComputeMethod.h"
#include "Tensor.h"
#include "CoilGroup.h"

#include <cmath>
#include <cstdio>
#include <vector>
#include <chrono>

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
    torusGroup.setDefaultPrecisionFactor(PrecisionFactor(3.0));

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
    int coilCount2 = 8 * threadCount;

    begin_time = high_resolution_clock::now();
    testCoilGroupMTD(coilCount1, numPoints, threadCount, false);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MT perf  : %.0f kPoints/s\n", 1e-3 * coilCount1 * numPoints / interval);

    begin_time = high_resolution_clock::now();
    testCoilGroupMTD(coilCount2, numPoints, threadCount, false);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MTD perf : %.0f kPoints/s", 1e-3 * coilCount2 * numPoints / interval);
}

void testCoilGroupMTDFields(int threadCount)
{
    using namespace std::chrono;

    high_resolution_clock::time_point begin_time;
    double interval;

    int numCoils = 8 * threadCount;
    const int numPoints = 1'000;

    CoilGroup coilGroup = CoilGroup();
    Coil referenceCoil = Coil(0.1, 0.1, 0.1, 10000);

    for (int i = 1; i <= numCoils; ++i)
    {
        Coil tempCoil = Coil(0.1, 0.1, 0.1, 10000);
        tempCoil.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.15*i),
                0.0, 0.0);
        coilGroup.addCoil(tempCoil);
    }

    std::vector<vec3::CoordVector3> referencePoints(numPoints);
    std::vector<vec3::FieldVector3> computedAPotential;
    std::vector<vec3::FieldVector3> computedBField;
    std::vector<vec3::Matrix3> computedBGradient;

    for (int i = 0; i < numPoints; ++i)
        referencePoints[i] = vec3::CoordVector3(vec3::CARTESIAN, 0.1, 1.0 * i / numPoints, -0.1);

    computedAPotential = coilGroup.computeAllAPotentialComponents(referencePoints, CPU_ST);
    printf("%.15g\n", computedAPotential[numPoints / 2].xComponent);
    computedAPotential = coilGroup.computeAllAPotentialComponents(referencePoints, CPU_MT);
    printf("%.15g\n", computedAPotential[numPoints / 2].xComponent);
    printf("\n");

    computedBField = coilGroup.computeAllBFieldComponents(referencePoints, CPU_ST);
    printf("%.15g\n", computedBField[numPoints / 2].xComponent);
    computedBField = coilGroup.computeAllBFieldComponents(referencePoints, CPU_MT);
    printf("%.15g\n", computedBField[numPoints / 2].xComponent);
    printf("\n");

    computedBGradient = coilGroup.computeAllBGradientTensors(referencePoints, CPU_ST);
    printf("%.15g\n", computedBGradient[numPoints / 2].xyElement);
    computedBGradient = coilGroup.computeAllBGradientTensors(referencePoints, CPU_MT);
    printf("%.15g\n", computedBGradient[numPoints / 2].xyElement);
    printf("\n");
}

void testCoilGroupMTDInductanceAndForce(int threadCount)
{
    using namespace std::chrono;

    high_resolution_clock::time_point begin_time;
    double interval;

    int numCoils = 8 * threadCount;
    CoilGroup coilGroup = CoilGroup();
    Coil referenceCoil = Coil(0.1, 0.1, 0.1, 10000);

    for (int i = 1; i <= numCoils; ++i)
    {
        Coil tempCoil = Coil(0.1, 0.1, 0.1, 10000);
        tempCoil.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.15*i),
                0.0, 0.0);
        coilGroup.addCoil(tempCoil);
    }

    begin_time = high_resolution_clock::now();
    printf("%.15g\n", coilGroup.computeMutualInductance(referenceCoil, PrecisionFactor(5.0), CPU_ST));
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("ST  perf  : %.0f coils/s\n", numCoils / interval);


    begin_time = high_resolution_clock::now();
    printf("%.15g\n", coilGroup.computeMutualInductance(referenceCoil, PrecisionFactor(5.0), CPU_MT));
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MTD perf  : %.0f coils/s\n", numCoils / interval);

    printf("\n");

    begin_time = high_resolution_clock::now();
    printf("%.15g\n",
           coilGroup.computeAmpereForce(referenceCoil, PrecisionFactor(5.0), CPU_ST).first.zComponent);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("ST  perf  : %.0f coils/s\n", numCoils / interval);


    begin_time = high_resolution_clock::now();
    printf("%.15g\n",
           coilGroup.computeAmpereForce(referenceCoil,PrecisionFactor(5.0), CPU_MT).first.zComponent);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MTD perf  : %.0f coils/s\n", numCoils / interval);
}
