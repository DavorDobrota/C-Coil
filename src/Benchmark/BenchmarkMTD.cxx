#include "Benchmark.h"
#include "Coil.h"
#include "ComputeMethod.h"
#include "Tensor.h"
#include "Compare.h"
#include "CoilGroup.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>
#include <vector>
#include <chrono>


void benchCoilGroupComputeAllFields(PrecisionFactor precisionFactor, int numCoils, int opCount, int threadCount)
{
    using namespace std::chrono;

    high_resolution_clock::time_point begin_time;
    double interval;

    printf("Benchmarking MT and GPU field compute performance for %d coils in %d points and precision factor %.1f\n\n",
           numCoils, opCount, precisionFactor.relativePrecision);

    double torusRadius = 1.0;

    CoilGroup torusGroupThick = CoilGroup();
    CoilGroup torusGroupFlat = CoilGroup();

    vec3::Vector3Array fieldPoints(opCount);
    vec3::Vector3Array computedAPotential;
    vec3::Vector3Array computedBField;
    vec3::Matrix3Array computedGradient;

    for (int i = 0; i < opCount; ++i)
        fieldPoints[i] = vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / opCount);

    for (int i = 0; i < numCoils; ++i)
    {
        Coil tempCoil = Coil(torusRadius / 10.0, torusRadius / 100.0, torusRadius / 100.0, 10000, 10);
        tempCoil.setPositionAndOrientation(
                vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / numCoils),
                M_PI_2, 2*M_PI * i / numCoils + M_PI_2);
        torusGroupThick.addCoil(tempCoil);
    }
    torusGroupThick.setThreadCount(threadCount);
    torusGroupThick.setDefaultPrecisionFactor(precisionFactor);

    for (int i = 0; i < numCoils; ++i)
    {
        Coil tempCoil = Coil(torusRadius / 10.0, torusRadius / 100.0, 0.0, 100, 10);
        tempCoil.setPositionAndOrientation(
                vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / numCoils),
                M_PI_2, 2*M_PI * i / numCoils + M_PI_2);
        torusGroupFlat.addCoil(tempCoil);
    }
    torusGroupFlat.setThreadCount(threadCount);
    torusGroupFlat.setDefaultPrecisionFactor(precisionFactor);

    // MT slow tests
    begin_time = high_resolution_clock::now();
    computedAPotential = torusGroupFlat.computeAllAPotentialVectors(fieldPoints, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_MT slow : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * numCoils / interval);

    begin_time = high_resolution_clock::now();
    computedBField = torusGroupFlat.computeAllBFieldVectors(fieldPoints, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Field     B CPU_MT slow : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * numCoils / interval);

    begin_time = high_resolution_clock::now();
    computedGradient = torusGroupFlat.computeAllBGradientMatrices(fieldPoints, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G CPU_MT slow : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * numCoils / interval);

    printf("\n");

    // MT fast tests
    begin_time = high_resolution_clock::now();
    computedAPotential = torusGroupThick.computeAllAPotentialVectors(fieldPoints, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_MT fast : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * numCoils / interval);

    begin_time = high_resolution_clock::now();
    computedBField = torusGroupThick.computeAllBFieldVectors(fieldPoints, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Field     B CPU_MT fast : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * numCoils / interval);

    begin_time = high_resolution_clock::now();
    computedGradient = torusGroupThick.computeAllBGradientMatrices(fieldPoints, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G CPU_MT fast : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * numCoils / interval);

    printf("\n");

    // GPU slow tests
    computedAPotential = torusGroupThick.computeAllAPotentialVectors(fieldPoints, GPU); // warmup

    begin_time = high_resolution_clock::now();
    computedAPotential = torusGroupFlat.computeAllAPotentialVectors(fieldPoints, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A GPU    slow : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * numCoils / interval);

    begin_time = high_resolution_clock::now();
    computedBField = torusGroupFlat.computeAllBFieldVectors(fieldPoints, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Field     B GPU    slow : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * numCoils / interval);

    begin_time = high_resolution_clock::now();
    computedGradient = torusGroupFlat.computeAllBGradientMatrices(fieldPoints, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G GPU    slow : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * numCoils / interval);

    printf("\n");

    // GPU fast tests

    begin_time = high_resolution_clock::now();
    computedAPotential = torusGroupThick.computeAllAPotentialVectors(fieldPoints, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A GPU    fast : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * numCoils / interval);

    begin_time = high_resolution_clock::now();
    computedBField = torusGroupThick.computeAllBFieldVectors(fieldPoints, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Field     B GPU    fast : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * numCoils / interval);

    begin_time = high_resolution_clock::now();
    computedGradient = torusGroupThick.computeAllBGradientMatrices(fieldPoints, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G GPU    fast : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * numCoils / interval);

    printf("\n");
}

void benchCoilGroupComputeAllFieldsGPU(int numCoils, int opCount)
{
    using namespace std::chrono;

    high_resolution_clock::time_point begin_time;
    double interval;

    printf("Benchmarking GPU field compute performance for %d coils in %d points\n\n", numCoils, opCount);

    double torusRadius = 1.0;

    CoilGroup torusGroupThick = CoilGroup();
    CoilGroup torusGroupFlat = CoilGroup();

    vec3::Vector3Array fieldPoints(opCount);
    vec3::Vector3Array computedAPotential;
    vec3::Vector3Array computedBField;
    vec3::Matrix3Array computedGradient;

    for (int i = 0; i < opCount; ++i)
        fieldPoints[i] = vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / opCount);

    for (int i = 0; i < numCoils; ++i)
    {
        Coil tempCoil = Coil(torusRadius / 10.0, torusRadius / 100.0, torusRadius / 100.0, 10000, 10);
        tempCoil.setPositionAndOrientation(
                vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / numCoils),
                M_PI_2, 2*M_PI * i / numCoils + M_PI_2);
        torusGroupThick.addCoil(tempCoil);
    }

    for (int i = 0; i < numCoils; ++i)
    {
        Coil tempCoil = Coil(torusRadius / 10.0, torusRadius / 100.0, 0.0, 100, 10);
        tempCoil.setPositionAndOrientation(
                vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / numCoils),
                M_PI_2, 2*M_PI * i / numCoils + M_PI_2);
        torusGroupFlat.addCoil(tempCoil);
    }

    computedAPotential = torusGroupThick.computeAllAPotentialVectors(fieldPoints, GPU); // warmup

    for (int i = 1; i <= 8; ++i)
    {
        torusGroupFlat.setDefaultPrecisionFactor(PrecisionFactor(double(i)));
        torusGroupThick.setDefaultPrecisionFactor(PrecisionFactor(double(i)));

        printf("Precision factor %.1f\n", double(i));

        begin_time = high_resolution_clock::now();
        computedAPotential = torusGroupFlat.computeAllAPotentialVectors(fieldPoints, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("Potential A  slow : %.1f kPoints/s | eff %.1f MPoints/s\n",
               0.001 * opCount / interval, 1e-6 * opCount * numCoils / interval);

        begin_time = high_resolution_clock::now();
        computedBField = torusGroupFlat.computeAllBFieldVectors(fieldPoints, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("Field     B  slow : %.1f kPoints/s | eff %.1f MPoints/s\n",
               0.001 * opCount / interval, 1e-6 * opCount * numCoils / interval);

        begin_time = high_resolution_clock::now();
        computedGradient = torusGroupFlat.computeAllBGradientMatrices(fieldPoints, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("Gradient  G  slow : %.1f kPoints/s | eff %.1f MPoints/s\n",
               0.001 * opCount / interval, 1e-6 * opCount * numCoils / interval);

        begin_time = high_resolution_clock::now();
        computedAPotential = torusGroupThick.computeAllAPotentialVectors(fieldPoints, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("Potential A  fast : %.1f kPoints/s | eff %.1f MPoints/s\n",
               0.001 * opCount / interval, 1e-6 * opCount * numCoils / interval);

        begin_time = high_resolution_clock::now();
        computedBField = torusGroupThick.computeAllBFieldVectors(fieldPoints, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("Field     B  fast : %.1f kPoints/s | eff %.1f MPoints/s\n",
               0.001 * opCount / interval, 1e-6 * opCount * numCoils / interval);

        begin_time = high_resolution_clock::now();
        computedGradient = torusGroupThick.computeAllBGradientMatrices(fieldPoints, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("Gradient  G  fast : %.1f kPoints/s | eff %.1f MPoints/s\n",
               0.001 * opCount / interval, 1e-6 * opCount * numCoils / interval);

        printf("\n");
    }
}

void benchCoilGroupMTvsMTD(int threadCount, int pointCount)
{
    using namespace std::chrono;

    high_resolution_clock::time_point begin_time;
    double interval;

    int coilCount1 = threadCount;
    int coilCount2 = 4 * threadCount;

    begin_time = high_resolution_clock::now();
    compCoilGroupMTD(coilCount1, pointCount, threadCount, false);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MT perf  : %.0f kPoints/s\n", 1e-3 * coilCount1 * pointCount / interval);

    begin_time = high_resolution_clock::now();
    compCoilGroupMTD(coilCount2, pointCount, threadCount, false);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MTD perf : %.0f kPoints/s\n", 1e-3 * coilCount2 * pointCount / interval);

    printf("\n");
}

void benchCoilGroupComputeAllFieldsMTD(int threadCount)
{
    using namespace std::chrono;

    high_resolution_clock::time_point begin_time;
    double interval;

    int numCoils = 8 * threadCount;
    const int pointCount = 1'000;

    CoilGroup coilGroup = CoilGroup();
    Coil referenceCoil = Coil(0.1, 0.1, 0.1, 10000);

    for (int i = 1; i <= numCoils; ++i)
    {
        Coil tempCoil = Coil(0.1, 0.1, 0.1, 10000);
        tempCoil.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.15*i),0.0, 0.0);
        coilGroup.addCoil(tempCoil);
    }

    vec3::Vector3Array referencePoints(pointCount);
    vec3::Vector3Array computedAPotential;
    vec3::Vector3Array computedBField;
    vec3::Matrix3Array computedBGradient;

    for (int i = 0; i < pointCount; ++i)
        referencePoints[i] = vec3::Vector3(0.1, 1.0 * i / pointCount, -0.1);

    computedAPotential = coilGroup.computeAllAPotentialVectors(referencePoints, CPU_ST);
    printf("%.15g\n", computedAPotential[pointCount / 2].x);
    computedAPotential = coilGroup.computeAllAPotentialVectors(referencePoints, CPU_MT);
    printf("%.15g\n", computedAPotential[pointCount / 2].x);
    printf("\n");

    computedBField = coilGroup.computeAllBFieldVectors(referencePoints, CPU_ST);
    printf("%.15g\n", computedBField[pointCount / 2].x);
    computedBField = coilGroup.computeAllBFieldVectors(referencePoints, CPU_MT);
    printf("%.15g\n", computedBField[pointCount / 2].x);
    printf("\n");

    computedBGradient = coilGroup.computeAllBGradientMatrices(referencePoints, CPU_ST);
    printf("%.15g\n", computedBGradient[pointCount / 2].xy);
    computedBGradient = coilGroup.computeAllBGradientMatrices(referencePoints, CPU_MT);
    printf("%.15g\n", computedBGradient[pointCount / 2].xy);
    printf("\n");
}

void benchCoilGroupMInductanceAndForceMTD(int threadCount)
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
        tempCoil.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.15*i),0.0, 0.0);
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
           coilGroup.computeAmpereForce(referenceCoil, PrecisionFactor(5.0), CPU_ST).first.z);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("ST  perf  : %.0f coils/s\n", numCoils / interval);


    begin_time = high_resolution_clock::now();
    printf("%.15g\n",
           coilGroup.computeAmpereForce(referenceCoil,PrecisionFactor(5.0), CPU_MT).first.z);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MTD perf  : %.0f coils/s\n", numCoils / interval);
}

void benchMInductanceAndForceComputeAll(PrecisionFactor precisionFactor, int threadCount)
{
    using namespace std::chrono;

    Coil prim = Coil(0.1, 0.1, 0.1, 10000);
    Coil sec = Coil(0.1, 0.1, 0.1, 10000);

    int numOps = threadCount * 32;

    vec3::Vector3Array primPositions(numOps);
    vec3::Vector3Array secPositions(numOps);
    std::vector<double> primYAxisAngle(numOps);
    std::vector<double> primZAxisAngle(numOps);
    std::vector<double> secYAxisAngle(numOps);
    std::vector<double> secZAxisAngle(numOps);

    std::vector<double> mutualInductanceMT(numOps);
    std::vector<double> mutualInductanceAll(numOps);
    std::vector<std::pair<vec3::Vector3, vec3::Vector3>> forceAndTorqueMT(numOps);
    std::vector<std::pair<vec3::Vector3, vec3::Vector3>> forceAndTorqueAll(numOps);

    for (int i = 0; i < numOps; ++i)
    {
        primPositions[i] = vec3::Vector3(0.0, 0.0, 0.0);
        secPositions[i] = vec3::Vector3(0.0, 0.1, 0.2 + double(i) * 0.005);
        primYAxisAngle[i] = 0.0;
        primZAxisAngle[i] = 0.0;
        secYAxisAngle[i] = 0.5;
        secZAxisAngle[i] = 0.5;
    }

    high_resolution_clock::time_point begin_time;
    double interval;

    printf("Mutual inductance:\n");

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < numOps; ++i)
    {
        prim.setPositionAndOrientation(primPositions[i], primYAxisAngle[i], primZAxisAngle[i]);
        sec.setPositionAndOrientation(secPositions[i], secYAxisAngle[i], secZAxisAngle[i]);
        mutualInductanceMT[i] = Coil::computeMutualInductance(prim, sec, precisionFactor, CPU_MT);
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MT perf  : %.1f Ops/s\n", numOps / interval);

    begin_time = high_resolution_clock::now();
    mutualInductanceAll = Coil::computeAllMutualInductanceArrangements(prim, sec, primPositions,secPositions,
                                                                       primYAxisAngle, primZAxisAngle,
                                                                       secYAxisAngle, secZAxisAngle,
                                                                       precisionFactor, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MTD perf : %.1f Ops/s\n", numOps / interval);

    printf("\nAmpere force:\n");

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < numOps; ++i)
    {
        prim.setPositionAndOrientation(primPositions[i], primYAxisAngle[i], primZAxisAngle[i]);
        sec.setPositionAndOrientation(secPositions[i], secYAxisAngle[i], secZAxisAngle[i]);
        forceAndTorqueMT[i] = Coil::computeAmpereForce(prim, sec, precisionFactor, CPU_MT);
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MT perf  : %.1f Ops/s\n", numOps / interval);

    begin_time = high_resolution_clock::now();
    forceAndTorqueAll = Coil::computeAllAmpereForceArrangements(prim, sec, primPositions,secPositions,
                                                                primYAxisAngle, primZAxisAngle,
                                                                secYAxisAngle, secZAxisAngle,
                                                                precisionFactor, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MTD perf : %.1f Ops/s\n", numOps / interval);

    printf("\n");
}
