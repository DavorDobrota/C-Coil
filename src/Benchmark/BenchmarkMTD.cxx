#include "Benchmark.h"
#include "Coil.h"
#include "Coil/EnumsAndConstants/ComputeMethod.h"
#include "Tensor.h"
#include "Compare.h"
#include "CoilGroup.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>
#include <vector>
#include <chrono>


void benchCoilGroupMTvsMTD(int threadCount, int pointCount)
{
    using namespace std::chrono;

    high_resolution_clock::time_point begin_time;
    double interval;

    int coilCount1 = threadCount;
    int coilCount2 = 3 * threadCount;

    printf("Quick performance benchmark for %d coils and %d points\n\n", threadCount, pointCount);

    begin_time = high_resolution_clock::now();
    compCoilGroupMTD(coilCount1, pointCount, 1, false);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("ST  perf : %.0f kPoints/s\n", 1e-3 * coilCount1 * pointCount / interval);

    begin_time = high_resolution_clock::now();
    compCoilGroupMTD(coilCount1, pointCount, threadCount, false);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MT  perf : %.0f kPoints/s\n", 1e-3 * coilCount1 * pointCount / interval);

    begin_time = high_resolution_clock::now();
    compCoilGroupMTD(coilCount2, pointCount, threadCount, false);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MTD perf : %.0f kPoints/s\n", 1e-3 * coilCount2 * pointCount / interval);

    printf("\n");
}

void benchCoilGroupComputeAllFields(PrecisionFactor precisionFactor, int coilCount, int opCount, int threadCount)
{
    using namespace std::chrono;

    high_resolution_clock::time_point begin_time;
    double interval;

    printf("Benchmarking MT and GPU field compute performance for %d coils in %d points and precision factor %.1f\n\n",
           coilCount, opCount, precisionFactor.relativePrecision);

    double torusRadius = 1.0;

    CoilGroup torusGroupThick = CoilGroup();
    CoilGroup torusGroupFlat = CoilGroup();

    vec3::Vector3Array fieldPoints(opCount);
    vec3::Vector3Array computedAPotential;
    vec3::Vector3Array computedBField;
    vec3::Matrix3Array computedGradient;

    for (int i = 0; i < opCount; ++i)
        fieldPoints[i] = vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / opCount);

    for (int i = 0; i < coilCount; ++i)
    {
        Coil tempCoil = Coil(torusRadius / 10.0, torusRadius / 100.0, torusRadius / 100.0, 10000, 10);
        tempCoil.setPositionAndOrientation(
                vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / coilCount),
                M_PI_2, 2*M_PI * i / coilCount + M_PI_2);
        torusGroupThick.addCoil(tempCoil);
    }
    torusGroupThick.setThreadCount(threadCount);
    torusGroupThick.setDefaultPrecisionFactor(precisionFactor);

    for (int i = 0; i < coilCount; ++i)
    {
        Coil tempCoil = Coil(torusRadius / 10.0, torusRadius / 100.0, 0.0, 100, 10);
        tempCoil.setPositionAndOrientation(
                vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / coilCount),
                M_PI_2, 2*M_PI * i / coilCount + M_PI_2);
        torusGroupFlat.addCoil(tempCoil);
    }
    torusGroupFlat.setThreadCount(threadCount);
    torusGroupFlat.setDefaultPrecisionFactor(precisionFactor);

    // MT slow tests
    begin_time = high_resolution_clock::now();
    computedAPotential = torusGroupFlat.computeAllAPotentialVectors(fieldPoints, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_MT slow : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    begin_time = high_resolution_clock::now();
    computedBField = torusGroupFlat.computeAllBFieldVectors(fieldPoints, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Field     B CPU_MT slow : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    begin_time = high_resolution_clock::now();
    computedGradient = torusGroupFlat.computeAllBGradientMatrices(fieldPoints, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G CPU_MT slow : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    printf("\n");

    // MT fast tests
    begin_time = high_resolution_clock::now();
    computedAPotential = torusGroupThick.computeAllAPotentialVectors(fieldPoints, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A CPU_MT fast : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    begin_time = high_resolution_clock::now();
    computedBField = torusGroupThick.computeAllBFieldVectors(fieldPoints, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Field     B CPU_MT fast : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    begin_time = high_resolution_clock::now();
    computedGradient = torusGroupThick.computeAllBGradientMatrices(fieldPoints, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G CPU_MT fast : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    printf("\n");

    // GPU slow tests
    computedAPotential = torusGroupThick.computeAllAPotentialVectors(fieldPoints, GPU); // warmup

    begin_time = high_resolution_clock::now();
    computedAPotential = torusGroupFlat.computeAllAPotentialVectors(fieldPoints, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A GPU    slow : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    begin_time = high_resolution_clock::now();
    computedBField = torusGroupFlat.computeAllBFieldVectors(fieldPoints, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Field     B GPU    slow : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    begin_time = high_resolution_clock::now();
    computedGradient = torusGroupFlat.computeAllBGradientMatrices(fieldPoints, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G GPU    slow : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    printf("\n");

    // GPU fast tests

    begin_time = high_resolution_clock::now();
    computedAPotential = torusGroupThick.computeAllAPotentialVectors(fieldPoints, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Potential A GPU    fast : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    begin_time = high_resolution_clock::now();
    computedBField = torusGroupThick.computeAllBFieldVectors(fieldPoints, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Field     B GPU    fast : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    begin_time = high_resolution_clock::now();
    computedGradient = torusGroupThick.computeAllBGradientMatrices(fieldPoints, GPU);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("Gradient  G GPU    fast : %.1f kPoints/s | eff %.1f kPoints/s\n",
           0.001 * opCount / interval, 0.001 * opCount * coilCount / interval);

    printf("\n");
}

void benchCoilGroupComputeAllFieldsGPU(int coilCount, int opCount)
{
    using namespace std::chrono;

    high_resolution_clock::time_point begin_time;
    double interval;

    printf("Benchmarking GPU field compute performance for %d coils in %d points\n\n", coilCount, opCount);

    double torusRadius = 1.0;

    CoilGroup torusGroupThick = CoilGroup();
    CoilGroup torusGroupFlat = CoilGroup();

    vec3::Vector3Array fieldPoints(opCount);
    vec3::Vector3Array computedAPotential;
    vec3::Vector3Array computedBField;
    vec3::Matrix3Array computedGradient;

    for (int i = 0; i < opCount; ++i)
        fieldPoints[i] = vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / opCount);

    for (int i = 0; i < coilCount; ++i)
    {
        Coil tempCoil = Coil(torusRadius / 10.0, torusRadius / 100.0, torusRadius / 100.0, 10000, 10);
        tempCoil.setPositionAndOrientation(
                vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / coilCount),
                M_PI_2, 2*M_PI * i / coilCount + M_PI_2);
        torusGroupThick.addCoil(tempCoil);
    }

    for (int i = 0; i < coilCount; ++i)
    {
        Coil tempCoil = Coil(torusRadius / 10.0, torusRadius / 100.0, 0.0, 100, 10);
        tempCoil.setPositionAndOrientation(
                vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / coilCount),
                M_PI_2, 2*M_PI * i / coilCount + M_PI_2);
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
               0.001 * opCount / interval, 1e-6 * opCount * coilCount / interval);

        begin_time = high_resolution_clock::now();
        computedBField = torusGroupFlat.computeAllBFieldVectors(fieldPoints, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("Field     B  slow : %.1f kPoints/s | eff %.1f MPoints/s\n",
               0.001 * opCount / interval, 1e-6 * opCount * coilCount / interval);

        begin_time = high_resolution_clock::now();
        computedGradient = torusGroupFlat.computeAllBGradientMatrices(fieldPoints, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("Gradient  G  slow : %.1f kPoints/s | eff %.1f MPoints/s\n",
               0.001 * opCount / interval, 1e-6 * opCount * coilCount / interval);

        begin_time = high_resolution_clock::now();
        computedAPotential = torusGroupThick.computeAllAPotentialVectors(fieldPoints, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("Potential A  fast : %.1f kPoints/s | eff %.1f MPoints/s\n",
               0.001 * opCount / interval, 1e-6 * opCount * coilCount / interval);

        begin_time = high_resolution_clock::now();
        computedBField = torusGroupThick.computeAllBFieldVectors(fieldPoints, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("Field     B  fast : %.1f kPoints/s | eff %.1f MPoints/s\n",
               0.001 * opCount / interval, 1e-6 * opCount * coilCount / interval);

        begin_time = high_resolution_clock::now();
        computedGradient = torusGroupThick.computeAllBGradientMatrices(fieldPoints, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("Gradient  G  fast : %.1f kPoints/s | eff %.1f MPoints/s\n",
               0.001 * opCount / interval, 1e-6 * opCount * coilCount / interval);

        printf("\n");
    }
}


void benchCoilGroupMInductanceAndForce(int opCount, int threadCount)
{
    using namespace std::chrono;

    high_resolution_clock::time_point begin_time;
    double interval;

    int coilCountMTD = 3 * threadCount;
    int coilCount = threadCount;

    CoilGroup coilGroupMTD = CoilGroup();
    CoilGroup coilGroup = CoilGroup();
    Coil referenceCoil = Coil(0.1, 0.1, 0.1, 10000);

    for (int i = 1; i <= coilCountMTD; ++i)
    {
        Coil tempCoil = Coil(0.1, 0.1, 0.1, 10000);
        tempCoil.setPositionAndOrientation(vec3::Vector3(1e-8, 0.0, 0.15*i),0.0, 0.0);
        coilGroupMTD.addCoil(tempCoil);
    }

    for (int i = 0; i <= coilCount; ++i)
    {
        Coil tempCoil = Coil(0.1, 0.1, 0.1, 10000);
        tempCoil.setPositionAndOrientation(vec3::Vector3(1e-8, 0.0, 0.15*i),0.0, 0.0);
        coilGroup.addCoil(tempCoil);
    }

    printf("Benchmarking mutual inductance and force performance for various compute methods\n\n");

    for (int i = 1; i <= 8; ++i)
    {
        auto precisionFactor = PrecisionFactor(double(i));
        double tempInductance;
        double perf;

        printf("Mutual inductance performance for precision factor %.1f\n", double(i));

        begin_time = high_resolution_clock::now();
        for (int j = 0; j < opCount; ++j)
            tempInductance = coilGroup.computeMutualInductance(referenceCoil, precisionFactor, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = opCount * coilCount / interval;
        printf("ST  : %6.3f ms/coil | %.0f coils/s\n", 1e3/ perf, perf);

        begin_time = high_resolution_clock::now();
        for (int j = 0; j < opCount; ++j)
            tempInductance = coilGroup.computeMutualInductance(referenceCoil, precisionFactor, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = opCount * coilCount / interval;
        printf("MT  : %6.3f ms/coil | %.0f coils/s\n", 1e3 / perf, perf);

        begin_time = high_resolution_clock::now();
        for (int j = 0; j < opCount; ++j)
            tempInductance = coilGroupMTD.computeMutualInductance(referenceCoil, precisionFactor, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = opCount * coilCountMTD / interval;
        printf("MTD : %6.3f ms/coil | %.0f coils/s\n", 1e3 / perf, perf);

        tempInductance = coilGroupMTD.computeMutualInductance(referenceCoil, precisionFactor, GPU); // GPU warmup
        begin_time = high_resolution_clock::now();
        for (int j = 0; j < opCount; ++j)
            tempInductance = coilGroupMTD.computeMutualInductance(referenceCoil, precisionFactor, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = opCount * coilCountMTD / interval;
        printf("GPU : %6.3f ms/coil | %.0f coils/s\n", 1e3 / perf, perf);

        printf("\n");
    }

    printf("\n");

    for (int i = 1; i <= 8; ++i)
    {
        auto precisionFactor = PrecisionFactor(double(i));
        std::pair<vec3::Vector3, vec3::Vector3> tempForce;
        double perf;

        printf("Ampere force performance for precision factor %.1f\n", double(i));

        begin_time = high_resolution_clock::now();
        for (int j = 0; j < opCount; ++j)
            tempForce = coilGroup.computeAmpereForce(referenceCoil, precisionFactor, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = opCount * coilCount / interval;
        printf("ST  : %6.3f ms/coil | %.0f coils/s\n", 1e3/ perf, perf);

        begin_time = high_resolution_clock::now();
        for (int j = 0; j < opCount; ++j)
            tempForce = coilGroup.computeAmpereForce(referenceCoil, precisionFactor, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = opCount * coilCount / interval;
        printf("MT  : %6.3f ms/coil | %.0f coils/s\n", 1e3 / perf, perf);

        begin_time = high_resolution_clock::now();
        for (int j = 0; j < opCount; ++j)
            tempForce = coilGroupMTD.computeAmpereForce(referenceCoil, precisionFactor, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = opCount * coilCountMTD / interval;
        printf("MTD : %6.3f ms/coil | %.0f coils/s\n", 1e3 / perf, perf);

        tempForce = coilGroupMTD.computeAmpereForce(referenceCoil, precisionFactor, GPU); // GPU warmup
        begin_time = high_resolution_clock::now();
        for (int j = 0; j < opCount; ++j)
            tempForce = coilGroupMTD.computeAmpereForce(referenceCoil, precisionFactor, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = opCount * coilCountMTD / interval;
        printf("GPU : %6.3f ms/coil | %.0f coils/s\n", 1e3 / perf, perf);

        printf("\n");
    }
}

void benchCoilGroupMInductanceAndForceAll(int coilCount, int opCount, int threadCount)
{
    using namespace std::chrono;

    double torusRadius = 1.0;

    CoilGroup torusGroup = CoilGroup();

    vec3::Vector3Array secPositions(opCount);
    std::vector<double> secYAxisAngle(opCount);
    std::vector<double> secZAxisAngle(opCount);

    for (int i = 0; i < coilCount; ++i)
    {
        Coil tempCoil = Coil(0.1, 0.1, 0.1, 10000, 10);
        tempCoil.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.2 * double(i)));

        torusGroup.addCoil(tempCoil);
    }

    torusGroup.setDefaultPrecisionFactor(PrecisionFactor(5.0));

    Coil secondary = Coil(0.1, 0.1, 0.1, 10000, 5);

    for (int i = 0; i < opCount; ++i)
    {
        secPositions[i] = vec3::Vector3(0.5, 0.0, 0.1 * double(i));
        secYAxisAngle[i] = 0.5;
        secZAxisAngle[i] = 0.6;
    }

    high_resolution_clock::time_point begin_time;
    double interval;
    double perf;

    std::vector<double> tempInductances = torusGroup.computeAllMutualInductanceArrangements(
            secondary, secPositions, secYAxisAngle, secZAxisAngle,
            PrecisionFactor(1.0), GPU
    ); // GPU warmup

    for (int i = 1; i <= 8; ++i)
    {
        printf("Mutual inductance performance for precision factor %.1f\n", double(i));

        double precision = ((double(i) - 1.0) / 2.0) + 1.0;
        torusGroup.setDefaultPrecisionFactor(PrecisionFactor(precision));

        begin_time = high_resolution_clock::now();
        tempInductances = torusGroup.computeAllMutualInductanceArrangements(
            secondary, secPositions, secYAxisAngle, secZAxisAngle,
            PrecisionFactor(precision), GPU
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = opCount * coilCount / interval;
        printf("GPU : %6.2f microseconds/Op | %.0f Ops/s\n", 1e6 / perf, perf);

        printf("\n");
    }

}

void benchMInductanceAndForceComputeAllMTvsMTD(PrecisionFactor precisionFactor, int threadCount)
{
    using namespace std::chrono;

    Coil prim = Coil(0.1, 0.1, 0.1, 10000);
    Coil sec = Coil(0.1, 0.1, 0.1, 10000);

    int opCount = threadCount * 32;

    vec3::Vector3Array primPositions(opCount);
    vec3::Vector3Array secPositions(opCount);
    std::vector<double> primYAxisAngle(opCount);
    std::vector<double> primZAxisAngle(opCount);
    std::vector<double> secYAxisAngle(opCount);
    std::vector<double> secZAxisAngle(opCount);

    std::vector<double> mutualInductanceMT(opCount);
    std::vector<double> mutualInductanceAll(opCount);
    std::vector<std::pair<vec3::Vector3, vec3::Vector3>> forceAndTorqueMT(opCount);
    std::vector<std::pair<vec3::Vector3, vec3::Vector3>> forceAndTorqueAll(opCount);

    for (int i = 0; i < opCount; ++i)
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
    for (int i = 0; i < opCount; ++i)
    {
        prim.setPositionAndOrientation(primPositions[i], primYAxisAngle[i], primZAxisAngle[i]);
        sec.setPositionAndOrientation(secPositions[i], secYAxisAngle[i], secZAxisAngle[i]);
        mutualInductanceMT[i] = Coil::computeMutualInductance(prim, sec, precisionFactor, CPU_MT);
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MT perf  : %.1f Ops/s\n", opCount / interval);

    begin_time = high_resolution_clock::now();
    mutualInductanceAll = Coil::computeAllMutualInductanceArrangements(prim, sec, primPositions,secPositions,
                                                                       primYAxisAngle, primZAxisAngle,
                                                                       secYAxisAngle, secZAxisAngle,
                                                                       precisionFactor, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MTD perf : %.1f Ops/s\n", opCount / interval);

    printf("\nAmpere force:\n");

    begin_time = high_resolution_clock::now();
    for (int i = 0; i < opCount; ++i)
    {
        prim.setPositionAndOrientation(primPositions[i], primYAxisAngle[i], primZAxisAngle[i]);
        sec.setPositionAndOrientation(secPositions[i], secYAxisAngle[i], secZAxisAngle[i]);
        forceAndTorqueMT[i] = Coil::computeAmpereForce(prim, sec, precisionFactor, CPU_MT);
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MT perf  : %.1f Ops/s\n", opCount / interval);

    begin_time = high_resolution_clock::now();
    forceAndTorqueAll = Coil::computeAllAmpereForceArrangements(prim, sec, primPositions,secPositions,
                                                                primYAxisAngle, primZAxisAngle,
                                                                secYAxisAngle, secZAxisAngle,
                                                                precisionFactor, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MTD perf : %.1f Ops/s\n", opCount / interval);

    printf("\n");
}

void benchMInductanceAndForceComputeAllGPU(int configCount)
{
    using namespace std::chrono;

    Coil prim = Coil(0.1, 0.1, 0.1, 10000);
    Coil sec = Coil(0.1, 0.1, 0.1, 10000);

    vec3::Vector3Array primPositions(configCount);
    vec3::Vector3Array secPositions(configCount);
    std::vector<double> primYAxisAngle(configCount);
    std::vector<double> primZAxisAngle(configCount);
    std::vector<double> secYAxisAngle(configCount);
    std::vector<double> secZAxisAngle(configCount);

    for (int i = 0; i < configCount; ++i)
    {
        primPositions[i] = vec3::Vector3(0.0, 0.0, 0.0);
        secPositions[i] = vec3::Vector3(0.0, 0.1, 0.2 + double(i) * 0.005);
        primYAxisAngle[i] = 0.5;
        primZAxisAngle[i] = 0.5;
        secYAxisAngle[i] = 0.6;
        secZAxisAngle[i] = 0.6;
    }

    std::vector<double> warmupOutput = Coil::computeAllMutualInductanceArrangements(
        prim, sec, primPositions,secPositions,
        primYAxisAngle, primZAxisAngle, secYAxisAngle, secZAxisAngle,
        PrecisionFactor(2.0), GPU
    ); // GPU warmup

    high_resolution_clock::time_point begin_time;
    double interval;

    printf("Benchmarking mutual inductance, force and torque on GPU for %d configurations:\n\n", configCount);

    for (int i = 1; i <= 9; ++i)
    {
        std::vector<double> mutualInductance(configCount);
        std::vector<std::pair<vec3::Vector3, vec3::Vector3>> forceTorque(configCount);
        auto precisionFactor = PrecisionFactor(double(i));
        double performance;

        printf("Performance for precision factor %.1f\n", double(i));

        begin_time = high_resolution_clock::now();
        mutualInductance = Coil::computeAllMutualInductanceArrangements(prim, sec, primPositions,secPositions,
                                                                        primYAxisAngle, primZAxisAngle,
                                                                        secYAxisAngle, secZAxisAngle,
                                                                        precisionFactor, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        performance = configCount / interval;
        printf("Mutual inductance : %6.2f microseconds/Op | %.0f Ops/s\n", 1e6 / performance, performance);

        begin_time = high_resolution_clock::now();
        forceTorque = Coil::computeAllAmpereForceArrangements(prim, sec, primPositions,secPositions,
                                                              primYAxisAngle, primZAxisAngle,
                                                              secYAxisAngle, secZAxisAngle,
                                                              precisionFactor, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        performance = configCount / interval;
        printf("Force and torque  : %6.2f microseconds/Op | %.0f Ops/s\n", 1e6 / performance, performance);


        printf("\n");
    }
}

void benchMInductanceAndForceComputeAll(int configCount, int threadCount)
{
    using namespace std::chrono;

    Coil prim = Coil(0.1, 0.1, 0.1, 10000);
    Coil sec = Coil(0.1, 0.1, 0.1, 10000);

    vec3::Vector3Array primPositions(configCount);
    vec3::Vector3Array secPositions(configCount);
    std::vector<double> primYAxisAngle(configCount);
    std::vector<double> primZAxisAngle(configCount);
    std::vector<double> secYAxisAngle(configCount);
    std::vector<double> secZAxisAngle(configCount);

    for (int i = 0; i < configCount; ++i)
    {
        primPositions[i] = vec3::Vector3(0.0, 0.0, 0.0);
        secPositions[i] = vec3::Vector3(0.0, 0.1, 0.2 + double(i) * 0.005);
        primYAxisAngle[i] = 0.5;
        primZAxisAngle[i] = 0.5;
        secYAxisAngle[i] = 0.6;
        secZAxisAngle[i] = 0.6;
    }

    high_resolution_clock::time_point begin_time;
    double interval;

    printf("Benchmarking mutual inductance for %d configurations:\n\n", configCount);

    for (int i = 1; i <= 9; ++i)
    {
        std::vector<double> mutualInductance(configCount);
        auto precisionFactor = PrecisionFactor(double(i));
        double performance;

        printf("Performance for precision factor %.1f\n", double(i));

        begin_time = high_resolution_clock::now();
        mutualInductance = Coil::computeAllMutualInductanceArrangements(prim, sec, primPositions,secPositions,
                                                                           primYAxisAngle, primZAxisAngle,
                                                                           secYAxisAngle, secZAxisAngle,
                                                                           precisionFactor, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        performance = configCount / interval;
        printf("MT  perf : %6.3f milliseconds/Op | %.1f Ops/s\n", 1e3 / performance, performance);

        mutualInductance = Coil::computeAllMutualInductanceArrangements(prim, sec, primPositions,secPositions,
                                                                        primYAxisAngle, primZAxisAngle,
                                                                        secYAxisAngle, secZAxisAngle,
                                                                        precisionFactor, GPU); // GPU warmup

        begin_time = high_resolution_clock::now();
        mutualInductance = Coil::computeAllMutualInductanceArrangements(prim, sec, primPositions,secPositions,
                                                                           primYAxisAngle, primZAxisAngle,
                                                                           secYAxisAngle, secZAxisAngle,
                                                                           precisionFactor, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        performance = configCount / interval;
        printf("GPU perf : %6.2f microseconds/Op | %.0f Ops/s\n", 1e6 / performance, performance);

        printf("\n");
    }

    printf("\n");
    printf("Benchmarking force and torque for %d configurations:\n\n", configCount);

    for (int i = 1; i <= 9; ++i)
    {
        std::vector<std::pair<vec3::Vector3, vec3::Vector3>> forceAndTorque(configCount);
        auto precisionFactor = PrecisionFactor(double(i));
        double performance;

        printf("Performance for precision factor %.1f\n", double(i));

        begin_time = high_resolution_clock::now();
        forceAndTorque = Coil::computeAllAmpereForceArrangements(prim, sec, primPositions,secPositions,
                                                                    primYAxisAngle, primZAxisAngle,
                                                                    secYAxisAngle, secZAxisAngle,
                                                                    precisionFactor, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        performance = configCount / interval;
        printf("MT  perf : %6.3f milliseconds/Op | %.1f Ops/s\n", 1e3 / performance, performance);

        begin_time = high_resolution_clock::now();
        forceAndTorque = Coil::computeAllAmpereForceArrangements(prim, sec, primPositions,secPositions,
                                                                    primYAxisAngle, primZAxisAngle,
                                                                    secYAxisAngle, secZAxisAngle,
                                                                    precisionFactor, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        performance = configCount / interval;
        printf("GPU perf : %6.2f microseconds/Op | %.0f Ops/s\n", 1e6 / performance, performance);

        printf("\n");
    }
}
