#include "Benchmark.h"

#include <cstdio>
#include <vector>
#include <chrono>


void Benchmark::coilGroupMInductanceAndForce(int repeatCount, int threadCount)
{
    using namespace std::chrono;

    high_resolution_clock::time_point begin_time;
    double interval;

    int coilCountMTD = 2 * threadCount;
    int coilCount = 2 * threadCount - 1;

    CoilGroup coilGroupMTD = CoilGroup();
    CoilGroup coilGroup = CoilGroup();
    Coil referenceCoil = Coil(0.1, 0.1, 0.1, 10000);

    for (int i = 1; i <= coilCountMTD; ++i)
    {
        coilGroupMTD.addCoil(
            0.1, 0.1, 0.1, 10000, 10,
            PrecisionFactor(), 8,
            vec3::Vector3(1e-8, 0.0, 0.15*double(i)),0.0, 0.0
        );
    }

    for (int i = 0; i <= coilCount; ++i)
    {
        coilGroup.addCoil(
            0.1, 0.1, 0.1, 10000, 10,
            PrecisionFactor(), 8,
            vec3::Vector3(1e-8, 0.0, 0.15*double(i)),0.0, 0.0
        );
    }

    printf("Benchmarking mutual inductance and force performance for various compute methods\n\n");

    for (int i = 1; i <= 8; ++i)
    {
        auto precisionFactor = PrecisionFactor(double(i));
        double tempInductance;
        double perf;

        printf("Mutual inductance performance for precision factor %.1f\n", double(i));

        begin_time = high_resolution_clock::now();
        for (int j = 0; j < repeatCount; ++j)
            tempInductance = coilGroup.computeMutualInductance(referenceCoil, precisionFactor, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = repeatCount * coilCount / interval;
        printf("ST  : %6.3f ms/coil | %.0f coils/s\n", 1e3/ perf, perf);

        begin_time = high_resolution_clock::now();
        for (int j = 0; j < repeatCount; ++j)
            tempInductance = coilGroup.computeMutualInductance(referenceCoil, precisionFactor, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = repeatCount * coilCount / interval;
        printf("MT  : %6.3f ms/coil | %.0f coils/s\n", 1e3 / perf, perf);

        begin_time = high_resolution_clock::now();
        for (int j = 0; j < repeatCount; ++j)
            tempInductance = coilGroupMTD.computeMutualInductance(referenceCoil, precisionFactor, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = repeatCount * coilCountMTD / interval;
        printf("MTD : %6.3f ms/coil | %.0f coils/s\n", 1e3 / perf, perf);

        tempInductance = coilGroupMTD.computeMutualInductance(referenceCoil, precisionFactor, GPU); // GPU warmup
        begin_time = high_resolution_clock::now();
        for (int j = 0; j < repeatCount; ++j)
            tempInductance = coilGroupMTD.computeMutualInductance(referenceCoil, precisionFactor, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = repeatCount * coilCountMTD / interval;
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
        for (int j = 0; j < repeatCount; ++j)
            tempForce = coilGroup.computeForceTorque(referenceCoil, precisionFactor, CPU_ST);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = repeatCount * coilCount / interval;
        printf("ST  : %6.3f ms/coil | %.0f coils/s\n", 1e3/ perf, perf);

        begin_time = high_resolution_clock::now();
        for (int j = 0; j < repeatCount; ++j)
            tempForce = coilGroup.computeForceTorque(referenceCoil, precisionFactor, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = repeatCount * coilCount / interval;
        printf("MT  : %6.3f ms/coil | %.0f coils/s\n", 1e3 / perf, perf);

        begin_time = high_resolution_clock::now();
        for (int j = 0; j < repeatCount; ++j)
            tempForce = coilGroupMTD.computeForceTorque(referenceCoil, precisionFactor, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = repeatCount * coilCountMTD / interval;
        printf("MTD : %6.3f ms/coil | %.0f coils/s\n", 1e3 / perf, perf);

        tempForce = coilGroupMTD.computeForceTorque(referenceCoil, precisionFactor, GPU); // GPU warmup
        begin_time = high_resolution_clock::now();
        for (int j = 0; j < repeatCount; ++j)
            tempForce = coilGroupMTD.computeForceTorque(referenceCoil, precisionFactor, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = repeatCount * coilCountMTD / interval;
        printf("GPU : %6.3f ms/coil | %.0f coils/s\n", 1e3 / perf, perf);

        printf("\n");
    }
}

void Benchmark::coilGroupMInductanceAndForceAll(int coilCount, int opCount, int threadCount)
{
    using namespace std::chrono;

    double torusRadius = 1.0;

    CoilGroup coilGroup = CoilGroup();

    vec3::Vector3Array secPositions(opCount);
    std::vector<double> secYAxisAngle(opCount);
    std::vector<double> secZAxisAngle(opCount);

    for (int i = 0; i < coilCount; ++i)
    {
        Coil tempCoil = Coil(0.1, 0.1, 0.1, 10000, 10);
        tempCoil.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.2 * double(i)));

        coilGroup.addCoil(
                0.1, 0.1, 0.1, 10000, 10,
                PrecisionFactor(), 8,
                vec3::Vector3(0.0, 0.0, 0.2 * double(i))
        );
    }

    coilGroup.setDefaultPrecisionFactor(PrecisionFactor(5.0));
    coilGroup.setThreadCount(threadCount);

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

    std::vector<double> tempInductances;
    std::vector<std::pair<vec3::Vector3, vec3::Vector3>> tempForceTorque;

    printf("Benchmarking CoilGroup Mutual Inductance and Force performance for %d coils and %d configurations\n\n",
           coilCount, opCount
    );

    tempInductances = coilGroup.computeAllMutualInductanceArrangements(
            secondary, secPositions, secYAxisAngle, secZAxisAngle,
            PrecisionFactor(1.0), GPU
    ); // GPU warmup

    printf("Mutual Inductance:\n\n");

    for (int i = 1; i <= 8; ++i)
    {
        printf("Precision factor %.1f\n", double(i));

        begin_time = high_resolution_clock::now();
        tempInductances = coilGroup.computeAllMutualInductanceArrangements(
                secondary, secPositions, secYAxisAngle, secZAxisAngle,
                PrecisionFactor(PrecisionFactor(i)), CPU_MT
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = opCount * coilCount / interval;
        printf("CPU_MT : %6.3f milliseconds/Op | %.1f Ops/s\n", 1e3 / perf, perf);

        begin_time = high_resolution_clock::now();
        tempInductances = coilGroup.computeAllMutualInductanceArrangements(
                secondary, secPositions, secYAxisAngle, secZAxisAngle,
                PrecisionFactor(PrecisionFactor(i)), GPU
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = opCount * coilCount / interval;
        printf("GPU    : %6.2f microseconds/Op | %.0f Ops/s\n", 1e6 / perf, perf);

        printf("\n");
    }

    printf("Force and torque:\n\n");

    for (int i = 1; i <= 8; ++i)
    {
        printf("Precision factor %.1f\n", double(i));

        begin_time = high_resolution_clock::now();
        tempForceTorque = coilGroup.computeAllForceTorqueArrangements(
                secondary, secPositions, secYAxisAngle, secZAxisAngle,
                PrecisionFactor(i), CPU_MT
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = opCount * coilCount / interval;
        printf("CPU_MT : %6.3f milliseconds/Op | %.1f Ops/s\n", 1e3 / perf, perf);

        begin_time = high_resolution_clock::now();
        tempForceTorque = coilGroup.computeAllForceTorqueArrangements(
                secondary, secPositions, secYAxisAngle, secZAxisAngle,
                PrecisionFactor(i), GPU
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = opCount * coilCount / interval;
        printf("GPU    : %6.2f microseconds/Op | %.0f Ops/s\n", 1e6 / perf, perf);

        printf("\n");
    }
}

void Benchmark::coilGroupMInductanceAndForceAllGPU(int coilCount, int opCount)
{
    using namespace std::chrono;

    double torusRadius = 1.0;

    CoilGroup coilGroup = CoilGroup();

    vec3::Vector3Array secPositions(opCount);
    std::vector<double> secYAxisAngle(opCount);
    std::vector<double> secZAxisAngle(opCount);

    for (int i = 0; i < coilCount; ++i)
    {
        Coil tempCoil = Coil(0.1, 0.1, 0.1, 10000, 10);
        tempCoil.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.2 * double(i)));

        coilGroup.addCoil(
                0.1, 0.1, 0.1, 10000, 10,
                PrecisionFactor(), 8,
                vec3::Vector3(0.0, 0.0, 0.2 * double(i))
        );
    }

    coilGroup.setDefaultPrecisionFactor(PrecisionFactor(5.0));

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

    std::vector<double> tempInductances;
    std::vector<std::pair<vec3::Vector3, vec3::Vector3>> tempForceTorque;

    printf("Benchmarking CoilGroup Mutual Inductance and Force GPU performance for %d coils and %d configurations\n",
           coilCount, opCount
    );

    tempInductances = coilGroup.computeAllMutualInductanceArrangements(
            secondary, secPositions, secYAxisAngle, secZAxisAngle,
            PrecisionFactor(1.0), GPU
    ); // GPU warmup

    for (int i = 1; i <= 8; ++i)
    {
        printf("Precision factor %.1f\n", double(i));

        begin_time = high_resolution_clock::now();
        tempInductances = coilGroup.computeAllMutualInductanceArrangements(
                secondary, secPositions, secYAxisAngle, secZAxisAngle,
                PrecisionFactor(i), GPU
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = opCount * coilCount / interval;
        printf("Mutual inductance : %6.2f microseconds/Op | %.0f Ops/s\n", 1e6 / perf, perf);

        begin_time = high_resolution_clock::now();
        tempForceTorque = coilGroup.computeAllForceTorqueArrangements(
                secondary, secPositions, secYAxisAngle, secZAxisAngle,
                PrecisionFactor(i), GPU
        );
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        perf = opCount * coilCount / interval;
        printf("Force and torque  : %6.2f microseconds/Op | %.0f Ops/s\n", 1e6 / perf, perf);

        printf("\n");
    }
}

