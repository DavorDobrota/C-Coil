#include "Benchmark.h"

#include <chrono>


void Benchmark::coilMInductanceAndForceComputeAll(int configCount, int threadCount)
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
        forceAndTorque = Coil::computeAllForceTorqueArrangements(prim, sec, primPositions, secPositions,
                                                                 primYAxisAngle, primZAxisAngle,
                                                                 secYAxisAngle, secZAxisAngle,
                                                                 precisionFactor, CPU_MT);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        performance = configCount / interval;
        printf("MT  perf : %6.3f milliseconds/Op | %.1f Ops/s\n", 1e3 / performance, performance);

        begin_time = high_resolution_clock::now();
        forceAndTorque = Coil::computeAllForceTorqueArrangements(prim, sec, primPositions, secPositions,
                                                                 primYAxisAngle, primZAxisAngle,
                                                                 secYAxisAngle, secZAxisAngle,
                                                                 precisionFactor, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        performance = configCount / interval;
        printf("GPU perf : %6.2f microseconds/Op | %.0f Ops/s\n", 1e6 / performance, performance);

        printf("\n");
    }
}

void Benchmark::coilMInductanceAndForceComputeAllMTvsMTD(PrecisionFactor precisionFactor, int threadCount)
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
        forceAndTorqueMT[i] = Coil::computeForceTorque(prim, sec, precisionFactor, CPU_MT);
    }
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MT perf  : %.1f Ops/s\n", opCount / interval);

    begin_time = high_resolution_clock::now();
    forceAndTorqueAll = Coil::computeAllForceTorqueArrangements(prim, sec, primPositions, secPositions,
                                                                primYAxisAngle, primZAxisAngle,
                                                                secYAxisAngle, secZAxisAngle,
                                                                precisionFactor, CPU_MT);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("MTD perf : %.1f Ops/s\n", opCount / interval);

    printf("\n");
}

void Benchmark::coilMInductanceAndForceComputeAllGPU(int configCount)
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

    printf("Benchmarking mutual inductance, force and torque GPU performance for %d configurations:\n\n", configCount);

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
        forceTorque = Coil::computeAllForceTorqueArrangements(prim, sec, primPositions, secPositions,
                                                              primYAxisAngle, primZAxisAngle,
                                                              secYAxisAngle, secZAxisAngle,
                                                              precisionFactor, GPU);
        interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        performance = configCount / interval;
        printf("Force and torque  : %6.2f microseconds/Op | %.0f Ops/s\n", 1e6 / performance, performance);


        printf("\n");
    }
}
