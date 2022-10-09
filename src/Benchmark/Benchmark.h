#ifndef GENERAL_COIL_PROGRAM_BENCHMARK_H
#define GENERAL_COIL_PROGRAM_BENCHMARK_H

#include "Coil.h"
#include "CoilGroup.h"


const int g_maxPot = 22;
const int g_maxPotGroup = 16;

/// @brief Contains functions that are used to benchmark performance of various methods concerning Coil and CoilGroup.
/// Performance data obtained from them is very useful for assessing implementation efficiency.
namespace Benchmark
{
    /// @brief Benchmarks performance of basic math functions such as +, cos, log, and atan2.
    void mathFunctions();

    /// @brief Benchmarks single thread (CPU_ST) performance of Coil field compute methods for a given number of points.
    void computeFieldsST(int opCount = 50'000);
    /// @brief Benchmarks performance of Coil field compute methods for different ComputeMethods, a given
    /// number of points, and number of threads. The number of repetitions is determined by repeatCount.
    void computeAllFields(PrecisionFactor precisionFactor = PrecisionFactor(),
                          int opCount = 20'000, int repeatCount = 1, int threadCount = g_defaultThreadCount);

    /// @brief A very comprehensive benchmark of precision factor and threadCount influence on performance
    /// of field calculations. Filaments, thin, flat, and rectangular coils are included.
    void computeAllFieldsEveryCoilType(int opCount = 100'000, int threadCount = g_defaultThreadCount);

    /// @brief Benchmarks CPU_MT performance scaling of Coil field compute methods with the number of points for
    /// {2^0, 2^1, 2^2, ..., 2^maxPointsLog2}, given PrecisionFactor, and number of threads.
    void computeAllFieldsWorkloadScalingMT(PrecisionFactor precisionFactor = PrecisionFactor(),
                                           int threadCount = g_defaultThreadCount, int maxPointsLog2 = g_maxPot);
    /// @brief Benchmarks GPU performance scaling of Coil field compute methods with the number of points from
    /// {2^0, 2^1, 2^2, ..., 2^maxPointsLog2}, and given PrecisionFactor,
    void computeAllFieldsWorkloadScalingGPU(PrecisionFactor precisionFactor = PrecisionFactor(),
                                            int maxPointsLog2 = g_maxPot);

    /// @brief Benchmarks Coil::computeMutualInductance for z-axis case with two thick coils.
    void mInductanceZAxis(ComputeMethod computeMethod = CPU_ST, int threadCount = g_defaultThreadCount);
    /// @brief Benchmarks Coil::computeMutualInductance for z-axis case with two thick coils, CPU_MT ComputeMethod,
    /// and number of threads from {1,..., maxThreadCount}.
    void mInductanceZAxisMTScaling(int maxThreadCount = g_defaultThreadCount);
    /// @brief Benchmarks Coil::computeAndSetSelfInductance for a thick coil and PrecisionFactors from {1, 2,..., 15}.
    void selfInductance();

    /// @brief Benchmarks Coil::computeMutualInductance for general case with two thick coils.
    void mInductanceGeneral(ComputeMethod computeMethod = CPU_ST, int threadCount = g_defaultThreadCount);
    /// @brief Benchmarks Coil::computeMutualInductance for general case with two thick coils, CPU_MT ComputeMethod,
    /// and number of threads from {1,..., maxThreadCount}.
    void mInductanceGeneralMTScaling(int maxThreadCount = g_defaultThreadCount);

    /// @brief Benchmarks Coil::computeForceTorque for z-axis case with two thick coils.
    void forceZAxis(ComputeMethod computeMethod = CPU_ST, int threadCount = g_defaultThreadCount);
    /// @brief Benchmarks Coil::computeForceTorque for z-axis case with two thick coils, CPU_MT ComputeMethod,
    /// and number of threads from {1,..., maxThreadCount}.
    void forceZAxisMTScaling(int maxThreadCount = g_defaultThreadCount);

    /// @brief Benchmarks Coil::computeForceTorque for general case with two thick coils.
    void forceGeneral(ComputeMethod computeMethod = CPU_ST, int threadCount = g_defaultThreadCount);
    /// @brief Benchmarks Coil::computeForceTorque for general case with two thick coils for CPU_MT ComputeMethod,
    /// and number of threads from {1,..., maxThreadCount}.
    void forceGeneralMTScaling(int maxThreadCount = g_defaultThreadCount);

    /// @brief Benchmarks Coil::computeAllMutualInductanceArrangements and Coil::computeAllForceTorqueArrangements
    /// for a given number of configurations.
    void coilMInductanceAndForceComputeAll(int configCount = 100, int threadCount = g_defaultThreadCount);
    /// @brief Benchmarks CPU_MT (MT vs MTD) performance of Coil::computeAllMutualInductanceArrangements and
    /// Coil::computeAllForceTorqueArrangements for given PrecisionFactor and threadCount.
    void coilMInductanceAndForceComputeAllMTvsMTD(PrecisionFactor precisionFactor = PrecisionFactor(),
                                                  int threadCount = g_defaultThreadCount);
    /// @brief Benchmarks GPU performance of Coil::computeAllMutualInductanceArrangements and
    /// Coil::computeAllForceTorqueArrangements for a given number of configurations.
    void coilMInductanceAndForceComputeAllGPU(int configCount = 10'000);

    /// @brief Benchmarks performance of CoilGroup field calculations with different compute methods
    /// with given PrecisionFactor, number of coils, and number of points.
    void coilGroupComputeAllFields(PrecisionFactor precisionFactor = PrecisionFactor(),
                                   int coilCount = 50, int opCount = 100'000, int threadCount = g_defaultThreadCount);
    /// @brief Benchmarks performance of CPU_MT ComputeMethod (MT vs MTD) for CoilGroup fieldCalculations.
    void coilGroupComputeAllFieldsMTvsMTD(int threadCount = g_defaultThreadCount, int pointCount = 20'000);
    /// @brief Benchmarks GPU performance of CoilGroup fieldCalculations for a given number of coils and points.
    void coilGroupComputeAllFieldsGPU(int coilCount = 100, int opCount = 131'072);

    /// @brief Benchmarks CPU_MT performance scaling of CoilGroup field compute methods with the number of points from
    /// {2^0, 2^1, 2^2, ..., 2^maxPointsLog2}, given PrecisionFactor, number of coils, and number of threads.
    void coilGroupComputeAllFieldsMTScaling(PrecisionFactor precisionFactor = PrecisionFactor(),
                                            int threadCount = g_defaultThreadCount,
                                            int coilCount = 100, int maxPointsLog2 = g_maxPotGroup);
    /// @brief Benchmarks GPU performance scaling of CoilGroup field compute methods with the number of points for
    /// {2^0, 2^1, 2^2, ..., 2^maxPointsLog2}, given PrecisionFactor, and number of coils.
    void coilGroupComputeAllFieldsGPUScaling(PrecisionFactor precisionFactor = PrecisionFactor(),
                                             int coilCount = 100, int maxPointsLog2 = g_maxPotGroup);

    /// @brief Benchmarks performance of CoilGroup::computeMInductance and CoilGroup::computeForceTorque for
    /// different ComputeMethods and a given number of threads. The number of repetitions is determined by repeatCount.
    void coilGroupMInductanceAndForce(int repeatCount = 2, int threadCount = g_defaultThreadCount);
    /// @brief Benchmarks performance of CoilGroup::computeAllMutualInductanceArrangements and
    /// CoilGroup::computeAllForceTorqueArrangements for different ComputeMethods and a given number of threads.
    /// The number of repetitions is determined by repeatCount.
    void coilGroupMInductanceAndForceAll(int coilCount = 50, int opCount = 10, int threadCount = g_defaultThreadCount);
    /// @brief Benchmarks GPU performance of CoilGroup::computeAllMutualInductanceArrangements and
    /// CoilGroup::computeAllForceTorqueArrangements and a given number of threads.
    /// The number of repetitions is determined by repeatCount.
    void coilGroupMInductanceAndForceAllGPU(int coilCount = 100, int opCount = 1000);
}


#endif //GENERAL_COIL_PROGRAM_BENCHMARK_H
