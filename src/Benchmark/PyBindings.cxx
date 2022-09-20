#include <Benchmark.h>
#include <PyMain.h>

#include <pybind11/pybind11.h>
namespace py = pybind11;


void initBenchmark(py::module_ &mainModule)
{
    py::module_ benchmarkModule = mainModule.def_submodule("benchmark");

    benchmarkModule.def("bench_math_functions", Benchmark::benchMathFunctions);

    benchmarkModule.def(
            "bench_compute_fields_ST", Benchmark::benchComputeFieldsST,
            py::arg("op_count") = 50'000)
        .def(
            "bench_compute_all_fields", Benchmark::benchComputeAllFields,
            py::arg("precision_factor") = PrecisionFactor(),
            py::arg("op_count") = 20'000, py::arg("repeat_count") = 1,
            py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "bench_compute_all_fields_every_coil_type", Benchmark::benchComputeAllFieldsEveryCoilType,
            py::arg("op_count") = 100'000, py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "bench_compute_all_fields_workload_scaling_MT", Benchmark::benchComputeAllFieldsWorkloadScalingMT,
            py::arg("precision_factor") = PrecisionFactor(),
            py::arg("thread_count") = g_defaultThreadCount,
            py::arg("max_points_log2") = g_maxPot)
        .def(
            "bench_compute_all_fields_workload_scaling_GPU", Benchmark::benchComputeAllFieldsWorkloadScalingGPU,
            py::arg("precision_factor") = PrecisionFactor(),
            py::arg("max_points_log2") = g_maxPot);

    benchmarkModule.def(
            "bench_m_inductance_Z_axis", Benchmark::benchMInductanceZAxis,
            py::arg("compute_method") = CPU_ST, py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "bench_m_inductance_Z_axis_MT_scaling", Benchmark::benchMInductanceZAxisMTScaling,
            py::arg("max_thread_count") = g_defaultThreadCount)
        .def("bench_self_inductance", Benchmark::benchSelfInductance);

    benchmarkModule.def(
            "bench_m_inductance_general", Benchmark::benchMInductanceGeneral,
            py::arg("compute_method") = CPU_ST, py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "bench_m_inductance_general_MT_scaling", Benchmark::benchMInductanceGeneralMTScaling,
            py::arg("max_thread_count") = g_defaultThreadCount);

    benchmarkModule.def(
            "bench_force_general", Benchmark::benchForceGeneral,
            py::arg("compute_method") = CPU_ST, py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "bench_force_general_MT_scaling", Benchmark::benchForceGeneralMTScaling,
            py::arg("max_thread_count") = g_defaultThreadCount);

    benchmarkModule.def(
            "bench_coil_m_inductance_and_force_compute_all_MT_vs_MTD",
            Benchmark::benchCoilMInductanceAndForceComputeAllMTvsMTD,
            py::arg("precision_factor") = PrecisionFactor(), py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "bench_coil_m_inductance_and_force_compute_all_GPU",
            Benchmark::benchCoilMInductanceAndForceComputeAllGPU, py::arg("config_count") = 10'000)
        .def(
            "bench_coil_m_inductance_and_force_compute_all", Benchmark::benchCoilMInductanceAndForceComputeAll,
            py::arg("config_count") = 10'000, py::arg("thread_count") = g_defaultThreadCount);

    benchmarkModule.def(
            "bench_coil_group_compute_all_fields_MT_vs_MTD", Benchmark::benchCoilGroupComputeAllFieldsMTvsMTD,
            py::arg("thread_count") = g_defaultThreadCount, py::arg("point_count") = 20'000)
        .def(
            "bench_coil_group_compute_all_fields", Benchmark::benchCoilGroupComputeAllFields,
            py::arg("precision_factor") = PrecisionFactor(), py::arg("num_coils") = 50,
            py::arg("op_count") = 100'000, py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "bench_coil_group_compute_all_fields_GPU", Benchmark::benchCoilGroupComputeAllFieldsGPU,
            py::arg("num_coils") = 100, py::arg("op_count") = 131'072)
        .def(
            "bench_coil_group_m_inductance_and_force", Benchmark::benchCoilGroupMInductanceAndForce,
            py::arg("op_count") = 2, py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "bench_coil_group_m_inductance_and_force_all", Benchmark::benchCoilGroupMInductanceAndForceAll,
            py::arg("coil_count") = 50, py::arg("op_count") = 10,
            py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "bench_coil_group_m_inductance_and_force_all_GPU", Benchmark::benchCoilGroupMInductanceAndForceAllGPU,
            py::arg("coil_count") = 50, py::arg("op_count") = 10);
}
