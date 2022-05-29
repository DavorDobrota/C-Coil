#include <Benchmark.h>
#include <PyMain.h>

#include <pybind11/pybind11.h>
namespace py = pybind11;


void initBenchmark(py::module_ &mainModule)
{
    py::module_ benchmarkModule = mainModule.def_submodule("benchmark");

    benchmarkModule.def("bench_math_functions", benchMathFunctions);

    benchmarkModule.def(
            "bench_compute_fields_ST", benchComputeFieldsST,
            py::arg("op_count") = 50'000)
        .def(
            "bench_compute_all_fields", benchComputeAllFields,
            py::arg("precision_factor") = PrecisionFactor(),
            py::arg("op_count") = 20'000, py::arg("repeat_count") = 1,
            py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "bench_compute_all_fields_every_coil_type", benchComputeAllFieldsEveryCoilType,
            py::arg("op_count") = 100'000, py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "bench_compute_all_fields_workload_scaling_MT", benchComputeAllFieldsWorkloadScalingMT,
            py::arg("precision_factor") = PrecisionFactor(),
            py::arg("thread_count") = g_defaultThreadCount,
            py::arg("max_points_log2") = g_maxPot);

    benchmarkModule.def(
            "bench_m_inductance_Z_axis", benchMInductanceZAxis,
            py::arg("compute_method") = CPU_ST, py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "bench_m_inductance_Z_axis_MT_scaling", benchMInductanceZAxisMTScaling,
            py::arg("max_thread_count") = g_defaultThreadCount)
        .def("bench_self_inductance", benchSelfInductance);

    benchmarkModule.def(
            "bench_m_inductance_general", benchMInductanceGeneral,
            py::arg("compute_method") = CPU_ST, py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "bench_m_inductance_general_MT_scaling", benchMInductanceGeneralMTScaling,
            py::arg("max_thread_count") = g_defaultThreadCount);

    benchmarkModule.def(
            "bench_force_Z_axis", benchForceZAxis,
            py::arg("compute_method") = CPU_ST, py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "bench_force_Z_axis_MT_scaling", benchForceZAxisMTScaling,
            py::arg("max_thread_count") = g_defaultThreadCount);

    benchmarkModule.def(
            "bench_force_general", benchForceGeneral,
            py::arg("compute_method") = CPU_ST, py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "bench_force_general_MT_scaling", benchForceGeneralMTScaling,
            py::arg("max_thread_count") = g_defaultThreadCount);

    benchmarkModule.def(
            "bench_coil_group_MT_vs_MTD", benchCoilGroupMTvsMTD,
            py::arg("thread_count") = g_defaultThreadCount, py::arg("point_count") = 20'000)
        .def(
            "bench_coil_group_compute_all_fields_MTD", benchCoilGroupComputeAllFieldsMTD,
            py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "bench_coil_group_m_inductance_and_force_MTD", benchCoilGroupMInductanceAndForceMTD,
            py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "bench_m_inductance_and_force_compute_all", benchMInductanceAndForceComputeAll,
            py::arg("precision_factor") = PrecisionFactor(), py::arg("thread_count") = g_defaultThreadCount);
}
