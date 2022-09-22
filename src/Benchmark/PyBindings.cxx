#include <Benchmark.h>
#include <PyMain.h>

#include <pybind11/pybind11.h>
namespace py = pybind11;


void initBenchmark(py::module_ &mainModule)
{
    py::module_ benchmarkModule = mainModule.def_submodule("benchmark");

    benchmarkModule.def("math_functions", Benchmark::mathFunctions);

    benchmarkModule.def(
            "compute_fields_ST", Benchmark::computeFieldsST,
            py::arg("op_count") = 50'000)
        .def(
            "compute_all_fields", Benchmark::computeAllFields,
            py::arg("precision_factor") = PrecisionFactor(),
            py::arg("op_count") = 20'000, py::arg("repeat_count") = 1,
            py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "compute_all_fields_every_coil_type", Benchmark::computeAllFieldsEveryCoilType,
            py::arg("op_count") = 100'000, py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "compute_all_fields_workload_scaling_MT", Benchmark::computeAllFieldsWorkloadScalingMT,
            py::arg("precision_factor") = PrecisionFactor(),
            py::arg("thread_count") = g_defaultThreadCount,
            py::arg("max_points_log2") = g_maxPot)
        .def(
            "compute_all_fields_workload_scaling_GPU", Benchmark::computeAllFieldsWorkloadScalingGPU,
            py::arg("precision_factor") = PrecisionFactor(),
            py::arg("max_points_log2") = g_maxPot);

    benchmarkModule.def(
            "m_inductance_Z_axis", Benchmark::mInductanceZAxis,
            py::arg("compute_method") = CPU_ST, py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "m_inductance_Z_axis_MT_scaling", Benchmark::mInductanceZAxisMTScaling,
            py::arg("max_thread_count") = g_defaultThreadCount)
        .def("self_inductance", Benchmark::selfInductance);

    benchmarkModule.def(
            "m_inductance_general", Benchmark::mInductanceGeneral,
            py::arg("compute_method") = CPU_ST, py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "m_inductance_general_MT_scaling", Benchmark::mInductanceGeneralMTScaling,
            py::arg("max_thread_count") = g_defaultThreadCount);

    benchmarkModule.def(
            "force_general", Benchmark::forceGeneral,
            py::arg("compute_method") = CPU_ST, py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "force_general_MT_scaling", Benchmark::forceGeneralMTScaling,
            py::arg("max_thread_count") = g_defaultThreadCount);

    benchmarkModule.def(
            "coil_m_inductance_and_force_compute_all_MT_vs_MTD",
            Benchmark::coilMInductanceAndForceComputeAllMTvsMTD,
            py::arg("precision_factor") = PrecisionFactor(), py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "coil_m_inductance_and_force_compute_all_GPU",
            Benchmark::coilMInductanceAndForceComputeAllGPU, py::arg("config_count") = 10'000)
        .def(
            "coil_m_inductance_and_force_compute_all", Benchmark::coilMInductanceAndForceComputeAll,
            py::arg("config_count") = 10'000, py::arg("thread_count") = g_defaultThreadCount);

    benchmarkModule.def(
            "coil_group_compute_all_fields_MT_vs_MTD", Benchmark::coilGroupComputeAllFieldsMTvsMTD,
            py::arg("thread_count") = g_defaultThreadCount, py::arg("point_count") = 20'000)
        .def(
            "coil_group_compute_all_fields", Benchmark::coilGroupComputeAllFields,
            py::arg("precision_factor") = PrecisionFactor(), py::arg("num_coils") = 50,
            py::arg("op_count") = 100'000, py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "coil_group_compute_all_fields_GPU", Benchmark::coilGroupComputeAllFieldsGPU,
            py::arg("num_coils") = 100, py::arg("op_count") = 131'072)
        .def(
            "coil_group_m_inductance_and_force", Benchmark::coilGroupMInductanceAndForce,
            py::arg("op_count") = 2, py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "coil_group_m_inductance_and_force_all", Benchmark::coilGroupMInductanceAndForceAll,
            py::arg("coil_count") = 50, py::arg("op_count") = 10,
            py::arg("thread_count") = g_defaultThreadCount)
        .def(
            "coil_group_m_inductance_and_force_all_GPU", Benchmark::coilGroupMInductanceAndForceAllGPU,
            py::arg("coil_count") = 50, py::arg("op_count") = 10);
}
