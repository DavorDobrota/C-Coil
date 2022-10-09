#include <Compare.h>
#include <PyMain.h>

#include <pybind11/pybind11.h>
namespace py = pybind11;


void initCompare(py::module_ &mainModule)
{
    py::module_ compareModule = mainModule.def_submodule("compare");

    compareModule.def("fields_precision_CPU_vs_GPU", Compare::fieldsPrecisionCPUvsGPU)
                 .def("mutual_inductance_and_force_torque_precision_CPU_vs_GPU",
                      Compare::mutualInductanceAndForceTorquePrecisionCPUvsGPU)
                 .def("mutual_inductance_special_case", Compare::mutualInductanceSpecialCase);

    compareModule.def("force_torque_filaments_Z_axis", Compare::forceTorqueFilamentsZAxis)
                 .def("force_torque_thick_coils_general", Compare::forceTorqueThickCoilsGeneral)
                 .def("force_torque_thin_coils_Z_axis", Compare::forceTorqueThinCoilsZAxis)
                 .def("force_torque_filaments_general", Compare::forceTorqueFilamentsGeneral)
                 .def("force_torque_Z_axis", Compare::forceTorqueZAxis)
                 .def("force_on_dipole_vs_force_torque", Compare::forceOnDipoleVsForceTorque);

    compareModule.def("mutual_inductance_misaligned_coils", Compare::mutualInductanceMisalignedCoils)
                 .def("mutual_inductance_parallel_axes", Compare::mutualInductanceParallelAxes)
                 .def("mutual_inductance_parallel_axes_graphs", Compare::mutualInductanceParallelAxesGraphs)
                 .def("mutual_inductance_general_case", Compare::mutualInductanceGeneralCase)
                 .def("mutual_inductance_general_graphs", Compare::mutualInductanceGeneralGraphs)
                 .def("mutual_inductance_general_edge_cases", Compare::mutualInductanceGeneralEdgeCases);

    compareModule.def("mutual_inductance_Z_axis", Compare::mutualInductanceZAxis)
                 .def("self_inductance", Compare::selfInductance);

    compareModule.def(
        "fields_coil_group_MTD",Compare::fieldsCoilGroupMTD,
        py::arg("num_coils") = 100, py::arg("num_points") = 10'000,
        py::arg("thread_count") = g_defaultThreadCount, py::arg("print") = true);
}
