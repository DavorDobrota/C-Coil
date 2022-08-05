#include <Compare.h>
#include <PyMain.h>

#include <pybind11/pybind11.h>
namespace py = pybind11;


void initCompare(py::module_ &mainModule)
{
    py::module_ compareModule = mainModule.def_submodule("compare");

    compareModule.def("comp_method_precision_CPU_vs_GPU", compMethodPrecisionCPUvsGPU)
                 .def("comp_m_inductance_for_special_case", compMInductanceForSpecialCase);

    compareModule.def("comp_ampere_force_filaments_Z_axis", compAmpereForceFilamentsZAxis)
                 .def("comp_ampere_force_thick_coils_general", compAmpereForceThickCoilsGeneral)
                 .def("comp_ampere_force_thin_coils_Z_axis", compAmpereForceThinCoilsZAxis)
                 .def("comp_ampere_force_filaments_general", compAmpereForceFilamentsGeneral)
                 .def("comp_ampere_force_Z_axis", compAmpereForceZAxis)
                 .def("comp_force_on_dipole_vs_ampere_force", compForceOnDipoleVsAmpereForce);

    compareModule.def("comp_m_inductance_general_misaligned_coils", compMInductanceGeneralMisalignedCoils)
                 .def("comp_m_inductance_general_parallel_axes", compMInductanceGeneralParallelAxes)
                 .def("comp_m_inductance_general_case", compMInductanceGeneralCase)
                 .def("comp_m_inductance_general_graphs", compMInductanceGeneralGraphs)
                 .def("comp_m_inductance_general_edge_cases", compMInductanceGeneralEdgeCases);

    compareModule.def("comp_mutual_inductance_Z_axis", compMutualInductanceZAxis)
                 .def("comp_precision_CPU_vs_GPU", compPrecisionCPUvsGPU)
                 .def("comp_self_inductance", compSelfInductance);

    compareModule.def(
        "comp_coil_group_MTD",
        compCoilGroupMTD,
        py::arg("num_coils") = 100, py::arg("num_points") = 10'000,
        py::arg("thread_count") = g_defaultThreadCount, py::arg("print") = true);
}
