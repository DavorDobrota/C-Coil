#include <sstream>

#include <Coil.h>
#include <CoilGroup.h>
#include <ComputeMethod.h>
#include <PyMain.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


void initCoilGroup(py::module_ &mainModule)
{
    py::module_ coilGroupModule = mainModule.def_submodule("coil_group");

    py::class_<CoilGroup> coilGroup(coilGroupModule, "CoilGroup");


    // CoilGroup

    coilGroup.def(
        py::init<std::vector<Coil>, PrecisionFactor, int>(),
        py::arg("member_coils") = std::vector<Coil>(),
        py::arg("precision_factor") = PrecisionFactor(),
        py::arg("thread_count") = defaultThreadCount);

    coilGroup.def("get_default_precision_factor", &CoilGroup::getDefaultPrecisionFactor)
        .def("get_thread_count", &CoilGroup::getThreadCount)
        .def("get_member_coils", &CoilGroup::getMemberCoils);

    coilGroup.def(
            "set_default_precision_factor", &CoilGroup::setDefaultPrecisionFactor,
            py::arg("precision_factor") = PrecisionFactor())
        .def("set_thread_count", &CoilGroup::setThreadCount, py::arg("thread_count"))
        .def("add_coil", &CoilGroup::addCoil, py::arg("coil"));

    coilGroup.def("compute_B_field_vector", &CoilGroup::computeBFieldVector, py::arg("point_vector"))
        .def("compute_A_potential_vector", &CoilGroup::computeAPotentialVector, py::arg("point_vector"))
        .def("compute_E_field_vector", &CoilGroup::computeEFieldVector, py::arg("point_vector"))
        .def("compute_B_gradient_tensor", &CoilGroup::computeBGradientTensor, py::arg("point_vector"));

    coilGroup.def(
            "compute_all_A_potential_components", &CoilGroup::computeAllAPotentialComponents,
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_A_potential_x", &CoilGroup::computeAllAPotentialX,
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_A_potential_y", &CoilGroup::computeAllAPotentialY,
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_A_potential_z", &CoilGroup::computeAllAPotentialZ,
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_A_potential_abs", &CoilGroup::computeAllAPotentialAbs,
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST);

    coilGroup.def(
            "compute_all_B_field_components", &CoilGroup::computeAllBFieldComponents,
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_B_field_x", &CoilGroup::computeAllBFieldX,
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_B_field_y", &CoilGroup::computeAllBFieldY,
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_B_field_z", &CoilGroup::computeAllBFieldZ,
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_B_field_abs", &CoilGroup::computeAllBFieldAbs,
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST);

    coilGroup.def(
            "compute_all_E_field_components", &CoilGroup::computeAllEFieldComponents,
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_E_field_x", &CoilGroup::computeAllEFieldX,
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_E_field_y", &CoilGroup::computeAllEFieldY,
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_E_field_z", &CoilGroup::computeAllEFieldZ,
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_E_field_abs", &CoilGroup::computeAllEFieldAbs,
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST);

    coilGroup.def(
            "compute_all_B_gradient_tensors", &CoilGroup::computeAllBGradientTensors,
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST);

    coilGroup.def(
            "compute_mutual_inductance", &CoilGroup::computeMutualInductance,
            py::arg("secondary"),
            py::arg("precision_factor") = PrecisionFactor(),
            py::arg("compute_method") = CPU_ST)
        .def(
            "compute_ampere_force", &CoilGroup::computeAmpereForce,
            py::arg("secondary"),
            py::arg("precision_factor") = PrecisionFactor(),
            py::arg("compute_method") = CPU_ST);

    coilGroup.def(
        "compute_force_on_dipole_moment", &CoilGroup::computeForceOnDipoleMoment,
        py::arg("point_vector"), py::arg("dipole_moment"));

    coilGroup.def("__repr__", &CoilGroup::operator std::string);
}