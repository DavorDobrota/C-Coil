#include <PyMain.h>
#include <Coil.h>
#include <CoilType.h>
#include <ComputeMethod.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


void initCoil(py::module_ &mainModule)
{
    py::module_ coilModule = mainModule.def_submodule("coil");

    py::enum_<ComputeMethod> computeMethod(coilModule, "ComputeMethod");
    py::enum_<CoilType> coilType(coilModule, "CoilType");
    py::class_<PrecisionFactor> precisionFactor(coilModule, "PrecisionFactor");
    py::class_<PrecisionArguments> precisionArguments(coilModule, "PrecisionArguments");
    py::class_<CoilPairArguments> coilPairArguments(coilModule, "CoilPairArguments");
    py::class_<Coil> coil(coilModule, "Coil");


    // ComputeMethod

    computeMethod.value("CPU_ST", ComputeMethod::CPU_ST)
        .value("CPU_MT", ComputeMethod::CPU_MT)
        .value("GPU", ComputeMethod::GPU)
        .export_values();


    // CoilType

    coilType.value("RECTANGULAR", CoilType::RECTANGULAR)
        .value("THIN", CoilType::THIN)
        .value("FLAT", CoilType::FLAT)
        .value("FILAMENT", CoilType::FILAMENT)
        .export_values();


    // PrecisionFactor

    precisionFactor.def(py::init<>())
        .def(py::init<double>(), py::arg("relative_precision"));

    precisionFactor.def_readwrite("relative_precision", &PrecisionFactor::relativePrecision);

    precisionFactor.def("__repr__", &PrecisionFactor::operator std::string);


    // PrecisionArguments

    precisionArguments.def(py::init<>())
        .def(
            py::init<int, int, int, int, int, int>(),
            py::arg("angular_blocks"), py::arg("thickness_blocks"), py::arg("length_blocks"),
            py::arg("angular_increments"), py::arg("thickness_increments"), py::arg("length_increments"));

    precisionArguments.def_readwrite("angular_block_count", &PrecisionArguments::angularBlockCount)
        .def_readwrite("thickness_block_count", &PrecisionArguments::thicknessBlockCount)
        .def_readwrite("length_block_count", &PrecisionArguments::lengthBlockCount)
        .def_readwrite("angular_increment_count", &PrecisionArguments::angularIncrementCount)
        .def_readwrite("thickness_increment_count", &PrecisionArguments::thicknessIncrementCount)
        .def_readwrite("length_increment_count", &PrecisionArguments::lengthIncrementCount);

    precisionArguments.def_static("get_coil_precision_arguments_CPU", &PrecisionArguments::getCoilPrecisionArgumentsCPU)
        .def_static("get_coil_precision_arguments_GPU", &PrecisionArguments::getCoilPrecisionArgumentsGPU);

    precisionArguments.def("__repr__", &PrecisionArguments::operator std::string);


    // CoilPairArguments

    coilPairArguments.def(py::init<>())
        .def(
            py::init<PrecisionArguments, PrecisionArguments>(),
            py::arg("primary_precision"), py::arg("secondary_precision"));

    coilPairArguments.def_readwrite("primary_precision", &CoilPairArguments::primaryPrecision)
        .def_readwrite("secondary_precision", &CoilPairArguments::secondaryPrecision);

    coilPairArguments.def_static(
        "get_appropriate_coil_pair_arguments", &CoilPairArguments::getAppropriateCoilPairArguments,
        py::arg("primary"), py::arg("secondary"),
        py::arg("precision_factor"),
        py::arg("computeMethod") = CPU_ST, py::arg("z_axis_case") = true);

    coilPairArguments.def("__repr__", &CoilPairArguments::operator std::string);


    // Coil

    coil.def(py::init<>())
        .def(
            py::init<double, double, double, int, double, double, double,
                     PrecisionFactor, int, vec3::CoordVector3, double, double>(),
            py::arg("inner_radius"), py::arg("thickness"), py::arg("length"), py::arg("num_of_turns"),
            py::arg("current"), py::arg("wire_resistivity"), py::arg("sine_frequency"),
            py::arg("precision_factor") = PrecisionFactor(), py::arg("thread_count") = defaultThreadCount,
            py::arg("coordinate_position") = vec3::CoordVector3(),
            py::arg("y_axis_angle") = 0.0, py::arg("z_axis_angle") = 0.0)
        .def(
            py::init<double, double, double, int, double, double, double,
                    PrecisionArguments, PrecisionArguments, int, vec3::CoordVector3, double, double>(),
            py::arg("inner_radius"), py::arg("thickness"), py::arg("length"), py::arg("num_of_turns"),
            py::arg("current"), py::arg("wire_resistivity"), py::arg("sine_frequency"),
            py::arg("precision_settings_CPU"), py::arg("precision_settings_GPU"),
            py::arg("thread_count") = defaultThreadCount, py::arg("coordinate_position") = vec3::CoordVector3(),
            py::arg("y_axis_angle") = 0.0, py::arg("z_axis_angle") = 0.0)
        .def(
            py::init<double, double, double, int, double, double, PrecisionFactor, int,
                    vec3::CoordVector3, double, double>(),
            py::arg("inner_radius"), py::arg("thickness"), py::arg("length"),
            py::arg("num_of_turns"), py::arg("current"), py::arg("sine_frequency"),
            py::arg("precision_factor") = PrecisionFactor(), py::arg("thread_count") = defaultThreadCount,
            py::arg("coordinate_position") = vec3::CoordVector3(),
            py::arg("y_axis_angle") = 0.0, py::arg("z_axis_angle") = 0.0)
        .def(
            py::init<double, double, double, int, double, double, PrecisionArguments, PrecisionArguments, int,
            vec3::CoordVector3, double, double>(),
            py::arg("inner_radius"), py::arg("thickness"), py::arg("length"), py::arg("num_of_turns"),
            py::arg("current"), py::arg("sine_frequency"),
            py::arg("precision_settings_CPU"), py::arg("precision_settings_GPU"),
            py::arg("thread_count") = defaultThreadCount, py::arg("coordinate_position") = vec3::CoordVector3(),
            py::arg("y_axis_angle") = 0.0, py::arg("z_axis_angle") = 0.0)
        .def(
            py::init<double, double, double, int, double, PrecisionFactor, int, vec3::CoordVector3, double, double>(),
            py::arg("inner_radius"), py::arg("thickness"), py::arg("length"),
            py::arg("num_of_turns"), py::arg("current"), py::arg("precision_factor") = PrecisionFactor(),
            py::arg("thread_count") = defaultThreadCount, py::arg("coordinate_position") = vec3::CoordVector3(),
            py::arg("y_axis_angle") = 0.0, py::arg("z_axis_angle") = 0.0)
        .def(
            py::init<double, double, double, int, double, PrecisionArguments, PrecisionArguments, int,
            vec3::CoordVector3, double, double>(),
            py::arg("inner_radius"), py::arg("thickness"), py::arg("length"), py::arg("num_of_turns"),
            py::arg("current"), py::arg("precision_settings_CPU"), py::arg("precision_settings_GPU"),
            py::arg("thread_count") = defaultThreadCount, py::arg("coordinate_position") = vec3::CoordVector3(),
            py::arg("y_axis_angle") = 0.0, py::arg("z_axis_angle") = 0.0)
        .def(
            py::init<double, double, double, int, PrecisionFactor, int, vec3::CoordVector3, double, double>(),
            py::arg("inner_radius"), py::arg("thickness"), py::arg("length"),
            py::arg("num_of_turns"), py::arg("precision_factor") = PrecisionFactor(),
            py::arg("thread_count") = defaultThreadCount, py::arg("coordinate_position") = vec3::CoordVector3(),
            py::arg("y_axis_angle") = 0.0, py::arg("z_axis_angle") = 0.0)
        .def(
            py::init<double, double, double, int, PrecisionArguments, PrecisionArguments, int,
            vec3::CoordVector3, double, double>(),
            py::arg("inner_radius"), py::arg("thickness"), py::arg("length"), py::arg("num_of_turns"),
            py::arg("precision_settings_CPU"), py::arg("precision_settings_GPU"),
            py::arg("thread_count") = defaultThreadCount, py::arg("coordinate_position") = vec3::CoordVector3(),
            py::arg("y_axis_angle") = 0.0, py::arg("z_axis_angle") = 0.0);

    coil.def("get_id", &Coil::getId)
        .def("get_inner_radius", &Coil::getInnerRadius)
        .def("get_thickness", &Coil::getThickness)
        .def("get_length", &Coil::getLength)
        .def("get_num_of_turns", &Coil::getNumOfTurns);

    coil.def("get_current_density", &Coil::getCurrentDensity)
        .def("get_current", &Coil::getCurrent);

    coil.def("get_wire_resistivity", &Coil::getWireResistivity)
        .def("is_sine_driven", &Coil::isSineDriven)
        .def("get_sine_frequency", &Coil::getSineFrequency);

    coil.def("get_magnetic_moment", &Coil::getMagneticMoment)
        .def("get_average_wire_thickness", &Coil::getAverageWireThickness);

    coil.def("get_self_inductance", &Coil::getSelfInductance)
        .def("get_resistance", &Coil::getResistance)
        .def("get_reactance", &Coil::getReactance)
        .def("get_impedance", &Coil::getImpedance);

    coil.def("get_precision_settings_CPU", &Coil::getPrecisionSettingsCPU)
        .def("get_precision_settings_GPU", &Coil::getPrecisionSettingsGPU)
        .def("get_thread_count", &Coil::getThreadCount)
        .def("is_using_fast_method", &Coil::isUsingFastMethod)
        .def("get_coil_type", &Coil::getCoilType);

    coil.def("get_position_vector", &Coil::getPositionVector)
        .def("get_rotation_angles", &Coil::getRotationAngles);

    coil.def("set_current_density", &Coil::setCurrentDensity, py::arg("current_density"))
        .def("set_current", &Coil::setCurrent, py::arg("current"))
        .def("set_wire_resistivity", &Coil::setWireResistivity, py::arg("wire_resistivity"))
        .def("set_sine_frequency", &Coil::setSineFrequency, py::arg("sine_frequency"))
        .def(
            "set_default_precision_GPU",
            static_cast<void (Coil::*)(const PrecisionArguments&)>(&Coil::setDefaultPrecisionGPU),
            py::arg("precision_settings"))
        .def(
            "set_default_precision_GPU",
            static_cast<void (Coil::*)(PrecisionFactor)>(&Coil::setDefaultPrecisionGPU),
            py::arg("precision_factor"))
        .def(
            "set_default_precision_CPU",
            static_cast<void (Coil::*)(const PrecisionArguments&)>(&Coil::setDefaultPrecisionCPU),
            py::arg("precision_settings"))
        .def(
            "set_default_precision_CPU",
            static_cast<void (Coil::*)(PrecisionFactor)>(&Coil::setDefaultPrecisionCPU),
            py::arg("precision_factor"))
        .def("set_thread_count", &Coil::setThreadCount, py::arg("thread_count"));

    coil.def("set_self_inductance", &Coil::setSelfInductance, py::arg("self_inductance"));

    coil.def(
            "set_position_and_orientation", &Coil::setPositionAndOrientation,
            py::arg("position_vector"), py::arg("y_axis_angle"), py::arg("z_axis_angle"));

    coil.def(
            "compute_B_field_x",
            static_cast<double (Coil::*)(vec3::CoordVector3) const>(&Coil::computeBFieldX),
            py::arg("point_vector"))
        .def(
            "compute_B_field_x",
            static_cast<double (Coil::*)(vec3::CoordVector3, const PrecisionArguments&) const>(&Coil::computeBFieldX),
            py::arg("point_vector"), py::arg("used_precision"));

    coil.def(
            "compute_B_field_y",
            static_cast<double (Coil::*)(vec3::CoordVector3) const>(&Coil::computeBFieldY),
            py::arg("point_vector"))
        .def(
            "compute_B_field_y",
            static_cast<double (Coil::*)(vec3::CoordVector3, const PrecisionArguments&) const>(&Coil::computeBFieldY),
            py::arg("point_vector"), py::arg("used_precision"));

    coil.def(
            "compute_B_field_z",
            static_cast<double (Coil::*)(vec3::CoordVector3) const>(&Coil::computeBFieldZ),
            py::arg("point_vector"))
        .def(
            "compute_B_field_z",
            static_cast<double (Coil::*)(vec3::CoordVector3, const PrecisionArguments&) const>(&Coil::computeBFieldZ),
            py::arg("point_vector"), py::arg("used_precision"));

    coil.def(
            "compute_B_field_abs",
            static_cast<double (Coil::*)(vec3::CoordVector3) const>(&Coil::computeBFieldAbs),
            py::arg("point_vector"))
        .def(
            "compute_B_field_abs",
            static_cast<double (Coil::*)(vec3::CoordVector3, const PrecisionArguments&) const>(&Coil::computeBFieldAbs),
            py::arg("point_vector"), py::arg("used_precision"));

    coil.def(
            "compute_B_field_vector",
            static_cast<vec3::Vector3 (Coil::*)(vec3::CoordVector3) const>(&Coil::computeBFieldVector),
            py::arg("point_vector"))
        .def(
            "compute_B_field_vector",
            static_cast<vec3::Vector3 (Coil::*)(vec3::CoordVector3, const PrecisionArguments&) const>(&Coil::computeBFieldVector),
            py::arg("point_vector"), py::arg("used_precision"));

    coil.def(
            "compute_A_potential_x",
            static_cast<double (Coil::*)(vec3::CoordVector3) const>(&Coil::computeAPotentialX),
            py::arg("point_vector"))
        .def(
            "compute_A_potential_x",
            static_cast<double (Coil::*)(vec3::CoordVector3, const PrecisionArguments&) const>(&Coil::computeAPotentialX),
            py::arg("point_vector"), py::arg("used_precision"));

    coil.def(
            "compute_A_potential_y",
            static_cast<double (Coil::*)(vec3::CoordVector3) const>(&Coil::computeAPotentialY),
            py::arg("point_vector"))
        .def(
            "compute_A_potential_y",
            static_cast<double (Coil::*)(vec3::CoordVector3, const PrecisionArguments&) const>(&Coil::computeAPotentialY),
            py::arg("point_vector"), py::arg("used_precision"));

    coil.def(
            "compute_A_potential_z",
            static_cast<double (Coil::*)(vec3::CoordVector3) const>(&Coil::computeAPotentialZ),
            py::arg("point_vector"))
        .def(
            "compute_A_potential_z",
            static_cast<double (Coil::*)(vec3::CoordVector3, const PrecisionArguments&) const>(&Coil::computeAPotentialZ),
            py::arg("point_vector"), py::arg("used_precision"));

    coil.def(
            "compute_A_potential_abs",
            static_cast<double (Coil::*)(vec3::CoordVector3) const>(&Coil::computeAPotentialAbs),
            py::arg("point_vector"))
        .def(
            "compute_A_potential_abs",
            static_cast<double (Coil::*)(vec3::CoordVector3, const PrecisionArguments&) const>(&Coil::computeAPotentialAbs),
            py::arg("point_vector"), py::arg("used_precision"));

    coil.def(
            "compute_A_potential_vector",
            static_cast<vec3::Vector3 (Coil::*)(vec3::CoordVector3) const>(&Coil::computeAPotentialVector),
            py::arg("point_vector"))
        .def(
            "compute_A_potential_vector",
            static_cast<vec3::Vector3 (Coil::*)(vec3::CoordVector3, const PrecisionArguments&) const>(&Coil::computeAPotentialVector),
            py::arg("point_vector"), py::arg("used_precision"));

    coil.def(
            "compute_E_field_x",
            static_cast<double (Coil::*)(vec3::CoordVector3) const>(&Coil::computeEFieldX),
            py::arg("point_vector"))
        .def(
            "compute_E_field_x",
            static_cast<double (Coil::*)(vec3::CoordVector3, const PrecisionArguments&) const>(&Coil::computeEFieldX),
            py::arg("point_vector"), py::arg("used_precision"));

    coil.def(
            "compute_E_field_y",
            static_cast<double (Coil::*)(vec3::CoordVector3) const>(&Coil::computeEFieldY),
            py::arg("point_vector"))
        .def(
            "compute_E_field_y",
            static_cast<double (Coil::*)(vec3::CoordVector3, const PrecisionArguments&) const>(&Coil::computeEFieldY),
            py::arg("point_vector"), py::arg("used_precision"));

    coil.def(
            "compute_E_field_z",
            static_cast<double (Coil::*)(vec3::CoordVector3) const>(&Coil::computeEFieldZ),
            py::arg("point_vector"))
        .def(
            "compute_E_field_z",
            static_cast<double (Coil::*)(vec3::CoordVector3, const PrecisionArguments&) const>(&Coil::computeEFieldZ),
            py::arg("point_vector"), py::arg("used_precision"));

    coil.def(
            "compute_E_field_abs",
            static_cast<double (Coil::*)(vec3::CoordVector3) const>(&Coil::computeEFieldAbs),
            py::arg("point_vector"))
        .def(
            "compute_E_field_abs",
            static_cast<double (Coil::*)(vec3::CoordVector3, const PrecisionArguments&) const>(&Coil::computeEFieldAbs),
            py::arg("point_vector"), py::arg("used_precision"));

    coil.def(
            "compute_E_field_vector",
            static_cast<vec3::Vector3 (Coil::*)(vec3::CoordVector3) const>(&Coil::computeEFieldVector),
            py::arg("point_vector"))
        .def(
            "compute_E_field_vector",
            static_cast<vec3::Vector3 (Coil::*)(vec3::CoordVector3, const PrecisionArguments&) const>(&Coil::computeEFieldVector),
            py::arg("point_vector"), py::arg("used_precision"));

    coil.def(
            "compute_B_gradient_tensor",
            static_cast<vec3::Matrix3 (Coil::*)(vec3::CoordVector3) const>(&Coil::computeBGradientTensor),
            py::arg("point_vector"))
        .def(
            "compute_B_gradient_tensor",
            static_cast<vec3::Matrix3 (Coil::*)(vec3::CoordVector3, const PrecisionArguments&) const>(&Coil::computeBGradientTensor),
            py::arg("point_vector"), py::arg("used_precision"));

    coil.def(
            "compute_all_B_field_x",
            static_cast<std::vector<double> (Coil::*)(const std::vector<vec3::CoordVector3>&, ComputeMethod) const>(&Coil::computeAllBFieldX),
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_B_field_x",
            static_cast<std::vector<double> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, const PrecisionArguments&, ComputeMethod
            ) const>(&Coil::computeAllBFieldX),
            py::arg("point_vectors"), py::arg("used_precision"), py::arg("compute_method") = CPU_ST);

    coil.def(
            "compute_all_B_field_y",
            static_cast<std::vector<double> (Coil::*)(const std::vector<vec3::CoordVector3>&, ComputeMethod) const>(&Coil::computeAllBFieldY),
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_B_field_y",
            static_cast<std::vector<double> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, const PrecisionArguments&, ComputeMethod
            ) const>(&Coil::computeAllBFieldY),
            py::arg("point_vectors"), py::arg("used_precision"), py::arg("compute_method") = CPU_ST);

    coil.def(
            "compute_all_B_field_z",
            static_cast<std::vector<double> (Coil::*)(const std::vector<vec3::CoordVector3>&, ComputeMethod) const>(&Coil::computeAllBFieldZ),
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_B_field_z",
            static_cast<std::vector<double> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, const PrecisionArguments&, ComputeMethod
            ) const>(&Coil::computeAllBFieldZ),
            py::arg("point_vectors"), py::arg("used_precision"), py::arg("compute_method") = CPU_ST);

    coil.def(
            "compute_all_B_field_abs",
            static_cast<std::vector<double> (Coil::*)(const std::vector<vec3::CoordVector3>&, ComputeMethod) const>(&Coil::computeAllBFieldAbs),
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_B_field_abs",
            static_cast<std::vector<double> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, const PrecisionArguments&, ComputeMethod
            ) const>(&Coil::computeAllBFieldAbs),
            py::arg("point_vectors"), py::arg("used_precision"), py::arg("compute_method") = CPU_ST);

    coil.def(
            "compute_all_B_field_components",
            static_cast<std::vector<vec3::Vector3> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, ComputeMethod
            ) const>(&Coil::computeAllBFieldComponents),
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_B_field_components",
            static_cast<std::vector<vec3::Vector3> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, const PrecisionArguments&, ComputeMethod
            ) const>(&Coil::computeAllBFieldComponents),
            py::arg("point_vectors"), py::arg("used_precision"), py::arg("compute_method") = CPU_ST);

    coil.def(
            "compute_all_A_potential_x",
            static_cast<std::vector<double> (Coil::*)(const std::vector<vec3::CoordVector3>&, ComputeMethod) const>(&Coil::computeAllAPotentialX),
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_A_potential_x",
            static_cast<std::vector<double> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, const PrecisionArguments&, ComputeMethod
            ) const>(&Coil::computeAllAPotentialX),
            py::arg("point_vectors"), py::arg("used_precision"), py::arg("compute_method") = CPU_ST);

    coil.def(
            "compute_all_A_potential_y",
            static_cast<std::vector<double> (Coil::*)(const std::vector<vec3::CoordVector3>&, ComputeMethod) const>(&Coil::computeAllAPotentialY),
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_A_potential_y",
            static_cast<std::vector<double> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, const PrecisionArguments&, ComputeMethod
            ) const>(&Coil::computeAllAPotentialY),
            py::arg("point_vectors"), py::arg("used_precision"), py::arg("compute_method") = CPU_ST);

    coil.def(
            "compute_all_A_potential_z",
            static_cast<std::vector<double> (Coil::*)(const std::vector<vec3::CoordVector3>&, ComputeMethod) const>(&Coil::computeAllAPotentialZ),
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_A_potential_z",
            static_cast<std::vector<double> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, const PrecisionArguments&, ComputeMethod
            ) const>(&Coil::computeAllAPotentialZ),
            py::arg("point_vectors"), py::arg("used_precision"), py::arg("compute_method") = CPU_ST);

    coil.def(
            "compute_all_A_potential_abs",
            static_cast<std::vector<double> (Coil::*)(const std::vector<vec3::CoordVector3>&, ComputeMethod) const>(&Coil::computeAllAPotentialAbs),
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_A_potential_abs",
            static_cast<std::vector<double> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, const PrecisionArguments&, ComputeMethod
            ) const>(&Coil::computeAllAPotentialAbs),
            py::arg("point_vectors"), py::arg("used_precision"), py::arg("compute_method") = CPU_ST);

    coil.def(
            "compute_all_A_potential_components",
            static_cast<std::vector<vec3::Vector3> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, ComputeMethod
            ) const>(&Coil::computeAllAPotentialComponents),
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_A_potential_components",
            static_cast<std::vector<vec3::Vector3> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, const PrecisionArguments&, ComputeMethod
            ) const>(&Coil::computeAllAPotentialComponents),
            py::arg("point_vectors"), py::arg("used_precision"), py::arg("compute_method") = CPU_ST);

    coil.def(
            "compute_all_E_field_x",
            static_cast<std::vector<double> (Coil::*)(const std::vector<vec3::CoordVector3>&, ComputeMethod) const>(&Coil::computeAllEFieldX),
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_E_field_x",
            static_cast<std::vector<double> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, const PrecisionArguments&, ComputeMethod
            ) const>(&Coil::computeAllEFieldX),
            py::arg("point_vectors"), py::arg("used_precision"), py::arg("compute_method") = CPU_ST);

    coil.def(
            "compute_all_E_field_y",
            static_cast<std::vector<double> (Coil::*)(const std::vector<vec3::CoordVector3>&, ComputeMethod) const>(&Coil::computeAllEFieldY),
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_E_field_y",
            static_cast<std::vector<double> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, const PrecisionArguments&, ComputeMethod
            ) const>(&Coil::computeAllEFieldY),
            py::arg("point_vectors"), py::arg("used_precision"), py::arg("compute_method") = CPU_ST);

    coil.def(
            "compute_all_E_field_z",
            static_cast<std::vector<double> (Coil::*)(const std::vector<vec3::CoordVector3>&, ComputeMethod) const>(&Coil::computeAllEFieldZ),
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_E_field_z",
            static_cast<std::vector<double> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, const PrecisionArguments&, ComputeMethod
            ) const>(&Coil::computeAllEFieldZ),
            py::arg("point_vectors"), py::arg("used_precision"), py::arg("compute_method") = CPU_ST);

    coil.def(
            "compute_all_E_field_abs",
            static_cast<std::vector<double> (Coil::*)(const std::vector<vec3::CoordVector3>&, ComputeMethod) const>(&Coil::computeAllEFieldAbs),
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_E_field_abs",
            static_cast<std::vector<double> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, const PrecisionArguments&, ComputeMethod
            ) const>(&Coil::computeAllEFieldAbs),
            py::arg("point_vectors"), py::arg("used_precision"), py::arg("compute_method") = CPU_ST);

    coil.def(
            "compute_all_E_field_components",
            static_cast<std::vector<vec3::Vector3> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, ComputeMethod
            ) const>(&Coil::computeAllEFieldComponents),
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_E_field_components",
            static_cast<std::vector<vec3::Vector3> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, const PrecisionArguments&, ComputeMethod
            ) const>(&Coil::computeAllEFieldComponents),
            py::arg("point_vectors"), py::arg("used_precision"), py::arg("compute_method") = CPU_ST);

    coil.def(
            "compute_all_B_gradient_tensors",
            static_cast<std::vector<vec3::Matrix3> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, ComputeMethod
            ) const>(&Coil::computeAllBGradientTensors),
            py::arg("point_vectors"), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_all_B_gradient_tensors",
            static_cast<std::vector<vec3::Matrix3> (Coil::*)(
                const std::vector<vec3::CoordVector3>&, const PrecisionArguments&, ComputeMethod
            ) const>(&Coil::computeAllBGradientTensors),
            py::arg("point_vectors"), py::arg("used_precision"), py::arg("compute_method") = CPU_ST);

    coil.def_static(
            "compute_mutual_inductance",
            static_cast<double (*)(const Coil&, const Coil&, PrecisionFactor, ComputeMethod)>(&Coil::computeMutualInductance),
            py::arg("primary"), py::arg("secondary"),
            py::arg("precision_factor") = PrecisionFactor(), py::arg("compute_method") = CPU_ST)
        .def_static(
            "compute_mutual_inductance",
            static_cast<double (*)(const Coil&, const Coil&, CoilPairArguments, ComputeMethod)>(&Coil::computeMutualInductance),
            py::arg("primary"), py::arg("secondary"), py::arg("inductance_arguments"), py::arg("compute_method") = CPU_ST);

    coil.def(
            "compute_secondary_induced_voltage",
            static_cast<double (Coil::*)(const Coil&, PrecisionFactor, ComputeMethod) const>(&Coil::computeSecondaryInducedVoltage),
            py::arg("secondary"), py::arg("precision_factor") = PrecisionFactor(), py::arg("compute_method") = CPU_ST)
        .def(
            "compute_secondary_induced_voltage",
            static_cast<double (Coil::*)(const Coil&, CoilPairArguments, ComputeMethod) const>(&Coil::computeSecondaryInducedVoltage),
            py::arg("secondary"), py::arg("inductance_arguments"), py::arg("compute_method") = CPU_ST);

    coil.def(
            "compute_and_set_self_inductance", &Coil::computeAndSetSelfInductance,
            py::arg("precision_factor") = PrecisionFactor(), py::arg("compute_method") = CPU_ST);

    coil.def_static(
            "compute_ampere_force",
            static_cast<std::pair<vec3::Vector3, vec3::Vector3> (*)(
                const Coil&, const Coil&, PrecisionFactor, ComputeMethod
            )>(&Coil::computeAmpereForce),
            py::arg("primary"), py::arg("secondary"),
            py::arg("precision_factor") = PrecisionFactor(), py::arg("compute_method") = CPU_ST)
        .def_static(
            "compute_ampere_force",
            static_cast<std::pair<vec3::Vector3, vec3::Vector3> (*)(
                    const Coil&, const Coil&, CoilPairArguments, ComputeMethod
            )>(&Coil::computeAmpereForce),
            py::arg("primary"), py::arg("secondary"), py::arg("force_arguments"), py::arg("compute_method") = CPU_ST);

    coil.def(
            "compute_force_on_dipole_moment",
            static_cast<std::pair<vec3::Vector3, vec3::Vector3> (Coil::*)(
                vec3::CoordVector3, vec3::Vector3
            ) const>(&Coil::computeForceOnDipoleMoment),
            py::arg("point_vector"), py::arg("dipole_moment"))
        .def(
            "compute_force_on_dipole_moment",
            static_cast<std::pair<vec3::Vector3, vec3::Vector3> (Coil::*)(
                    vec3::CoordVector3, vec3::Vector3, const PrecisionArguments&
            ) const>(&Coil::computeForceOnDipoleMoment),
            py::arg("point_vector"), py::arg("dipole_moment"), py::arg("used_precision"));

    coil.def_static(
            "compute_all_mutual_inductance_arrangements", &Coil::computeAllMutualInductanceArrangements,
            py::arg("primary"), py::arg("secondary"), py::arg("primary_positions"), py::arg("secondary_positions"),
            py::arg("primary_Y_angles"), py::arg("primary_Z_angles"),
            py::arg("secondary_Y_angles"), py::arg("secondary_Z_angles"),
            py::arg("precision_factor") = PrecisionFactor(), py::arg("compute_method") = CPU_ST)
        .def_static(
            "compute_all_ampere_force_arrangements", &Coil::computeAllAmpereForceArrangements,
            py::arg("primary"), py::arg("secondary"), py::arg("primary_positions"), py::arg("secondary_positions"),
            py::arg("primary_Y_angles"), py::arg("primary_Z_angles"),
            py::arg("secondary_Y_angles"), py::arg("secondary_Z_angles"),
            py::arg("precision_factor") = PrecisionFactor(), py::arg("compute_method") = CPU_ST);

    coil.def("__repr__", &Coil::operator std::string);
}