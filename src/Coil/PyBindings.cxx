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
        py::arg("method") = CPU_ST, py::arg("z_axis_case") = true);

    coilPairArguments.def("__repr__", &CoilPairArguments::operator std::string);


    // Coil

    // TODO: Finish Coil
}