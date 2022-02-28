#include <Coil.h>
#include <CoilData.h>
#include <PyMain.h>

#include <pybind11/pybind11.h>
namespace py = pybind11;


void initCoilData(py::module_ &mainModule)
{
    py::module_ coilDataModule = mainModule.def_submodule("coil_data");

    py::class_<CoilData> coilData(coilDataModule, "CoilData");


    // CoilData

    coilData.def(py::init<const Coil&>());

    coilData.def_readwrite("num_of_turns", &CoilData::numOfTurns)
            .def_readwrite("current_density", &CoilData::currentDensity)
            .def_readwrite("inner_radius", &CoilData::innerRadius)
            .def_readwrite("thickness", &CoilData::thickness)
            .def_readwrite("length", &CoilData::length)
            .def_readwrite("angular_iterations", &CoilData::angularIterations)
            .def_readwrite("length_iterations", &CoilData::lengthIterations)
            .def_readwrite("thickness_increments", &CoilData::thicknessIncrements)
            .def_readonly("position_array", &CoilData::positionArray)
            .def_readonly("weight_array", &CoilData::weightArray);
}
