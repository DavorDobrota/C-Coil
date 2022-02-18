#ifndef GENERAL_COIL_PROGRAM_PYMAIN_H
#define GENERAL_COIL_PROGRAM_PYMAIN_H

#include<pybind11/pybind11.h>
namespace py = pybind11;

void initTensor(py::module_ &mainModule);

#endif //GENERAL_COIL_PROGRAM_PYMAIN_H
