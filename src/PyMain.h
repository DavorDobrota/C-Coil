#ifndef GENERAL_COIL_PROGRAM_PYMAIN_H
#define GENERAL_COIL_PROGRAM_PYMAIN_H

#include <pybind11/pybind11.h>

void initTensor(pybind11::module_ &mainModule);
void initCoil(pybind11::module_ &mainModule);
void initCoilData(pybind11::module_ &mainModule);
void initCoilGroup(pybind11::module_ &mainModule);
void initBenchmark(pybind11::module_ &mainModule);
void initCompare(pybind11::module_ &mainModule);

#endif //GENERAL_COIL_PROGRAM_PYMAIN_H
