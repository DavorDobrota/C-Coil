#ifndef COPPER_HARDWARE_ACCELERATED_FUNCTIONS
#define COPPER_HARDWARE_ACCELERATED_FUNCTIONS

#include "CUDAConstants.h"
#include "CoilData.h"


void Calculate_hardware_accelerated_a(long long numOps, CoilData coil,
                                      const DataVector *posArr,
                                      DataVector *resArr = nullptr);

void Calculate_hardware_accelerated_b(long long numOps, CoilData coil,
                                      const DataVector *posArr,
                                      DataVector *resArr = nullptr);

void Calculate_hardware_accelerated_g(long long numOps, CoilData coil,
                                      const DataVector *posArr,
                                      DataMatrix *resArr = nullptr);

void Calculate_hardware_accelerated_a_group(long long numCoils, long long numOps,
                                            const CoilData *coilArr,
                                            const DataVector *posArr,
                                            DataVector *resArr = nullptr);

void Calculate_hardware_accelerated_b_group(long long numCoils, long long numOps,
                                            const CoilData *coilArr,
                                            const DataVector *posArr,
                                            DataVector *resArr = nullptr);

void Calculate_hardware_accelerated_g_group(long long numCoils, long long numOps,
                                            const CoilData *coilArr,
                                            const DataVector *posArr,
                                            DataMatrix *resArr = nullptr);

#endif // COPPER_HARDWARE_ACCELERATED_FUNCTIONS
