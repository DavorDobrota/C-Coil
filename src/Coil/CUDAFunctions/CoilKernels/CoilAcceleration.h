#ifndef COPPER_HARDWARE_ACCELERATED_FUNCTIONS
#define COPPER_HARDWARE_ACCELERATED_FUNCTIONS

#include "CUDAFunctions/ConstantsAndStructs/CUDAConstants.h"
#include "CUDAFunctions/ConstantsAndStructs/CoilDataStructs.h"


void Calculate_hardware_accelerated_a(long long opCount, CoilData coil,
                                      const DataVector *posArr,
                                      DataVector *resArr = nullptr);

void Calculate_hardware_accelerated_b(long long opCount, CoilData coil,
                                      const DataVector *posArr,
                                      DataVector *resArr = nullptr);

void Calculate_hardware_accelerated_g(long long opCount, CoilData coil,
                                      const DataVector *posArr,
                                      DataMatrix *resArr = nullptr);


void Calculate_mutual_inductance_configurations(long long configCount, long long pointCount,
                                                CoilPairArgumentsData coilPair,
                                                const CoilPairPositionData *configArr,
                                                TYPE *inductanceArr = nullptr);

void Calculate_force_and_torque_configurations(long long configCount, long long pointCount,
                                               CoilPairArgumentsData coilPair,
                                               const CoilPairPositionData *configArr,
                                               ForceTorqueData *forceTorqueArr = nullptr);

#endif // COPPER_HARDWARE_ACCELERATED_FUNCTIONS
