#ifndef COPPER_HARDWARE_ACCELERATED_FUNCTIONS
#define COPPER_HARDWARE_ACCELERATED_FUNCTIONS

#include "CUDAUtils/ConstantsAndStructs/CUDAConstants.h"
#include "CUDAUtils/ConstantsAndStructs/CoilDataStructs.h"


void Calculate_hardware_accelerated_a(long long opCount, CoilData coil,
                                      const VectorData *posArr,
                                      VectorData *resArr = nullptr);

void Calculate_hardware_accelerated_b(long long opCount, CoilData coil,
                                      const VectorData *posArr,
                                      VectorData *resArr = nullptr);

void Calculate_hardware_accelerated_g(long long opCount, CoilData coil,
                                      const VectorData *posArr,
                                      MatrixData *resArr = nullptr);


void Calculate_mutual_inductance_configurations(long long configCount, long long pointCount,
                                                const CoilPairArgumentsData *coilPair,
                                                const CoilPairPositionData *configArr,
                                                TYPE *inductanceArr = nullptr);

void Calculate_force_and_torque_configurations(long long configCount, long long pointCount,
                                               const CoilPairArgumentsData *coilPair,
                                               const CoilPairPositionData *configArr,
                                               ForceTorqueData *forceTorqueArr = nullptr);

#endif // COPPER_HARDWARE_ACCELERATED_FUNCTIONS
