#ifndef COIL_EVOLUTION_COILGROUPACCELERATION_H
#define COIL_EVOLUTION_COILGROUPACCELERATION_H

#include "CUDAFunctions/ConstantsAndStructs/CUDAConstants.h"
#include "CUDAFunctions/ConstantsAndStructs/CoilDataStructs.h"

void Calculate_hardware_accelerated_a_group(long long coilCount, long long opCount,
                                            const CoilData *coilArr,
                                            const DataVector *posArr,
                                            DataVector *resArr = nullptr);

void Calculate_hardware_accelerated_b_group(long long coilCount, long long opCount,
                                            const CoilData *coilArr,
                                            const DataVector *posArr,
                                            DataVector *resArr = nullptr);

void Calculate_hardware_accelerated_g_group(long long coilCount, long long opCount,
                                            const CoilData *coilArr,
                                            const DataVector *posArr,
                                            DataMatrix *resArr = nullptr);

void Calculate_mutual_inductance_configurations_group(long long coilCount, long long configCount, long long pointCount,
                                                      SecondaryCoilData secondaryCoil,
                                                      const CoilData *coils,
                                                      const SecondaryCoilPositionData *secondaryPositions,
                                                      TYPE *inductanceArr = nullptr);

#endif //COIL_EVOLUTION_COILGROUPACCELERATION_H
