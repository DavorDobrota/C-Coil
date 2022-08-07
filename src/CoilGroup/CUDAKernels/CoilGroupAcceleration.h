#ifndef COIL_EVOLUTION_COILGROUPACCELERATION_H
#define COIL_EVOLUTION_COILGROUPACCELERATION_H

#include "CUDAUtils/ConstantsAndStructs/CUDAConstants.h"
#include "CUDAUtils/ConstantsAndStructs/CoilDataStructs.h"

void Calculate_hardware_accelerated_a_group(long long coilCount, long long opCount,
                                            const CoilData *coilArr,
                                            const DataVector *posArr,
                                            DataVector *resArr = nullptr);

void Calculate_hardware_accelerated_b_group(long long coilCount, long long opCount,
                                            const CoilData *coilArr,
                                            const DataVector *posArr,
                                            DataVector *resArr = nullptr);

void Calculate_hardware_accelerated_e_group(long long coilCount, long long opCount,
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

void Calculate_force_and_torque_configurations_group(long long coilCount, long long configCount, long long pointCount,
                                                     SecondaryCoilData secondaryCoil,
                                                     const CoilData *coils,
                                                     const SecondaryCoilPositionData *secondaryPositions,
                                                     ForceTorqueData *forceTorqueArr = nullptr);

#endif //COIL_EVOLUTION_COILGROUPACCELERATION_H
