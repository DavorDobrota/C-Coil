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

void Calculate_hardware_accelerated_g(long long num_ops,
                                      const TYPE *z_ar,
                                      const TYPE *r_ar,
                                      TYPE current_density,
                                      TYPE innerRadius,
                                      TYPE length,
                                      TYPE thickness,
                                      int lengthIncrements,
                                      int thicknessIncrements,
                                      int angularIncrements,
                                      TYPE *gradientRPArr = nullptr,
                                      TYPE *gradientRRArr = nullptr,
                                      TYPE *gradientRZArr = nullptr,
                                      TYPE *gradientZZArr = nullptr);

#endif // COPPER_HARDWARE_ACCELERATED_FUNCTIONS
