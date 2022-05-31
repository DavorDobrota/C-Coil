#ifndef COPPER_HARDWARE_ACCELERATED_FUNCTIONS
#define COPPER_HARDWARE_ACCELERATED_FUNCTIONS

#include "constants.h"


void Calculate_hardware_accelerated_a
(
    long long num_ops,
    const TYPE *z_ar,
    const TYPE *r_ar,
    TYPE current_density,
    TYPE innerRadius,
    TYPE length,
    TYPE thickness,
    int lengthIncrements,
    int thicknessIncrements,
    int angularIncrements,
    TYPE *potentialArray = nullptr
);

void Calculate_hardware_accelerated_b
(
    long long num_ops,
    const TYPE *z_ar,
    const TYPE *r_ar,
    TYPE current_density,
    TYPE innerRadius,
    TYPE length,
    TYPE thickness,
    int lengthIncrements,
    int thicknessIncrements,
    int angularIncrements,
    TYPE *fieldHArray = nullptr,
    TYPE *fieldZArray = nullptr
);

void Calculate_hardware_accelerated_g
(
    long long num_ops,
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
    TYPE *gradientZZArr = nullptr
);

#endif // COPPER_HARDWARE_ACCELERATED_FUNCTIONS
