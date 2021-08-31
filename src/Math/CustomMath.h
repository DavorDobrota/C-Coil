#ifndef GENERAL_COIL_PROGRAM_CUSTOMMATH_H
#define GENERAL_COIL_PROGRAM_CUSTOMMATH_H

namespace customMath
{
    const extern double taylorWeightsLn[30];
    const extern double taylorWeightsCos[9];

    double ln(double x);
    double cos(double x);
}

#endif //GENERAL_COIL_PROGRAM_CUSTOMMATH_H
