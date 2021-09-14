#ifndef GENERAL_COIL_PROGRAM_CUSTOMMATH_H
#define GENERAL_COIL_PROGRAM_CUSTOMMATH_H

namespace customMath
{
    const extern double taylorWeightsLn[30];
    const extern double taylorWeightsCos[9];

    const extern double taylorTableLn[64][8];

    double ln(double x);
    double cos(double x);

    float lnf(float x);
    float cosf(float x);  // TODO: implement
}

#endif //GENERAL_COIL_PROGRAM_CUSTOMMATH_H
