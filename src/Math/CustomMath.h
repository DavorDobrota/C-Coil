#ifndef GENERAL_COIL_PROGRAM_CUSTOMMATH_H
#define GENERAL_COIL_PROGRAM_CUSTOMMATH_H

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) || defined(WIN64) || defined(_WIN64) || defined(__WIN64)
    #define LN customMath::ln
    #define COS customMath::cos

    #define LNF customMath::lnf
    #define COSF customMath::cosf
#else
    #define LN std::log
    #define COS std::cos

    #define LNF std::log
    #define COSF std::cos
#endif

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
