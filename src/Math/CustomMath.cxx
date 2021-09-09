#include "CustomMath.h"

#include <cmath>
#include <cstdint>
#include <cstring>

#include <immintrin.h>

const double
customMath::taylorWeightsLn[] = {1.09861228866810978211, 0.33333333333333331483, -0.055555555555555552472,
                                 0.012345679012345678327, -0.0030864197530864195818, 0.00082304526748971191738,
                                 -0.00022862368541380886293, 6.53210529753739589e-05, -1.905197378448407304e-05,
                                 5.6450292694767621626e-06, -1.6935087808430286488e-06, 5.1318447904334197628e-07,
                                 -1.5680636859657673262e-07, 4.8248113414331303872e-08, -1.4933939866340640099e-08,
                                 4.6461146250837549549e-09, -1.45191082033867332e-09, 4.5550143383174067577e-10,
                                 -1.4339859953962206747e-10, 4.5283768275670128609e-11, -1.433985995396220707e-11,
                                 4.5523364933213352186e-12, -1.4484707024204249708e-12, 4.6183123845288911381e-13,
                                 -1.4752942339467290434e-13, 4.72094154862953304e-14, -1.5131222912274144035e-14,
                                 4.8569357496188614282e-15, -1.5611579195203482035e-15, 5.0244162927091664506e-16};

const double
customMath::taylorWeightsCos[] = {0.16666666666666665741, -0.0083333333333333332177, 0.00019841269841269841253,
                                  -2.7557319223985892511e-06, 2.5052108385441720224e-08, -1.6059043836821613341e-10,
                                  7.6471637318198164055e-13, -2.8114572543455205981e-15, 8.2206352466243294955e-18};

double customMath::ln(double x)
{
    uint64_t xBits;
    double x1;
    int exponent;
    // basic principle, ln(a) = ln(b * 2^n) = ln(b) + n * ln(2), b element [2.0, 4.0], a is positive
    // converting double to an int so bitwise operations can be performed
    std::memcpy(&xBits, &x, sizeof(x));
    // extracting the exponent bits and manipulating them to their final form
    unsigned long long exp = xBits & 0x7ff0000000000000;
    exp >>= 52;
    exp -= 1024;
    exponent = (int) exp;
    // changing the exponent so all numbers land in [2.0, 4.0]
    xBits &= 0x400fffffffffffff;
    xBits |= 0x4000000000000000;
    // converting the altered number back to double
    std::memcpy(&x1, &xBits, sizeof(xBits));

    double output = exponent * 0.69314718055994528623;
    x1 -= 3;

    double x2 = x1 * x1;
    double x3 = x1 * x2;
    double x4 = x2 * x2;
    double x5 = x2 * x3;
    double x6 = x2 * x4;
    double x7 = x1 * x6;
    double x8 = x4 * x4;
    double x10 = x5 * x5;

    __m256d inputVector = _mm256_set_pd(x4, x5, x6, x10);
    __m256d outputVector = _mm256_mul_pd(inputVector, inputVector);
    
    double outAr[4];
    _mm256_storeu_pd(outAr, outputVector);

    double x12 = outAr[0];
    double x16 = outAr[1];
    double x20 = outAr[2];
    double x24 = outAr[3];

    output += taylorWeightsLn[0] + taylorWeightsLn[1] * x1;
    output += taylorWeightsLn[2] * x2 + taylorWeightsLn[3] * x3;
    output += taylorWeightsLn[4] * x4 + taylorWeightsLn[5] * x5;
    output += taylorWeightsLn[6] * x6 + taylorWeightsLn[7] * x7;
    output += taylorWeightsLn[8] * x8 + taylorWeightsLn[9] * x1 * x8;
    output += taylorWeightsLn[10] * x10 + taylorWeightsLn[11] * x1 * x10;
    output += taylorWeightsLn[12] * x12 + + taylorWeightsLn[13] * x1 * x12;
    output += taylorWeightsLn[14] * x2 * x12 + taylorWeightsLn[15] * x3 * x12;
    output += taylorWeightsLn[16] * x16 + taylorWeightsLn[17] * x1 * x16;
    output += taylorWeightsLn[18] * x2 * x16 + taylorWeightsLn[19] * x3 * x16;
    output += taylorWeightsLn[20] * x20 + taylorWeightsLn[21] * x1 * x20;
    output += taylorWeightsLn[22] * x2 * x20 + taylorWeightsLn[23] * x3 * x20;
    output += taylorWeightsLn[24] * x24 + taylorWeightsLn[25] * x1 * x24;
    output += taylorWeightsLn[26] * x2 * x24 + taylorWeightsLn[27] * x3 * x24;
    output += taylorWeightsLn[28] * x4 * x24 + taylorWeightsLn[29] * x5 * x24;

    return output;
}

double customMath::cos(double x)
{
    double x1 = x - M_PI_2;
    double output = 0.0;

    double x2 = x1 * x1;
    double x3 = x1 * x2;
    double x5 = x3 * x2;
    double x7 = x5 * x2;
    double x9 = x7 * x2;
    double x11 = x9 * x2;
    double x13 = x11 * x2;
    double x15 = x13 * x2;
    double x17 = x15 * x2;

    output += -x1 + taylorWeightsCos[0] * x3;
    output += taylorWeightsCos[1] * x5 + taylorWeightsCos[2] * x7;
    output += taylorWeightsCos[3] * x9 + taylorWeightsCos[4] * x11;
    output += taylorWeightsCos[5] * x13 + taylorWeightsCos[6] * x15;
    output += taylorWeightsCos[7] * x17 + taylorWeightsCos[8] * x17 * x2;

    return output;
}

float customMath::lnf(float x)
{
    auto xd = (double) x;
    uint64_t xBits;
    double x1;
    int exponent;
    // basic principle, ln(a) = ln(b * 2^n) = ln(b) + n * ln(2), b element [2.0, 4.0], a is positive
    // converting double to an int so bitwise operations can be performed
    std::memcpy(&xBits, &xd, sizeof(xd));
    // extracting the exponent bits and manipulating them to their final form
    unsigned long long exp = xBits & 0x7ff0000000000000;
    exp >>= 52;
    exp -= 1024;
    exponent = (int) exp;
    // changing the exponent so all numbers land in [2.0, 4.0]
    xBits &= 0x400fffffffffffff;
    xBits |= 0x4000000000000000;
    // converting the altered number back to double
    std::memcpy(&x1, &xBits, sizeof(xBits));

    double output = exponent * 0.69314718055994528623;
    x1 -= 3;

    double x2 = x1 * x1;
    double x4 = x2 * x2;
    double x6 = x2 * x4;
    double x8 = x4 * x4;
    double x10 = x6 * x4;
    double x12 = x6 * x6;

    output += taylorWeightsLn[0] + taylorWeightsLn[1] * x1;
    output += taylorWeightsLn[2] * x2 + taylorWeightsLn[3] * x2 * x1;
    output += taylorWeightsLn[4] * x4 + taylorWeightsLn[5] * x4 * x1;
    output += taylorWeightsLn[6] * x6 + taylorWeightsLn[7] * x6 * x1;
    output += taylorWeightsLn[8] * x8 + taylorWeightsLn[9] * x1 * x8;
    output += taylorWeightsLn[10] * x10 + taylorWeightsLn[11] * x1 * x10;
    output += taylorWeightsLn[12] * x12 + taylorWeightsLn[13] * x1 * x12;

    return (float) output;
}