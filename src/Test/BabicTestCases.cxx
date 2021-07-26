#include "Test.h"
#include "Coil.h"

#include <cstdio>
#include <cmath>

void testCoilMutualInductanceGeneralThinCoilAndFilament()
{
    Coil primaryGeneral = Coil(0.06, 1e-15, 0.12, 120);
    Coil secondaryGeneral = Coil(0.05, 1e-15, 1e-15, 1);

    auto precision = PrecisionFactor(6.0);

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                              0.00, 0.0,
                                                              acos(i * 0.1), precision));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                              0.03, 0.0,
                                                              acos(i * 0.1), precision));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                              0.06, 0.0,
                                                              acos(i * 0.1), precision));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                              0.12, 0.0,
                                                              acos(i * 0.1), precision));
    }
    printf("\n");
}

void testCoilMutualInductanceGeneralThinCoilAndThinCoil()
{
    Coil primaryGeneral = Coil(0.06, 1e-18, 0.12, 120);
    Coil secondaryGeneral = Coil(0.05, 1e-18, 0.04, 60);

    auto precision = PrecisionFactor(6.0);

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                              0.0, 0.0,
                                                              acos(i * 0.1), precision));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
        //   printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                              0.03, 0.0,
                                                              acos(i * 0.1)), precision);
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                              0.06, 0.0,
                                                              acos(i * 0.1)), precision);
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                              0.12, 0.0,
                                                              acos(i * 0.1), precision));
    }
    printf("\n");
}

void testCoilMutualInductanceGeneralPancakeAndPancake()
{
    Coil primaryGeneral = Coil(0.04, 0.02, 1e-15, 200);
    Coil secondaryGeneral = Coil(0.015, 0.01, 1e-15, 100);

    auto precision = PrecisionFactor(6.0);

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                              0.05, 0.0,
                                                              acos(i * 0.1), precision));
    }
    printf("\n");
}

void testCoilMutualInductanceGeneralRectangularAndFilament()
{
    Coil primaryGeneral = Coil(0.04, 0.02, 0.01, 100);
    Coil secondaryGeneral = Coil(0.02, 1e-15, 1e-15, 1);

    auto precision = PrecisionFactor(6.0);

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primaryGeneral, secondaryGeneral,
                                                              0.0, 0.0,
                                                              acos(i * 0.1), precision));
    }
    printf("\n");
}

