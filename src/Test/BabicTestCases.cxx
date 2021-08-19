#include "Test.h"
#include "Coil.h"

#include <cstdio>
#include <cmath>

void testCoilMutualInductanceGeneralThinCoilAndFilament()
{
    Coil primary = Coil(0.06, 0.0, 0.12, 120);
    Coil secondary = Coil(0.05, 0.0, 0.0, 1);

    auto precision = PrecisionFactor(7.0);

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(vec3::CoordVector3(), std::acos(i * 0.1));
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primary, secondary, precision));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.03), std::acos(i * 0.1));
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primary, secondary, precision));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.06), std::acos(i * 0.1));
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primary, secondary, precision));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.12), std::acos(i * 0.1));
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primary, secondary, precision));
    }
    printf("\n");
}

void testCoilMutualInductanceGeneralThinCoilAndThinCoil()
{
    Coil primary = Coil(0.06, 0.0, 0.12, 120);
    Coil secondary = Coil(0.05, 0.0, 0.04, 60);

    auto precision = PrecisionFactor(7.0);

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(vec3::CoordVector3(), std::acos(i * 0.1));
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primary, secondary, precision));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
        //   printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.03), std::acos(i * 0.1));
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primary, secondary, precision));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.06), std::acos(i * 0.1));
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primary, secondary, precision));
    }
    printf("\n");

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.12), std::acos(i * 0.1));
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primary, secondary, precision));
    }
    printf("\n");
}

void testCoilMutualInductanceGeneralPancakeAndPancake()
{
    Coil primary = Coil(0.04, 0.02, 0.0, 200);
    Coil secondary = Coil(0.015, 0.01, 0.0, 100);

    auto precision = PrecisionFactor(6.0);

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.05), std::acos(i * 0.1));
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primary, secondary, precision));
    }
    printf("\n");
}

void testCoilMutualInductanceGeneralRectangularAndFilament()
{
    Coil primary = Coil(0.04, 0.02, 0.01, 100);
    Coil secondary = Coil(0.02, 0.0, 0.0, 1);

    auto precision = PrecisionFactor(6.0);

    for (int i = 10; i >= 0; --i){
        //    printf("cos(alpha) = %.1f: ", i * 0.1);
        secondary.setPositionAndOrientation(vec3::CoordVector3(), std::acos(i * 0.1));
        printf("%.15f\n", 1e6 * Coil::computeMutualInductance(primary, secondary, precision));
    }
    printf("\n");
}

