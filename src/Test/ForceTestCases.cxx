#include "Test.h"
#include "Coil.h"

#include <cstdio>
#include <cmath>

void testAmpereForceForFilamentsZAxis()
{
    Coil coil1 = Coil(0.5, 0.0, 0.0, 1, 100);
    Coil coil2 = Coil(0.3, 0.0, 0.0, 1, 200);
    coil2.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.4));

    for (int i = 1; i <= 8; i++)
    {
        printf("%.16g\n", 1000 * Coil::computeAmpereForce(coil1, coil2, PrecisionFactor(i)).first.zComponent);
    }

    printf("\n");
}

void testAmpereForceGeneralCase()
{
    Coil coil1 = Coil(0.0602, 0.0728, 0.5292, 126, 16500);
    Coil coil2 = Coil(0.168, 0.0285, 0.552, 1890, 725);
    Coil coil3 = Coil(0.1965, 0.04365, 0.552, 3792, 725);

    printf("%.15g T\n",
           coil1.computeBFieldZ(vec3::CoordVector3()) +
           coil2.computeBFieldZ(vec3::CoordVector3()) +
           coil3.computeBFieldZ(vec3::CoordVector3()));

    printf("%.15g MJ\n", 1e-6 *
           (0.5 * coil1.computeAndSetSelfInductance(PrecisionFactor(12)) * coil1.getCurrent() * coil1.getCurrent() +
           0.5 * coil2.computeAndSetSelfInductance(PrecisionFactor(12)) * coil2.getCurrent() * coil2.getCurrent() +
           0.5 * coil3.computeAndSetSelfInductance(PrecisionFactor(12)) * coil3.getCurrent() * coil3.getCurrent()));

    std::pair<vec3::FieldVector3, vec3::FieldVector3> forcePair1, forcePair2;
    auto precision = PrecisionFactor(9.0);

    printf("Force and Torque in displacement\n");
    for (double dz = 0.0; dz <= 0.004; dz += 0.001)
    {
        for (double dr = 0.0; dr <= 0.004; dr += 0.002)
        {
            coil1.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, dr, 0.0, dz));
            forcePair1 = Coil::computeAmpereForce(coil3, coil1, precision, CPU_ST);
            forcePair2 = Coil::computeAmpereForce(coil2, coil1, precision, CPU_ST);
            printf("%16.10g %16.10g %16.10g %16.10g %16.10g %16.10g\n",
                   forcePair1.first.xComponent + forcePair2.first.xComponent,
                   forcePair1.first.yComponent + forcePair2.first.yComponent,
                   forcePair1.first.zComponent + forcePair2.first.zComponent,
                   forcePair1.second.xComponent + forcePair2.second.xComponent,
                   forcePair1.second.yComponent + forcePair2.second.yComponent,
                   forcePair1.second.zComponent + forcePair2.second.zComponent);
        }
    }
    printf("Force and Torque in rotation\n");
    for (int i = 0; i <= 10; ++i)
    {
        coil1.setPositionAndOrientation(vec3::CoordVector3(), M_PI/360 * i);
        forcePair1 = Coil::computeAmpereForce(coil3, coil1, precision, CPU_ST);
        forcePair2 = Coil::computeAmpereForce(coil2, coil1, precision, CPU_ST);
        printf("%16.10g %16.10g %16.10g %16.10g %16.10g %16.10g\n",
               forcePair1.first.xComponent + forcePair2.first.xComponent,
               forcePair1.first.yComponent + forcePair2.first.yComponent,
               forcePair1.first.zComponent + forcePair2.first.zComponent,
               forcePair1.second.xComponent + forcePair2.second.xComponent,
               forcePair1.second.yComponent + forcePair2.second.yComponent,
               forcePair1.second.zComponent + forcePair2.second.zComponent);
    }
    printf("\n");

    coil1.setPositionAndOrientation(vec3::CoordVector3(), M_PI/36 + 1e-7);
    double M1 = coil1.getCurrent() * coil2.getCurrent() * Coil::computeMutualInductance(coil1, coil2, precision) +
                coil1.getCurrent() * coil3.getCurrent() * Coil::computeMutualInductance(coil1, coil3, precision);
    coil1.setPositionAndOrientation(vec3::CoordVector3(), M_PI/36 - 1e-7);
    double M2 = coil1.getCurrent() * coil2.getCurrent() * Coil::computeMutualInductance(coil1, coil2, precision) +
            coil1.getCurrent() * coil3.getCurrent() * Coil::computeMutualInductance(coil1, coil3, precision);
    printf("By mutual inductance gradient : %.10g", (M1 - M2) / 2e-7);
}

void testAmpereForceThinCoils()
{
    ComputeMethod method = CPU_ST;
    std::pair<vec3::FieldVector3, vec3::FieldVector3> tempForce;

    Coil coil1 = Coil(0.5, 0.0, 0.3, 200);
    Coil coil2 = Coil(0.5, 0.0, 0.2, 100);
    coil2.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.4));
    for (int i = 1; i <= 10; ++i)
    {
        tempForce = Coil::computeAmpereForce(coil1, coil2, PrecisionFactor(i), method);
        printf("%.16g mN\n", 1000 * tempForce.first.zComponent);
    }
    printf("\n");

    Coil coil3 = Coil(0.1, 0.0, 0.15, 225);
    Coil coil4 = Coil(0.09, 0.0, 0.15, 300);
    coil4.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.05));
    for (int i = 1; i <= 12; ++i)
    {
        tempForce = Coil::computeAmpereForce(coil3, coil4, PrecisionFactor(i), method);
        printf("%.16g mN\n", 1000 * tempForce.first.zComponent);
    }
    printf("\n");

    Coil coil5 = Coil(0.25, 0.0, 0.06, 40);
    Coil coil6 = Coil(0.2, 0.0, 0.04, 120);
    coil6.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.3));
    for (int i = 1; i <= 10; ++i)
    {
        tempForce = Coil::computeAmpereForce(coil5, coil6, PrecisionFactor(i), method);
        printf("%.16g mN\n", 1000 * tempForce.first.zComponent);
    }
    printf("\n");

    Coil coil7 = Coil(0.0762, 0.0832, 0.0, 516, 1.42);
    Coil coil8 = Coil(0.0762, 0.0832, 0.0, 516, 1.42);
    coil8.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.0468));
    for (int i = 1; i <= 10; ++i)
    {
        tempForce = Coil::computeAmpereForce(coil7, coil8, PrecisionFactor(i), method);
        printf("%.16g N\n", tempForce.first.zComponent);
    }
    printf("\n");

    Coil coil9 = Coil(0.16, 0.12, 0.0, 100, 10);
    Coil coil10 = Coil(0.11, 0.15, 0.0, 100, 10);
    for (int i = 1; i <= 10; ++i)
    {
        coil10.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.05));
        tempForce = Coil::computeAmpereForce(coil9, coil10, PrecisionFactor(i), method);
        printf("%.16g N\n", tempForce.first.zComponent);
    }
    printf("\n");

    Coil coil11 = Coil(0.12, 0.11, 0.0, 100, 10);
    Coil coil12 = Coil(0.12, 0.11, 0.0, 100, 10);
    for (int i = 1; i <= 10; ++i)
    {
        coil12.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.02));
        tempForce = Coil::computeAmpereForce(coil11, coil12, PrecisionFactor(i), method);
        printf("%.16g N\n", tempForce.first.zComponent);
    }
    printf("\n");

    Coil coil13 = Coil(0.1, 0.0, 0.2, 100);
    Coil coil14 = Coil(0.2, 0.2, 0.0, 100);
    coil14.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.6));
    for (int i = 1; i <= 10; ++i)
    {
        tempForce = Coil::computeAmpereForce(coil13, coil14, PrecisionFactor(i), method);
        printf("%.16g mN\n", 1000 * tempForce.first.zComponent);
    }
    printf("\n");
}

void testAmpereForceFilamentsGeneral()
{
    std::pair<vec3::FieldVector3, vec3::FieldVector3> tempForce;
    double a, b, c;
    printf("Singular cases (different precisions):\n\n");

    Coil coil1 = Coil(0.2, 0.0, 0.0, 1);
    Coil coil2 = Coil(0.1, 0.0, 0.0, 1);
    a = 1.0; b = 1.0; c = 1.0;
    coil2.setPositionAndOrientation(
            vec3::CoordVector3(vec3::CARTESIAN, 0.1, 0.1, 0.1),
            std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    printf("%g %g\n", std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    for (int i = 1; i <= 9; ++i)
    {
        tempForce = Coil::computeAmpereForce(coil1, coil2, PrecisionFactor(i));
        printf("(%20.16g, %20.16g, %20.16g) microN\n",
               1e6 * tempForce.first.xComponent, 1e6 * tempForce.first.yComponent, 1e6 * tempForce.first.zComponent);
    }
    printf("\n");

    Coil coil3 = Coil(0.4, 0.0, 0.0, 1);
    Coil coil4 = Coil(0.05, 0.0, 0.0, 1);
    a = 3.0; b = 2.0; c = 1.0;
    coil4.setPositionAndOrientation(
            vec3::CoordVector3(vec3::CARTESIAN, 0.1, 0.15, 0.0),
            std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    printf("%g %g\n", std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    for (int i = 1; i <= 9; ++i)
    {
        tempForce = Coil::computeAmpereForce(coil3, coil4, PrecisionFactor(i));
        printf("(%20.16g, %20.16g, %20.16g) nanoN\n",
               1e9 * tempForce.first.xComponent, 1e9 * tempForce.first.yComponent, 1e9 * tempForce.first.zComponent);
    }
    printf("\n");

    Coil coil5 = Coil(0.9, 0.0, 0.0, 1);
    Coil coil6 = Coil(0.6, 0.0, 0.0, 1);
    a = 1.0; b = 1.0; c = 1.0;
    coil6.setPositionAndOrientation(
            vec3::CoordVector3(vec3::CARTESIAN, 0.3, 0.2, 0.5),
            std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    printf("%g %g\n", std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    for (int i = 1; i <= 9; ++i)
    {
        tempForce = Coil::computeAmpereForce(coil5, coil6, PrecisionFactor(i));
        printf("(%20.16g, %20.16g, %20.16g) microN\n",
               1e6 * tempForce.first.xComponent, 1e6 * tempForce.first.yComponent, 1e6 * tempForce.first.zComponent);
    }
    printf("\n");

    Coil coil7 = Coil(0.005, 0.0, 0.0, 1);
    Coil coil8 = Coil(0.001, 0.0, 0.0, 1);
    a = 3.0; b = 1.0; c = 2.0;
    coil8.setPositionAndOrientation(
            vec3::CoordVector3(vec3::CARTESIAN, 0.003, 0.001, 0.0005),
            std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    printf("%g %g\n", std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    for (int i = 1; i <= 9; ++i)
    {
        tempForce = Coil::computeAmpereForce(coil7, coil8, PrecisionFactor(i));
        printf("(%20.16g, %20.16g, %20.16g) microN\n",
               1e6 * tempForce.first.xComponent, 1e6 * tempForce.first.yComponent, 1e6 * tempForce.first.zComponent);
    }
    printf("\n");

    Coil coil9 = Coil(0.3, 0.0, 0.0, 1, 1.0);
    Coil coil10 = Coil(0.3, 0.0, 0.0, 1, -1.0);
    a = 1.0; b = -2.0; c = 1.0;
    coil10.setPositionAndOrientation(
            vec3::CoordVector3(vec3::CARTESIAN, 0.1, -0.3, 0.2),
            std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    printf("%g %g\n", std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    for (int i = 1; i <= 9; ++i)
    {
        tempForce = Coil::computeAmpereForce(coil9, coil10, PrecisionFactor(i));
        printf("(%20.16g, %20.16g, %20.16g) microN\n",
               1e6 * tempForce.first.xComponent, 1e6 * tempForce.first.yComponent, 1e6 * tempForce.first.zComponent);
    }
    printf("\n");

    printf("Varying distance cases (fixed precision):\n\n");

    Coil coil11 = Coil(0.0425, 0.0, 0.0, 1);
    Coil coil12 = Coil(0.02, 0.0, 0.0, 1);
    for (int i = 0; i <= 11; ++i)
    {
        coil12.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, 0.003, 0.0, 0.001 * i));
        tempForce = Coil::computeAmpereForce(coil11, coil12, PrecisionFactor(8.0));
        printf("%20.16g\t %20.16g\n",
               1e7 * tempForce.first.xComponent, 1e7 * tempForce.first.zComponent);
    }
    printf("\n");

    double rArr[] = {0.0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.4999, 0.5, 0.5001, 0.6, 0.7, 0.8, 0.9, 1.0,
                     1.1, 1.2, 1.3, 1.4, 1.4999, 1.5, 1.5001, 1.6, 1.8, 1.9, 2.0, 3.0};
    Coil coil13 = Coil(1.0, 0.0, 0.0, 1);
    Coil coil14 = Coil(0.5, 0.0, 0.0, 1);
    for (double i : rArr)
    {
        coil14.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CARTESIAN, i, 0.0, 1e-18));
        tempForce = Coil::computeAmpereForce(coil13, coil14, PrecisionFactor(10.0));
        printf("%6g %20.16g\t %20.16g\n",i, 1e6 * tempForce.first.xComponent, 1e6 * tempForce.first.zComponent);
    }
    printf("\n");

    printf("Some more singular cases (different precisions):\n\n");

    Coil coil15 = Coil(1.0, 0.0, 0.0, 1);
    Coil coil16 = Coil(0.5, 0.0, 0.0, 1);
    coil16.setPositionAndOrientation(
            vec3::CoordVector3(vec3::CARTESIAN, 2.0, 2.0, 2.0));
    for (int i = 1; i <= 9; ++i)
    {
        tempForce = Coil::computeAmpereForce(coil15, coil16, PrecisionFactor(i));
        printf("(%20.16g, %20.16g, %20.16g) nanoN\n",
               1e9 * tempForce.first.xComponent, 1e9 * tempForce.first.yComponent, 1e9 * tempForce.first.zComponent);
    }
    printf("\n");

    Coil coil17 = Coil(1.0, 0.0, 0.0, 1);
    Coil coil18 = Coil(0.5, 0.0, 0.0, 1);
    coil18.setPositionAndOrientation(
            vec3::CoordVector3(vec3::CARTESIAN, 1.0, 2.0, 3.0), M_PI_2);
    for (int i = 1; i <= 9; ++i)
    {
        tempForce = Coil::computeAmpereForce(coil17, coil18, PrecisionFactor(i));
        printf("(%20.16g, %20.16g, %20.16g) nanoN\n",
               1e9 * tempForce.first.xComponent, 1e9 * tempForce.first.yComponent, 1e9 * tempForce.first.zComponent);
    }
    printf("\n");

    Coil coil19 = Coil(1.0, 0.0, 0.0, 1);
    Coil coil20 = Coil(0.5, 0.0, 0.0, 1);
    coil20.setPositionAndOrientation(
            vec3::CoordVector3(vec3::CARTESIAN, 2.0, 2.0, 2.0), M_PI_2, M_PI_2);
    for (int i = 1; i <= 9; ++i)
    {
        tempForce = Coil::computeAmpereForce(coil19, coil20, PrecisionFactor(i));
        printf("(%20.16g, %20.16g, %20.16g) nanoN\n",
               1e9 * tempForce.first.xComponent, 1e9 * tempForce.first.yComponent, 1e9 * tempForce.first.zComponent);
    }
    printf("\n");
}