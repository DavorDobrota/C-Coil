#include "Compare.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>


void Compare::forceTorqueFilamentsZAxis()
{
    Coil coil1 = Coil(0.5, 0.0, 0.0, 1, 100);
    Coil coil2 = Coil(0.3, 0.0, 0.0, 1, 200);
    coil2.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.4));

    for (int i = 1; i <= 8; i++)
    {
        printf("%.16g\n", 1000 * Coil::computeForceTorque(coil1, coil2, PrecisionFactor(i)).first.z);
    }

    printf("\n");
}

void Compare::forceTorqueThickCoilsGeneral()
{
    CoilGroup group = CoilGroup();
    group.addCoil(0.0602, 0.0728, 0.5292, 126, 16500);
    group.addCoil(0.168, 0.0285, 0.552, 1890, 725);
    group.addCoil(0.1965, 0.04365, 0.552, 3792, 725);

    printf("%.15g T\n", group.computeBFieldVector(vec3::Vector3()).z);

    printf("%.15g MJ\n", 1e-6 * (
           0.5 * group[0].computeAndSetSelfInductance(PrecisionFactor(12.0)) *
           group[0].getCurrent() * group[0].getCurrent() +
           0.5 * group[1].computeAndSetSelfInductance(PrecisionFactor(12.0)) *
           group[1].getCurrent() * group[1].getCurrent() +
           0.5 * group[2].computeAndSetSelfInductance(PrecisionFactor(12.0)) *
           group[2].getCurrent() * group[2].getCurrent() +
           Coil::computeMutualInductance(group[0], group[1]) * group[0].getCurrent() * group[1].getCurrent() +
           Coil::computeMutualInductance(group[1], group[2]) * group[1].getCurrent() * group[2].getCurrent() +
           Coil::computeMutualInductance(group[0], group[2]) * group[0].getCurrent() * group[2].getCurrent())
    );

    std::pair<vec3::Vector3, vec3::Vector3> forcePair;
    auto precision = PrecisionFactor(8.0);

    printf("Force and Torque in displacement\n");
    for (int i = 0; i <= 4; ++i)
    {
        for (int j = 0; j <= 2; ++j)
        {
            group[0].setPositionAndOrientation(vec3::Vector3(0.0, 0.002 * j, 0.001 * i));
            forcePair = group.computeForceTorque(group[0], precision, CPU_MT);
            printf("%21.15g %21.15g %21.15g\n%21.15g %21.15g %21.15g\n\n",
                   forcePair.first.x, forcePair.first.y, forcePair.first.z,
                   forcePair.second.x, forcePair.second.y, forcePair.second.z);
        }
    }
    printf("Force and Torque in rotation\n");
    for (int i = 0; i <= 10; ++i)
    {
        group[0].setPositionAndOrientation(vec3::Vector3(), M_PI/360 * i, M_PI_2);
        forcePair = group.computeForceTorque(group[0], precision, CPU_MT);
        printf("%21.15g %21.15g %21.15g\n%21.15g %21.15g %21.15g\n\n",
               forcePair.first.x, forcePair.first.y, forcePair.first.z,
               forcePair.second.x, forcePair.second.y, forcePair.second.z);
    }
    printf("\n");

    group[0].setPositionAndOrientation(vec3::Vector3(), M_PI/36 + 1e-7);
    double M1 = group[0].getCurrent() * group[1].getCurrent() *
            Coil::computeMutualInductance(group[0], group[1], precision, CPU_MT) +
            group[0].getCurrent() * group[1].getCurrent() * Coil::computeMutualInductance(group[0], group[2], precision, CPU_MT);

    group[0].setPositionAndOrientation(vec3::Vector3(), M_PI/36 - 1e-7);
    double M2 = group[0].getCurrent() * group[1].getCurrent() *
            Coil::computeMutualInductance(group[0], group[1], precision, CPU_MT) +
            group[0].getCurrent() * group[2].getCurrent() * Coil::computeMutualInductance(group[0], group[2], precision, CPU_MT);

    printf("By mutual inductance gradient : %.15g\n\n", (M1 - M2) / 2e-7);

    group[0].setPositionAndOrientation(
            vec3::Vector3(0.001, 0.0, 0.0), M_PI/180, 3 * M_PI_2);
    forcePair = group.computeForceTorque(group[0], precision, CPU_MT);
    printf("%21.15g %21.15g %21.15g\n%21.15g %21.15g %21.15g\n\n",
           forcePair.first.x, forcePair.first.y, forcePair.first.z,
           forcePair.second.x, forcePair.second.y, forcePair.second.z);
    group[0].setPositionAndOrientation(
            vec3::Vector3(0.001, 0.0, 0.001), M_PI/180, 3 * M_PI_2);
    forcePair = group.computeForceTorque(group[0], precision, CPU_MT);
    printf("%21.15g %21.15g %21.15g\n%21.15g %21.15g %21.15g\n\n",
           forcePair.first.x, forcePair.first.y, forcePair.first.z,
           forcePair.second.x, forcePair.second.y, forcePair.second.z);

    group[0].setPositionAndOrientation(
            vec3::Vector3(0.001, 0.001, 0.001), M_PI/180, 3 * M_PI_2);
    forcePair = group.computeForceTorque(group[0], precision, CPU_MT);
    printf("%21.15g %21.15g %21.15g\n%21.15g %21.15g %21.15g\n\n",
           forcePair.first.x, forcePair.first.y, forcePair.first.z,
           forcePair.second.x, forcePair.second.y, forcePair.second.z);
}

void Compare::forceTorqueThinCoilsZAxis()
{
    ComputeMethod computeMethod = CPU_ST;
    std::pair<vec3::Vector3, vec3::Vector3> tempForce;

    Coil coil1 = Coil(0.5, 0.0, 0.3, 200);
    Coil coil2 = Coil(0.5, 0.0, 0.2, 100);
    coil2.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.4));
    for (int i = 1; i <= 10; ++i)
    {
        tempForce = Coil::computeForceTorque(coil1, coil2, PrecisionFactor(i), computeMethod);
        printf("%.16g mN\n", 1000 * tempForce.first.z);
    }
    printf("\n");

    Coil coil3 = Coil(0.1, 0.0, 0.15, 225);
    Coil coil4 = Coil(0.09, 0.0, 0.15, 300);
    coil4.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.05));
    for (int i = 1; i <= 12; ++i)
    {
        tempForce = Coil::computeForceTorque(coil3, coil4, PrecisionFactor(i), computeMethod);
        printf("%.16g mN\n", 1000 * tempForce.first.z);
    }
    printf("\n");

    Coil coil5 = Coil(0.25, 0.0, 0.06, 40);
    Coil coil6 = Coil(0.2, 0.0, 0.04, 120);
    coil6.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.3));
    for (int i = 1; i <= 10; ++i)
    {
        tempForce = Coil::computeForceTorque(coil5, coil6, PrecisionFactor(i), computeMethod);
        printf("%.16g mN\n", 1000 * tempForce.first.z);
    }
    printf("\n");

    Coil coil7 = Coil(0.0762, 0.0832, 0.0, 516, 1.42);
    Coil coil8 = Coil(0.0762, 0.0832, 0.0, 516, 1.42);
    coil8.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.0468));
    for (int i = 1; i <= 10; ++i)
    {
        tempForce = Coil::computeForceTorque(coil7, coil8, PrecisionFactor(i), computeMethod);
        printf("%.16g N\n", tempForce.first.z);
    }
    printf("\n");

    Coil coil9 = Coil(0.16, 0.12, 0.0, 100, 10);
    Coil coil10 = Coil(0.11, 0.15, 0.0, 100, 10);
    for (int i = 1; i <= 10; ++i)
    {
        coil10.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.05));
        tempForce = Coil::computeForceTorque(coil9, coil10, PrecisionFactor(i), computeMethod);
        printf("%.16g N\n", tempForce.first.z);
    }
    printf("\n");

    Coil coil11 = Coil(0.12, 0.11, 0.0, 100, 10);
    Coil coil12 = Coil(0.12, 0.11, 0.0, 100, 10);
    for (int i = 1; i <= 10; ++i)
    {
        coil12.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.02));
        tempForce = Coil::computeForceTorque(coil11, coil12, PrecisionFactor(i), computeMethod);
        printf("%.16g N\n", tempForce.first.z);
    }
    printf("\n");

    Coil coil13 = Coil(0.1, 0.0, 0.2, 100);
    Coil coil14 = Coil(0.2, 0.2, 0.0, 100);
    coil14.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.6));
    for (int i = 1; i <= 10; ++i)
    {
        tempForce = Coil::computeForceTorque(coil13, coil14, PrecisionFactor(i), computeMethod);
        printf("%.16g mN\n", 1000 * tempForce.first.z);
    }
    printf("\n");
}

void Compare::forceTorqueFilamentsGeneral()
{
    std::pair<vec3::Vector3, vec3::Vector3> tempForce;
    double a, b, c;
    printf("Singular cases (different precisions):\n\n");

    Coil coil1 = Coil(0.2, 0.0, 0.0, 1);
    Coil coil2 = Coil(0.1, 0.0, 0.0, 1);
    a = 1.0; b = 1.0; c = 1.0;
    coil2.setPositionAndOrientation(
            vec3::Vector3(0.1, 0.1, 0.1),
            std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    printf("%g %g\n", std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    for (int i = 1; i <= 9; ++i)
    {
        tempForce = Coil::computeForceTorque(coil1, coil2, PrecisionFactor(i));
        printf("(%20.16g, %20.16g, %20.16g) microN\n",
               1e6 * tempForce.first.x, 1e6 * tempForce.first.y, 1e6 * tempForce.first.z);
    }
    printf("\n");

    Coil coil3 = Coil(0.4, 0.0, 0.0, 1);
    Coil coil4 = Coil(0.05, 0.0, 0.0, 1);
    a = 3.0; b = 2.0; c = 1.0;
    coil4.setPositionAndOrientation(
            vec3::Vector3(0.1, 0.15, 0.0),
            std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    printf("%g %g\n", std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    for (int i = 1; i <= 9; ++i)
    {
        tempForce = Coil::computeForceTorque(coil3, coil4, PrecisionFactor(i));
        printf("(%20.16g, %20.16g, %20.16g) nanoN\n",
               1e9 * tempForce.first.x, 1e9 * tempForce.first.y, 1e9 * tempForce.first.z);
    }
    printf("\n");

    Coil coil5 = Coil(0.9, 0.0, 0.0, 1);
    Coil coil6 = Coil(0.6, 0.0, 0.0, 1);
    a = 1.0; b = 1.0; c = 1.0;
    coil6.setPositionAndOrientation(
            vec3::Vector3(0.3, 0.2, 0.5),
            std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    printf("%g %g\n", std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    for (int i = 1; i <= 9; ++i)
    {
        tempForce = Coil::computeForceTorque(coil5, coil6, PrecisionFactor(i));
        printf("(%20.16g, %20.16g, %20.16g) microN\n",
               1e6 * tempForce.first.x, 1e6 * tempForce.first.y, 1e6 * tempForce.first.z);
    }
    printf("\n");

    Coil coil7 = Coil(0.005, 0.0, 0.0, 1);
    Coil coil8 = Coil(0.001, 0.0, 0.0, 1);
    a = 3.0; b = 1.0; c = 2.0;
    coil8.setPositionAndOrientation(
            vec3::Vector3(0.003, 0.001, 0.0005),
            std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    printf("%g %g\n", std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    for (int i = 1; i <= 9; ++i)
    {
        tempForce = Coil::computeForceTorque(coil7, coil8, PrecisionFactor(i));
        printf("(%20.16g, %20.16g, %20.16g) microN\n",
               1e6 * tempForce.first.x, 1e6 * tempForce.first.y, 1e6 * tempForce.first.z);
    }
    printf("\n");

    Coil coil9 = Coil(0.3, 0.0, 0.0, 1, 1.0);
    Coil coil10 = Coil(0.3, 0.0, 0.0, 1, -1.0);
    a = 1.0; b = -2.0; c = 1.0;
    coil10.setPositionAndOrientation(
            vec3::Vector3(0.1, -0.3, 0.2),
            std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    printf("%g %g\n", std::atan2(std::sqrt(a*a + b*b), c), std::atan2(b, a));
    for (int i = 1; i <= 9; ++i)
    {
        tempForce = Coil::computeForceTorque(coil9, coil10, PrecisionFactor(i));
        printf("(%20.16g, %20.16g, %20.16g) microN\n",
               1e6 * tempForce.first.x, 1e6 * tempForce.first.y, 1e6 * tempForce.first.z);
    }
    printf("\n");

    printf("Varying distance cases (fixed precision):\n\n");

    Coil coil11 = Coil(0.0425, 0.0, 0.0, 1);
    Coil coil12 = Coil(0.02, 0.0, 0.0, 1);
    for (int i = 0; i <= 11; ++i)
    {
        coil12.setPositionAndOrientation(
                vec3::Vector3(0.003, 0.0, 0.001 * i));
        tempForce = Coil::computeForceTorque(coil11, coil12, PrecisionFactor(8.0));
        printf("%20.16g\t %20.16g\n",
               1e7 * tempForce.first.x, 1e7 * tempForce.first.z);
    }
    printf("\n");

    double rArr[] = {0.0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.4999, 0.5, 0.5001, 0.6, 0.7, 0.8, 0.9, 1.0,
                     1.1, 1.2, 1.3, 1.4, 1.4999, 1.5, 1.5001, 1.6, 1.8, 1.9, 2.0, 3.0};
    Coil coil13 = Coil(1.0, 0.0, 0.0, 1);
    Coil coil14 = Coil(0.5, 0.0, 0.0, 1);
    for (double i : rArr)
    {
        coil14.setPositionAndOrientation(
                vec3::Vector3(i, 0.0, 1e-18));
        tempForce = Coil::computeForceTorque(coil13, coil14, PrecisionFactor(10.0));
        printf("%6g %20.16g\t %20.16g\n",i, 1e6 * tempForce.first.x, 1e6 * tempForce.first.z);
    }
    printf("\n");

    printf("Some more singular cases (different precisions):\n\n");

    Coil coil15 = Coil(1.0, 0.0, 0.0, 1);
    Coil coil16 = Coil(0.5, 0.0, 0.0, 1);
    coil16.setPositionAndOrientation(
            vec3::Vector3(2.0, 2.0, 2.0));
    for (int i = 1; i <= 9; ++i)
    {
        tempForce = Coil::computeForceTorque(coil15, coil16, PrecisionFactor(i));
        printf("(%20.16g, %20.16g, %20.16g) nanoN\n",
               1e9 * tempForce.first.x, 1e9 * tempForce.first.y, 1e9 * tempForce.first.z);
    }
    printf("\n");

    Coil coil17 = Coil(1.0, 0.0, 0.0, 1);
    Coil coil18 = Coil(0.5, 0.0, 0.0, 1);
    coil18.setPositionAndOrientation(
            vec3::Vector3(1.0, 2.0, 3.0), M_PI_2);
    for (int i = 1; i <= 9; ++i)
    {
        tempForce = Coil::computeForceTorque(coil17, coil18, PrecisionFactor(i));
        printf("(%20.16g, %20.16g, %20.16g) nanoN\n",
               1e9 * tempForce.first.x, 1e9 * tempForce.first.y, 1e9 * tempForce.first.z);
    }
    printf("\n");

    Coil coil19 = Coil(1.0, 0.0, 0.0, 1);
    Coil coil20 = Coil(0.5, 0.0, 0.0, 1);
    coil20.setPositionAndOrientation(
            vec3::Vector3(2.0, 2.0, 2.0), M_PI_2, M_PI_2);
    for (int i = 1; i <= 9; ++i)
    {
        tempForce = Coil::computeForceTorque(coil19, coil20, PrecisionFactor(i));
        printf("(%20.16g, %20.16g, %20.16g) nanoN\n",
               1e9 * tempForce.first.x, 1e9 * tempForce.first.y, 1e9 * tempForce.first.z);
    }
    printf("\n");
}

void Compare::forceTorqueZAxis()
{
    Coil prim1 = Coil(0.03, 0.03, 0.12, 3600, PrecisionFactor(6.0), 16);
    Coil sec1 = Coil(0.02, 0.025, 0.04, 1000, PrecisionFactor(6.0), 16);

    for (int i = 0; i < 500; ++i)
    {
        sec1.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.08 + i*0.001));
        printf("%.15f\n",
               Coil::computeForceTorque(prim1, sec1, PrecisionFactor(8.0), CPU_MT).first.z);
    }
    printf("\n");
}

void Compare::forceOnDipoleVsForceTorque()
{
    printf("This comparison shows the application of dipole approximation for force calculation\n");
    printf("The first column represents the force | and the second represents torque\n\n");
    printf("The first part involves a large coil that generates a field and a small one that feels the force\n");
    printf("The small coil is slightly off the center of the large coil and rotated\n\n");

    Coil coil1 = Coil(0.5, 0.1, 0.1, 100, 1000);
    Coil coil2 = Coil(0.005, 0.001, 0.006, 60, 100);
    coil2.setPositionAndOrientation(vec3::Vector3(0.01, 0.01, 0.02), 0.5, 1.0);

    std::pair amperePair = Coil::computeForceTorque(coil1, coil2);
    std::pair dipolePair = coil1.computeForceOnDipoleMoment(coil2.getPositionVector(), coil2.getMagneticMoment());

    printf("Ampere: %14.7g %14.7g %14.7g | %14.7g %14.7g %14.7g\n",
           amperePair.first.x, amperePair.first.y, amperePair.first.z,
           amperePair.second.x, amperePair.second.y, amperePair.second.z);

    printf("Dipole: %14.7g %14.7g %14.7g | %14.7g %14.7g %14.7g\n\n",
           dipolePair.first.x, dipolePair.first.y, dipolePair.first.z,
           dipolePair.second.x, dipolePair.second.y, dipolePair.second.z);

    printf("The second part involves two duplicate coils\n");
    printf("One coil is placed far away from the other and rotated\n\n");

    Coil coil3 = Coil(0.2, 0.1, 0.1, 100, 100);
    Coil coil4 = Coil(0.2, 0.1, 0.1, 100, 100);
    coil4.setPositionAndOrientation(vec3::Vector3(1.2, 1.6, 3.0), 2.0, 2.5);

    amperePair = Coil::computeForceTorque(coil3, coil4);
    dipolePair = coil3.computeForceOnDipoleMoment(coil4.getPositionVector(), coil4.getMagneticMoment());

    printf("Ampere: %14.7g %14.7g %14.7g | %14.7g %14.7g %14.7g\n",
           amperePair.first.x, amperePair.first.y, amperePair.first.z,
           amperePair.second.x, amperePair.second.y, amperePair.second.z);

    printf("Dipole: %14.7g %14.7g %14.7g | %14.7g %14.7g %14.7g\n\n",
           dipolePair.first.x, dipolePair.first.y, dipolePair.first.z,
           dipolePair.second.x, dipolePair.second.y, dipolePair.second.z);
}
