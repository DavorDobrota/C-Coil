#include "Test.h"
#include "Coil.h"

#include <cstdio>
#include <cmath>

void testCoilAmpereForceForFilamentsZAxis()
{
    Coil coil1 = Coil(0.5, 1e-15, 1e-15, 1, 100);
    Coil coil2 = Coil(0.3, 1e-15, 1e-15, 1, 200);

    for (int i = 1; i <= 8; i++)
        printf("%.16g\n", 1000 * Coil::computeAmpereForceZAxis(coil1, coil2, 0.4, PrecisionFactor(i)));
    printf("\n");
}

void testCoilAmpereForceGeneralCase()
{
    Coil coil1 = Coil(0.1204, 0.1456, 0.5292, 126, 16500.0);
    Coil coil2 = Coil(0.336, 0.057, 0.552, 1890, 725.0);
    Coil coil3 = Coil(0.393, 0.0873, 0.552, 3792, 725.0);

    coil1.setThreadCount(16);
    coil2.setThreadCount(16);
    coil3.setThreadCount(16);

    std::vector<double> tempVector1, tempVector2;

    printf("Force and Torque in displacement\n");
    for (double dz = 0.0; dz <= 0.004; dz += 0.001)
    {
        for (double dr = 0.0; dr <= 0.004; dr += 0.002)
        {
            tempVector1 = Coil::computeAmpereForceGeneral(coil2, coil1, 0.0 + dz, 0.0 + dr,
                                                          PrecisionFactor(7.0), CPU_MT);
            tempVector2 = Coil::computeAmpereForceGeneral(coil3, coil1, 0.0 + dz, 0.0 + dr,
                                                          PrecisionFactor(7.0), CPU_MT);
            printf("%.10f %.10f %.10f %.10f %.10f %.10f\n",
                   tempVector1[0] + tempVector2[0], tempVector1[1] + tempVector2[1], tempVector1[2] + tempVector2[2],
                   tempVector1[3] + tempVector2[3], tempVector1[4] + tempVector2[4], tempVector1[5] + tempVector2[5]);
        }
    }
    printf("Force and Torque in rotation\n");
    for (int i = 0; i <= 10; ++i)
    {
        tempVector1 = Coil::computeAmpereForceGeneral(coil2, coil1, 0.0, 0.0, M_PI/360 * i,
                                                      PrecisionFactor(7.0), CPU_MT);
        tempVector2 = Coil::computeAmpereForceGeneral(coil3, coil1, 0.0, 0.0, M_PI/360 * i,
                                                      PrecisionFactor(7.0), CPU_MT);
        printf("%.10f %.10f %.10f %.10f %.10f %.10f\n",
               tempVector1[0] + tempVector2[0], tempVector1[1] + tempVector2[1], tempVector1[2] + tempVector2[2],
               tempVector1[3] + tempVector2[3], tempVector1[4] + tempVector2[4], tempVector1[5] + tempVector2[5]);
    }
    printf("\n");
}

void testCoilAmpereForceThinCoils()
{
    for (int i = 1; i <= 10; ++i)
    {
        Coil coil1 = Coil(0.5, 1e-15, 0.2, 100);
        Coil coil2 = Coil(0.5, 1e-15, 0.3, 200);
        printf("%.16g mN\n", 1000 * Coil::computeAmpereForceZAxis(coil1, coil2, 0.4, PrecisionFactor(i), CPU_ST));
    }
    printf("\n");

    for (int i = 1; i <= 10; ++i)
    {
        Coil coil1 = Coil(0.1, 1e-15, 0.15, 225);
        Coil coil2 = Coil(0.09, 1e-15, 0.15, 300);
        printf("%.16g mN\n", 1000 * Coil::computeAmpereForceZAxis(coil1, coil2, 0.05, PrecisionFactor(i), CPU_ST));
    }
    printf("\n");

    for (int i = 1; i <= 10; ++i)
    {
        Coil coil1 = Coil(0.25, 1e-15, 0.06, 40);
        Coil coil2 = Coil(0.2, 1e-15, 0.04, 120);
        printf("%.16g mN\n", 1000 * Coil::computeAmpereForceZAxis(coil1, coil2, 0.3, PrecisionFactor(i), CPU_ST));
    }
    printf("\n");

    for (int i = 1; i <= 10; ++i)
    {
        Coil coil1 = Coil(0.0762, 0.0832, 1e-15, 516, 1.42);
        Coil coil2 = Coil(0.0762, 0.0832, 1e-15, 516, 1.42);
        printf("%.16g N\n", Coil::computeAmpereForceZAxis(coil1, coil2, 0.0468, PrecisionFactor(i), CPU_ST));
    }
    printf("\n");

    for (int i = 1; i <= 10; ++i)
    {
        Coil coil1 = Coil(0.16, 0.12, 1e-15, 100, 10);
        Coil coil2 = Coil(0.11, 0.15, 1e-15, 100, 10);
        printf("%.16g N\n", Coil::computeAmpereForceZAxis(coil1, coil2, 0.05, PrecisionFactor(i), CPU_ST));
    }
    printf("\n");

    for (int i = 1; i <= 10; ++i)
    {
        Coil coil1 = Coil(0.1, 1e-15, 0.2, 100);
        Coil coil2 = Coil(0.2, 0.2, 1e-15, 100);
        printf("%.16g mN\n", 1000 * Coil::computeAmpereForceZAxis(coil1, coil2, 0.6, PrecisionFactor(i), CPU_ST));
    }
    printf("\n");
}