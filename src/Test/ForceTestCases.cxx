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
    Coil coil1 = Coil(0.1204, 0.1456, 0.5292, 126, 16500);
    Coil coil2 = Coil(0.336, 0.057, 0.552, 1890, 725);
    Coil coil3 = Coil(0.393, 0.0873, 0.552, 3792, 725);


    coil1.setThreadCount(16);
    coil2.setThreadCount(16);
    coil3.setThreadCount(16);

    std::pair<vec3::FieldVector3, vec3::FieldVector3> forcePair1, forcePair2;

    printf("Force and Torque in displacement\n");
    for (double dz = 0.0; dz <= 0.004; dz += 0.001)
    {
        for (double dr = 0.0; dr <= 0.004; dr += 0.002)
        {
            forcePair1 = Coil::computeAmpereForceGeneral(coil2, coil1, 0.0 + dz, 0.0 + dr,
                                                          PrecisionFactor(8.0), CPU_MT);
            forcePair2 = Coil::computeAmpereForceGeneral(coil3, coil1, 0.0 + dz, 0.0 + dr,
                                                          PrecisionFactor(8.0), CPU_MT);
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
        forcePair1 = Coil::computeAmpereForceGeneral(coil2, coil1, 0.0, 0.0, M_PI/360 * i,
                                                      PrecisionFactor(7.0), CPU_MT);
        forcePair1 = Coil::computeAmpereForceGeneral(coil3, coil1, 0.0, 0.0, M_PI/360 * i,
                                                      PrecisionFactor(7.0), CPU_MT);
        printf("%.10f %.10f %.10f %.10f %.10f %.10f\n",
               forcePair1.first.xComponent + forcePair2.first.xComponent,
               forcePair1.first.yComponent + forcePair2.first.yComponent,
               forcePair1.first.zComponent + forcePair2.first.zComponent,
               forcePair1.second.xComponent + forcePair2.second.xComponent,
               forcePair1.second.yComponent + forcePair2.second.yComponent,
               forcePair1.second.zComponent + forcePair2.second.zComponent);
    }
    printf("\n");
}

void testCoilAmpereForceThinCoils()
{
    ComputeMethod method = CPU_ST;

    for (int i = 1; i <= 10; ++i)
    {
        Coil coil1 = Coil(0.5, 1e-15, 0.2, 100);
        Coil coil2 = Coil(0.5, 1e-15, 0.3, 200);
        printf("%.16g mN\n", 1000 * Coil::computeAmpereForceZAxis(coil1, coil2, 0.4, PrecisionFactor(i), method));
    }
    printf("\n");

    for (int i = 1; i <= 10; ++i)
    {
        Coil coil3 = Coil(0.1, 1e-15, 0.15, 225);
        Coil coil4 = Coil(0.09, 1e-15, 0.15, 300);
        printf("%.16g mN\n", 1000 * Coil::computeAmpereForceZAxis(coil3, coil4, 0.05, PrecisionFactor(i), method));
    }
    printf("\n");

    for (int i = 1; i <= 10; ++i)
    {
        Coil coil5 = Coil(0.25, 1e-15, 0.06, 40);
        Coil coil6 = Coil(0.2, 1e-15, 0.04, 120);
        printf("%.16g mN\n", 1000 * Coil::computeAmpereForceZAxis(coil5, coil6, 0.3, PrecisionFactor(i), method));
    }
    printf("\n");

    for (int i = 1; i <= 10; ++i)
    {
        Coil coil7 = Coil(0.0762, 0.0832, 1e-15, 516, 1.42);
        Coil coil8 = Coil(0.0762, 0.0832, 1e-15, 516, 1.42);
        printf("%.16g N\n", Coil::computeAmpereForceZAxis(coil7, coil8, 0.0468, PrecisionFactor(i), method));
    }
    printf("\n");

    for (int i = 1; i <= 10; ++i)
    {
        Coil coil9 = Coil(0.16, 0.12, 1e-15, 100, 10);
        Coil coil10 = Coil(0.11, 0.15, 1e-15, 100, 10);
        printf("%.16g N\n", Coil::computeAmpereForceZAxis(coil9, coil10, 0.05, PrecisionFactor(i), method));
    }
    printf("\n");

    for (int i = 1; i <= 10; ++i)
    {
        Coil coil11 = Coil(0.1, 1e-15, 0.2, 100);
        Coil coil12 = Coil(0.2, 0.2, 1e-15, 100);
        printf("%.16g mN\n", 1000 * Coil::computeAmpereForceZAxis(coil11, coil12, 0.6, PrecisionFactor(i), method));
    }
    printf("\n");
}