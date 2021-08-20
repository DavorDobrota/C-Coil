#include "Test.h"
#include "Coil.h"

#include <cstdio>
#include <cmath>

void testCoilAmpereForceForFilamentsZAxis()
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

void testCoilAmpereForceGeneralCase()
{
    Coil coil1 = Coil(0.1204, 0.1456, 0.5292, 126, 16500);
    Coil coil2 = Coil(0.336, 0.057, 0.552, 1890, 725);
    Coil coil3 = Coil(0.393, 0.0873, 0.552, 3792, 725);


    coil1.setThreadCount(16);
    coil2.setThreadCount(16);
    coil3.setThreadCount(16);

    std::pair<vec3::FieldVector3, vec3::FieldVector3> forcePair1, forcePair2;
    auto precision = PrecisionFactor(8.0);

    printf("Force and Torque in displacement\n");
    for (double dz = 0.0; dz <= 0.004; dz += 0.001)
    {
        for (double dr = 0.0; dr <= 0.004; dr += 0.002)
        {
            coil1.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, dr, 0.0, dz));
            forcePair1 = Coil::computeAmpereForce(coil2, coil1, precision, CPU_MT);
            forcePair2 = Coil::computeAmpereForce(coil3, coil1, precision, CPU_MT);
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
        forcePair1 = Coil::computeAmpereForce(coil2, coil1, precision, CPU_MT);
        forcePair1 = Coil::computeAmpereForce(coil3, coil1, precision, CPU_MT);
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

    std::pair<vec3::FieldVector3, vec3::FieldVector3> tempForce;

    for (int i = 1; i <= 10; ++i)
    {
        Coil coil1 = Coil(0.5, 0.0, 0.2, 100);
        Coil coil2 = Coil(0.5, 0.0, 0.3, 200);
        coil2.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.4));
        tempForce = Coil::computeAmpereForce(coil1, coil2, PrecisionFactor(i), method);
        printf("%.16g mN\n", 1000 * tempForce.first.zComponent);
    }
    printf("\n");

    for (int i = 1; i <= 10; ++i)
    {
        Coil coil3 = Coil(0.1, 0.0, 0.15, 225);
        Coil coil4 = Coil(0.09, 0.0, 0.15, 300);
        coil4.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.05));
        tempForce = Coil::computeAmpereForce(coil3, coil4, PrecisionFactor(i), method);
        printf("%.16g mN\n", 1000 * tempForce.first.zComponent);
    }
    printf("\n");

    for (int i = 1; i <= 10; ++i)
    {
        Coil coil5 = Coil(0.25, 0.0, 0.06, 40);
        Coil coil6 = Coil(0.2, 0.0, 0.04, 120);
        coil6.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.3));
        tempForce = Coil::computeAmpereForce(coil5, coil6, PrecisionFactor(i), method);
        printf("%.16g mN\n", 1000 * tempForce.first.zComponent);
    }
    printf("\n");

    for (int i = 1; i <= 10; ++i)
    {
        Coil coil7 = Coil(0.0762, 0.0832, 0.0, 516, 1.42);
        Coil coil8 = Coil(0.0762, 0.0832, 0.0, 516, 1.42);
        coil8.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.3));
        tempForce = Coil::computeAmpereForce(coil7, coil8, PrecisionFactor(i), method);
        printf("%.16g N\n", tempForce.first.zComponent);
    }
    printf("\n");

    for (int i = 1; i <= 10; ++i)
    {
        Coil coil9 = Coil(0.16, 0.12, 0.0, 100, 10);
        Coil coil10 = Coil(0.11, 0.15, 0.0, 100, 10);
        coil10.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.3));
        tempForce = Coil::computeAmpereForce(coil9, coil10, PrecisionFactor(i), method);
        printf("%.16g N\n", tempForce.first.zComponent);
    }
    printf("\n");

    for (int i = 1; i <= 10; ++i)
    {
        Coil coil11 = Coil(0.1, 0.0, 0.2, 100);
        Coil coil12 = Coil(0.2, 0.2, 0.0, 100);
        coil12.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.3));
        tempForce = Coil::computeAmpereForce(coil11, coil12, PrecisionFactor(i), method);
        printf("%.16g mN\n", 1000 * tempForce.first.zComponent);
    }
    printf("\n");
}