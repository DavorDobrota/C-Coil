#include "Test.h"
#include "Coil.h"
#include "Tensor.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>


void compAmpereForceZAxis()
{
    Coil prim1 = Coil(0.03, 0.03, 0.12, 3600, PrecisionFactor(6.0), 16);
    Coil sec1 = Coil(0.02, 0.025, 0.04, 1000, PrecisionFactor(6.0), 16);

    for (int i = 0; i < 500; ++i)
    {
        sec1.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.08 + i*0.001));
        printf("%.15f\n",
               Coil::computeAmpereForce(prim1, sec1, PrecisionFactor(8.0), CPU_MT).first.z);
    }
    printf("\n");
}

void testAmpereForceGeneralForZAxis()
{
    Coil prim = Coil(0.03, 0.03, 0.12, 3600, PrecisionFactor(6.0), 16);
    Coil sec = Coil(0.02, 0.025, 0.04, 1000, PrecisionFactor(6.0), 16);
    prim.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.1, 0.0, 0.0));


    std::pair<vec3::FieldVector3, vec3::FieldVector3> forcePair;

    for (int i = 0; i < 100; ++i)
    {
        sec.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.1, 0.0, 0.08 + i * 0.001));
        forcePair = Coil::computeAmpereForce(prim, sec, PrecisionFactor(8.0), CPU_MT);
        printf("%.15f %.15f %.15f %.15f %.15f %.15f\n",
               forcePair.first.x, forcePair.first.y, forcePair.first.z,
               forcePair.second.x, forcePair.second.y, forcePair.second.z);
    }
    printf("\n");
}

void testGradientTensor()
{
    Coil coil = Coil(1.0, 0.0, 1e-10, 1);

    vec3::Matrix3 tensor;

    printf("Z Axis test\n");
    for (int i = 0; i < 1000; ++i)
    {
        tensor = coil.computeBGradientTensor(vec3::CoordVector3(vec3::CYLINDRICAL, 0.001 * i, 0.0, 0.0));
        printf("%.15f %.15f %.15f\n", tensor.xx / (1e-7), tensor.yy / (1e-7), tensor.zz / (1e-7));
    }
    printf("\n");

    printf("Off axis test\n");
    for (int i = 0; i < 1000; ++i)
    {
        tensor = coil.computeBGradientTensor(vec3::CoordVector3(vec3::CYLINDRICAL, i * 0.001, 0.5, M_PI / 4));
        printf("%.8f %.8f %.8f | %.8f %.8f %.8f | %.8f %.8f %.8f\n",
               tensor.xx / (1e-7), tensor.xy / (1e-7), tensor.xz / (1e-7),
               tensor.yx / (1e-7), tensor.yy / (1e-7), tensor.yz / (1e-7),
               tensor.zx / (1e-7), tensor.zy / (1e-7), tensor.zz / (1e-7));
    }
    printf("\n");
}

void compForceOnDipoleVsAmpereForce()
{
    printf("This comparison shows the application of dipole approximation for force calculation\n");
    printf("The first column represents the force | and the second represents torque\n\n");
    printf("The first part involves a large coil that generates a field and a small one that feels the force\n");
    printf("The small coil is slightly off the center of the large coil and rotated\n\n");

    Coil coil1 = Coil(0.5, 0.1, 0.1, 100, 1000);
    Coil coil2 = Coil(0.005, 0.001, 0.006, 60, 100);
    coil2.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.01, 0.01, 0.02), 0.5, 1.0);

    std::pair amperePair = Coil::computeAmpereForce(coil1, coil2);
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
    coil4.setPositionAndOrientation(
            vec3::CoordVector3(vec3::CARTESIAN, 1.2, 1.6, 3.0), 2.0, 2.5);

    amperePair = Coil::computeAmpereForce(coil3, coil4);
    dipolePair = coil3.computeForceOnDipoleMoment(coil4.getPositionVector(), coil4.getMagneticMoment());

    printf("Ampere: %14.7g %14.7g %14.7g | %14.7g %14.7g %14.7g\n",
           amperePair.first.x, amperePair.first.y, amperePair.first.z,
           amperePair.second.x, amperePair.second.y, amperePair.second.z);

    printf("Dipole: %14.7g %14.7g %14.7g | %14.7g %14.7g %14.7g\n\n",
           dipolePair.first.x, dipolePair.first.y, dipolePair.first.z,
           dipolePair.second.x, dipolePair.second.y, dipolePair.second.z);
}
