#include "Test.h"
#include "Coil.h"
#include "Tensor.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>


void testAmpereForceGeneralForZAxis()
{
    Coil prim = Coil(0.03, 0.03, 0.12, 3600, PrecisionFactor(6.0), 16);
    Coil sec = Coil(0.02, 0.025, 0.04, 1000, PrecisionFactor(6.0), 16);
    prim.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.1, 0.0, 0.0));


    std::pair<vec3::Vector3, vec3::Vector3> forcePair;

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
