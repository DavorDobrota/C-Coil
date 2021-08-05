#include "Test.h"
#include "Coil.h"
#include "OldCoil.h"
#include "Tensor.h"

#include <cstdio>
#include <cmath>


void testCoilAmpereForceZAxis()
{
    Coil prim1 = Coil(0.03, 0.03, 0.12, 3600, PrecisionFactor(6.0), 16);
    Coil sec1 = Coil(0.02, 0.025, 0.04, 1000, PrecisionFactor(6.0), 16);

    for (int i = 0; i < 500; ++i)
    {
        printf("%.15f\n",
               Coil::computeAmpereForceZAxis(prim1, sec1, 0.08 + i*0.001, PrecisionFactor(8.0), CPU_MT));
    }
    printf("\n");
}

void testCoilAmpereForceGeneralForZAxis()
{
    Coil prim1 = Coil(0.03, 0.03, 0.12, 3600, PrecisionFactor(6.0), 16);
    Coil sec1 = Coil(0.02, 0.025, 0.04, 1000, PrecisionFactor(6.0), 16);

    std::pair<vec3::FieldVector3, vec3::FieldVector3> forcePair;

    for (int i = 0; i < 100; ++i)
    {
        forcePair = Coil::computeAmpereForceGeneral(prim1, sec1, 0.08 + i*0.001, 1e-18,
                                                     PrecisionFactor(8.0), CPU_MT);
        printf("%.15f %.15f %.15f %.15f %.15f %.15f\n",
               forcePair.first.xComponent, forcePair.first.yComponent, forcePair.first.zComponent,
               forcePair.second.xComponent, forcePair.second.yComponent, forcePair.second.zComponent);
    }
    printf("\n");
}

void testCoilGradientTensor()
{
    Coil coil = Coil(1.0, 1e-15, 1e-15, 1);

    vec3::Matrix3 tensor;

    printf("Z Axis test\n");
    for (int i = 0; i < 1000; ++i)
    {
        tensor = coil.computeBGradientTensor(vec3::CoordVector3(vec3::CYLINDRICAL, 0.001 * i, 0.0, 0.0));
        printf("%.15f %.15f %.15f\n", tensor.xxElement / (1e-7), tensor.yyElement / (1e-7), tensor.zzElement / (1e-7));
    }
    printf("\n");

    printf("Off axis test\n");
    for (int i = 0; i < 1000; ++i)
    {
        tensor = coil.computeBGradientTensor(vec3::CoordVector3(vec3::CYLINDRICAL, i * 0.001, 0.5, M_PI / 4));
        printf("%.8f %.8f %.8f | %.8f %.8f %.8f | %.8f %.8f %.8f\n",
               tensor.xxElement / (1e-7), tensor.xyElement / (1e-7), tensor.xzElement / (1e-7),
               tensor.yxElement / (1e-7), tensor.yyElement / (1e-7), tensor.yzElement / (1e-7),
               tensor.zxElement / (1e-7), tensor.zyElement / (1e-7), tensor.zzElement / (1e-7));
    }
    printf("\n");
}