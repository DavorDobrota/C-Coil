#include "Test.h"
#include "Coil.h"
#include "Tensor.h"

#include <cstdio>
#include <cmath>
#include <chrono>


void testAmpereForceZAxis()
{
    Coil prim1 = Coil(0.03, 0.03, 0.12, 3600, PrecisionFactor(6.0), 16);
    Coil sec1 = Coil(0.02, 0.025, 0.04, 1000, PrecisionFactor(6.0), 16);

    for (int i = 0; i < 500; ++i)
    {
        sec1.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.08 + i*0.001));
        printf("%.15f\n",
               Coil::computeAmpereForce(prim1, sec1, PrecisionFactor(8.0), CPU_MT).first.zComponent);
    }
    printf("\n");
}

void testAmpereForceZAxisPerformance(ComputeMethod method, int nThreads)
{
    using namespace std::chrono;

    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);
    primary.setThreadCount(nThreads);
    secondary.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.2));

    int nOps = 8192;
    double temp;

    printf("Expected execution time for one Ampere force z-axis calculation of specified precision\n");

    for (int i = 1; i <= 9; ++i)
    {
        int currentOperations = nOps / (int) pow(2, i);

        high_resolution_clock::time_point begin_time = high_resolution_clock::now();
        for (int j = 0; j < currentOperations; ++j)
            temp = Coil::computeAmpereForce(primary, secondary, PrecisionFactor(i), method).first.zComponent;
        double interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("precisionFactor(%.1f) : %6.4f ms/op\n", (double) i, 1'000.0 * interval / currentOperations);
    }
}

void testAmpereForceZAxisMTScaling(int maxThreads)
{
    printf("Performance comparison between different numbers of threads:\n");

    printf(" -> single thread:\n");
    testAmpereForceZAxisPerformance(CPU_ST);
    printf("\n");

    for (int i = 2; i <= maxThreads; ++i)
    {
        printf(" -> %2d threads:\n", i);
        testAmpereForceZAxisPerformance(CPU_MT, i);
        printf("\n");
    }
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
               forcePair.first.xComponent, forcePair.first.yComponent, forcePair.first.zComponent,
               forcePair.second.xComponent, forcePair.second.yComponent, forcePair.second.zComponent);
    }
    printf("\n");
}

void testAmpereForceGeneralPerformance(ComputeMethod method, int nThreads)
{
    using namespace std::chrono;

    Coil primary = Coil(0.1, 0.1, 0.1, 100);
    Coil secondary = Coil(0.3, 0.1, 0.1, 100);
    primary.setThreadCount(nThreads);

    primary.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.1, 0.0, 0.0));
    secondary.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.1, 0.0, 0.2));

    int nOps = 1024;
    std::pair<vec3::FieldVector3, vec3::FieldVector3> temp;

    printf("Expected execution time for one Ampere force general case calculation of specified precision\n");

    for (int i = 1; i <= 9; ++i)
    {
        int currentOperations = nOps / (int) pow(2, i);

        high_resolution_clock::time_point begin_time = high_resolution_clock::now();
        for (int j = 0; j < currentOperations; ++j)
            temp = Coil::computeAmpereForce(primary, secondary, PrecisionFactor(i), method);
        double interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
        printf("precisionFactor(%.1f) : %6.3f ms/op\n", (double) i, 1'000.0 * interval / currentOperations);
    }
}

void testAmpereForceGeneralMTScaling(int maxThreads)
{
    printf("Performance comparison between different numbers of threads:\n");

    printf(" -> single thread:\n");
    testAmpereForceGeneralPerformance(CPU_ST);
    printf("\n");

    for (int i = 2; i <= maxThreads; ++i)
    {
        printf(" -> %2d threads:\n", i);
        testAmpereForceGeneralPerformance(CPU_MT, i);
        printf("\n");
    }
}

void testGradientTensor()
{
    Coil coil = Coil(1.0, 0.0, 1e-10, 1);

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

void testForceOnDipole()
{
    Coil coil1 = Coil(0.5, 0.1, 0.1, 100, 1000);
    Coil coil2 = Coil(0.005, 0.001, 0.006, 60, 100);
    coil2.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.01, 0.01, 0.02), 0.5, 1.0);

    std::pair amperePair = Coil::computeAmpereForce(coil1, coil2);
    std::pair momentPair = coil1.computeForceOnDipoleMoment(coil2.getPositionVector(), coil2.getMagneticMoment());

    printf("%21.15g %21.15g %21.15g\n%21.15g %21.15g %21.15g\n\n",
           amperePair.first.xComponent + amperePair.first.xComponent,
           amperePair.first.yComponent + amperePair.first.yComponent,
           amperePair.first.zComponent + amperePair.first.zComponent,
           amperePair.second.xComponent + amperePair.second.xComponent,
           amperePair.second.yComponent + amperePair.second.yComponent,
           amperePair.second.zComponent + amperePair.second.zComponent);

    printf("%21.15g %21.15g %21.15g\n%21.15g %21.15g %21.15g\n\n",
           momentPair.first.xComponent + momentPair.first.xComponent,
           momentPair.first.yComponent + momentPair.first.yComponent,
           momentPair.first.zComponent + momentPair.first.zComponent,
           momentPair.second.xComponent + momentPair.second.xComponent,
           momentPair.second.yComponent + momentPair.second.yComponent,
           momentPair.second.zComponent + momentPair.second.zComponent);

    Coil coil3 = Coil(0.2, 0.1, 0.1, 100, 100);
    Coil coil4 = Coil(0.2, 0.1, 0.1, 100, 100);
    coil4.setPositionAndOrientation(
            vec3::CoordVector3(vec3::CARTESIAN, 1.2, 1.6, 3.0), 2.0, 2.5);

    amperePair = Coil::computeAmpereForce(coil3, coil4);
    momentPair = coil3.computeForceOnDipoleMoment(coil4.getPositionVector(), coil4.getMagneticMoment());

    printf("%21.15g %21.15g %21.15g\n%21.15g %21.15g %21.15g\n\n",
           amperePair.first.xComponent + amperePair.first.xComponent,
           amperePair.first.yComponent + amperePair.first.yComponent,
           amperePair.first.zComponent + amperePair.first.zComponent,
           amperePair.second.xComponent + amperePair.second.xComponent,
           amperePair.second.yComponent + amperePair.second.yComponent,
           amperePair.second.zComponent + amperePair.second.zComponent);

    printf("%21.15g %21.15g %21.15g\n%21.15g %21.15g %21.15g\n\n",
           momentPair.first.xComponent + momentPair.first.xComponent,
           momentPair.first.yComponent + momentPair.first.yComponent,
           momentPair.first.zComponent + momentPair.first.zComponent,
           momentPair.second.xComponent + momentPair.second.xComponent,
           momentPair.second.yComponent + momentPair.second.yComponent,
           momentPair.second.zComponent + momentPair.second.zComponent);
}
