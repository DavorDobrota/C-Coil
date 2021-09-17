#include <cstdio>
#include <cmath>
#include <string>

#include "Coil.h"
#include "Test.h"


int main()
{

//    testAmpereForceThinCoils();
//    testAmpereForceFilamentsGeneral();

//    testFunctionPerformance();
//
//    testPerformanceForComputeAll(PrecisionFactor(1.0), 50000, 5, 2);
//    testPerformanceForComputeAll(PrecisionFactor(2.0), 50000, 5, 3);
//    testPerformanceForComputeAll(PrecisionFactor(3.0), 50000, 5, 4);
//    testPerformanceForComputeAll(PrecisionFactor(4.0), 50000, 5, 5);
//    testPerformanceForComputeAll(PrecisionFactor(5.0), 50000, 5, 6);
//    testPerformanceForComputeAll(PrecisionFactor(6.0), 50000, 5, 7);
//    testPerformanceForComputeAll(PrecisionFactor(7.0), 50000, 5, 8);
//
//    testMutualInductanceZAxisMTScaling(12);
//    testMutualInductanceGeneralMTScaling(12);
//
//    testAmpereForceZAxisMTScaling(12);
//    testAmpereForceGeneralMTScaling(12);

//    testMutualInductanceZAxis();
//    testMutualInductanceGeneralParallelAxes();
//    testMutualInductanceGeneralConway();
//    testAmpereForceGeneralCase();

//    testPerformanceForVariousCoilTypes(700'001);

double R1 = 0.1;
double a1 = 0.06;
double b1 = 0.1;

double R2 = 0.1;
double a2 = 0.1;
double b2 = 0.1;

Coil coil1 = Coil(R1, a1, b1, 100);
Coil coil2 = Coil(R2, a2, b2, 100);
coil2.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN,0.0, 0.0, (b1 + b2) * 0.5 + (b1 + b2) / 40.0));

for (int i = 1; i <= 15; ++i)
    printf("%.15g\n", Coil::computeMutualInductance(coil1, coil2, PrecisionFactor(i), CPU_MT));
printf("\n");

for (int i = 1; i <= 15; ++i)
    printf("%.15g\n", Coil::computeMutualInductance(coil2, coil1, PrecisionFactor(i), CPU_MT));
printf("\n");


    return 0;
}
