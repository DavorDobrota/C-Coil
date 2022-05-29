#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>

#include "Benchmark.h"
#include "Coil.h"
#include "CoilGroup.h"
#include "Test.h"
#include "Compare.h"


#pragma clang diagnostic push
#pragma ide diagnostic ignored "Simplify"
int main()
{

//    const int xDim = 2501;
//    const int yDim = 2501;
//    const double xSize = 0.25;
//    const double ySize = 0.25;
//    const double zPos = 0.0;
//
//    FILE *output = fopen("output.txt", "w");
//
//    Coil coil = Coil(0.03, 0.03, 0.12, 3600, 10, PrecisionFactor(1.0), 12);
//    coil.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.0, 0.0, 0.0), M_PI_2, 0.0);
//
//    printf("1\n");
//    std::vector<vec3::CoordVector3> positions(xDim * yDim);
//    printf("1\n");
//    std::vector<double> fields(xDim * yDim);
//    printf("1\n");
//
//    for (int i = 0; i < yDim; ++i)
//        for (int j = 0; j < xDim; ++j)
//        {
//            double yPos = ySize / (yDim-1) * (i - (yDim-1) / 2);
//            double xPos = xSize / (xDim-1) * (j - (xDim-1) / 2);
//            positions[i * xDim + j] = vec3::CoordVector3(vec3::CARTESIAN, xPos, yPos, zPos);
//        }
//    printf("1\n");
//    fields = coil.computeAllBFieldAbs(positions, CPU_MT);
//    printf("3\n");
//    for (int i = 0; i < yDim; ++i)
//    {
//        for (int j = 0; j < xDim; ++j)
//            fprintf(output, "%.15g\t", fields[i*xDim + j]);
//
//        fprintf(output, "\n");
//    }

//    benchComputeAllFieldsEveryCoilType(200'003, 8);
//    benchMathFunctions();

//    benchMInductanceZAxisMTScaling(16);
//    benchMInductanceGeneralMTScaling(16);
//    benchSelfInductance();
//    benchForceZAxisMTScaling(16);
//    benchForceGeneralMTScaling(16);

    benchComputeAllFieldsWorkloadScalingMT(PrecisionFactor(3.0), 16, 26);

    return 0;
}
#pragma clang diagnostic pop
