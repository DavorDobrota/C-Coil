#include "Compare.h"
#include "Coil.h"
#include "Coil/EnumsAndConstants/ComputeMethod.h"
#include "Tensor.h"
#include "CoilGroup.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>
#include <vector>


void compCoilGroupMTD(int coilCount, int pointCount, int threadCount, bool print)
{
    double torusRadius = 1.0;

    CoilGroup torusGroup = CoilGroup();
    vec3::Vector3Array fieldPoints(pointCount);
    vec3::Vector3Array computedBField;

    for (int i = 0; i < pointCount; ++i)
        fieldPoints[i] = vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / pointCount);

    for (int i = 0; i < coilCount; ++i)
    {
        torusGroup.addCoil(
            torusRadius / 10.0, torusRadius / 100.0, torusRadius / 100.0,
            10000, 10, PrecisionFactor(), 8,
            vec3::Vector3::getFromCylindricalCoords(0.0, torusRadius, 2*M_PI * i / coilCount),
            M_PI_2, 2*M_PI * i / coilCount + M_PI_2
        );
    }

    torusGroup.setThreadCount(threadCount);
    torusGroup.setDefaultPrecisionFactor(PrecisionFactor(3.0));

    computedBField = torusGroup.computeAllBFieldVectors(fieldPoints, CPU_MT);

    if (print)
        for (int i = 0; i < pointCount; ++i)
            printf("%.15g\n",
                   std::sqrt(computedBField[i].x * computedBField[i].x +
                             computedBField[i].y * computedBField[i].y));
}
