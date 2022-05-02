#include "Compare.h"
#include "Coil.h"
#include "ComputeMethod.h"
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
    std::vector<vec3::CoordVector3> fieldPoints(pointCount);
    std::vector<vec3::FieldVector3> computedBField;

    for (int i = 0; i < pointCount; ++i)
        fieldPoints[i] = vec3::CoordVector3(vec3::CYLINDRICAL, 0.0, torusRadius, 2*M_PI * i / pointCount);

    for (int i = 0; i < coilCount; ++i)
    {
        Coil tempCoil = Coil(torusRadius / 10.0, torusRadius / 100.0, torusRadius / 100.0, 10000, 10);
        tempCoil.setPositionAndOrientation(
                vec3::CoordVector3(vec3::CYLINDRICAL, 0.0, torusRadius, 2*M_PI * i / coilCount),
                M_PI_2, 2*M_PI * i / coilCount + M_PI_2);
        torusGroup.addCoil(tempCoil);
    }

    torusGroup.setThreadCount(threadCount);
    torusGroup.setDefaultPrecisionFactor(PrecisionFactor(3.0));

    computedBField = torusGroup.computeAllBFieldComponents(fieldPoints, CPU_MT);

    if (print)
        for (int i = 0; i < pointCount; ++i)
            printf("%.15g\n",
                   std::sqrt(computedBField[i].x * computedBField[i].x +
                             computedBField[i].y * computedBField[i].y));
}