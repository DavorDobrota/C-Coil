#include "Compare.h"
#include "Coil.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>


void compMethodPrecisionCPUvsGPU()
{
    Coil testCoil = Coil(0.03, 0.03, 0.12, 3600, PrecisionFactor(6.0), 12);

    int pointCount = 2000;
    double radius = 0.075;

    std::vector<vec3::CoordVector3> positionValues(pointCount);

    for (int i = 0; i < pointCount; i++)
        positionValues[i] = vec3::CoordVector3(vec3::SPHERICAL, radius, i * M_PI / pointCount, 0.0);

    std::vector<double> cpuPotential;
    std::vector<double> gpuPotential;

    std::vector<vec3::FieldVector3> cpuFieldVectors;
    std::vector<vec3::FieldVector3> gpuFieldVectors;

    cpuPotential = testCoil.computeAllAPotentialAbs(positionValues, CPU_ST);
    cpuFieldVectors = testCoil.computeAllBFieldComponents(positionValues, CPU_ST);

    gpuPotential = testCoil.computeAllAPotentialAbs(positionValues, GPU);
    gpuFieldVectors = testCoil.computeAllBFieldComponents(positionValues, GPU);

    FILE *output = fopen("output.txt", "w");

    for (int i = 0; i < pointCount; ++i)
    {
        fprintf(output, "%.20f\t%.20f\t%.20f\t%.20f\t%.20f\t%.20f\n",
                cpuFieldVectors[i].x, gpuFieldVectors[i].x,
                cpuFieldVectors[i].z, gpuFieldVectors[i].z,
                cpuPotential[i], gpuPotential[i]);
    }
}

void compPrecisionCPUvsGPU()
{
    printf("Comparing calculation precision in GPU and CPU cases (CPU | GPU)\n\n");

    Coil prim = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec = Coil(0.03, 0.03, 0.12, 3600);

    prim.setPositionAndOrientation(vec3::CoordVector3(), 0.2, 0.2);
    sec.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.1, 0.0, 0.2));

    printf("M inductance: %.15g | %.15g\n",
           Coil::computeMutualInductance(prim, sec, PrecisionFactor(7.0), CPU_MT),
           Coil::computeMutualInductance(prim, sec, PrecisionFactor(7.0), GPU));

    printf("Ampere force: %.15g | %.15g\n",
           Coil::computeAmpereForce(prim, sec, PrecisionFactor(7.0), CPU_MT).first.z,
           Coil::computeAmpereForce(prim, sec, PrecisionFactor(7.0), GPU).first.z);
}

void compMInductanceForSpecialCase()
{
    Coil primary = Coil(0.071335, 0.01397, 0.142748, 1142);
    Coil secondary = Coil(0.096945, 0.041529, 0.02413, 516);

    secondary.setPositionAndOrientation(vec3::CoordVector3(vec3::CARTESIAN, 0.30988, 0.0, 0.07366));

    printf("%.15f\n", Coil::computeMutualInductance(primary, secondary));
}
