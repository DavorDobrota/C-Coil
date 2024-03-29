#include "Compare.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>


void Compare::fieldsPrecisionCPUvsGPU()
{
    printf("Compares potential A and field B values for specific points, values are written in output.txt\n\n");

    Coil testCoil = Coil(0.03, 0.03, 0.12, 3600, PrecisionFactor(6.0), 12);

    int pointCount = 2000;
    double radius = 0.075;

    vec3::Vector3Array positionValues(pointCount);

    for (int i = 0; i < pointCount; i++)
        positionValues[i] = vec3::Vector3::getFromSphericalCoords(radius, i * M_PI / pointCount, 0.0);

    std::vector<double> cpuPotential;
    std::vector<double> gpuPotential;

    vec3::Vector3Array cpuFieldVectors;
    vec3::Vector3Array gpuFieldVectors;

    cpuPotential = testCoil.computeAllAPotentialVectors(positionValues, CPU_ST).abs();
    cpuFieldVectors = testCoil.computeAllBFieldVectors(positionValues, CPU_ST);

    gpuPotential = testCoil.computeAllAPotentialVectors(positionValues, GPU).abs();
    gpuFieldVectors = testCoil.computeAllBFieldVectors(positionValues, GPU);

    FILE *output = fopen("output.txt", "w");

    for (int i = 0; i < pointCount; ++i)
    {
        fprintf(output, "%.20f\t%.20f\t%.20f\t%.20f\t%.20f\t%.20f\n",
                cpuFieldVectors[i].x, gpuFieldVectors[i].x,
                cpuFieldVectors[i].z, gpuFieldVectors[i].z,
                cpuPotential[i], gpuPotential[i]);
        printf("%9.4g %9.4g %9.4g\n",
               std::abs(cpuFieldVectors[i].x - gpuFieldVectors[i].x) / std::abs(cpuFieldVectors[i].x),
               std::abs(cpuFieldVectors[i].z - gpuFieldVectors[i].z) / std::abs(cpuFieldVectors[i].z),
               std::abs(cpuPotential[i] - gpuPotential[i]) / std::abs(cpuPotential[i])
        );
    }

    fclose(output);
}

void Compare::mutualInductanceAndForceTorquePrecisionCPUvsGPU()
{
    printf("Comparing calculation precision in GPU and CPU cases (CPU | GPU) fast\n\n");

    Coil prim1 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec1 = Coil(0.03, 0.03, 0.12, 3600);

    prim1.setPositionAndOrientation(vec3::Vector3(-0.1, -0.1, -0.1), 0.2, 0.2);
    sec1.setPositionAndOrientation(vec3::Vector3(0.1, 0.0, 0.2));

    printf("M inductance: %.15g | %.15g\n",
           Coil::computeMutualInductance(prim1, sec1, PrecisionFactor(7.0), CPU_MT),
           Coil::computeMutualInductance(prim1, sec1, PrecisionFactor(7.0), GPU));

    printf("Ampere force: %.15g | %.15g\n",
           Coil::computeForceTorque(prim1, sec1, PrecisionFactor(7.0), CPU_MT).first.z,
           Coil::computeForceTorque(prim1, sec1, PrecisionFactor(7.0), GPU).first.z);

    vec3::Vector3Array coordArr;
    coordArr.append(vec3::Vector3());
    coordArr.append(vec3::Vector3());

    vec3::Matrix3Array gradientArr = prim1.computeAllBGradientMatrices(coordArr, GPU);

    printf("Gradient xx : %.15g | %.15g\n", prim1.computeBGradientMatrix(vec3::Vector3()).xx, gradientArr[0].xx);
    printf("Gradient yy : %.15g | %.15g\n", prim1.computeBGradientMatrix(vec3::Vector3()).yy, gradientArr[0].yy);
    printf("Gradient zz : %.15g | %.15g\n", prim1.computeBGradientMatrix(vec3::Vector3()).zz, gradientArr[0].zz);

    printf("\n");
    printf("Comparing calculation precision in GPU and CPU cases (CPU | GPU) slow\n\n");

    Coil prim2 = Coil(0.03, 0.06, 0.0, 3600);
    Coil sec2 = Coil(0.03, 0.03, 0.0, 3600);

    prim2.setPositionAndOrientation(vec3::Vector3(-0.1, -0.1, -0.1), 0.2, 0.2);
    sec2.setPositionAndOrientation(vec3::Vector3(0.1, 0.0, 0.2));

    printf("M inductance: %.15g | %.15g\n",
           Coil::computeMutualInductance(prim2, sec2, PrecisionFactor(7.0), CPU_MT),
           Coil::computeMutualInductance(prim2, sec2, PrecisionFactor(7.0), GPU));

    printf("Ampere force: %.15g | %.15g\n",
           Coil::computeForceTorque(prim2, sec2, PrecisionFactor(7.0), CPU_MT).first.z,
           Coil::computeForceTorque(prim2, sec2, PrecisionFactor(7.0), GPU).first.z);

    gradientArr = prim2.computeAllBGradientMatrices(coordArr, GPU);

    printf("Gradient xx : %.15g | %.15g\n", prim2.computeBGradientMatrix(vec3::Vector3()).xx, gradientArr[0].xx);
    printf("Gradient yy : %.15g | %.15g\n", prim2.computeBGradientMatrix(vec3::Vector3()).yy, gradientArr[0].yy);
    printf("Gradient zz : %.15g | %.15g\n", prim2.computeBGradientMatrix(vec3::Vector3()).zz, gradientArr[0].zz);
}

void Compare::mutualInductanceSpecialCase()
{
    Coil primary = Coil(0.071335, 0.01397, 0.142748, 1142);
    Coil secondary = Coil(0.096945, 0.041529, 0.02413, 516);

    secondary.setPositionAndOrientation(vec3::Vector3(0.30988, 0.0, 0.07366));

    printf("%.15f\n", Coil::computeMutualInductance(primary, secondary));
}
