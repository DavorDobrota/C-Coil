#include "Coil.h"

#define _USE_MATH_DEFINES
#include <math.h>


std::vector<vec3::Vector3> Coil::computeAllBFieldVectors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                         const PrecisionArguments &usedPrecision,
                                                         ComputeMethod computeMethod) const
{
    if (computeMethod == CPU_MT) {

        return calculateAllBFieldMT(pointVectors, usedPrecision);
    }
    else if (computeMethod == GPU) {

        return calculateAllBFieldGPU(pointVectors, usedPrecision);
    }
    else
    {
        std::vector<vec3::Vector3> computedFieldArr(pointVectors.size());

        for (int i = 0; i < pointVectors.size(); ++i)
            computedFieldArr[i] = computeBFieldVector(pointVectors[i], usedPrecision);

        return computedFieldArr;
    }
}

std::vector<vec3::Vector3> Coil::computeAllBFieldVectors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                         ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllBFieldVectors(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllBFieldVectors(pointVectors, defaultPrecisionCPU, computeMethod);
}


std::vector<vec3::Vector3> Coil::computeAllAPotentialVectors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                             const PrecisionArguments &usedPrecision,
                                                             ComputeMethod computeMethod) const
{
    if (computeMethod == CPU_MT)
    {
        return calculateAllAPotentialMT(pointVectors, usedPrecision);
    }
    else if (computeMethod == GPU)
    {
        return calculateAllAPotentialGPU(pointVectors, usedPrecision);
    }
    else
    {
        std::vector<vec3::Vector3> computedPotentialArr(pointVectors.size());

        for (int i = 0; i < pointVectors.size(); ++i)
            computedPotentialArr[i] = computeAPotentialVector(pointVectors[i], usedPrecision);

        return computedPotentialArr;
    }
}

std::vector<vec3::Vector3> Coil::computeAllAPotentialVectors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                             ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllAPotentialVectors(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllAPotentialVectors(pointVectors, defaultPrecisionCPU, computeMethod);
}


std::vector<vec3::Vector3> Coil::computeAllEFieldVectors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                         const PrecisionArguments &usedPrecision,
                                                         ComputeMethod computeMethod) const
{
    std::vector<vec3::Vector3> output = computeAllAPotentialVectors(pointVectors, usedPrecision, computeMethod);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (auto & i : output)
        i *= frequencyFactor;

    return output;
}

std::vector<vec3::Vector3> Coil::computeAllEFieldVectors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                         ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllEFieldVectors(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllEFieldVectors(pointVectors, defaultPrecisionCPU, computeMethod);
}


std::vector<vec3::Matrix3> Coil::computeAllBGradientMatrices(const std::vector<vec3::CoordVector3> &pointVectors,
                                                             const PrecisionArguments &usedPrecision,
                                                             ComputeMethod computeMethod) const
{
    if (computeMethod == CPU_MT)
    {
        return calculateAllBGradientMT(pointVectors, usedPrecision);
    }
    else if (computeMethod == GPU)
    {
        return calculateAllBGradientGPU(pointVectors, usedPrecision);
    }
    else
    {
        std::vector<vec3::Matrix3> computedGradientArr(pointVectors.size());

        for (int i = 0; i < pointVectors.size(); ++i)
            computedGradientArr[i] = computeBGradientMatrix(pointVectors[i], usedPrecision);

        return computedGradientArr;
    }
}

std::vector<vec3::Matrix3> Coil::computeAllBGradientMatrices(const std::vector<vec3::CoordVector3> &pointVectors,
                                                             ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllBGradientMatrices(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllBGradientMatrices(pointVectors, defaultPrecisionCPU, computeMethod);
}