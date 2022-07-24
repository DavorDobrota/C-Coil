#include "Coil.h"

#define _USE_MATH_DEFINES
#include <math.h>


vec3::Vector3Array Coil::computeAllBFieldVectors(const std::vector<vec3::CoordVector3> &pointVectors,
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
        vec3::Vector3Array computedFieldArr = vec3::Vector3Array();
        computedFieldArr.reserve(pointVectors.size());

        for (const auto & pointVector : pointVectors)
            computedFieldArr += computeBFieldVector(pointVector, usedPrecision);

        return computedFieldArr;
    }
}

vec3::Vector3Array Coil::computeAllBFieldVectors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                 ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllBFieldVectors(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllBFieldVectors(pointVectors, defaultPrecisionCPU, computeMethod);
}


vec3::Vector3Array Coil::computeAllAPotentialVectors(const std::vector<vec3::CoordVector3> &pointVectors,
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
        vec3::Vector3Array computedPotentialArr = vec3::Vector3Array();
        computedPotentialArr.reserve(pointVectors.size());

        for (const auto & pointVector : pointVectors)
            computedPotentialArr += computeAPotentialVector(pointVector, usedPrecision);

        return computedPotentialArr;
    }
}

vec3::Vector3Array Coil::computeAllAPotentialVectors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                     ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllAPotentialVectors(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllAPotentialVectors(pointVectors, defaultPrecisionCPU, computeMethod);
}


vec3::Vector3Array Coil::computeAllEFieldVectors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                 const PrecisionArguments &usedPrecision,
                                                 ComputeMethod computeMethod) const
{
    vec3::Vector3Array output = computeAllAPotentialVectors(pointVectors, usedPrecision, computeMethod);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (int i = 0; i < output.size(); ++i)
        output[i] *= frequencyFactor;

    return output;
}

vec3::Vector3Array Coil::computeAllEFieldVectors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                 ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllEFieldVectors(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllEFieldVectors(pointVectors, defaultPrecisionCPU, computeMethod);
}


vec3::Matrix3Array Coil::computeAllBGradientMatrices(const std::vector<vec3::CoordVector3> &pointVectors,
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
        vec3::Matrix3Array computedGradientArr = vec3::Matrix3Array();
        computedGradientArr.reserve(pointVectors.size());

        for (const auto & pointVector : pointVectors)
            computedGradientArr += computeBGradientMatrix(pointVector, usedPrecision);

        return computedGradientArr;
    }
}

vec3::Matrix3Array Coil::computeAllBGradientMatrices(const std::vector<vec3::CoordVector3> &pointVectors,
                                                     ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllBGradientMatrices(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllBGradientMatrices(pointVectors, defaultPrecisionCPU, computeMethod);
}