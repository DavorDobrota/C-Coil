#include "Coil.h"

#define _USE_MATH_DEFINES
#include <math.h>


vec3::Vector3Array Coil::computeAllBFieldVectors(const vec3::Vector3Array &pointVectors,
                                                 const PrecisionArguments &usedPrecision,
                                                 ComputeMethod computeMethod) const
{
    if (computeMethod == CPU_MT)
    {
        return calculateAllBFieldMT(pointVectors, usedPrecision);
    }
    else if (computeMethod == GPU)
    {
        return calculateAllBFieldGPU(pointVectors, usedPrecision);
    }
    else
    {
        vec3::Vector3Array computedFieldArr = vec3::Vector3Array();
        computedFieldArr.reserve(pointVectors.size());

        for (int i = 0; i < pointVectors.size(); ++i)
            computedFieldArr += computeBFieldVector(pointVectors[i], usedPrecision);

        return computedFieldArr;
    }
}

vec3::Vector3Array Coil::computeAllBFieldVectors(const vec3::Vector3Array &pointVectors,
                                                 ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllBFieldVectors(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllBFieldVectors(pointVectors, defaultPrecisionCPU, computeMethod);
}


vec3::Vector3Array Coil::computeAllAPotentialVectors(const vec3::Vector3Array &pointVectors,
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
        vec3::Vector3Array computedPotentialArr;
        computedPotentialArr.reserve(pointVectors.size());

        for (int i = 0; i < pointVectors.size(); ++i)
            computedPotentialArr += computeAPotentialVector(pointVectors[i], usedPrecision);

        return computedPotentialArr;
    }
}

vec3::Vector3Array Coil::computeAllAPotentialVectors(const vec3::Vector3Array &pointVectors,
                                                     ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllAPotentialVectors(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllAPotentialVectors(pointVectors, defaultPrecisionCPU, computeMethod);
}


vec3::Vector3Array Coil::computeAllEFieldVectors(const vec3::Vector3Array &pointVectors,
                                                 const PrecisionArguments &usedPrecision,
                                                 ComputeMethod computeMethod) const
{
    vec3::Vector3Array output = computeAllAPotentialVectors(pointVectors, usedPrecision, computeMethod);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (int i = 0; i < output.size(); ++i)
        output[i] *= frequencyFactor;

    return output;
}

vec3::Vector3Array Coil::computeAllEFieldVectors(const vec3::Vector3Array &pointVectors,
                                                 ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllEFieldVectors(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllEFieldVectors(pointVectors, defaultPrecisionCPU, computeMethod);
}


vec3::Matrix3Array Coil::computeAllBGradientMatrices(const vec3::Vector3Array &pointVectors,
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
        vec3::Matrix3Array computedGradientArr;
        computedGradientArr.reserve(pointVectors.size());

        for (int i = 0; i < pointVectors.size(); ++i)
            computedGradientArr += computeBGradientMatrix(pointVectors[i], usedPrecision);

        return computedGradientArr;
    }
}

vec3::Matrix3Array Coil::computeAllBGradientMatrices(const vec3::Vector3Array &pointVectors,
                                                     ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllBGradientMatrices(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllBGradientMatrices(pointVectors, defaultPrecisionCPU, computeMethod);
}