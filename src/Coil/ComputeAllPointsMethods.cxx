#include "Coil.h"

#define _USE_MATH_DEFINES
#include <math.h>


std::vector<vec3::FieldVector3> Coil::computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &pointVectors,
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
        std::vector<vec3::FieldVector3> computedFieldArr(pointVectors.size());

        for (int i = 0; i < pointVectors.size(); ++i)
            computedFieldArr[i] = computeBFieldVector(pointVectors[i], usedPrecision);

        return computedFieldArr;
    }
}

std::vector<vec3::FieldVector3> Coil::computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                 ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllBFieldComponents(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllBFieldComponents(pointVectors, defaultPrecisionCPU, computeMethod);
}

std::vector<double>Coil::computeAllBFieldX(const std::vector<vec3::CoordVector3> &pointVectors,
                                           const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(pointVectors, usedPrecision, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].x;

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldX(const std::vector<vec3::CoordVector3> &pointVectors,
                                            ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllBFieldX(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllBFieldX(pointVectors, defaultPrecisionCPU, computeMethod);
}

std::vector<double>Coil::computeAllBFieldY(const std::vector<vec3::CoordVector3> &pointVectors,
                                           const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(pointVectors, usedPrecision, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].y;

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldY(const std::vector<vec3::CoordVector3> &pointVectors,
                                           ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllBFieldY(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllBFieldY(pointVectors, defaultPrecisionCPU, computeMethod);
}

std::vector<double>Coil::computeAllBFieldZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                           const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(pointVectors, usedPrecision, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].z;

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                            ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllBFieldZ(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllBFieldZ(pointVectors, defaultPrecisionCPU, computeMethod);
}

std::vector<double>Coil::computeAllBFieldAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                             const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(pointVectors, usedPrecision, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = std::sqrt(computedFieldArr[i].x * computedFieldArr[i].x +
                                 computedFieldArr[i].y * computedFieldArr[i].y +
                                 computedFieldArr[i].z * computedFieldArr[i].z);

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                              ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllBFieldAbs(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllBFieldAbs(pointVectors, defaultPrecisionCPU, computeMethod);
}


std::vector<vec3::FieldVector3> Coil::computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &pointVectors,
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
        std::vector<vec3::FieldVector3> computedPotentialArr(pointVectors.size());

        for (int i = 0; i < pointVectors.size(); ++i)
            computedPotentialArr[i] = computeAPotentialVector(pointVectors[i], usedPrecision);

        return computedPotentialArr;
    }
}

std::vector<vec3::FieldVector3> Coil::computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                     ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllAPotentialComponents(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllAPotentialComponents(pointVectors, defaultPrecisionCPU, computeMethod);
}

std::vector<double> Coil::computeAllAPotentialX(const std::vector<vec3::CoordVector3> &pointVectors,
                                                const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllAPotentialComponents(pointVectors, usedPrecision, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].x;

    return outputArr;
}

std::vector<double> Coil::computeAllAPotentialX(const std::vector<vec3::CoordVector3> &pointVectors,
                                                ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllAPotentialX(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllAPotentialX(pointVectors, defaultPrecisionCPU, computeMethod);
}

std::vector<double> Coil::computeAllAPotentialY(const std::vector<vec3::CoordVector3> &pointVectors,
                                                const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllAPotentialComponents(pointVectors, usedPrecision, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].y;

    return outputArr;
}

std::vector<double> Coil::computeAllAPotentialY(const std::vector<vec3::CoordVector3> &pointVectors,
                                                ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllAPotentialY(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllAPotentialY(pointVectors, defaultPrecisionCPU, computeMethod);
}

std::vector<double> Coil::computeAllAPotentialZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                                const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllAPotentialComponents(pointVectors, usedPrecision, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].z;

    return outputArr;
}

std::vector<double> Coil::computeAllAPotentialZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                                ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllAPotentialZ(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllAPotentialZ(pointVectors, defaultPrecisionCPU, computeMethod);
}

std::vector<double> Coil::computeAllAPotentialAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                                  const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllAPotentialComponents(pointVectors, usedPrecision, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = std::sqrt(computedFieldArr[i].x * computedFieldArr[i].x +
                                 computedFieldArr[i].y * computedFieldArr[i].y +
                                 computedFieldArr[i].z * computedFieldArr[i].z);

    return outputArr;
}

std::vector<double> Coil::computeAllAPotentialAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                                  ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllAPotentialAbs(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllAPotentialAbs(pointVectors, defaultPrecisionCPU, computeMethod);
}


std::vector<vec3::FieldVector3> Coil::computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                 const PrecisionArguments &usedPrecision,
                                                                 ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> output = computeAllAPotentialComponents(pointVectors, usedPrecision, computeMethod);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (auto & i : output)
        i *= frequencyFactor;

    return output;
}

std::vector<vec3::FieldVector3> Coil::computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                 ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllEFieldComponents(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllEFieldComponents(pointVectors, defaultPrecisionCPU, computeMethod);
}

std::vector<double> Coil::computeAllEFieldX(const std::vector<vec3::CoordVector3> &pointVectors,
                                            const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<double> output = computeAllAPotentialX(pointVectors, usedPrecision, computeMethod);

    for (double & i : output)
        i *= (2 * M_PI * sineFrequency);

    return output;
}

std::vector<double> Coil::computeAllEFieldX(const std::vector<vec3::CoordVector3> &pointVectors,
                                            ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllEFieldX(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllEFieldX(pointVectors, defaultPrecisionCPU, computeMethod);
}

std::vector<double> Coil::computeAllEFieldY(const std::vector<vec3::CoordVector3> &pointVectors,
                                            const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<double> output = computeAllAPotentialY(pointVectors, usedPrecision, computeMethod);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (double & i : output)
        i *= frequencyFactor;

    return output;
}

std::vector<double> Coil::computeAllEFieldY(const std::vector<vec3::CoordVector3> &pointVectors,
                                            ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllEFieldY(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllEFieldY(pointVectors, defaultPrecisionCPU, computeMethod);
}

std::vector<double> Coil::computeAllEFieldZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                            const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<double> output = computeAllAPotentialZ(pointVectors, usedPrecision, computeMethod);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (double & i : output)
        i *= frequencyFactor;

    return output;
}

std::vector<double> Coil::computeAllEFieldZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                            ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllEFieldZ(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllEFieldZ(pointVectors, defaultPrecisionCPU, computeMethod);
}

std::vector<double> Coil::computeAllEFieldAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                              const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<double> output = computeAllAPotentialAbs(pointVectors, usedPrecision, computeMethod);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (double & i : output)
        i *= frequencyFactor;

    return output;
}

std::vector<double> Coil::computeAllEFieldAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                              ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllEFieldAbs(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllEFieldAbs(pointVectors, defaultPrecisionCPU, computeMethod);
}


std::vector<vec3::Matrix3> Coil::computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &pointVectors,
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
            computedGradientArr[i] = computeBGradientTensor(pointVectors[i], usedPrecision);

        return computedGradientArr;
    }
}

std::vector<vec3::Matrix3> Coil::computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                            ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllBGradientTensors(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllBGradientTensors(pointVectors, defaultPrecisionCPU, computeMethod);
}