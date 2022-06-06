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

        return calculateAllBFieldGPU(pointVectors);
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
    return computeAllBFieldComponents(pointVectors, defaultPrecision, computeMethod);
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
    return computeAllBFieldX(pointVectors, defaultPrecision, computeMethod);
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
    return computeAllBFieldY(pointVectors, defaultPrecision, computeMethod);
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
    return computeAllBFieldZ(pointVectors, defaultPrecision, computeMethod);
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
    return computeAllBFieldAbs(pointVectors, defaultPrecision, computeMethod);
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
        return calculateAllAPotentialGPU(pointVectors);
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
    return computeAllAPotentialComponents(pointVectors, defaultPrecision, computeMethod);
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
    return computeAllAPotentialX(pointVectors, defaultPrecision, computeMethod);
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
    return computeAllAPotentialY(pointVectors, defaultPrecision, computeMethod);
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
    return computeAllAPotentialZ(pointVectors, defaultPrecision, computeMethod);
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
    return computeAllAPotentialAbs(pointVectors, defaultPrecision, computeMethod);
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
    return computeAllEFieldX(pointVectors, defaultPrecision, computeMethod);
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
    return computeAllEFieldY(pointVectors, defaultPrecision, computeMethod);
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
    return computeAllEFieldZ(pointVectors, defaultPrecision, computeMethod);
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
    return computeAllEFieldAbs(pointVectors, defaultPrecision, computeMethod);
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
    return computeAllEFieldComponents(pointVectors, defaultPrecision, computeMethod);
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
        return calculateAllBGradientGPU(pointVectors);
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
    return computeAllBGradientTensors(pointVectors, defaultPrecision, computeMethod);
}