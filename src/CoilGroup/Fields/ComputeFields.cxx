#include "CoilGroup.h"


vec3::Vector3 CoilGroup::computeAPotentialVector(vec3::Vector3 pointVector) const
{
    vec3::Vector3 totalField = vec3::Vector3();

    for (const auto& memberCoil : memberCoils)
        totalField += memberCoil->computeAPotentialVector(pointVector);

    return totalField;
}

vec3::Vector3 CoilGroup::computeBFieldVector(vec3::Vector3 pointVector) const
{
    vec3::Vector3 totalField = vec3::Vector3();

    for (const auto& memberCoil : memberCoils)
        totalField += memberCoil->computeBFieldVector(pointVector);

    return totalField;
}

vec3::Vector3 CoilGroup::computeEFieldVector(vec3::Vector3 pointVector) const
{
    vec3::Vector3 totalField = vec3::Vector3();

    for (const auto& memberCoil : memberCoils)
        totalField += memberCoil->computeEFieldVector(pointVector);

    return totalField;
}

vec3::Matrix3 CoilGroup::computeBGradientMatrix(vec3::Vector3 pointVector) const
{
    vec3::Matrix3 totalGradient = vec3::Matrix3();

    for (const auto& memberCoil : memberCoils)
        totalGradient += memberCoil->computeBGradientMatrix(pointVector);

    return totalGradient;
}


vec3::Vector3Array CoilGroup::computeAllAPotentialVectors(const vec3::Vector3Array &pointVectors,
                                                          ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
    {
        return calculateAllAPotentialGPU(pointVectors);
    }
    else if (memberCoils.size() >= threadCount && computeMethod == CPU_MT)
    {
        return calculateAllAPotentialMTD(pointVectors);
    }
    else
    {
        vec3::Vector3Array tempArr(pointVectors.size());
        vec3::Vector3Array outputArr(pointVectors.size());

        for (const auto& memberCoil : memberCoils)
        {
            tempArr = memberCoil->computeAllAPotentialVectors(pointVectors, computeMethod);
            for (int i = 0; i < pointVectors.size(); ++i)
                outputArr[i] += tempArr[i];
        }
        return outputArr;
    }
}

vec3::Vector3Array CoilGroup::computeAllBFieldVectors(const vec3::Vector3Array &pointVectors,
                                                      ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
    {
        return calculateAllBFieldGPU(pointVectors);
    }
    else if (memberCoils.size() >= threadCount && computeMethod == CPU_MT)
    {
        return calculateAllBFieldMTD(pointVectors);
    }
    else
    {
        vec3::Vector3Array tempArr(pointVectors.size());
        vec3::Vector3Array outputArr(pointVectors.size());

        for (const auto& memberCoil : memberCoils)
        {
            tempArr = memberCoil->computeAllBFieldVectors(pointVectors, computeMethod);
            for (int i = 0; i < pointVectors.size(); ++i)
                outputArr[i] += tempArr[i];
        }
        return outputArr;
    }
}

vec3::Vector3Array CoilGroup::computeAllEFieldVectors(const vec3::Vector3Array &pointVectors,
                                                      ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
    {
        return calculateAllEFieldGPU(pointVectors);
    }
    else if (memberCoils.size() >= threadCount && computeMethod == CPU_MT)
    {
        return calculateAllEFieldMTD(pointVectors);
    }
    else
    {
        vec3::Vector3Array tempArr(pointVectors.size());
        vec3::Vector3Array outputArr(pointVectors.size());

        for (const auto& memberCoil : memberCoils)
        {
            tempArr = memberCoil->computeAllEFieldVectors(pointVectors, computeMethod);
            for (int i = 0; i < pointVectors.size(); ++i)
                outputArr[i] += tempArr[i];
        }
        return outputArr;
    }
}

vec3::Matrix3Array CoilGroup::computeAllBGradientMatrices(const vec3::Vector3Array &pointVectors,
                                                          ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
    {
        return calculateAllBGradientGPU(pointVectors);
    }
    else if (memberCoils.size() >= threadCount && computeMethod == CPU_MT)
    {
        return calculateAllBGradientMTD(pointVectors);
    }
    else
    {
        vec3::Matrix3Array tempArr(pointVectors.size());
        vec3::Matrix3Array outputArr(pointVectors.size());

        for (const auto& memberCoil : memberCoils)
        {
            tempArr = memberCoil->computeAllBGradientMatrices(pointVectors, computeMethod);
            for (int i = 0; i < pointVectors.size(); ++i)
                outputArr[i] += tempArr[i];
        }
        return outputArr;
    }
}
