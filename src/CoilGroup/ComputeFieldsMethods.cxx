#include "CoilGroup.h"

#include <math.h>


std::vector<vec3::FieldVector3>
CoilGroup::computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                          ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return calculateAllAPotentialComponentsGPU(pointVectors);

    else if (memberCoils.size() < 2 * threadCount || computeMethod != CPU_MT)
    {
        std::vector<vec3::FieldVector3> tempArr(pointVectors.size());
        std::vector<vec3::FieldVector3> outputArr(pointVectors.size());

        for (const auto & memberCoil : memberCoils)
        {
            tempArr = memberCoil.computeAllAPotentialComponents(pointVectors, computeMethod);
            for (int i = 0; i < pointVectors.size(); ++i)
                outputArr[i] += tempArr[i];
        }
        return outputArr;
    }
    else
        return calculateAllAPotentialComponentsMTD(pointVectors);
}

std::vector<double> CoilGroup::computeAllAPotentialX(const std::vector<vec3::CoordVector3> &pointVectors,
                                                     ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedPotentialArr = computeAllAPotentialComponents(pointVectors, computeMethod);
    std::vector<double> outputArr(computedPotentialArr.size());

    for (int i = 0; i < computedPotentialArr.size(); ++i)
        outputArr[i] = computedPotentialArr[i].x;

    return outputArr;
}

std::vector<double> CoilGroup::computeAllAPotentialY(const std::vector<vec3::CoordVector3> &pointVectors,
                                                     ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedPotentialArr = computeAllAPotentialComponents(pointVectors, computeMethod);
    std::vector<double> outputArr(computedPotentialArr.size());

    for (int i = 0; i < computedPotentialArr.size(); ++i)
        outputArr[i] = computedPotentialArr[i].y;

    return outputArr;
}

std::vector<double> CoilGroup::computeAllAPotentialZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                                     ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedPotentialArr = computeAllAPotentialComponents(pointVectors, computeMethod);
    std::vector<double> outputArr(computedPotentialArr.size());

    for (int i = 0; i < computedPotentialArr.size(); ++i)
        outputArr[i] = computedPotentialArr[i].z;

    return outputArr;
}

std::vector<double> CoilGroup::computeAllAPotentialAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                                       ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedPotentialArr = computeAllAPotentialComponents(pointVectors, computeMethod);
    std::vector<double> outputArr(computedPotentialArr.size());

    for (int i = 0; i < computedPotentialArr.size(); ++i)
        outputArr[i] = std::sqrt(computedPotentialArr[i].x * computedPotentialArr[i].x +
                                 computedPotentialArr[i].y * computedPotentialArr[i].y +
                                 computedPotentialArr[i].z * computedPotentialArr[i].z);

    return outputArr;
}


std::vector<vec3::FieldVector3>
CoilGroup::computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                      ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return calculateAllBFieldComponentsGPU(pointVectors);

    else if (memberCoils.size() < 2 * threadCount || computeMethod != CPU_MT)
    {
        std::vector<vec3::FieldVector3> tempArr(pointVectors.size());
        std::vector<vec3::FieldVector3> outputArr(pointVectors.size());

        for (const auto & memberCoil : memberCoils)
        {
            tempArr = memberCoil.computeAllBFieldComponents(pointVectors, computeMethod);
            for (int i = 0; i < pointVectors.size(); ++i)
                outputArr[i] += tempArr[i];
        }
        return outputArr;
    }
    else
        return calculateAllBFieldComponentsMTD(pointVectors);
}

std::vector<double> CoilGroup::computeAllBFieldX(const std::vector<vec3::CoordVector3> &pointVectors,
                                                 ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(pointVectors, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].x;

    return outputArr;
}

std::vector<double> CoilGroup::computeAllBFieldY(const std::vector<vec3::CoordVector3> &pointVectors,
                                                 ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(pointVectors, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].y;

    return outputArr;
}

std::vector<double> CoilGroup::computeAllBFieldZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                                 ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(pointVectors, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].z;

    return outputArr;
}

std::vector<double> CoilGroup::computeAllBFieldAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                                   ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(pointVectors, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = std::sqrt(computedFieldArr[i].x * computedFieldArr[i].x +
                                 computedFieldArr[i].y * computedFieldArr[i].y +
                                 computedFieldArr[i].z * computedFieldArr[i].z);

    return outputArr;
}


std::vector<vec3::FieldVector3>
CoilGroup::computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                      ComputeMethod computeMethod) const
{
    if (memberCoils.size() < 2 * threadCount || computeMethod != CPU_MT)
    {
        std::vector<vec3::FieldVector3> tempArr(pointVectors.size());
        std::vector<vec3::FieldVector3> outputArr(pointVectors.size());

        for (const auto & memberCoil : memberCoils)
        {
            tempArr = memberCoil.computeAllEFieldComponents(pointVectors, computeMethod);
            for (int i = 0; i < pointVectors.size(); ++i)
                outputArr[i] += tempArr[i];
        }
        return outputArr;
    }
    else
        return calculateAllAPotentialComponentsMTD(pointVectors);
}

std::vector<double> CoilGroup::computeAllEFieldX(const std::vector<vec3::CoordVector3> &pointVectors,
                                                 ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllEFieldComponents(pointVectors, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].x;

    return outputArr;
}

std::vector<double> CoilGroup::computeAllEFieldY(const std::vector<vec3::CoordVector3> &pointVectors,
                                                 ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllEFieldComponents(pointVectors, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].y;

    return outputArr;
}

std::vector<double> CoilGroup::computeAllEFieldZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                                 ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllEFieldComponents(pointVectors, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].z;

    return outputArr;
}

std::vector<double> CoilGroup::computeAllEFieldAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                                   ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllEFieldComponents(pointVectors, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = std::sqrt(computedFieldArr[i].x * computedFieldArr[i].x +
                                 computedFieldArr[i].y * computedFieldArr[i].y +
                                 computedFieldArr[i].z * computedFieldArr[i].z);

    return outputArr;
}


std::vector<vec3::Matrix3> CoilGroup::computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                 ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return calculateAllBGradientTensorsGPU(pointVectors);

    else if (memberCoils.size() < 2 * threadCount || computeMethod != CPU_MT)
    {
        std::vector<vec3::Matrix3> tempArr(pointVectors.size());
        std::vector<vec3::Matrix3> outputArr(pointVectors.size());

        for (const auto & memberCoil : memberCoils)
        {
            tempArr = memberCoil.computeAllBGradientTensors(pointVectors, computeMethod);
            for (int i = 0; i < pointVectors.size(); ++i)
                outputArr[i] += tempArr[i];
        }
        return outputArr;
    }
    else
        return calculateAllBGradientTensorsMTD(pointVectors);
}