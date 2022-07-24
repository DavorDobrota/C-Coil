#include "CoilGroup.h"

#include <math.h>


vec3::Vector3Array CoilGroup::computeAllAPotentialVectors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                          ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return calculateAllAPotentialGPU(pointVectors);

    else if (memberCoils.size() < 2 * threadCount || computeMethod != CPU_MT)
    {
        vec3::Vector3Array tempArr(pointVectors.size());
        vec3::Vector3Array outputArr(pointVectors.size());

        for (const auto & memberCoil : memberCoils)
        {
            tempArr = memberCoil.computeAllAPotentialVectors(pointVectors, computeMethod);
            for (int i = 0; i < pointVectors.size(); ++i)
                outputArr[i] += tempArr[i];
        }
        return outputArr;
    }
    else
        return calculateAllAPotentialMTD(pointVectors);
}


vec3::Vector3Array CoilGroup::computeAllBFieldVectors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                      ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return calculateAllBFieldGPU(pointVectors);

    else if (memberCoils.size() < 2 * threadCount || computeMethod != CPU_MT)
    {
        vec3::Vector3Array tempArr(pointVectors.size());
        vec3::Vector3Array outputArr(pointVectors.size());

        for (const auto & memberCoil : memberCoils)
        {
            tempArr = memberCoil.computeAllBFieldVectors(pointVectors, computeMethod);
            for (int i = 0; i < pointVectors.size(); ++i)
                outputArr[i] += tempArr[i];
        }
        return outputArr;
    }
    else
        return calculateAllBFieldMTD(pointVectors);
}


vec3::Vector3Array CoilGroup::computeAllEFieldVectors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                      ComputeMethod computeMethod) const
{
    if (memberCoils.size() < 2 * threadCount || computeMethod != CPU_MT)
    {
        vec3::Vector3Array tempArr(pointVectors.size());
        vec3::Vector3Array outputArr(pointVectors.size());

        for (const auto & memberCoil : memberCoils)
        {
            tempArr = memberCoil.computeAllEFieldVectors(pointVectors, computeMethod);
            for (int i = 0; i < pointVectors.size(); ++i)
                outputArr[i] += tempArr[i];
        }
        return outputArr;
    }
    else
        return calculateAllAPotentialMTD(pointVectors);
}


vec3::Matrix3Array CoilGroup::computeAllBGradientMatrices(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                  ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return calculateAllBGradientGPU(pointVectors);

    else if (memberCoils.size() < 2 * threadCount || computeMethod != CPU_MT)
    {
        vec3::Matrix3Array tempArr(pointVectors.size());
        vec3::Matrix3Array outputArr(pointVectors.size());

        for (const auto & memberCoil : memberCoils)
        {
            tempArr = memberCoil.computeAllBGradientMatrices(pointVectors, computeMethod);
            for (int i = 0; i < pointVectors.size(); ++i)
                outputArr[i] += tempArr[i];
        }
        return outputArr;
    }
    else
        return calculateAllBGradientMTD(pointVectors);
}