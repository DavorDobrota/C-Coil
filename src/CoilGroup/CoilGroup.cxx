#include <utility>

#include "CoilGroup.h"


CoilGroup::CoilGroup(std::vector<Coil> memberCoils, PrecisionFactor precisionFactor, int threadCount) :
        memberCoils(std::move(memberCoils)), threadCount(threadCount)
{

}


PrecisionArguments CoilGroup::getDefaultPrecision() const { return defaultPrecision; }

int CoilGroup::getThreadCount() const { return threadCount; }

const std::vector<Coil> &CoilGroup::getMemberCoils() const { return memberCoils; }


void CoilGroup::setDefaultPrecision(const PrecisionArguments &defaultPrecision)
{
    this->defaultPrecision = defaultPrecision;
}

void CoilGroup::setDefaultPrecision(PrecisionFactor precisionFactor, ComputeMethod method)
{
    if (method == GPU)
    {
        for (auto & memberCoil : memberCoils)
            memberCoil.setDefaultPrecision(PrecisionArguments::getCoilPrecisionArgumentsGPU(memberCoil, precisionFactor));
    }
    else
        for (auto & memberCoil : memberCoils)
            memberCoil.setDefaultPrecision(PrecisionArguments::getCoilPrecisionArgumentsCPU(memberCoil, precisionFactor));
}

void CoilGroup::setThreadCount(int threadCount)
{
    this->threadCount = threadCount;
}

void CoilGroup::addCoil(Coil coil)
{
    this->memberCoils.push_back(coil);
}


std::vector<vec3::FieldVector3>
CoilGroup::computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                      const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    if (memberCoils.size() < 2 * threadCount || method != CPU_MT)
    {
        std::vector<vec3::FieldVector3> tempArr(pointVectorArr.size());
        std::vector<vec3::FieldVector3> outputArr(pointVectorArr.size());

        for (const auto & memberCoil : memberCoils)
        {
            tempArr = memberCoil.computeAllBFieldComponents(pointVectorArr, usedPrecision, method);
            for (int i = 0; i < pointVectorArr.size(); ++i)
                outputArr[i] += tempArr[i];
        }
        return outputArr;
    }
    else
        return calculateAllBFieldComponentsMTD(pointVectorArr, usedPrecision);
}

std::vector<vec3::FieldVector3>
CoilGroup::computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                    ComputeMethod method) const
{
    return computeAllBFieldComponents(pointVectorArr, defaultPrecision, method);
}

std::vector<vec3::FieldVector3>
CoilGroup::calculateAllBFieldComponentsMTD(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                           const PrecisionArguments &usedPrecision) const
{
    // TODO implement a different kind of multithreading where each thread receives ONE coil to compute
    return std::vector<vec3::FieldVector3>();
}
