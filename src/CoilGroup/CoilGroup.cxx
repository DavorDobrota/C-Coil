#include "CoilGroup.h"
#include "ThreadPool.h"

#include <sstream>

#define _USE_MATH_DEFINES
#include <math.h>
#include <utility>
#include <functional>

using namespace std::placeholders;


namespace
{
    threadPool::ThreadPoolControl g_threadPool;
}


CoilGroup::CoilGroup(std::vector<Coil> memberCoils, PrecisionFactor precisionFactor, int threadCount) :
            memberCoils(std::move(memberCoils)), defaultPrecisionFactor(precisionFactor), threadCount(threadCount) {}


PrecisionFactor CoilGroup::getDefaultPrecisionFactor() const { return defaultPrecisionFactor; }

int CoilGroup::getThreadCount() const { return threadCount; }

const std::vector<Coil> &CoilGroup::getMemberCoils() const { return memberCoils; }


void CoilGroup::setDefaultPrecisionFactor(PrecisionFactor precisionFactor)
{
    defaultPrecisionFactor = precisionFactor;

    for (auto & memberCoil : memberCoils)
        memberCoil.setDefaultPrecision(precisionFactor);
}

void CoilGroup::setThreadCount(int threadCount)
{
    this->threadCount = threadCount;

    for (auto & memberCoil : memberCoils)
        memberCoil.setThreadCount(threadCount);
}

void CoilGroup::addCoil(Coil coil)
{
    this->memberCoils.push_back(coil);
}


vec3::Vector3 CoilGroup::computeAPotentialVector(vec3::CoordVector3 pointVector) const
{
    vec3::Vector3 totalField = vec3::Vector3();

    for (const auto & memberCoil : memberCoils)
        totalField += memberCoil.computeAPotentialVector(pointVector);

    return totalField;
}

vec3::Vector3 CoilGroup::computeBFieldVector(vec3::CoordVector3 pointVector) const
{
    vec3::Vector3 totalField = vec3::Vector3();

    for (const auto & memberCoil : memberCoils)
        totalField += memberCoil.computeBFieldVector(pointVector);

    return totalField;
}

vec3::Vector3 CoilGroup::computeEFieldVector(vec3::CoordVector3 pointVector) const
{
    vec3::Vector3 totalField = vec3::Vector3();

    for (const auto & memberCoil : memberCoils)
        totalField += memberCoil.computeEFieldVector(pointVector);

    return totalField;
}

vec3::Matrix3 CoilGroup::computeBGradientMatrix(vec3::CoordVector3 pointVector) const
{
    vec3::Matrix3 totalGradient = vec3::Matrix3();

    for (const auto & memberCoil : memberCoils)
        totalGradient += memberCoil.computeBGradientMatrix(pointVector);

    return totalGradient;
}


double CoilGroup::computeMutualInductance(const Coil &secondary, PrecisionFactor precisionFactor, ComputeMethod computeMethod) const
{
    double totalMutualInductance = 0.0;

    if (memberCoils.size() < 2 * threadCount || computeMethod != CPU_MT)
    {
        for (const auto &memberCoil: memberCoils) {
            if (memberCoil.getId() != secondary.getId())
                totalMutualInductance += Coil::computeMutualInductance(memberCoil, secondary, precisionFactor, computeMethod);
        }
        return totalMutualInductance;
    }
    else
        return calculateMutualInductanceMTD(secondary, precisionFactor);
}

std::pair<vec3::Vector3, vec3::Vector3>
CoilGroup::computeAmpereForce(const Coil &secondary, PrecisionFactor precisionFactor, ComputeMethod computeMethod) const
{
    vec3::Vector3 totalForce{};
    vec3::Vector3 totalTorque{};
    std::pair<vec3::Vector3, vec3::Vector3> tempPair;

    if (memberCoils.size() < 2 * threadCount || computeMethod != CPU_MT)
    {
        for (const auto &memberCoil: memberCoils) {
            if (memberCoil.getId() != secondary.getId()) {
                tempPair = Coil::computeAmpereForce(memberCoil, secondary, precisionFactor, computeMethod);
                totalForce += tempPair.first;
                totalTorque += tempPair.second;
            }
        }
        return {totalForce, totalTorque};
    }
    else
        return calculateAmpereForceMTD(secondary, precisionFactor);
}


std::pair<vec3::Vector3, vec3::Vector3>
CoilGroup::computeForceOnDipoleMoment(vec3::CoordVector3 pointVector, vec3::Vector3 dipoleMoment) const
{
    vec3::Vector3 totalForce = vec3::Vector3();
    vec3::Vector3 totalTorque = vec3::Vector3();
    std::pair<vec3::Vector3, vec3::Vector3> tempPair;

    for (const auto & memberCoil : memberCoils)
    {
        tempPair = memberCoil.computeForceOnDipoleMoment(pointVector, dipoleMoment);
        totalForce += tempPair.first;
        totalTorque += tempPair.second;
    }
    return {totalForce, totalTorque};
}

CoilGroup::operator std::string() const
{
    std::stringstream output;

    auto stringifyVector = [](auto &ar) -> std::string
    {
        std::stringstream output;

        output << "[";

        for(int i = 0; i < ar.size(); i++)
        {
            if(i != 0)
                output << ", ";
            output << std::string(ar[i]);
        }

        output << "]";

        return output.str();
    };

    output << "CoilGroup("
           << "member_coils=" << stringifyVector(memberCoils)
           << ", default_precision_factor=" << std::string(defaultPrecisionFactor)
           << ", thread_count=" << threadCount
           << ")";

    return output.str();
}
