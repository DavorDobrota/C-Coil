#include "CoilGroup.h"
#include "ThreadPool.h"

#include <sstream>

#define _USE_MATH_DEFINES
#include <math.h>
#include <numeric>
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


void CoilGroup::setDefaultPrecisionFactor(PrecisionFactor precisionFactor, ComputeMethod computeMethod)
{
    defaultPrecisionFactor = precisionFactor;

    if (computeMethod == GPU)
    {
        for (auto & memberCoil : memberCoils)
            memberCoil.setDefaultPrecisionCPU(
                    PrecisionArguments::getCoilPrecisionArgumentsGPU(memberCoil, precisionFactor));
    }
    else
        for (auto & memberCoil : memberCoils)
            memberCoil.setDefaultPrecisionCPU(
                    PrecisionArguments::getCoilPrecisionArgumentsCPU(memberCoil, precisionFactor));
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


vec3::FieldVector3 CoilGroup::computeBFieldVector(vec3::CoordVector3 pointVector) const
{
    vec3::FieldVector3 totalField = vec3::FieldVector3();

    for (const auto & memberCoil : memberCoils)
        totalField += memberCoil.computeBFieldVector(pointVector);

    return totalField;
}

vec3::FieldVector3 CoilGroup::computeAPotentialVector(vec3::CoordVector3 pointVector) const
{
    vec3::FieldVector3 totalField = vec3::FieldVector3();

    for (const auto & memberCoil : memberCoils)
        totalField += memberCoil.computeAPotentialVector(pointVector);

    return totalField;
}

vec3::FieldVector3 CoilGroup::computeEFieldVector(vec3::CoordVector3 pointVector) const
{
    vec3::FieldVector3 totalField = vec3::FieldVector3();

    for (const auto & memberCoil : memberCoils)
        totalField += memberCoil.computeEFieldVector(pointVector);

    return totalField;
}

vec3::Matrix3 CoilGroup::computeBGradientTensor(vec3::CoordVector3 pointVector) const
{
    vec3::Matrix3 totalGradient = vec3::Matrix3();

    for (const auto & memberCoil : memberCoils)
        totalGradient += memberCoil.computeBGradientTensor(pointVector);

    return totalGradient;
}

std::vector<vec3::FieldVector3>
CoilGroup::computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                    ComputeMethod computeMethod) const
{
    if (memberCoils.size() < 4 * threadCount || computeMethod != CPU_MT)
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

std::vector<vec3::FieldVector3>
CoilGroup::computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                          ComputeMethod computeMethod) const
{
    if (memberCoils.size() < 4 * threadCount || computeMethod != CPU_MT)
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

std::vector<vec3::FieldVector3>
CoilGroup::computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                      ComputeMethod computeMethod) const
{
    if (memberCoils.size() < 4 * threadCount || computeMethod != CPU_MT)
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

std::vector<vec3::Matrix3> CoilGroup::computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                 ComputeMethod computeMethod) const
{
    if (memberCoils.size() < 4 * threadCount || computeMethod != CPU_MT)
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

std::vector<vec3::FieldVector3>
CoilGroup::calculateAllBFieldComponentsMTD(const std::vector<vec3::CoordVector3> &pointVectors) const
{
    g_threadPool.setTaskCount(memberCoils.size());
    g_threadPool.getCompletedTasks().store(0ull);

    std::vector<std::vector<vec3::FieldVector3>> intermediateValues(memberCoils.size());
    for(auto vec : intermediateValues)
        vec.resize(pointVectors.size());

    auto calcThread = []
    (
        int idx,
        const Coil &coil,
        const std::vector<vec3::CoordVector3> &pointVectors,
        std::vector<vec3::FieldVector3> &outputVector
    )
    {
        outputVector = coil.computeAllBFieldComponents(pointVectors);

        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (int i = 0; i < memberCoils.size(); i++)
    {
        g_threadPool.push
        (
            calcThread,
            std::ref(memberCoils[i]),
            std::ref(pointVectors),
            std::ref(intermediateValues[i])
        );
    }

    g_threadPool.synchronizeThreads();

    std::vector<vec3::FieldVector3> ret(pointVectors.size());

    for(auto values : intermediateValues)
    {
        for(int i = 0; i < values.size(); i++)
        {
            ret[i] += values[i];
        }
    }

    return ret;
}

std::vector<vec3::FieldVector3>
CoilGroup::calculateAllAPotentialComponentsMTD(const std::vector<vec3::CoordVector3> &pointVectors) const
{
    g_threadPool.setTaskCount(memberCoils.size());
    g_threadPool.getCompletedTasks().store(0ull);

    std::vector<std::vector<vec3::FieldVector3>> intermediateValues(memberCoils.size());
    for(auto vec : intermediateValues)
        vec.resize(pointVectors.size());

    auto calcThread = []
    (
        int idx,
        const Coil &coil,
        const std::vector<vec3::CoordVector3> &pointVectors,
        std::vector<vec3::FieldVector3> &outputVector
    )
    {
        outputVector = coil.computeAllAPotentialComponents(pointVectors);

        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (int i = 0; i < memberCoils.size(); i++)
    {
        g_threadPool.push
        (
            calcThread,
            std::ref(memberCoils[i]),
            std::ref(pointVectors),
            std::ref(intermediateValues[i])
        );
    }

    g_threadPool.synchronizeThreads();

    std::vector<vec3::FieldVector3> ret(pointVectors.size());

    for(auto values : intermediateValues)
    {
        for(int i = 0; i < values.size(); i++)
        {
            ret[i] += values[i];
        }
    }

    return ret;
}

std::vector<vec3::FieldVector3>
CoilGroup::calculateAllEFieldComponentsMTD(const std::vector<vec3::CoordVector3> &pointVectors) const
{
    g_threadPool.setTaskCount(memberCoils.size());
    g_threadPool.getCompletedTasks().store(0ull);

    std::vector<std::vector<vec3::FieldVector3>> intermediateValues(memberCoils.size());
    for(auto vec : intermediateValues)
        vec.resize(pointVectors.size());

    auto calcThread = []
    (
        int idx,
        const Coil &coil,
        const std::vector<vec3::CoordVector3> &pointVectors,
        std::vector<vec3::FieldVector3> &outputVector
    )
    {
        outputVector = coil.computeAllEFieldComponents(pointVectors);

        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (int i = 0; i < memberCoils.size(); i++)
    {
        g_threadPool.push
        (
            calcThread,
            std::ref(memberCoils[i]),
            std::ref(pointVectors),
            std::ref(intermediateValues[i])
        );
    }

    g_threadPool.synchronizeThreads();

    std::vector<vec3::FieldVector3> ret(pointVectors.size());

    for(auto values : intermediateValues)
    {
        for(int i = 0; i < values.size(); i++)
        {
            ret[i] += values[i];
        }
    }

    return ret;
}

std::vector<vec3::Matrix3>
CoilGroup::calculateAllBGradientTensorsMTD(const std::vector<vec3::CoordVector3> &pointVectors) const
{
    g_threadPool.setTaskCount(memberCoils.size());
    g_threadPool.getCompletedTasks().store(0ull);

    std::vector<std::vector<vec3::Matrix3>> intermediateValues(memberCoils.size());
    for(auto vec : intermediateValues)
        vec.resize(pointVectors.size());

    auto calcThread = []
    (
        int idx,
        const Coil &coil,
        const std::vector<vec3::CoordVector3> &pointVectors,
        std::vector<vec3::Matrix3> &outputVector
    )
    {
        outputVector = coil.computeAllBGradientTensors(pointVectors);

        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (int i = 0; i < memberCoils.size(); i++)
    {
        g_threadPool.push
        (
            calcThread,
            std::ref(memberCoils[i]),
            std::ref(pointVectors),
            std::ref(intermediateValues[i])
        );
    }

    g_threadPool.synchronizeThreads();

    std::vector<vec3::Matrix3> ret(pointVectors.size());

    for(auto values : intermediateValues)
    {
        for(int i = 0; i < values.size(); i++)
        {
            ret[i] += values[i];
        }
    }

    return ret;
}


double CoilGroup::computeMutualInductance(const Coil &secondary, PrecisionFactor precisionFactor, ComputeMethod computeMethod) const
{
    double totalMutualInductance = 0.0;

    if (memberCoils.size() < 4 * threadCount || computeMethod != CPU_MT)
    {
        for (const auto &memberCoil: memberCoils) {
            if (memberCoil.getId() != secondary.getId())
                totalMutualInductance += Coil::computeMutualInductance(memberCoil, secondary, precisionFactor, computeMethod);
        }
        return totalMutualInductance;
    }
    else
        return computeMutualInductanceMTD(secondary, precisionFactor);
}

std::pair<vec3::FieldVector3, vec3::FieldVector3>
CoilGroup::computeAmpereForce(const Coil &secondary, PrecisionFactor precisionFactor, ComputeMethod computeMethod) const
{
    vec3::FieldVector3 totalForce{};
    vec3::FieldVector3 totalTorque{};
    std::pair<vec3::FieldVector3, vec3::FieldVector3> tempPair;

    if (memberCoils.size() < 4 * threadCount || computeMethod != CPU_MT)
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
        return computeAmpereForceMTD(secondary, precisionFactor);
}

double CoilGroup::computeMutualInductanceMTD(const Coil &secondary, PrecisionFactor precisionFactor) const
{
    g_threadPool.setTaskCount(memberCoils.size());
    g_threadPool.getCompletedTasks().store(0ull);

    std::vector<double> intermediateValues(memberCoils.size());

    auto calcThread = []
    (
        int idx,
        const Coil &coil,
        const Coil &secondary,
        PrecisionFactor precisionFactor,
        double &mutualInductance
    )
    {
        mutualInductance = Coil::computeMutualInductance(coil, secondary, precisionFactor);

        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (int i = 0; i < memberCoils.size(); i++)
    {
        if (memberCoils[i].getId() != secondary.getId())
            g_threadPool.push
            (
                calcThread,
                std::ref(memberCoils[i]),
                std::ref(secondary),
                precisionFactor,
                std::ref(intermediateValues[i])
            );
    }

    g_threadPool.synchronizeThreads();

    double mutualInductance = std::accumulate(intermediateValues.begin(), intermediateValues.end(), 0.0);

    return mutualInductance;
}

std::pair<vec3::FieldVector3, vec3::FieldVector3>
CoilGroup::computeAmpereForceMTD(const Coil &secondary, PrecisionFactor precisionFactor) const
{
    g_threadPool.setTaskCount(memberCoils.size());
    g_threadPool.getCompletedTasks().store(0ull);

    std::vector<std::pair<vec3::FieldVector3, vec3::FieldVector3>> intermediateValues(memberCoils.size());

    auto calcThread = []
    (
        int idx,
        const Coil &coil,
        const Coil &secondary,
        PrecisionFactor precisionFactor,
        std::pair<vec3::FieldVector3, vec3::FieldVector3> &ampereForce
    )
    {
        ampereForce = Coil::computeAmpereForce(coil, secondary, precisionFactor);

        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (int i = 0; i < memberCoils.size(); i++)
    {
        if (memberCoils[i].getId() != secondary.getId())
            g_threadPool.push
            (
                calcThread,
                std::ref(memberCoils[i]),
                std::ref(secondary),
                precisionFactor,
                std::ref(intermediateValues[i])
            );
    }

    g_threadPool.synchronizeThreads();

    vec3::FieldVector3 force{}, torque{};

    for(auto value : intermediateValues)
    {
        force += value.first;
        torque += value.second;
    }

    return {force, torque};
}

std::pair<vec3::FieldVector3, vec3::FieldVector3>
CoilGroup::computeForceOnDipoleMoment(vec3::CoordVector3 pointVector, vec3::FieldVector3 dipoleMoment) const
{
    vec3::FieldVector3 totalForce = vec3::FieldVector3();
    vec3::FieldVector3 totalTorque = vec3::FieldVector3();
    std::pair<vec3::FieldVector3, vec3::FieldVector3> tempPair;

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
