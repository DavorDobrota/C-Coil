#include "CoilGroup.h"

#include <utility>
#include <cmath>
#include <functional>

#include "ThreadPool.h"
using namespace std::placeholders;


namespace
{
    threadPool::ThreadPoolControl g_threadPool;
}

template <class C, class...Args>
void calcThread
(
    int idx,
    std::function<C(const std::vector<vec3::CoordVector3>&)> func,
    const std::vector<vec3::CoordVector3> &pointVectorArr,
    C &outputVector
)
{
    C tempArr;
    tempArr = func(pointVectorArr);
    for (int i = 0; i < pointVectorArr.size(); ++i)
        outputVector[i] += tempArr[i];

    g_threadPool.getCompletedTasks().fetch_add(1ull);
}

CoilGroup::CoilGroup(std::vector<Coil> memberCoils, PrecisionFactor precisionFactor, int threadCount) :
            memberCoils(std::move(memberCoils)), defaultPrecisionFactor(precisionFactor), threadCount(threadCount) {}


PrecisionFactor CoilGroup::getDefaultPrecisionFactor() const { return defaultPrecisionFactor; }

int CoilGroup::getThreadCount() const { return threadCount; }

const std::vector<Coil> &CoilGroup::getMemberCoils() const { return memberCoils; }


void CoilGroup::setDefaultPrecisionFactor(PrecisionFactor precisionFactor, ComputeMethod method)
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
CoilGroup::computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                    ComputeMethod method) const
{
    // true is temporary, until to-do is completed
    if (memberCoils.size() < 2 * threadCount || method != CPU_MT || true)
    {
        std::vector<vec3::FieldVector3> tempArr(pointVectorArr.size());
        std::vector<vec3::FieldVector3> outputArr(pointVectorArr.size());

        for (const auto & memberCoil : memberCoils)
        {
            tempArr = memberCoil.computeAllBFieldComponents(pointVectorArr, method);
            for (int i = 0; i < pointVectorArr.size(); ++i)
                outputArr[i] += tempArr[i];
        }
        return outputArr;
    }
    else
        return calculateAllBFieldComponentsMTD(pointVectorArr);
}

std::vector<vec3::FieldVector3>
CoilGroup::computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                          ComputeMethod method) const
{
    // true is temporary, until to-do is completed
    if (memberCoils.size() < 2 * threadCount || method != CPU_MT || true)
    {
        std::vector<vec3::FieldVector3> tempArr(pointVectorArr.size());
        std::vector<vec3::FieldVector3> outputArr(pointVectorArr.size());

        for (const auto & memberCoil : memberCoils)
        {
            tempArr = memberCoil.computeAllAPotentialComponents(pointVectorArr, method);
            for (int i = 0; i < pointVectorArr.size(); ++i)
                outputArr[i] += tempArr[i];
        }
        return outputArr;
    }
    else
        return calculateAllAPotentialComponentsMTD(pointVectorArr);
}

std::vector<vec3::FieldVector3>
CoilGroup::computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                      ComputeMethod method) const
{
    // true is temporary, until to-do is completed
    if (memberCoils.size() < 2 * threadCount || method != CPU_MT || true)
    {
        std::vector<vec3::FieldVector3> tempArr(pointVectorArr.size());
        std::vector<vec3::FieldVector3> outputArr(pointVectorArr.size());

        for (const auto & memberCoil : memberCoils)
        {
            tempArr = memberCoil.computeAllEFieldComponents(pointVectorArr, method);
            for (int i = 0; i < pointVectorArr.size(); ++i)
                outputArr[i] += tempArr[i];
        }
        return outputArr;
    }
    else
        return calculateAllAPotentialComponentsMTD(pointVectorArr);
}

std::vector<vec3::Matrix3> CoilGroup::computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                 ComputeMethod method) const
{
    // true is temporary, until to-do is completed
    if (memberCoils.size() < 2 * threadCount || method != CPU_MT || true)
    {
        std::vector<vec3::Matrix3> tempArr(pointVectorArr.size());
        std::vector<vec3::Matrix3> outputArr(pointVectorArr.size());

        for (const auto & memberCoil : memberCoils)
        {
            tempArr = memberCoil.computeAllBGradientTensors(pointVectorArr, method);
            for (int i = 0; i < pointVectorArr.size(); ++i)
                outputArr[i] += tempArr[i];
        }
        return outputArr;
    }
    else
        return calculateAllBGradientTensorsMTD(pointVectorArr);
}


std::vector<vec3::FieldVector3>
CoilGroup::calculateAllBFieldComponentsMTD(const std::vector<vec3::CoordVector3> &pointVectorArr, bool async) const {
    g_threadPool.setTaskCount(memberCoils.size());
    g_threadPool.getCompletedTasks().store(0ull);

    std::vector<vec3::FieldVector3> outputVector(pointVectorArr.size());

    for (Coil coil: memberCoils)
    {
//        auto func = std::bind
//        (
//            static_cast<std::vector<vec3::FieldVector3>(Coil::*)(const std::vector<vec3::CoordVector3>&, ComputeMethod) const>
//            (&Coil::computeAllBFieldComponents), &coil, _1, _2
//        );

        auto func = [&coil](const std::vector<vec3::CoordVector3>& pointVectorArr)
        {
            return coil.computeAllBFieldComponents(pointVectorArr);
        };

        g_threadPool.push
        (
            calcThread,
            func,
            std::ref(pointVectorArr),
            std::ref(outputVector)
        );
    }

    if(!async)
        g_threadPool.synchronizeThreads();
}

std::vector<vec3::FieldVector3>
CoilGroup::calculateAllAPotentialComponentsMTD(const std::vector<vec3::CoordVector3> &pointVectorArr, bool async) const
{
    // TODO implement a different kind of multithreading where each thread receives ONE coil to compute
    return std::vector<vec3::FieldVector3>();
}

std::vector<vec3::FieldVector3>
CoilGroup::calculateAllEFieldComponentsMTD(const std::vector<vec3::CoordVector3> &pointVectorArr, bool async) const
{
    // TODO implement a different kind of multithreading where each thread receives ONE coil to compute
    return std::vector<vec3::FieldVector3>();
}

std::vector<vec3::Matrix3>
CoilGroup::calculateAllBGradientTensorsMTD(const std::vector<vec3::CoordVector3> &pointVectorArr, bool async) const
{
    // TODO implement a different kind of multithreading where each thread receives ONE coil to compute
    return std::vector<vec3::Matrix3>();
}


double CoilGroup::computeMutualInductance(const Coil &secondary, PrecisionFactor precisionFactor, ComputeMethod method) const
{
    double totalMutualInductance = 0.0;

    for (const auto & memberCoil : memberCoils)
    {
        if (memberCoil.getId() != secondary.getId())
            totalMutualInductance += Coil::computeMutualInductance(memberCoil, secondary, precisionFactor, method);
    }
    return totalMutualInductance;
}

std::pair<vec3::FieldVector3, vec3::FieldVector3>
CoilGroup::computeAmpereForce(const Coil &secondary, PrecisionFactor precisionFactor, ComputeMethod method) const
{
    vec3::FieldVector3 totalForce = vec3::FieldVector3();
    vec3::FieldVector3 totalTorque = vec3::FieldVector3();
    std::pair<vec3::FieldVector3, vec3::FieldVector3> tempPair;

    for (const auto & memberCoil : memberCoils)
    {
        if (memberCoil.getId() != secondary.getId())
        {
            tempPair = Coil::computeAmpereForce(memberCoil, secondary, precisionFactor, method);
            totalForce += tempPair.first;
            totalTorque += tempPair.second;
        }
    }
    return std::make_pair(totalForce, totalTorque);
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
    return std::make_pair(totalForce, totalTorque);
}