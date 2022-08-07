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

    for (auto& memberCoil : memberCoils)
        memberCoil.setDefaultPrecision(precisionFactor);
}

void CoilGroup::setThreadCount(int threadCount)
{
    this->threadCount = threadCount;

    for (auto& memberCoil : memberCoils)
        memberCoil.setThreadCount(threadCount);
}

void CoilGroup::addCoil(Coil coil)
{
    this->memberCoils.push_back(coil);
}


vec3::Vector3 CoilGroup::computeAPotentialVector(vec3::Vector3 pointVector) const
{
    vec3::Vector3 totalField = vec3::Vector3();

    for (const auto& memberCoil : memberCoils)
        totalField += memberCoil.computeAPotentialVector(pointVector);

    return totalField;
}

vec3::Vector3 CoilGroup::computeBFieldVector(vec3::Vector3 pointVector) const
{
    vec3::Vector3 totalField = vec3::Vector3();

    for (const auto& memberCoil : memberCoils)
        totalField += memberCoil.computeBFieldVector(pointVector);

    return totalField;
}

vec3::Vector3 CoilGroup::computeEFieldVector(vec3::Vector3 pointVector) const
{
    vec3::Vector3 totalField = vec3::Vector3();

    for (const auto& memberCoil : memberCoils)
        totalField += memberCoil.computeEFieldVector(pointVector);

    return totalField;
}

vec3::Matrix3 CoilGroup::computeBGradientMatrix(vec3::Vector3 pointVector) const
{
    vec3::Matrix3 totalGradient = vec3::Matrix3();

    for (const auto& memberCoil : memberCoils)
        totalGradient += memberCoil.computeBGradientMatrix(pointVector);

    return totalGradient;
}


vec3::Vector3Array CoilGroup::computeAllAPotentialVectors(const vec3::Vector3Array &pointVectors,
                                                          ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
    {
        return calculateAllAPotentialGPU(pointVectors);
    }
    else if (memberCoils.size() >= 2 * threadCount && computeMethod == CPU_MT)
    {
        return calculateAllAPotentialMTD(pointVectors);
    }
    else
    {
        vec3::Vector3Array tempArr(pointVectors.size());
        vec3::Vector3Array outputArr(pointVectors.size());

        for (const auto& memberCoil : memberCoils)
        {
            tempArr = memberCoil.computeAllAPotentialVectors(pointVectors, computeMethod);
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
    else if (memberCoils.size() < 2 * threadCount || computeMethod != CPU_MT)
    {
        return calculateAllBFieldMTD(pointVectors);
    }
    else
    {
        vec3::Vector3Array tempArr(pointVectors.size());
        vec3::Vector3Array outputArr(pointVectors.size());

        for (const auto& memberCoil : memberCoils)
        {
            tempArr = memberCoil.computeAllBFieldVectors(pointVectors, computeMethod);
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
    else if (memberCoils.size() < 2 * threadCount && computeMethod == CPU_MT)
    {
        return calculateAllEFieldMTD(pointVectors);
    }
    else
    {
        vec3::Vector3Array tempArr(pointVectors.size());
        vec3::Vector3Array outputArr(pointVectors.size());

        for (const auto& memberCoil : memberCoils)
        {
            tempArr = memberCoil.computeAllEFieldVectors(pointVectors, computeMethod);
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
    else if (memberCoils.size() >= 2 * threadCount && computeMethod == CPU_MT)
    {
        return calculateAllBGradientMTD(pointVectors);
    }
    else
    {
        vec3::Matrix3Array tempArr(pointVectors.size());
        vec3::Matrix3Array outputArr(pointVectors.size());

        for (const auto& memberCoil : memberCoils)
        {
            tempArr = memberCoil.computeAllBGradientMatrices(pointVectors, computeMethod);
            for (int i = 0; i < pointVectors.size(); ++i)
                outputArr[i] += tempArr[i];
        }
        return outputArr;
    }
}


double CoilGroup::computeMutualInductance(const Coil &secondary, PrecisionFactor precisionFactor, ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
    {
        return calculateMutualInductanceGPU(secondary, precisionFactor);
    }
    else if (memberCoils.size() >= 2 * threadCount && computeMethod == CPU_MT)
    {
        return calculateMutualInductanceMTD(secondary, precisionFactor);
    }
    else
    {
        double totalMutualInductance = 0.0;

        for (const auto &memberCoil: memberCoils)
            if (memberCoil.getId() != secondary.getId())
                totalMutualInductance += Coil::computeMutualInductance(memberCoil, secondary, precisionFactor, computeMethod);

        return totalMutualInductance;
    }
}

std::pair<vec3::Vector3, vec3::Vector3>
CoilGroup::computeAmpereForce(const Coil &secondary, PrecisionFactor precisionFactor, ComputeMethod computeMethod) const
{


    if (computeMethod == GPU)
    {
        return calculateAmpereForceGPU(secondary, precisionFactor);
    }
    else if (memberCoils.size() >= 2 * threadCount && computeMethod == CPU_MT)
    {
        return calculateAmpereForceMTD(secondary, precisionFactor);
    }
    else
    {
        vec3::Vector3 totalForce;
        vec3::Vector3 totalTorque;
        std::pair<vec3::Vector3, vec3::Vector3> tempPair;

        for (const auto &memberCoil: memberCoils)
            if (memberCoil.getId() != secondary.getId())
            {
                tempPair = Coil::computeAmpereForce(memberCoil, secondary, precisionFactor, computeMethod);
                totalForce += tempPair.first;
                totalTorque += tempPair.second;
            }

        return {totalForce, totalTorque};
    }
}


std::pair<vec3::Vector3, vec3::Vector3>
CoilGroup::computeForceOnDipoleMoment(vec3::Vector3 pointVector, vec3::Vector3 dipoleMoment) const
{
    vec3::Vector3 totalForce = vec3::Vector3();
    vec3::Vector3 totalTorque = vec3::Vector3();
    std::pair<vec3::Vector3, vec3::Vector3> tempPair;

    for (const auto& memberCoil : memberCoils)
    {
        tempPair = memberCoil.computeForceOnDipoleMoment(pointVector, dipoleMoment);
        totalForce += tempPair.first;
        totalTorque += tempPair.second;
    }
    return {totalForce, totalTorque};
}


std::vector<double> CoilGroup::computeAllMutualInductanceArrangements(Coil secondary,
                                                                      const vec3::Vector3Array &secondaryPositions,
                                                                      const std::vector<double> &secondaryYAngles,
                                                                      const std::vector<double> &secondaryZAngles,
                                                                      PrecisionFactor precisionFactor,
                                                                      ComputeMethod computeMethod) const
{
    size_t arrangementCount = secondaryPositions.size();

    if (arrangementCount == secondaryPositions.size() &&
        arrangementCount == secondaryYAngles.size() &&
        arrangementCount == secondaryZAngles.size()) {
        if (computeMethod == GPU)
        {
            return calculateAllMutualInductanceArrangementsGPU(
                secondary, secondaryPositions, secondaryYAngles, secondaryZAngles,precisionFactor
            );
        }
        else if (arrangementCount >= 2 * this->threadCount && computeMethod == CPU_MT)
        {
            return calculateAllMutualInductanceArrangementsMTD(
                secondary, secondaryPositions, secondaryYAngles, secondaryZAngles, precisionFactor
            );
        }
        else
        {
            std::vector<double> outputMInductances;
            outputMInductances.reserve(arrangementCount);

            for (int i = 0; i < arrangementCount; ++i)
            {
                secondary.setPositionAndOrientation(secondaryPositions[i], secondaryYAngles[i], secondaryZAngles[i]);
                outputMInductances.emplace_back(computeMutualInductance(secondary, precisionFactor, computeMethod));
            }
            return outputMInductances;
        }
    }
    else
        throw std::logic_error("Array sizes do not match!");
}

std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
CoilGroup::computeAllAmpereForceArrangements(Coil secondary, const vec3::Vector3Array &secondaryPositions,
                                             const std::vector<double> &secondaryYAngles,
                                             const std::vector<double> &secondaryZAngles,
                                             PrecisionFactor precisionFactor, ComputeMethod computeMethod) const
{
    size_t arrangementCount = secondaryPositions.size();

    if (arrangementCount == secondaryPositions.size() &&
        arrangementCount == secondaryYAngles.size() &&
        arrangementCount == secondaryZAngles.size()) {
        if (computeMethod == GPU)
        {
            return calculateAllAmpereForceArrangementsGPU(
                secondary, secondaryPositions, secondaryYAngles, secondaryZAngles, precisionFactor
            );
        }
        else if (arrangementCount >= 2 * this->threadCount && computeMethod == CPU_MT)
        {
            return calculateAllAmpereForceArrangementsMTD(
                secondary, secondaryPositions, secondaryYAngles, secondaryZAngles, precisionFactor
            );
        }
        else
        {
            std::vector<std::pair<vec3::Vector3, vec3::Vector3>> outputMInductances;
            outputMInductances.reserve(arrangementCount);

            for (int i = 0; i < arrangementCount; ++i)
            {
                secondary.setPositionAndOrientation(secondaryPositions[i], secondaryYAngles[i], secondaryZAngles[i]);
                outputMInductances.emplace_back(computeAmpereForce(secondary, precisionFactor, computeMethod));
            }
            return outputMInductances;
        }
    }
    else
        throw std::logic_error("Array sizes do not match!");
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
