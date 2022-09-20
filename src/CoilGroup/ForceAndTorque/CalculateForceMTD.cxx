#include "CoilGroup.h"
#include "ThreadPool.h"

#include <utility>
#include <functional>


namespace
{
    threadPool::ThreadPoolControl g_threadPool;
}


std::pair<vec3::Vector3, vec3::Vector3>
CoilGroup::calculateForceTorqueMTD(const Coil &secondary, PrecisionFactor precisionFactor) const
{
    g_threadPool.setTaskCount(memberCoils.size());
    g_threadPool.getCompletedTasks().store(0ull);

    std::vector<std::pair<vec3::Vector3, vec3::Vector3>> intermediateValues(memberCoils.size());

    auto calcThread = []
    (
            int idx,
            const Coil &coil,
            const Coil &secondary,
            PrecisionFactor precisionFactor,
            std::pair<vec3::Vector3, vec3::Vector3> &ampereForce
    ){
        ampereForce = Coil::computeForceTorque(coil, secondary, precisionFactor);

        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (int i = 0; i < memberCoils.size(); i++)
    {
        if (memberCoils[i]->getId() != secondary.getId())
        {
            g_threadPool.push
            (
                calcThread,
                std::ref(*memberCoils[i]),
                std::ref(secondary),
                precisionFactor,
                std::ref(intermediateValues[i])
            );
        }
    }
    g_threadPool.synchronizeThreads();

    vec3::Vector3 force{}, torque{};

    for(auto value : intermediateValues)
    {
        force += value.first;
        torque += value.second;
    }

    return {force, torque};
}

std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
CoilGroup::calculateAllForceTorqueArrangementsMTD(const Coil &secondary, const vec3::Vector3Array &secondaryPositions,
                                                  const std::vector<double> &secondaryYAngles,
                                                  const std::vector<double> &secondaryZAngles,
                                                  PrecisionFactor precisionFactor) const
{
    size_t arrangementCount = secondaryPositions.size();

    std::vector<std::pair<vec3::Vector3, vec3::Vector3>> outputMInductances(arrangementCount);

    g_threadPool.setTaskCount(arrangementCount);
    g_threadPool.getCompletedTasks().store(0ull);

    auto calcThread = []
    (
            int idx,
            const std::vector<std::shared_ptr<Coil>> &group,
            const Coil &secondary,
            vec3::Vector3 secondaryPosition,
            double secondaryYAngle,
            double secondaryZAngle,
            PrecisionFactor precisionFactor,
            std::pair<vec3::Vector3, vec3::Vector3> &forceTorqueArr
    ){
        std::pair<vec3::Vector3, vec3::Vector3> forceTorque;

        Coil sec = Coil(secondary);
        sec.setPositionAndOrientation(secondaryPosition, secondaryYAngle, secondaryZAngle);

        for (const auto& memberCoil : group)
            if (memberCoil->getId() != secondary.getId())
            {
                std::pair<vec3::Vector3, vec3::Vector3> tempPair;
                tempPair = Coil::computeForceTorque(*memberCoil, sec, precisionFactor, CPU_ST);
                forceTorque.first += tempPair.first;
                forceTorque.second += tempPair.second;
            }

        forceTorqueArr = forceTorque;

        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (int i = 0; i < arrangementCount; ++i)
    {
        g_threadPool.push
        (
            calcThread,
            std::ref(this->memberCoils),
            std::ref(secondary),
            secondaryPositions[i],
            secondaryYAngles[i],
            secondaryZAngles[i],
            precisionFactor,
            std::ref(outputMInductances[i])
        );
    }
    g_threadPool.synchronizeThreads();

    return outputMInductances;
}