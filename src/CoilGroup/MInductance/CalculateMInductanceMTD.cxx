#include "CoilGroup.h"
#include "ThreadPool.h"

#include <numeric>
#include <functional>


namespace
{
    threadPool::ThreadPoolControl g_threadPool;
}


double CoilGroup::calculateMutualInductanceMTD(const Coil &secondary, PrecisionFactor precisionFactor) const
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
    ){
        mutualInductance = Coil::computeMutualInductance(coil, secondary, precisionFactor);

        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (int i = 0; i < memberCoils.size(); i++)
    {
        if (memberCoils[i].getId() != secondary.getId()) {
            g_threadPool.push
            (
                calcThread,
                std::ref(memberCoils[i]),
                std::ref(secondary),
                precisionFactor,
                std::ref(intermediateValues[i])
            );
        }
    }

    g_threadPool.synchronizeThreads();

    double mutualInductance = std::accumulate(intermediateValues.begin(), intermediateValues.end(), 0.0);

    return mutualInductance;
}


std::vector<double> CoilGroup::calculateAllMutualInductanceArrangementsMTD(Coil secondary,
                                                                           const vec3::Vector3Array &secondaryPositions,
                                                                           const std::vector<double> &secondaryYAngles,
                                                                           const std::vector<double> &secondaryZAngles,
                                                                           PrecisionFactor precisionFactor) const
{
    size_t arrangementCount = secondaryPositions.size();

    std::vector<double> outputMInductances(arrangementCount);

    g_threadPool.setTaskCount(arrangementCount);
    g_threadPool.getCompletedTasks().store(0ull);

    auto calcThread = []
    (
            int idx,
            const std::vector<Coil> &group,
            Coil secondary,
            vec3::Vector3 secondaryPosition,
            double secondaryYAngle,
            double secondaryZAngle,
            PrecisionFactor precisionFactor,
            double &mutualInductance
    ){
        double inductance;

        secondary.setPositionAndOrientation(secondaryPosition, secondaryYAngle, secondaryZAngle);

        for (const auto & memberCoil : group)
            if (memberCoil.getId() != secondary.getId())
                inductance += Coil::computeMutualInductance(memberCoil, secondary, precisionFactor, CPU_ST);

        mutualInductance = inductance;

        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (int i = 0; i < arrangementCount; ++i)
    {
        g_threadPool.push
        (
            calcThread,
            std::ref(this->memberCoils),
            secondary,
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
