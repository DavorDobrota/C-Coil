#include "CoilGroup.h"
#include "ThreadPool.h"

#include <numeric>
#include <utility>
#include <functional>


namespace
{
    threadPool::ThreadPoolControl g_threadPool;
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
CoilGroup::calculateAmpereForceMTD(const Coil &secondary, PrecisionFactor precisionFactor) const
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