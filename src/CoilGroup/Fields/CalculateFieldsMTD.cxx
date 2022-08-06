#include "CoilGroup.h"
#include "ThreadPool.h"

#include <functional>


namespace
{
    threadPool::ThreadPoolControl g_threadPool;
}


vec3::Vector3Array CoilGroup::calculateAllBFieldMTD(const vec3::Vector3Array &pointVectors) const
{
    g_threadPool.setTaskCount(memberCoils.size());
    g_threadPool.getCompletedTasks().store(0ull);

    std::vector<std::vector<vec3::Vector3>> intermediateValues(memberCoils.size());
    for(auto vec : intermediateValues)
        vec.resize(pointVectors.size());

    auto calcThread = []
    (
            int idx,
            const Coil &coil,
            const vec3::Vector3Array &pointVectors,
            std::vector<vec3::Vector3> &outputVector
    ){
        outputVector = coil.computeAllBFieldVectors(pointVectors).getItems();

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

    vec3::Vector3Array output(pointVectors.size());

    for(auto values : intermediateValues)
        for(int i = 0; i < values.size(); i++)
            output[i] += values[i];

    return output;
}

vec3::Vector3Array CoilGroup::calculateAllAPotentialMTD(const vec3::Vector3Array &pointVectors) const
{
    g_threadPool.setTaskCount(memberCoils.size());
    g_threadPool.getCompletedTasks().store(0ull);

    std::vector<std::vector<vec3::Vector3>> intermediateValues(memberCoils.size());
    for(auto vec : intermediateValues)
        vec.resize(pointVectors.size());

    auto calcThread = []
    (
            int idx,
            const Coil &coil,
            const vec3::Vector3Array &pointVectors,
            std::vector<vec3::Vector3> &outputVector
    ){
        outputVector = coil.computeAllAPotentialVectors(pointVectors).getItems();

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

    vec3::Vector3Array output(pointVectors.size());

    for(auto values : intermediateValues)
        for(int i = 0; i < values.size(); i++)
            output[i] += values[i];

    return output;
}

vec3::Vector3Array CoilGroup::calculateAllEFieldMTD(const vec3::Vector3Array &pointVectors) const
{
    g_threadPool.setTaskCount(memberCoils.size());
    g_threadPool.getCompletedTasks().store(0ull);

    std::vector<std::vector<vec3::Vector3>> intermediateValues(memberCoils.size());
    for(auto vec : intermediateValues)
        vec.resize(pointVectors.size());

    auto calcThread = []
    (
            int idx,
            const Coil &coil,
            const vec3::Vector3Array &pointVectors,
            std::vector<vec3::Vector3> &outputVector
    ){
        outputVector = coil.computeAllEFieldVectors(pointVectors).getItems();

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

    vec3::Vector3Array output(pointVectors.size());

    for(auto values : intermediateValues)
        for(int i = 0; i < values.size(); i++)
            output[i] += values[i];

    return output;
}

vec3::Matrix3Array CoilGroup::calculateAllBGradientMTD(const vec3::Vector3Array &pointVectors) const
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
            const vec3::Vector3Array &pointVectors,
            std::vector<vec3::Matrix3> &outputVector
    ){
        outputVector = coil.computeAllBGradientMatrices(pointVectors).getItems();

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

    vec3::Matrix3Array output(pointVectors.size());

    for(auto values : intermediateValues)
        for(int i = 0; i < values.size(); i++)
            output[i] += values[i];

    return output;
}