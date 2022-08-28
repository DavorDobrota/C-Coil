#include "CoilGroup.h"
#include "ThreadPool.h"

#include <functional>


namespace
{
    threadPool::ThreadPoolControl g_threadPool;
}


vec3::Vector3Array CoilGroup::calculateAllAPotentialMT(const vec3::Vector3Array &pointVectors) const
{
    std::vector<size_t> blockPositions = Coil::calculateChunkSize(pointVectors.size(), threadCount);

    g_threadPool.setTaskCount(threadCount);
    g_threadPool.getCompletedTasks().store(0ull);

    vec3::Vector3Array outputArr(pointVectors.size());

    auto calcThread = []
    (
            int idx,
            const std::vector<std::shared_ptr<Coil>> &coils,
            const vec3::Vector3Array &pointVectors,
            vec3::Vector3Array &resultArr,
            size_t startIdx, size_t stopIdx
    ){

        for (const auto& memberCoil : coils)
            for(size_t i = startIdx; i < stopIdx; i++)
                resultArr[i] += memberCoil->computeAPotentialVector(pointVectors[i]);

        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (int i = 0; i < threadCount; ++i)
    {
        g_threadPool.push
        (
            calcThread,
            std::ref(memberCoils),
            std::ref(pointVectors),
            std::ref(outputArr),
            blockPositions[i], blockPositions[i + 1]
        );
    }
    g_threadPool.synchronizeThreads();

    return outputArr;
}

vec3::Vector3Array CoilGroup::calculateAllBFieldMT(const vec3::Vector3Array &pointVectors) const
{
    std::vector<size_t> blockPositions = Coil::calculateChunkSize(pointVectors.size(), threadCount);

    g_threadPool.setTaskCount(threadCount);
    g_threadPool.getCompletedTasks().store(0ull);

    vec3::Vector3Array outputArr(pointVectors.size());

    auto calcThread = []
    (
            int idx,
            const std::vector<std::shared_ptr<Coil>> &coils,
            const vec3::Vector3Array &pointVectors,
            vec3::Vector3Array &resultArr,
            size_t startIdx, size_t stopIdx
    ){

        for (const auto& memberCoil : coils)
            for(size_t i = startIdx; i < stopIdx; i++)
                resultArr[i] += memberCoil->computeBFieldVector(pointVectors[i]);

        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (int i = 0; i < threadCount; ++i)
    {
        g_threadPool.push
        (
            calcThread,
            std::ref(memberCoils),
            std::ref(pointVectors),
            std::ref(outputArr),
            blockPositions[i], blockPositions[i + 1]
        );
    }
    g_threadPool.synchronizeThreads();

    return outputArr;
}

vec3::Vector3Array CoilGroup::calculateAllEFieldMT(const vec3::Vector3Array &pointVectors) const
{
    std::vector<size_t> blockPositions = Coil::calculateChunkSize(pointVectors.size(), threadCount);

    g_threadPool.setTaskCount(threadCount);
    g_threadPool.getCompletedTasks().store(0ull);

    vec3::Vector3Array outputArr(pointVectors.size());

    auto calcThread = []
    (
            int idx,
            const std::vector<std::shared_ptr<Coil>> &coils,
            const vec3::Vector3Array &pointVectors,
            vec3::Vector3Array &resultArr,
            size_t startIdx, size_t stopIdx
    ){

        for (const auto& memberCoil : coils)
            for(size_t i = startIdx; i < stopIdx; i++)
                resultArr[i] += memberCoil->computeEFieldVector(pointVectors[i]);

        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (int i = 0; i < threadCount; ++i)
    {
        g_threadPool.push
        (
            calcThread,
            std::ref(memberCoils),
            std::ref(pointVectors),
            std::ref(outputArr),
            blockPositions[i], blockPositions[i + 1]
        );
    }
    g_threadPool.synchronizeThreads();

    return outputArr;
}

vec3::Matrix3Array CoilGroup::calculateAllBGradientMT(const vec3::Vector3Array &pointVectors) const
{
    std::vector<size_t> blockPositions = Coil::calculateChunkSize(pointVectors.size(), threadCount);

    g_threadPool.setTaskCount(threadCount);
    g_threadPool.getCompletedTasks().store(0ull);

    vec3::Matrix3Array outputArr(pointVectors.size());

    auto calcThread = []
    (
            int idx,
            const std::vector<std::shared_ptr<Coil>> &coils,
            const vec3::Vector3Array &pointVectors,
            vec3::Matrix3Array &resultArr,
            size_t startIdx, size_t stopIdx
    ){

        for (const auto& memberCoil : coils)
            for(size_t i = startIdx; i < stopIdx; i++)
                resultArr[i] += memberCoil->computeBGradientMatrix(pointVectors[i]);

        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (int i = 0; i < threadCount; ++i)
    {
        g_threadPool.push
        (
            calcThread,
            std::ref(memberCoils),
            std::ref(pointVectors),
            std::ref(outputArr),
            blockPositions[i], blockPositions[i + 1]
        );
    }
    g_threadPool.synchronizeThreads();

    return outputArr;
}