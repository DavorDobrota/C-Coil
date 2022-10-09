#include "Coil.h"
#include "ThreadPool.h"

#define _USE_MATH_DEFINES
#include <math.h>


namespace
{
    threadPool::ThreadPoolControl g_threadPool;
}

// TODO - fix variable so it is external and setter returned to Coil.cxx


vec3::Vector3Array Coil::calculateAllAPotentialMT(const vec3::Vector3Array &pointVectors,
                                                  const PrecisionArguments &usedPrecision) const
{
    std::vector<size_t> blockPositions = calculateChunkSize(pointVectors.size(), threadCount);

    vec3::Vector3Array computedPotentials(pointVectors.size());
    std::vector<vec3::Vector3> &tempRef = computedPotentials.getItems();

    g_threadPool.setTaskCount(threadCount);
    g_threadPool.getCompletedTasks().store(0ull);

    auto calcThread = [](
            int idx,
            const Coil &coil,
            const PrecisionArguments &usedPrecision,
            const vec3::Vector3Array &pointVectors,
            std::vector<vec3::Vector3> &computedPotentials,
            size_t startIdx, size_t stopIdx
    ) -> void
    {
        for(size_t i = startIdx; i < stopIdx; i++)
            computedPotentials[i] = coil.calculateAPotential(pointVectors[i], usedPrecision);

        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for(size_t i = 0; i < threadCount; i++)
    {
        g_threadPool.push(
            calcThread,
            std::ref(*this),
            std::ref(usedPrecision),
            std::ref(pointVectors),
            std::ref(tempRef),
            blockPositions[i], blockPositions[i + 1]
        );
    }

    g_threadPool.synchronizeThreads();

    return computedPotentials;
}

vec3::Vector3Array Coil::calculateAllBFieldMT(const vec3::Vector3Array &pointVectors,
                                              const PrecisionArguments &usedPrecision) const
{
    std::vector<size_t> blockPositions = calculateChunkSize(pointVectors.size(), threadCount);

    vec3::Vector3Array computedFields(pointVectors.size());
    std::vector<vec3::Vector3> &tempRef = computedFields.getItems();

    g_threadPool.setTaskCount(threadCount);
    g_threadPool.getCompletedTasks().store(0ull);

    auto calcThread = [](
            int idx,
            const Coil &coil,
            const PrecisionArguments &usedPrecision,
            const vec3::Vector3Array &pointVectors,
            std::vector<vec3::Vector3> &computedFields,
            size_t startIdx, size_t stopIdx
    ) -> void
    {
        for(size_t i = startIdx; i < stopIdx; i++)
            computedFields[i] = coil.calculateBField(pointVectors[i], usedPrecision);

        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for(size_t i = 0; i < threadCount; i++)
    {
        g_threadPool.push(
            calcThread,
            std::ref(*this),
            std::ref(usedPrecision),
            std::ref(pointVectors),
            std::ref(tempRef),
            blockPositions[i], blockPositions[i + 1]
        );
    }
    g_threadPool.synchronizeThreads();

    return computedFields;
}

vec3::Matrix3Array Coil::calculateAllBGradientMT(const vec3::Vector3Array &pointVectors,
                                                 const PrecisionArguments &usedPrecision) const
{
    std::vector<size_t> blockPositions = calculateChunkSize(pointVectors.size(), threadCount);

    vec3::Matrix3Array computedGradients(pointVectors.size());
    std::vector<vec3::Matrix3> &tempRef = computedGradients.getItems();

    g_threadPool.setTaskCount(threadCount);
    g_threadPool.getCompletedTasks().store(0ull);

    auto calcThread = [](
            int idx,
            const Coil &coil,
            const PrecisionArguments &usedPrecision,
            const vec3::Vector3Array &pointVectors,
            std::vector<vec3::Matrix3> &computedGradients,
            size_t startIdx, size_t stopIdx
    ) -> void
    {
        for(size_t i = startIdx; i < stopIdx; i++)
            computedGradients[i] = coil.calculateBGradient(pointVectors[i], usedPrecision);

        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for(size_t i = 0; i < threadCount; i++)
    {
        g_threadPool.push(
            calcThread,
            std::ref(*this),
            std::ref(usedPrecision),
            std::ref(pointVectors),
            std::ref(tempRef),
            blockPositions[i], blockPositions[i + 1]
        );
    }
    g_threadPool.synchronizeThreads();

    return computedGradients;
}


std::vector<size_t> Coil::calculateChunkSize(size_t opCount, int threads)
{
    size_t average = std::floor((double)opCount / (double)threads);
    std::vector<size_t> ends(threads + 1);
    size_t remaining = opCount;
    ends[0] = 0;

    for(int i = 0; i < threads; i++)
    {
        size_t temp = (remaining % (threads - i) == 0 ? average : average + 1);
        ends[i+1] = (opCount - remaining) + temp;
        remaining -= temp;
    }

    return ends;
}
