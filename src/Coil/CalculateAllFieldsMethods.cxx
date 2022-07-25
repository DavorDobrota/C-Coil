#include "Coil.h"

#include "hardware_acceleration.h"
#include "ThreadPool.h"
#include "hardware_acceleration.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <numeric>


namespace
{
    threadPool::ThreadPoolControl g_threadPool;
}

// TODO - fix variable so it is external and setter returned to Coil.cxx
void Coil::setThreadCount(int threadCount)
{
    Coil::threadCount = threadCount;
    g_threadPool.setSize(threadCount);
}


vec3::Vector3Array Coil::calculateAllAPotentialMT(const vec3::Vector3Array &pointVectors,
                                                  const PrecisionArguments &usedPrecision) const
{
    std::vector<size_t> blockPositions = calculateChunkSize(pointVectors.size());

    vec3::Vector3Array computedPotentials(pointVectors.size());
    std::vector<vec3::Vector3> &tempRef = computedPotentials.getStdVectorRef();

    g_threadPool.setTaskCount(pointVectors.size());
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
        {
            auto result = coil.computeAPotentialVector(pointVectors[i], usedPrecision);
            computedPotentials[i] = result;

            g_threadPool.getCompletedTasks().fetch_add(1ull);
        }
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
    std::vector<size_t> blockPositions = calculateChunkSize(pointVectors.size());

    vec3::Vector3Array computedFields(pointVectors.size());
    std::vector<vec3::Vector3> &tempRef = computedFields.getStdVectorRef();

    g_threadPool.setTaskCount(pointVectors.size());
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
        {
            auto result = coil.computeBFieldVector(pointVectors[i], usedPrecision);
            computedFields[i] = result;

            g_threadPool.getCompletedTasks().fetch_add(1ull);
        }
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
    std::vector<size_t> blockPositions = calculateChunkSize(pointVectors.size());

    vec3::Matrix3Array computedGradients(pointVectors.size());
    std::vector<vec3::Matrix3> &tempRef = computedGradients.getStdVectorRef();

    g_threadPool.setTaskCount(pointVectors.size());
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
        {
            auto result = coil.computeBGradientMatrix(pointVectors[i], usedPrecision);
            computedGradients[i] = result;

            g_threadPool.getCompletedTasks().fetch_add(1ull);
        }
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

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"
vec3::Vector3Array Coil::calculateAllAPotentialGPU(const vec3::Vector3Array &pointVectors,
                                                   const PrecisionArguments &usedPrecision) const
{
    long long size = pointVectors.size();

    auto *coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *resultArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));

    if(!coordinateArr || !resultArr)
        throw std::bad_alloc();

    for (long long i = 0; i < size; ++i)
    {
        vec3::Vector3 tempVec = pointVectors[i];

        coordinateArr[i].x = tempVec.x;
        coordinateArr[i].y = tempVec.y;
        coordinateArr[i].z = tempVec.z;
    }

    CoilData coilData;
    generateCoilData(coilData, usedPrecision);

    #if USE_GPU == 1
        Calculate_hardware_accelerated_a(size, coilData, coordinateArr, resultArr);
    #else
        free(coordinateArr);
        free(resultArr);
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    free(coordinateArr);

    vec3::Vector3Array computedPotentialArr;
    computedPotentialArr.reserve(size);

    for (long long i = 0; i < pointVectors.size(); ++i)
        computedPotentialArr.append(resultArr[i].x, resultArr[i].y, resultArr[i].z);

    free(resultArr);

    return computedPotentialArr;
}
#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"
vec3::Vector3Array Coil::calculateAllBFieldGPU(const vec3::Vector3Array &pointVectors,
                                               const PrecisionArguments &usedPrecision) const
{
    long long size = pointVectors.size();

    auto *coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *resultArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));

    if(!coordinateArr || !resultArr)
        throw std::bad_alloc();

    for (long long i = 0; i < size; ++i)
    {
        vec3::Vector3 tempVec = pointVectors[i];

        coordinateArr[i].x = tempVec.x;
        coordinateArr[i].y = tempVec.y;
        coordinateArr[i].z = tempVec.z;
    }

    CoilData coilData;
    generateCoilData(coilData, usedPrecision);

    #if USE_GPU == 1
        Calculate_hardware_accelerated_b(size, coilData, coordinateArr, resultArr);
    #else
        free(coordinateArr);
        free(resultArr);
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    free(coordinateArr);

    vec3::Vector3Array computedFieldArr;
    computedFieldArr.reserve(size);

    for (long long i = 0; i < pointVectors.size(); ++i)
        computedFieldArr.append(resultArr[i].x, resultArr[i].y, resultArr[i].z);

    free(resultArr);

    return computedFieldArr;
}
#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"
vec3::Matrix3Array Coil::calculateAllBGradientGPU(const vec3::Vector3Array &pointVectors,
                                                  const PrecisionArguments &usedPrecision) const
{
    long long size = pointVectors.size();

    auto *coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *resultArr = static_cast<DataMatrix *>(calloc(size, sizeof(DataMatrix)));

    if(!coordinateArr || !resultArr)
        throw std::bad_alloc();

    for (long long i = 0; i < size; ++i)
    {
        vec3::Vector3 tempVec = pointVectors[i];

        coordinateArr[i].x = tempVec.x;
        coordinateArr[i].y = tempVec.y;
        coordinateArr[i].z = tempVec.z;
    }

    CoilData coilData;
    generateCoilData(coilData, usedPrecision);

    #if USE_GPU == 1
        Calculate_hardware_accelerated_g(size, coilData, coordinateArr, resultArr);
    #else
        free(coordinateArr);
        free(resultArr);
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    free(coordinateArr);

    vec3::Matrix3Array computedGradientArr;
    computedGradientArr.reserve(size);

    for (long long i = 0; i < pointVectors.size(); ++i)
        computedGradientArr.append(resultArr[i].xx, resultArr[i].xy, resultArr[i].xz,
                                   resultArr[i].yx, resultArr[i].yy, resultArr[i].yz,
                                   resultArr[i].zx, resultArr[i].zy, resultArr[i].zz);

    free(resultArr);

    return computedGradientArr;
}
#pragma clang diagnostic pop

std::vector<size_t> Coil::calculateChunkSize(size_t numOps) const
{
    size_t average = std::floor((double)numOps / (double)threadCount);
    std::vector<size_t> ends(threadCount + 1);
    size_t remaining = numOps;
    ends[0] = 0;

    for(int i = 0; i < threadCount; i++)
    {
        size_t temp = (remaining % (threadCount - i) == 0 ? average : average + 1);
        ends[i+1] = (numOps - remaining) + temp;
        remaining -= temp;
    }

    return ends;
}
