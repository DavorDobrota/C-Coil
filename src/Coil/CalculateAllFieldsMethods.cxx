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


std::vector<vec3::FieldVector3> Coil::calculateAllAPotentialMT(const std::vector<vec3::CoordVector3> &pointVectors,
                                                               const PrecisionArguments &usedPrecision) const
{
    std::vector<int> blockPositions = calculateChunkSize(pointVectors.size());

    std::vector<vec3::FieldVector3> computedPotentials(pointVectors.size());

    g_threadPool.setTaskCount(pointVectors.size());
    g_threadPool.getCompletedTasks().store(0ull);

    auto calcThread = [](
            int idx,
            const Coil &coil,
            const PrecisionArguments &usedPrecision,
            const std::vector<vec3::CoordVector3> &pointVectors,
            std::vector<vec3::FieldVector3> &computedPotentials,
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
                std::ref(computedPotentials),
                blockPositions[i], blockPositions[i + 1]
        );
    }

    g_threadPool.synchronizeThreads();

    return computedPotentials;
}

std::vector<vec3::FieldVector3> Coil::calculateAllBFieldMT(const std::vector<vec3::CoordVector3> &pointVectors,
                                                           const PrecisionArguments &usedPrecision) const
{
    std::vector<int> blockPositions = calculateChunkSize(pointVectors.size());

    std::vector<vec3::FieldVector3> computedFields(pointVectors.size());

    g_threadPool.setTaskCount(pointVectors.size());
    g_threadPool.getCompletedTasks().store(0ull);

    auto calcThread = [](
            int idx,
            const Coil &coil,
            const PrecisionArguments &usedPrecision,
            const std::vector<vec3::CoordVector3> &pointVectors,
            std::vector<vec3::FieldVector3> &computedFields,
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
                std::ref(computedFields),
                blockPositions[i], blockPositions[i + 1]
        );
    }
    g_threadPool.synchronizeThreads();

    return computedFields;
}

std::vector<vec3::Matrix3> Coil::calculateAllBGradientMT(const std::vector<vec3::CoordVector3> &pointVectors,
                                                         const PrecisionArguments &usedPrecision) const
{
    std::vector<int> blockPositions = calculateChunkSize(pointVectors.size());

    std::vector<vec3::Matrix3> computedGradients(pointVectors.size());

    g_threadPool.setTaskCount(pointVectors.size());
    g_threadPool.getCompletedTasks().store(0ull);

    auto calcThread = [](
            int idx,
            const Coil &coil,
            const PrecisionArguments &usedPrecision,
            const std::vector<vec3::CoordVector3> &pointVectors,
            std::vector<vec3::Matrix3> &computedGradients,
            size_t startIdx, size_t stopIdx
    ) -> void
    {
        for(size_t i = startIdx; i < stopIdx; i++)
        {
            auto result = coil.computeBGradientTensor(pointVectors[i], usedPrecision);
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
                std::ref(computedGradients),
                blockPositions[i], blockPositions[i + 1]
        );
    }
    g_threadPool.synchronizeThreads();

    return computedGradients;
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"
std::vector<vec3::FieldVector3> Coil::calculateAllAPotentialGPU(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                const PrecisionArguments &usedPrecision) const
{
    long long size = pointVectors.size();

    DataVector *coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    DataVector *resultArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));

    if(!coordinateArr || !resultArr)
        throw std::bad_alloc();

    for (long long i = 0; i < size; ++i)
    {
        vec3::CoordVector3 vector = pointVectors[i];
        vector.convertToCartesian();

        coordinateArr[i].x = vector.comp1;
        coordinateArr[i].y = vector.comp2;
        coordinateArr[i].z = vector.comp3;
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
    std::vector<vec3::FieldVector3> computedFieldArr;
    computedFieldArr.reserve(size);

    for (long long i = 0; i < pointVectors.size(); ++i)
        computedFieldArr.emplace_back(resultArr[i].x, resultArr[i].y, resultArr[i].z);

    free(resultArr);

    return computedFieldArr;
}
#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"
std::vector<vec3::FieldVector3> Coil::calculateAllBFieldGPU(const std::vector<vec3::CoordVector3> &pointVectors,
                                                            const PrecisionArguments &usedPrecision) const
{
    long long size = pointVectors.size();

    DataVector *coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    DataVector *resultArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));

    if(!coordinateArr || !resultArr)
        throw std::bad_alloc();

    for (long long i = 0; i < size; ++i)
    {
        vec3::CoordVector3 vector = pointVectors[i];
        vector.convertToCartesian();

        coordinateArr[i].x = vector.comp1;
        coordinateArr[i].y = vector.comp2;
        coordinateArr[i].z = vector.comp3;
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
    std::vector<vec3::FieldVector3> computedFieldArr;
    computedFieldArr.reserve(size);

    for (long long i = 0; i < pointVectors.size(); ++i)
        computedFieldArr.emplace_back(resultArr[i].x, resultArr[i].y, resultArr[i].z);

    free(resultArr);

    return computedFieldArr;
}
#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"
std::vector<vec3::Matrix3> Coil::calculateAllBGradientGPU(const std::vector<vec3::CoordVector3> &pointVectors,
                                                          const PrecisionArguments &usedPrecision) const
{
    long long size = pointVectors.size();

    DataVector *coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    DataMatrix *resultArr = static_cast<DataMatrix *>(calloc(size, sizeof(DataMatrix)));

    if(!coordinateArr || !resultArr)
        throw std::bad_alloc();

    for (long long i = 0; i < size; ++i)
    {
        vec3::CoordVector3 vector = pointVectors[i];
        vector.convertToCartesian();

        coordinateArr[i].x = vector.comp1;
        coordinateArr[i].y = vector.comp2;
        coordinateArr[i].z = vector.comp3;
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
    std::vector<vec3::Matrix3> computedFieldArr;
    computedFieldArr.reserve(size);

    for (long long i = 0; i < pointVectors.size(); ++i)
        computedFieldArr.emplace_back(resultArr[i].xx, resultArr[i].xy, resultArr[i].xz,
                                      resultArr[i].yx, resultArr[i].yy, resultArr[i].yz,
                                      resultArr[i].zx, resultArr[i].zy, resultArr[i].zz);

    free(resultArr);

    return computedFieldArr;
}
#pragma clang diagnostic pop

std::vector<int> Coil::calculateChunkSize(int numOps) const
{
    int average = std::floor((double)numOps / (double)threadCount);
    std::vector<int> ends(threadCount + 1);
    int remaining = numOps;
    ends[0] = 0;

    for(int i = 0; i < threadCount; i++)
    {
        int temp = (remaining % (threadCount - i) == 0 ? average : average + 1);
        ends[i+1] = (numOps - remaining) + temp;
        remaining -= temp;
    }

    return ends;
}
