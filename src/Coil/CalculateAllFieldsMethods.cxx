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
std::vector<vec3::FieldVector3> Coil::calculateAllAPotentialGPU(const std::vector<vec3::CoordVector3> &pointVectors) const
{
    long long size = pointVectors.size();

    std::vector<DataVector> coordinateArr(size);
    std::vector<DataVector> resultArr(size);

    for (int i = 0; i < size; ++i)
    {
        vec3::CoordVector3 vector = pointVectors[i];
        vector.convertToCartesian();

        coordinateArr[i].x = vector.comp1;
        coordinateArr[i].y = vector.comp2;
        coordinateArr[i].z = vector.comp3;
    }

    CoilData coilData;
    generateCoilData(coilData);

    Calculate_hardware_accelerated_a(size, coilData, &coordinateArr[0], &resultArr[0]);

    std::vector<vec3::FieldVector3> computedFieldArr(size);

    for (int i = 0; i < pointVectors.size(); ++i)
        computedFieldArr[i] = vec3::FieldVector3(resultArr[i].x, resultArr[i].y, resultArr[i].z);

    return computedFieldArr;
}
#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"
std::vector<vec3::FieldVector3> Coil::calculateAllBFieldGPU(const std::vector<vec3::CoordVector3> &pointVectors) const
{
    long long size = pointVectors.size();

    std::vector<DataVector> coordinateArr(size);
    std::vector<DataVector> resultArr(size);

    for (int i = 0; i < size; ++i)
    {
        vec3::CoordVector3 vector = pointVectors[i];
        vector.convertToCartesian();

        coordinateArr[i].x = vector.comp1;
        coordinateArr[i].y = vector.comp2;
        coordinateArr[i].z = vector.comp3;
    }

    CoilData coilData;
    generateCoilData(coilData);

    Calculate_hardware_accelerated_b(size, coilData, &coordinateArr[0], &resultArr[0]);

    std::vector<vec3::FieldVector3> computedFieldArr(size);

    for (int i = 0; i < pointVectors.size(); ++i)
        computedFieldArr[i] = vec3::FieldVector3(resultArr[i].x, resultArr[i].y, resultArr[i].z);

    return computedFieldArr;
}
#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"
void Coil::calculateAllBGradientGPU(const std::vector<double> &cylindricalZArr,
                                    const std::vector<double> &cylindricalRArr,
                                    std::vector<double> &computedGradientRPhi,
                                    std::vector<double> &computedGradientRR,
                                    std::vector<double> &computedGradientRZ,
                                    std::vector<double> &computedGradientZZ,
                                    const PrecisionArguments &usedPrecision) const
{
    computedGradientRPhi.resize(cylindricalZArr.size());
    computedGradientRR.resize(cylindricalZArr.size());
    computedGradientRZ.resize(cylindricalZArr.size());
    computedGradientZZ.resize(cylindricalZArr.size());

    std::vector<float> z_arr(cylindricalZArr.size());
    std::vector<float> r_arr(cylindricalZArr.size());
    std::vector<float> gradRP_arr(cylindricalZArr.size());
    std::vector<float> gradRR_arr(cylindricalZArr.size());
    std::vector<float> gradRZ_arr(cylindricalZArr.size());
    std::vector<float> gradZZ_arr(cylindricalZArr.size());

    for (int i = 0; i < cylindricalZArr.size(); ++i)
    {
        z_arr[i] = cylindricalZArr[i];
        r_arr[i] = cylindricalRArr[i];
    }

    #if USE_GPU == 1
        Calculate_hardware_accelerated_g(z_arr.size(), &z_arr[0], &r_arr[0],
                                         currentDensity, innerRadius, length, thickness,
                                         thickness/16, length/16, M_PI/48,
                                         &gradRP_arr[0], &gradRR_arr[0],
                                         &gradRZ_arr[0], &gradZZ_arr[0]);
    #else
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    for (int i = 0; i < gradRP_arr.size(); ++i)
    {
        computedGradientRPhi[i] = gradRP_arr[i];
        computedGradientRR[i] = gradRR_arr[i];
        computedGradientRZ[i] = gradRZ_arr[i];
        computedGradientZZ[i] = gradZZ_arr[i];
    }
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
