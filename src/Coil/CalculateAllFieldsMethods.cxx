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
    int numOps = pointVectors.size();
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
                ends[i], ends[i + 1]
        );
    }

    g_threadPool.synchronizeThreads();

    return computedPotentials;
}

std::vector<vec3::FieldVector3> Coil::calculateAllBFieldMT(const std::vector<vec3::CoordVector3> &pointVectors,
                                                           const PrecisionArguments &usedPrecision) const
{
    int numOps = pointVectors.size();
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
            ends[i], ends[i + 1]
        );
    }
    g_threadPool.synchronizeThreads();

    return computedFields;
}

std::vector<vec3::Matrix3> Coil::calculateAllBGradientMT(const std::vector<vec3::CoordVector3> &pointVectors,
                                                         const PrecisionArguments &usedPrecision) const
{
    int numOps = pointVectors.size();
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
                ends[i], ends[i + 1]
        );
    }
    g_threadPool.synchronizeThreads();

    return computedGradients;
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"
void Coil::calculateAllAPotentialGPU(const std::vector<double> &cylindricalZArr,
                                     const std::vector<double> &cylindricalRArr,
                                     std::vector<double> &computedPotentialArr,
                                     const PrecisionArguments &usedPrecision) const
{
    computedPotentialArr.resize(cylindricalZArr.size());

    std::vector<float> z_arr(cylindricalZArr.size());
    std::vector<float> r_arr(cylindricalZArr.size());
    std::vector<float> potential_arr(cylindricalZArr.size());

    for (int i = 0; i < cylindricalZArr.size(); ++i)
    {
        z_arr[i] = cylindricalZArr[i];
        r_arr[i] = cylindricalRArr[i];
    }

    #if USE_GPU == 1
        Calculate_hardware_accelerated_a(z_arr.size(), &z_arr[0], &r_arr[0],
                                         currentDensity, innerRadius, length, thickness,
                                         16, 16, 16,
                                         &potential_arr[0]);
    #else
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    for (int i = 0; i < potential_arr.size(); ++i)
        computedPotentialArr[i] = potential_arr[i];
}
#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"
void Coil::calculateAllBFieldGPU(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 std::vector<double> &computedFieldHArr,
                                 std::vector<double> &computedFieldZArr,
                                 const PrecisionArguments &usedPrecision) const
{
    computedFieldHArr.resize(cylindricalZArr.size());
    computedFieldZArr.resize(cylindricalZArr.size());

    std::vector<float> z_arr(cylindricalZArr.size());
    std::vector<float> r_arr(cylindricalZArr.size());
    std::vector<float> fieldH_arr(cylindricalZArr.size());
    std::vector<float> fieldZ_arr(cylindricalZArr.size());

    for (int i = 0; i < cylindricalZArr.size(); ++i)
    {
        z_arr[i] = cylindricalZArr[i];
        r_arr[i] = cylindricalRArr[i];
    }

    #if USE_GPU == 1
        Calculate_hardware_accelerated_b(z_arr.size(), &z_arr[0], &r_arr[0],
                                         currentDensity, innerRadius, length, thickness,
                                         thickness/16, length/16, M_PI/48,
                                         &fieldH_arr[0], &fieldZ_arr[0]);
    #else
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    for (int i = 0; i < fieldH_arr.size(); ++i)
    {
        computedFieldHArr[i] = fieldH_arr[i];
        computedFieldZArr[i] = fieldZ_arr[i];
    }
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

int Coil::calculateChunkSize(int numOps) const
{
    if (numOps < threadCount)
        return 1;
    else if (numOps % threadCount == 0)
        return numOps / threadCount;
    else
    {
        std::vector<double> fitnessArray;
        std::vector<int> chunkArray;
        int chunkCandidate, leftover;

        int modifier = 1;
        if (numOps > 10)
            modifier = std::floor(std::log10(numOps));

        for (int i = 1; i <= std::ceil(std::log2(numOps)); ++i)
        {
            chunkCandidate = numOps / (i * modifier * threadCount);
            leftover = numOps % (i * modifier * threadCount);

            fitnessArray.push_back((double) leftover / (chunkCandidate * i));
            chunkArray.push_back(chunkCandidate);
        }
        int chunkSize = chunkArray[0];
        double chunkFitness = fitnessArray[0];

        for (int i = 1; i < chunkArray.size(); ++i)
        {
            if (fitnessArray[i] < chunkFitness)
            {
                chunkSize = chunkArray[i];
                chunkFitness = fitnessArray[i];
            }
        }
        return chunkSize;
    }
}
