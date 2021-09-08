#include "Coil.h"
#include "ctpl.h"
#include "hardware_acceleration.h"

#include <cmath>
#include <algorithm>
#include <cstdio>

namespace
{
    ctpl::thread_pool g_threadPool;
}

// TODO - fix variable so it is external and setter returned to Coil.cxx
void Coil::setThreadCount(int threadCount)
{
    Coil::threadCount = threadCount;
    g_threadPool.resize(threadCount);
}

void Coil::calculateAllBFieldST(const std::vector<double> &cylindricalZArr,
                                const std::vector<double> &cylindricalRArr,
                                std::vector<double> &computedFieldHArr,
                                std::vector<double> &computedFieldZArr,
                                const PrecisionArguments &usedPrecision) const
{
    computedFieldHArr.resize(cylindricalZArr.size());
    computedFieldZArr.resize(cylindricalZArr.size());

    for (int i = 0; i < cylindricalZArr.size(); ++i)
    {
        std::pair values = calculateBField(cylindricalZArr[i], cylindricalRArr[i], usedPrecision);
        computedFieldHArr[i] = values.first;
        computedFieldZArr[i] = values.second;
    }
}

void Coil::calculateAllAPotentialST(const std::vector<double> &cylindricalZArr,
                                    const std::vector<double> &cylindricalRArr,
                                    std::vector<double> &computedPotentialArr,
                                    const PrecisionArguments &usedPrecision) const
{
    computedPotentialArr.resize(cylindricalZArr.size());

    for (int i = 0; i < cylindricalZArr.size(); ++i)
        computedPotentialArr[i] = calculateAPotential(cylindricalZArr[i], cylindricalRArr[i], usedPrecision);
}

void Coil::calculateAllBGradientST(const std::vector<double> &cylindricalZArr,
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

    for (int i = 0; i < cylindricalZArr.size(); ++i)
    {
        std::vector<double> gradient = calculateBGradient(cylindricalZArr[i], cylindricalRArr[i], usedPrecision);
        computedGradientRPhi[i] = gradient[0];
        computedGradientRR[i] = gradient[1];
        computedGradientRZ[i] = gradient[2];
        computedGradientZZ[i] = gradient[3];
    }
}

void Coil::calculateAllBFieldMT(const std::vector<double> &cylindricalZArr,
                                const std::vector<double> &cylindricalRArr,
                                std::vector<double> &computedFieldHArr,
                                std::vector<double> &computedFieldZArr,
                                const PrecisionArguments &usedPrecision,
                                int chunkSize,
                                bool async) const
{
    computedFieldHArr.resize(cylindricalZArr.size());
    computedFieldZArr.resize(cylindricalZArr.size());

    auto calcThread = [](
            int idx,
            const Coil &coil,
            const PrecisionArguments &usedPrecision,
            const std::vector<double> &cylindricalZ,
            const std::vector<double> &cylindricalR,
            std::vector<double> &computedFieldH,
            std::vector<double> &computedFieldZ,
            size_t startIdx, size_t stopIdx
    ) -> void
    {
        for(size_t i = startIdx; i < stopIdx; i++)
        {
            auto result = coil.calculateBField(cylindricalZ[i], cylindricalR[i], usedPrecision);
            computedFieldH[i] = result.first;
            computedFieldZ[i] = result.second;
        }
    };

    for(size_t i = 0; i < (size_t)std::ceil((double)cylindricalZArr.size() / (double)chunkSize); i++)
    {
        g_threadPool.push(
            calcThread,
            *this,
            usedPrecision,
            cylindricalZArr, cylindricalRArr,
            computedFieldHArr, computedFieldZArr,
            i * chunkSize, std::min((i + 1) * chunkSize, cylindricalZArr.size())
        );
    }

    if(!async)
        synchronizeThreads();
}

void Coil::calculateAllAPotentialMT(const std::vector<double> &cylindricalZArr,
                                    const std::vector<double> &cylindricalRArr,
                                    std::vector<double> &computedPotentialArr,
                                    const PrecisionArguments &usedPrecision,
                                    int chunkSize, bool async) const
{
    computedPotentialArr.resize(cylindricalZArr.size());

    auto calcThread = [](
            int idx,
            const Coil &coil,
            const PrecisionArguments &usedPrecision,
            const std::vector<double> &cylindricalZ,
            const std::vector<double> &cylindricalR,
            std::vector<double> &computedPotential,
            size_t startIdx, size_t stopIdx
    ) -> void
    {
        for(size_t i = startIdx; i < stopIdx; i++)
        {
            auto result = coil.calculateAPotential(cylindricalZ[i], cylindricalR[i], usedPrecision);
            computedPotential[i] = result;
        }
    };

    for(size_t i = 0; i < (size_t)std::ceil((double)cylindricalZArr.size() / (double)chunkSize); i++)
    {
        g_threadPool.push(
            calcThread,
            *this,
            usedPrecision,
            cylindricalZArr, cylindricalRArr,
            computedPotentialArr,
            size_t(i * chunkSize), size_t(std::min((i + 1) * chunkSize, cylindricalZArr.size()))
        );
    }

    if(!async)
        synchronizeThreads();
}

void Coil::calculateAllBGradientMT(const std::vector<double> &cylindricalZArr,
                                   const std::vector<double> &cylindricalRArr,
                                   std::vector<double> &computedGradientRPhiArr,
                                   std::vector<double> &computedGradientRRArr,
                                   std::vector<double> &computedGradientRZArr,
                                   std::vector<double> &computedGradientZZArr,
                                   const PrecisionArguments &usedPrecision,
                                   int chunkSize, bool async) const
{
    computedGradientRPhiArr.resize(cylindricalZArr.size());
    computedGradientRRArr.resize(cylindricalZArr.size());
    computedGradientRZArr.resize(cylindricalZArr.size());
    computedGradientZZArr.resize(cylindricalZArr.size());

    auto calcThread = [] (
        int idx,
        const Coil &coil,
        const PrecisionArguments &usedPrecision,
        const std::vector<double> &cylindricalZ,
        const std::vector<double> &cylindricalR,
        std::vector<double> &computedGradientRPhi,
        std::vector<double> &computedGradientRR,
        std::vector<double> &computedGradientRZ,
        std::vector<double> &computedGradientZZ,
        size_t startIdx, size_t stopIdx
    )
    {
        for(size_t i = startIdx; i < stopIdx; i++)
        {
            auto result = coil.calculateBGradient(cylindricalZ[i], cylindricalR[i], usedPrecision);
            computedGradientRPhi[i] = result[0];
            computedGradientRR[i] = result[1];
            computedGradientRZ[i] = result[2];
            computedGradientZZ[i] = result[3];
        }
    };

    for(size_t i = 0; i < (size_t)std::ceil((double)cylindricalZArr.size() / (double)chunkSize); i++)
    {
        g_threadPool.push (
            calcThread,
            *this,
            usedPrecision,
            cylindricalZArr, cylindricalRArr,
            computedGradientRPhiArr, computedGradientRRArr,
            computedGradientRZArr, computedGradientZZArr,
            i * chunkSize, std::min((i + 1) * chunkSize, cylindricalZArr.size())
        );
    }

    if(!async)
        synchronizeThreads();
}

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

    std::vector<float> polarR(cylindricalZArr.size());
    std::vector<float> polarTheta(cylindricalZArr.size());

    for (int i = 0; i < cylindricalZArr.size(); ++i)
    {
        polarR[i] = std::sqrt(cylindricalZArr[i] * cylindricalZArr[i] + cylindricalRArr[i] * cylindricalRArr[i]);
        polarTheta[i] = atan2(cylindricalRArr[i], cylindricalZArr[i]);
    }

    std::vector<float> fieldHArr(polarR.size());
    std::vector<float> fieldZArr(polarR.size());

    Calculate_hardware_accelerated_b(polarR.size(), &polarTheta[0], &polarR[0],
                                     currentDensity, innerRadius, length, thickness,
                                     thickness/16, length/16, M_PI/48,
                                     &fieldHArr[0], &fieldZArr[0]);

    for (int i = 0; i < polarR.size(); ++i)
    {
        computedFieldHArr[i] = fieldHArr[i];
        computedFieldZArr[i] = fieldZArr[i];
    }
}
#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"
void Coil::calculateAllAPotentialGPU(const std::vector<double> &cylindricalZArr,
                                     const std::vector<double> &cylindricalRArr,
                                     std::vector<double> &computedPotentialArr,
                                     const PrecisionArguments &usedPrecision) const
{
    computedPotentialArr.resize(cylindricalZArr.size());

    std::vector<float> polarR(cylindricalZArr.size());
    std::vector<float> polarTheta(cylindricalZArr.size());

    for (int i = 0; i < cylindricalZArr.size(); ++i)
    {
        polarR[i] = std::sqrt(cylindricalZArr[i] * cylindricalZArr[i] + cylindricalRArr[i] * cylindricalRArr[i]);
        polarTheta[i] = std::atan2(cylindricalRArr[i], cylindricalZArr[i]);
    }
    std::vector<float> potentialArr(polarR.size());

    Calculate_hardware_accelerated_a(polarR.size(), &polarTheta[0], &polarR[0],
                                     currentDensity, innerRadius, length, thickness,
                                     thickness / 16, length / 16, M_PI / 48,
                                     nullptr, nullptr, &potentialArr[0]);

    // TODO - fix frequency in GPU potential calculation, current temporary fix
    for (int i = 0; i < polarR.size(); ++i)
        computedPotentialArr[i] = potentialArr[i] / (2 * M_PI);
}
#pragma clang diagnostic pop

void Coil::calculateAllBGradientGPU(const std::vector<double> &cylindricalZArr,
                                    const std::vector<double> &cylindricalRArr,
                                    std::vector<double> &computedGradientRPhi,
                                    std::vector<double> &computedGradientRR,
                                    std::vector<double> &computedGradientRZ,
                                    std::vector<double> &computedGradientZZ,
                                    const PrecisionArguments &usedPrecision) const
{
    // TODO - sometime in the distant future, this may be implemented
}

void Coil::calculateAllBFieldSwitch(const std::vector<double> &cylindricalZArr,
                                    const std::vector<double> &cylindricalRArr,
                                    std::vector<double> &computedFieldHArr,
                                    std::vector<double> &computedFieldZArr,
                                    const PrecisionArguments &usedPrecision,
                                    ComputeMethod method) const
{
    int chunkSize = calculateChunkSize(cylindricalZArr.size());

    switch (method)
    {
        case GPU:
            calculateAllBFieldGPU(cylindricalZArr, cylindricalRArr,
                                  computedFieldHArr, computedFieldZArr, usedPrecision);
            break;
        case CPU_MT:
            calculateAllBFieldMT(cylindricalZArr, cylindricalRArr,
                                 computedFieldHArr, computedFieldZArr, usedPrecision, chunkSize);
            break;
        default:
            calculateAllBFieldST(cylindricalZArr, cylindricalRArr,
                                 computedFieldHArr, computedFieldZArr, usedPrecision);
    }
}

void Coil::calculateAllAPotentialSwitch(const std::vector<double> &cylindricalZArr,
                                        const std::vector<double> &cylindricalRArr,
                                        std::vector<double> &computedPotentialArr,
                                        const PrecisionArguments &usedPrecision,
                                        ComputeMethod method) const
{
    int chunkSize = calculateChunkSize(cylindricalZArr.size());
//    printf("%d\n", chunkSize);
    switch (method)
    {
        case GPU:

            calculateAllAPotentialGPU(cylindricalZArr, cylindricalRArr, computedPotentialArr, usedPrecision);
            break;
        case CPU_MT:
            calculateAllAPotentialMT(cylindricalZArr, cylindricalRArr, computedPotentialArr, usedPrecision, chunkSize);
            break;
        default:
            calculateAllAPotentialST(cylindricalZArr, cylindricalRArr, computedPotentialArr, usedPrecision);
    }
}

void Coil::calculateAllBGradientSwitch(const std::vector<double> &cylindricalZArr,
                                       const std::vector<double> &cylindricalRArr,
                                       std::vector<double> &computedGradientRPhi,
                                       std::vector<double> &computedGradientRR,
                                       std::vector<double> &computedGradientRZ,
                                       std::vector<double> &computedGradientZZ,
                                       const PrecisionArguments &usedPrecision,
                                       ComputeMethod method) const
{
    int chunkSize = calculateChunkSize(cylindricalZArr.size());

    switch (method)
    {
        case GPU:
            calculateAllBGradientGPU(cylindricalZArr, cylindricalRArr,
                                     computedGradientRPhi, computedGradientRR, computedGradientRZ, computedGradientZZ,
                                     usedPrecision);
            break;
        case CPU_MT:
            calculateAllBGradientMT(cylindricalZArr, cylindricalRArr,
                                     computedGradientRPhi, computedGradientRR, computedGradientRZ, computedGradientZZ,
                                     usedPrecision, chunkSize);
            break;
        default:
            calculateAllBGradientST(cylindricalZArr, cylindricalRArr,
                                    computedGradientRPhi, computedGradientRR, computedGradientRZ, computedGradientZZ,
                                    usedPrecision);
    }
}

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
        int modifier = std::floor(std::log10(numOps));

        for (int i = 1; i <= std::ceil(std::log2(numOps)); ++i)
        {
            chunkCandidate = numOps / (i * modifier * threadCount);
            leftover = numOps % (i * modifier * threadCount);

            fitnessArray.push_back((double) leftover / (chunkCandidate * i));
            chunkArray.push_back(chunkCandidate);
            //printf("%d : %d %d | %g\n", i, chunkCandidate, leftover, (double) leftover / (chunkCandidate * i));
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
        printf("%d\n", chunkSize);
        return chunkSize;
    }
}

void Coil::synchronizeThreads() const { while(g_threadPool.n_idle() < threadCount){} }
