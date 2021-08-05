#include "Coil.h"
#include "ctpl.h"
#include "hardware_acceleration.h"

#include <cmath>

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
        std::pair<double, double> values = calculateBField(cylindricalZArr[i], cylindricalRArr[i], usedPrecision);
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
                                const PrecisionArguments &usedPrecision) const
{
    computedFieldHArr.resize(cylindricalZArr.size());
    computedFieldZArr.resize(cylindricalZArr.size());

    auto calcThread = [this, &usedPrecision](int idx, double cylindricalZ, double cylindricalR, double &computedFieldH, double &computedFieldZ) -> void
    {
        auto result = calculateBField(cylindricalZ, cylindricalR, usedPrecision);
        computedFieldH = result.first;
        computedFieldZ = result.second;
    };

    for(int i = 0; i < cylindricalZArr.size(); i++)
    {
        g_threadPool.push(calcThread, cylindricalZArr[i], cylindricalRArr[i], std::ref(computedFieldHArr[i]),
                          std::ref(computedFieldZArr[i]));
    }

    while(g_threadPool.n_idle() < threadCount);
}

void Coil::calculateAllAPotentialMT(const std::vector<double> &cylindricalZArr,
                                    const std::vector<double> &cylindricalRArr,
                                    std::vector<double> &computedPotentialArr,
                                    const PrecisionArguments &usedPrecision) const
{
    computedPotentialArr.resize(cylindricalZArr.size());

    auto calcThread = [this, &usedPrecision](int idx, double cylindricalZ, double cylindricalR, double &computedPotential) -> void
    {
        auto result = calculateAPotential(cylindricalZ, cylindricalR, usedPrecision);
        computedPotential = result;
    };

    for(int i = 0; i < cylindricalZArr.size(); i++)
    {
        g_threadPool.push(calcThread, cylindricalZArr[i], cylindricalRArr[i], std::ref(computedPotentialArr[i]));
    }

    while(g_threadPool.n_idle() < threadCount);
}

void Coil::calculateAllBGradientMT(const std::vector<double> &cylindricalZArr,
                                   const std::vector<double> &cylindricalRArr,
                                   std::vector<double> &computedGradientRPhiArr,
                                   std::vector<double> &computedGradientRRArr,
                                   std::vector<double> &computedGradientRZArr,
                                   std::vector<double> &computedGradientZZArr,
                                   const PrecisionArguments &usedPrecision) const
{
    computedGradientRPhiArr.resize(cylindricalZArr.size());
    computedGradientRRArr.resize(cylindricalZArr.size());
    computedGradientRZArr.resize(cylindricalZArr.size());
    computedGradientZZArr.resize(cylindricalZArr.size());

    auto calcThread = [this, &usedPrecision] (
        int idx, double cylindricalZ, double cylindricalR, double computedGradientRPhi, double computedGradientRR,
        double computedGradientRZ, double computedGradientZZ
    )
    {
        auto result = calculateBGradient(cylindricalZ, cylindricalR, usedPrecision);
        computedGradientRPhi = result[0];
        computedGradientRR = result[1];
        computedGradientRZ = result[2];
        computedGradientZZ = result[3];
    };

    for(int i = 0; i < cylindricalZArr.size(); i++)
    {
        g_threadPool.push (
            calcThread, cylindricalZArr[i], cylindricalRArr[i],
            std::ref(computedGradientRPhiArr[i]), std::ref(computedGradientRRArr[i]),
            std::ref(computedGradientRZArr[i]), std::ref(computedGradientZZArr[i])
        );
    }

    while(g_threadPool.n_idle() < threadCount);
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
        polarR[i] = std::std::sqrt(cylindricalZArr[i] * cylindricalZArr[i] + cylindricalRArr[i] * cylindricalRArr[i]);
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
    switch (method)
    {
        case GPU:
            calculateAllBFieldGPU(cylindricalZArr, cylindricalRArr,
                                  computedFieldHArr, computedFieldZArr, usedPrecision);
            break;
        case CPU_MT:
            calculateAllBFieldMT(cylindricalZArr, cylindricalRArr,
                                 computedFieldHArr, computedFieldZArr, usedPrecision);
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
    switch (method)
    {
        case GPU:
            calculateAllAPotentialGPU(cylindricalZArr, cylindricalRArr, computedPotentialArr, usedPrecision);
            break;
        case CPU_MT:
            calculateAllAPotentialMT(cylindricalZArr, cylindricalRArr, computedPotentialArr, usedPrecision);
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
                                     usedPrecision);
            break;
        default:
            calculateAllBGradientST(cylindricalZArr, cylindricalRArr,
                                    computedGradientRPhi, computedGradientRR, computedGradientRZ, computedGradientZZ,
                                    usedPrecision);
    }
}
