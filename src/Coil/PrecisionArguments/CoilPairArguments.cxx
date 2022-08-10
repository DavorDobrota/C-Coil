#include "Coil.h"
#include "LegendreMatrix.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>
#include <sstream>


CoilPairArguments::CoilPairArguments() : CoilPairArguments(PrecisionArguments(), PrecisionArguments()) {}

CoilPairArguments::CoilPairArguments(const PrecisionArguments &primaryPrecision,
                                     const PrecisionArguments &secondaryPrecision)
{
    this->primaryPrecision = primaryPrecision;
    this->secondaryPrecision = secondaryPrecision;
}


CoilPairArguments CoilPairArguments::getAppropriateCoilPairArguments(const Coil &primary, const Coil &secondary,
                                                                     PrecisionFactor precisionFactor,
                                                                     ComputeMethod computeMethod, bool zAxisCase,
                                                                     bool pureGPU)
{
    if (computeMethod == GPU)
        return CoilPairArguments::calculateCoilPairArgumentsGPU(primary, secondary, precisionFactor, zAxisCase, pureGPU);

    else
        return CoilPairArguments::calculateCoilPairArgumentsCPU(primary, secondary, precisionFactor, zAxisCase);
}


std::vector<std::pair<int, int>>
CoilPairArguments::balanceIncrements(int totalIncrements, const std::vector<std::pair<int, double>> &components)
{
    std::vector<std::pair<int, int>> increments;
    increments.reserve(components.size());

    if (components.size() < 2) // trivial output
    {
        increments.emplace_back(components[0].first, totalIncrements);
        return increments;
    }

    double effectiveMeasure = 1.0;

    for (const auto & component : components)
        effectiveMeasure *= component.second;

    effectiveMeasure /= double(totalIncrements);
    double linearMeasure = std::pow(effectiveMeasure, 1.0 / double(components.size()));

    for (const auto & component : components)
    {
        int incrementCount = int(std::round(component.second / linearMeasure));
        int divisor = (incrementCount - 1) / Legendre::maxLegendreOrder + 1;

        if(divisor > 1)
            incrementCount = (incrementCount / divisor + 1) * divisor;

        increments.emplace_back(component.first, incrementCount);
    }
    return increments;
}


CoilPairArguments CoilPairArguments::calculateCoilPairArgumentsCPU(const Coil &primary, const Coil &secondary,
                                                                   PrecisionFactor precisionFactor, bool zAxisCase)
{
    int layerIncrements[6] = {1, 1, 1, 1, 1, 1};
    int totalIncrements = 1;

    double primAngularRoot = g_primAngularWeightModifier *
                             std::sqrt(M_PI * (primary.getInnerRadius() + 0.5 * primary.getThickness()));
    double primThicknessRoot = g_primLinearWeightModifier * std::sqrt(primary.getThickness());
    double primLengthRoot = g_primLinearWeightModifier * std::sqrt(primary.getLength());

    double secAngularRoot = g_secAngularWeightModifier *
                            std::sqrt(2 * M_PI * (secondary.getInnerRadius() + 0.5 * secondary.getThickness()));
    double secThicknessRoot = g_secLinearWeightModifier * std::sqrt(secondary.getThickness());
    double secLengthRoot = g_secLinearWeightModifier * std::sqrt(secondary.getLength());


    std::vector<std::pair<int, double>> precisionComponents;

    switch (primary.getCoilType())
    {
        case CoilType::RECTANGULAR:
        {
            totalIncrements *= g_baseLayerIncrementsCPU * g_baseLayerIncrementsCPU;
            precisionComponents.emplace_back(1, primThicknessRoot);
            precisionComponents.emplace_back(2, primAngularRoot);
            break;
        }
        case CoilType::THIN:
        {
            totalIncrements *= g_baseLayerIncrementsCPU;
            precisionComponents.emplace_back(2, primAngularRoot);
            break;
        }
        case CoilType::FLAT:
        {
            totalIncrements *= g_baseLayerIncrementsCPU * g_baseLayerIncrementsCPU;
            precisionComponents.emplace_back(1, primThicknessRoot);
            precisionComponents.emplace_back(2, primAngularRoot);
            break;
        }
        case CoilType::FILAMENT:
        {
            totalIncrements *= g_baseLayerIncrementsCPU;
            precisionComponents.emplace_back(2, primAngularRoot);
            break;
        }
    }

    if (zAxisCase)
    {
        switch (secondary.getCoilType())
        {
            case CoilType::RECTANGULAR:
            {
                totalIncrements *= g_baseLayerIncrementsCPU;
                precisionComponents.emplace_back(4, secThicknessRoot);
                break;
            }
            case CoilType::THIN:
            {
                totalIncrements *= 1;
                break;
            }
            case CoilType::FLAT:
            {
                totalIncrements *= g_baseLayerIncrementsCPU;
                precisionComponents.emplace_back(4, secThicknessRoot);
                break;
            }
            case CoilType::FILAMENT:
            {
                break;
            }
        }
    } else
    {
        switch (secondary.getCoilType())
        {
            case CoilType::RECTANGULAR:
            {
                totalIncrements *= g_baseLayerIncrementsCPU * g_baseLayerIncrementsCPU * g_baseLayerIncrementsCPU;
                precisionComponents.emplace_back(3, secLengthRoot);
                precisionComponents.emplace_back(4, secThicknessRoot);
                precisionComponents.emplace_back(5, secAngularRoot);
                break;
            }
            case CoilType::THIN:
            {
                totalIncrements *= g_baseLayerIncrementsCPU * g_baseLayerIncrementsCPU;
                precisionComponents.emplace_back(3, secLengthRoot);
                precisionComponents.emplace_back(5, secAngularRoot);
                break;
            }
            case CoilType::FLAT:
            {
                totalIncrements *= g_baseLayerIncrementsCPU * g_baseLayerIncrementsCPU;
                precisionComponents.emplace_back(4, secThicknessRoot);
                precisionComponents.emplace_back(5, secAngularRoot);
                break;
            }
            case CoilType::FILAMENT:
            {
                totalIncrements *= g_baseLayerIncrementsCPU;
                precisionComponents.emplace_back(5, secAngularRoot);
                break;
            }
        }
    }

    totalIncrements = int(double(totalIncrements) * std::pow(2, precisionFactor.relativePrecision - 1.0));

    std::vector<std::pair<int, int>> incrementCounts = balanceIncrements(totalIncrements, precisionComponents);

    for (auto & incrementCount : incrementCounts)
        layerIncrements[incrementCount.first] = incrementCount.second;

    int numBlocks[6];
    int numIncrements[6];

    for (int i = 0; i < 6; ++i)
    {
        numBlocks[i] = (layerIncrements[i] - 1) / Legendre::maxLegendreOrder + 1;
        numIncrements[i] = layerIncrements[i] / numBlocks[i];
    }

    auto primaryPrecision = PrecisionArguments(numBlocks[2],
                                               numBlocks[1],
                                               numBlocks[0],
                                               numIncrements[2],
                                               numIncrements[1],
                                               numIncrements[0]);

    auto secondaryPrecision = PrecisionArguments(numBlocks[5],
                                                 numBlocks[4],
                                                 numBlocks[3],
                                                 numIncrements[5],
                                                 numIncrements[4],
                                                 numIncrements[3]);

    #if PRINT_ENABLED
        int currentIncrements = layerIncrements[0] * layerIncrements[1] * layerIncrements[2] *
                                layerIncrements[3] * layerIncrements[4] * layerIncrements[5];
        printf("%d : %d %d %d | %d %d %d\n", currentIncrements,
               layerIncrements[0], layerIncrements[1], layerIncrements[2], layerIncrements[3], layerIncrements[4], layerIncrements[5]
        );

        printf("%.6g %.6g %.6g | %.6g %.6g %.6g\n",
               primLengthRoot / layerIncrements[0], primThicknessRoot / layerIncrements[1], primAngularRoot / layerIncrements[2],
               secLengthRoot / layerIncrements[3], secThicknessRoot / layerIncrements[4], secAngularRoot / layerIncrements[5]
        );
    #endif // PRINT_ENABLED

    return CoilPairArguments(primaryPrecision, secondaryPrecision);
}


CoilPairArguments CoilPairArguments::calculateCoilPairArgumentsGPU(const Coil &primary, const Coil &secondary,
                                                                   PrecisionFactor precisionFactor,
                                                                   bool zAxisCase, bool pureGPU)
{
    int layerIncrements[6] = {1, 1, 1, 1, 1, 1};
    int totalIncrements = 1;

    double primAngularRoot = g_primAngularWeightModifier *
                             std::sqrt(M_PI * (primary.getInnerRadius() + 0.5 * primary.getThickness()));
    double primThicknessRoot = g_primLinearWeightModifier * std::sqrt(primary.getThickness());
    double primLengthRoot = g_primLinearWeightModifier * std::sqrt(primary.getLength());

    double secAngularRoot = g_secAngularWeightModifier *
                            std::sqrt(2 * M_PI * (secondary.getInnerRadius() + 0.5 * secondary.getThickness()));
    double secThicknessRoot = g_secLinearWeightModifier * std::sqrt(secondary.getThickness());
    double secLengthRoot = g_secLinearWeightModifier * std::sqrt(secondary.getLength());

    int secBaseIncrements = g_baseLayerIncrementsCPU;

    if (pureGPU)
        secBaseIncrements = g_baseLayerIncrementsGPU;

    std::vector<std::pair<int, double>> precisionComponents;

    switch (primary.getCoilType())
    {
        case CoilType::RECTANGULAR:
        {
            totalIncrements *= g_baseLayerIncrementsGPU * g_baseLayerIncrementsGPU;
            precisionComponents.emplace_back(1, primThicknessRoot);
            precisionComponents.emplace_back(2, primAngularRoot);
            break;
        }
        case CoilType::THIN:
        {
            totalIncrements *= g_baseLayerIncrementsGPU;
            precisionComponents.emplace_back(2, primAngularRoot);
            break;
        }
        case CoilType::FLAT:
        {
            totalIncrements *= g_baseLayerIncrementsGPU * g_baseLayerIncrementsGPU;
            precisionComponents.emplace_back(1, primThicknessRoot);
            precisionComponents.emplace_back(2, primAngularRoot);
            break;
        }
        case CoilType::FILAMENT:
        {
            totalIncrements *= g_baseLayerIncrementsGPU;
            precisionComponents.emplace_back(2, primAngularRoot);
            break;
        }
    }

    if (zAxisCase)
    {
        switch (secondary.getCoilType())
        {
            case CoilType::RECTANGULAR:
            {
                totalIncrements *= secBaseIncrements;
                precisionComponents.emplace_back(4, secThicknessRoot);
                break;
            }
            case CoilType::THIN:
            {
                totalIncrements *= 1;
                break;
            }
            case CoilType::FLAT:
            {
                totalIncrements *= secBaseIncrements;
                precisionComponents.emplace_back(4, secThicknessRoot);
                break;
            }
            case CoilType::FILAMENT:
            {
                break;
            }
        }
    } else
    {
        switch (secondary.getCoilType())
        {
            case CoilType::RECTANGULAR:
            {
                totalIncrements *= secBaseIncrements * secBaseIncrements * secBaseIncrements;
                precisionComponents.emplace_back(3, secLengthRoot);
                precisionComponents.emplace_back(4, secThicknessRoot);
                precisionComponents.emplace_back(5, secAngularRoot);
                break;
            }
            case CoilType::THIN:
            {
                totalIncrements *= secBaseIncrements * secBaseIncrements;
                precisionComponents.emplace_back(3, secLengthRoot);
                precisionComponents.emplace_back(5, secAngularRoot);
                break;
            }
            case CoilType::FLAT:
            {
                totalIncrements *= secBaseIncrements * secBaseIncrements;
                precisionComponents.emplace_back(4, secThicknessRoot);
                precisionComponents.emplace_back(5, secAngularRoot);
                break;
            }
            case CoilType::FILAMENT:
            {
                totalIncrements *= secBaseIncrements;
                precisionComponents.emplace_back(5, secAngularRoot);
                break;
            }
        }
    }

    totalIncrements = int(double(totalIncrements) * std::pow(2, precisionFactor.relativePrecision - 1.0));

    std::vector<std::pair<int, int>> incrementCounts = {};
    bool incrementSpill;

    do
    {
        incrementCounts = balanceIncrements(totalIncrements, precisionComponents);
        incrementSpill = false;

        for (int i = 0; i < incrementCounts.size(); ++i)
        {
            if (incrementCounts[i].second >= GPU_INCREMENTS && (incrementCounts[i].first < 3 || pureGPU))
            {
                layerIncrements[incrementCounts[i].first] = int(GPU_INCREMENTS);
                precisionComponents.erase(precisionComponents.begin() + i);
                incrementCounts.erase(incrementCounts.begin() + i);

                totalIncrements = int(std::ceil(double(totalIncrements) / double(GPU_INCREMENTS)));

                incrementSpill = true;
                break;
            }
        }
    } while (incrementSpill && !precisionComponents.empty());

    for (auto & incrementCount : incrementCounts)
        layerIncrements[incrementCount.first] = incrementCount.second;

    int numBlocks[6];
    int numIncrements[6];

    for (int i = 0; i < 6; ++i)
    {
        numBlocks[i] = (layerIncrements[i] - 1) / Legendre::maxLegendreOrder + 1;
        numIncrements[i] = layerIncrements[i] / numBlocks[i];
    }

    auto primaryPrecision = PrecisionArguments(numBlocks[2],
                                               numBlocks[1],
                                               numBlocks[0],
                                               numIncrements[2],
                                               numIncrements[1],
                                               numIncrements[0]);

    auto secondaryPrecision = PrecisionArguments(numBlocks[5],
                                                 numBlocks[4],
                                                 numBlocks[3],
                                                 numIncrements[5],
                                                 numIncrements[4],
                                                 numIncrements[3]);

    #if PRINT_ENABLED
        int currentIncrements = layerIncrements[0] * layerIncrements[1] * layerIncrements[2] *
                                layerIncrements[3] * layerIncrements[4] * layerIncrements[5];
        printf("%d : %d %d %d | %d %d %d\n", currentIncrements,
               layerIncrements[0], layerIncrements[1], layerIncrements[2], layerIncrements[3], layerIncrements[4], layerIncrements[5]
        );

        printf("%.6g %.6g %.6g | %.6g %.6g %.6g\n",
               primLengthRoot / layerIncrements[0], primThicknessRoot / layerIncrements[1], primAngularRoot / layerIncrements[2],
               secLengthRoot / layerIncrements[3], secThicknessRoot / layerIncrements[4], secAngularRoot / layerIncrements[5]
        );
    #endif // PRINT_ENABLED

    return CoilPairArguments(primaryPrecision, secondaryPrecision);
}


CoilPairArguments::operator std::string() const
{
    std::stringstream output;

    output << "CoilPairArguments("
           << "primary_precision=" << std::string(primaryPrecision)
           << ", secondary_precision=" << std::string(secondaryPrecision)
           << ")";

    return output.str();
}
