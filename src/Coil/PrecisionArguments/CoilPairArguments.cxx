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
        increments.emplace_back(component.first, incrementCount);
    }
    return increments;
}


CoilPairArguments CoilPairArguments::getAppropriateCoilPairArguments(const Coil &primary, const Coil &secondary,
                                                                     PrecisionFactor precisionFactor,
                                                                     ComputeMethod computeMethod, bool zAxisCase,
                                                                     bool pureGPU)
{
    int layerIncrements[6] = {1, 1, 1, 1, 1, 1};
    double measures[6] = {};
    int totalIncrements = 1;

    measures[0] = g_primAngularWeightModifier * std::sqrt(M_PI * (primary.getInnerRadius() + 0.5*primary.getThickness()));
    measures[1] = g_primLinearWeightModifier * std::sqrt(primary.getThickness());
    measures[2] = g_primLinearWeightModifier * std::sqrt(primary.getLength());

    measures[3] = g_secAngularWeightModifier * std::sqrt(2*M_PI * (secondary.getInnerRadius() + 0.5*secondary.getThickness()));
    measures[4] = g_secLinearWeightModifier * std::sqrt(secondary.getThickness());
    measures[5] = g_secLinearWeightModifier * std::sqrt(secondary.getLength());

    int primBaseIncrements = g_baseLayerIncrementsCPU;
    int secBaseIncrements = g_baseLayerIncrementsCPU;

    if (computeMethod == ComputeMethod::GPU)
    {
        primBaseIncrements = g_baseLayerIncrementsGPU;

        if (pureGPU)
            secBaseIncrements = g_baseLayerIncrementsGPU;
    }

    std::vector<std::pair<int, double>> precisionComponents;
    precisionComponents.reserve(6);

    precisionComponents.emplace_back(0, measures[0]);
    totalIncrements *= primBaseIncrements;

    switch (primary.getCoilType())
    {
        case CoilType::RECTANGULAR:
        {
            totalIncrements *= primBaseIncrements;
            precisionComponents.emplace_back(1, measures[1]);
            break;
        }
        case CoilType::THIN:
        {
            break;
        }
        case CoilType::FLAT:
        {
            totalIncrements *= primBaseIncrements;
            precisionComponents.emplace_back(1, measures[1]);
            break;
        }
        case CoilType::FILAMENT:
        {
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
                precisionComponents.emplace_back(4, measures[4]);
                break;
            }
            case CoilType::THIN:
            {
                break;
            }
            case CoilType::FLAT:
            {
                totalIncrements *= secBaseIncrements;
                precisionComponents.emplace_back(4, measures[4]);
                break;
            }
            case CoilType::FILAMENT:
            {
                break;
            }
        }
    } else
    {
        precisionComponents.emplace_back(3, measures[3]);
        totalIncrements *= secBaseIncrements;

        switch (secondary.getCoilType())
        {
            case CoilType::RECTANGULAR:
            {
                totalIncrements *= secBaseIncrements * secBaseIncrements;
                precisionComponents.emplace_back(4, measures[4]);
                precisionComponents.emplace_back(5, measures[5]);
                break;
            }
            case CoilType::THIN:
            {
                totalIncrements *= secBaseIncrements;
                precisionComponents.emplace_back(5, measures[5]);
                break;
            }
            case CoilType::FLAT:
            {
                totalIncrements *= secBaseIncrements;
                precisionComponents.emplace_back(4, measures[4]);
                break;
            }
            case CoilType::FILAMENT:
            {
                break;
            }
        }
    }

    totalIncrements = int(double(totalIncrements) * std::pow(2, precisionFactor.relativePrecision - 1.0));

    std::vector<std::pair<int, int>> incrementCounts;
    bool incrementSpill;

    do
    {
        incrementCounts = balanceIncrements(totalIncrements, precisionComponents);
        incrementSpill = false;

        if (computeMethod == GPU)
        {
            for (int i = 0; i < incrementCounts.size(); ++i)
            {
                if (incrementCounts[i].second > GPU_INCREMENTS && (incrementCounts[i].first < 3 || pureGPU))
                {
                    layerIncrements[incrementCounts[i].first] = int(GPU_INCREMENTS);

                    precisionComponents.erase(precisionComponents.begin() + i);
                    incrementCounts.erase(incrementCounts.begin() + i);

                    totalIncrements = int(std::ceil(double(totalIncrements) / double(GPU_INCREMENTS)));
                    incrementSpill = true;
                    break;
                }
            }
        }
    } while (incrementSpill && !precisionComponents.empty());

    for (auto & incrementCount : incrementCounts)
        layerIncrements[incrementCount.first] = incrementCount.second;

    #if PRINT_PRECISION
        int currentIncrements = layerIncrements[0] * layerIncrements[1] * layerIncrements[2] *
                                layerIncrements[3] * layerIncrements[4] * layerIncrements[5];
        printf("%d : %d %d %d | %d %d %d\n", currentIncrements,
               layerIncrements[0], layerIncrements[1], layerIncrements[2], layerIncrements[3], layerIncrements[4], layerIncrements[5]
        );

        printf("%.6g %.6g %.6g | %.6g %.6g %.6g\n",
               measures[0] / layerIncrements[0], measures[1] / layerIncrements[1], measures[2] / layerIncrements[2],
               measures[3] / layerIncrements[3], measures[4] / layerIncrements[4], measures[5] / layerIncrements[5]
        );
    #endif // PRINT_PRECISION

    for (int i = 1; i < 6; ++i)
    {
        if (layerIncrements[i] == 1)
            layerIncrements[i] = int(std::round(double(layerIncrements[0]) * measures[i] / measures[0]));

        if (layerIncrements[i] == 0)
            layerIncrements[i] = 1;

        if (computeMethod == ComputeMethod::GPU && layerIncrements[i] >= GPU_INCREMENTS)
            layerIncrements[i] = GPU_INCREMENTS;
    }

    int numBlocks[6];
    int numIncrements[6];

    for (int i = 0; i < 6; ++i)
    {
        numBlocks[i] = (layerIncrements[i] - 1) / Legendre::maxLegendreOrder + 1;

        if (numBlocks[i] > 1)
            numIncrements[i] = layerIncrements[i] / numBlocks[i] + 1;
        else
            numIncrements[i] = layerIncrements[i];
    }

    auto primaryPrecision = PrecisionArguments
    (
        numBlocks[0],
        numBlocks[1],
        numBlocks[2],
        numIncrements[0],
        numIncrements[1],
        numIncrements[2]
    );
    auto secondaryPrecision = PrecisionArguments
    (
        numBlocks[3],
        numBlocks[4],
        numBlocks[5],
        numIncrements[3],
        numIncrements[4],
        numIncrements[5]
    );

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
