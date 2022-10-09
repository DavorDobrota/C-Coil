#include "Coil.h"
#include "LegendreMatrix.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <cstdio>
#include <sstream>


namespace
{
    const int g_defaultLegendreOrder = 20;
    const int g_defaultBlockCount = 1;
}


PrecisionArguments::PrecisionArguments() :
    PrecisionArguments
    (
        g_defaultBlockCount,
        g_defaultBlockCount,
        g_defaultBlockCount,
        g_defaultLegendreOrder,
        g_defaultLegendreOrder,
        g_defaultLegendreOrder
    ) {}

PrecisionArguments::PrecisionArguments(
        int angularBlocks, int thicknessBlocks, int lengthBlocks,
        int angularIncrements, int thicknessIncrements, int lengthIncrements) :
        angularBlocks(angularBlocks), thicknessBlocks(thicknessBlocks),
        lengthBlocks(lengthBlocks), angularIncrements(angularIncrements),
        thicknessIncrements(thicknessIncrements), lengthIncrements(lengthIncrements)
{
    if (angularIncrements > Legendre::maxLegendreOrder || angularIncrements < 1)
        PrecisionArguments::angularIncrements = g_defaultLegendreOrder;

    if (thicknessIncrements > Legendre::maxLegendreOrder || thicknessIncrements < 1)
        PrecisionArguments::thicknessIncrements = g_defaultLegendreOrder;

    if (lengthIncrements > Legendre::maxLegendreOrder || lengthIncrements < 1)
        PrecisionArguments::lengthIncrements = g_defaultLegendreOrder;


    if (angularBlocks < 1)
        PrecisionArguments::angularBlocks = g_defaultBlockCount;

    if (thicknessBlocks < 1)
        PrecisionArguments::thicknessBlocks = g_defaultBlockCount;

    if (lengthIncrements < 1)
        PrecisionArguments::lengthBlocks = g_defaultBlockCount;
}


PrecisionArguments PrecisionArguments::getCoilPrecisionArgumentsCPU(const Coil &coil, PrecisionFactor precisionFactor)
{
    return calculatePrecisionArguments(coil, precisionFactor, false);
}


PrecisionArguments PrecisionArguments::getCoilPrecisionArgumentsGPU(const Coil &coil, PrecisionFactor precisionFactor)
{
    return calculatePrecisionArguments(coil, precisionFactor, true);
}


PrecisionArguments PrecisionArguments::calculatePrecisionArguments(const Coil &coil, PrecisionFactor precisionFactor,
                                                                   bool useGPU)
{
    int layerIncrements[3] = {1, 1, 1};
    double measures[3] = {};
    int totalIncrements = 1;

    measures[0] = g_primAngularWeightModifier * std::sqrt(M_PI * (coil.getInnerRadius() + 0.5* coil.getThickness()));
    measures[1] = g_primLinearWeightModifier * std::sqrt(coil.getThickness());
    measures[2] = g_primLinearWeightModifier * std::sqrt(coil.getLength());

    std::vector<std::pair<int, double>> precisionComponents;
    precisionComponents.reserve(3);

    int baseIncrements = g_baseLayerIncrementsCPU;
    if (useGPU)
         baseIncrements = g_baseLayerIncrementsGPU;

    precisionComponents.emplace_back(0, measures[0]);
    totalIncrements *= baseIncrements;

    switch (coil.getCoilType())
    {
        case CoilType::RECTANGULAR:
        {
            totalIncrements *= baseIncrements;
            precisionComponents.emplace_back(1, measures[1]);
            break;
        }
        case CoilType::THIN:
        {
            break;
        }
        case CoilType::FLAT:
        {
            totalIncrements *= baseIncrements;
            precisionComponents.emplace_back(1, measures[1]);
            break;
        }
        case CoilType::FILAMENT:
        {
            break;
        }
    }

    totalIncrements = int(double(totalIncrements) * std::pow(2, precisionFactor.relativePrecision - 1.0));

    std::vector<std::pair<int, int>> incrementCounts;
    bool incrementSpill;

    do
    {
        incrementCounts = CoilPairArguments::balanceIncrements(totalIncrements, precisionComponents);
        incrementSpill = false;

        if (useGPU)
        {
            for (int i = 0; i < incrementCounts.size(); ++i)
            {
                if (incrementCounts[i].second > GPU_INCREMENTS)
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
        int currentIncrements = layerIncrements[0] * layerIncrements[1] * layerIncrements[2];
        printf("%d : %d %d %d | %.6g %.6g %.6g\n",
               currentIncrements, layerIncrements[0], layerIncrements[1], layerIncrements[2],
               measures[0] / layerIncrements[0], measures[1] / layerIncrements[1], measures[2] / layerIncrements[2]
        );
    #endif // PRINT_PRECISION

    if (layerIncrements[2] == 1)
        layerIncrements[2] = int(std::round(double(layerIncrements[0]) * measures[2] / measures[0]));

    if (layerIncrements[2] == 0)
        layerIncrements[2] = 1;

    if (useGPU && layerIncrements[2] > GPU_INCREMENTS)
        layerIncrements[2] = GPU_INCREMENTS;

    int numBlocks[3];
    int numIncrements[3];

    for (int i = 0; i < 3; ++i)
    {
        numBlocks[i] = (layerIncrements[i] - 1) / Legendre::maxLegendreOrder + 1;

        if (numBlocks[i] > 1)
            numIncrements[i] = layerIncrements[i] / numBlocks[i] + 1;
        else
            numIncrements[i] = layerIncrements[i];
    }

    return PrecisionArguments
    (
        numBlocks[0],
        numBlocks[1],
        numBlocks[2],
        numIncrements[0],
        numIncrements[1],
        numIncrements[2]
    );
}


PrecisionArguments PrecisionArguments::getSecondaryCoilPrecisionArgumentsGPU(const Coil &coil,
                                                                             PrecisionFactor precisionFactor)
{
    int layerIncrements[3] = {1, 1, 1};
    double measures[3] = {};
    int totalIncrements = 1;

    measures[0] = g_primAngularWeightModifier * std::sqrt(M_PI * (coil.getInnerRadius() + 0.5* coil.getThickness()));
    measures[1] = g_primLinearWeightModifier * std::sqrt(coil.getThickness());
    measures[2] = g_primLinearWeightModifier * std::sqrt(coil.getLength());

    std::vector<std::pair<int, double>> precisionComponents;
    precisionComponents.reserve(3);

    precisionComponents.emplace_back(0, measures[0]);
    totalIncrements *= g_baseLayerIncrementsGPU;

    switch (coil.getCoilType())
    {
        case CoilType::RECTANGULAR:
        {
            totalIncrements *= g_baseLayerIncrementsGPU * g_baseLayerIncrementsGPU;
            precisionComponents.emplace_back(1, measures[1]);
            precisionComponents.emplace_back(2, measures[2]);
            break;
        }
        case CoilType::THIN:
        {
            totalIncrements *= g_baseLayerIncrementsGPU;
            precisionComponents.emplace_back(2, measures[2]);
            break;
        }
        case CoilType::FLAT:
        {
            totalIncrements *= g_baseLayerIncrementsGPU;
            precisionComponents.emplace_back(1, measures[1]);
            break;
        }
        case CoilType::FILAMENT:
        {
            break;
        }
    }

    totalIncrements = int(double(totalIncrements) * std::pow(2, precisionFactor.relativePrecision - 1.0));

    std::vector<std::pair<int, int>> incrementCounts;
    bool incrementSpill;

    do
    {
        incrementCounts = CoilPairArguments::balanceIncrements(totalIncrements, precisionComponents);
        incrementSpill = false;

        for (int i = 0; i < incrementCounts.size(); ++i)
        {
            if (incrementCounts[i].second > GPU_INCREMENTS)
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

    #if PRINT_PRECISION
        int currentIncrements = layerIncrements[0] * layerIncrements[1] * layerIncrements[2];
        printf("%d : %d %d %d | %.6g %.6g %.6g\n",
               currentIncrements, layerIncrements[0], layerIncrements[1], layerIncrements[2],
               measures[0] / layerIncrements[0], measures[1] / layerIncrements[1], measures[2] / layerIncrements[2]
        );
    #endif // PRINT_PRECISION

    if (layerIncrements[2] == 1)
        layerIncrements[2] = int(std::round(double(layerIncrements[0]) * measures[2] / measures[0]));

    if (layerIncrements[2] == 0)
        layerIncrements[2] = 1;

    if (layerIncrements[2] > GPU_INCREMENTS)
        layerIncrements[2] = GPU_INCREMENTS;

    return PrecisionArguments
    (
        1,
        1,
        1,
        layerIncrements[0],
        layerIncrements[1],
        layerIncrements[2]
    );
}


PrecisionArguments::operator std::string() const
{
    std::stringstream output;

    output << "PrecisionArguments("
           << "angular_block_count=" << angularBlocks
           << ", thickness_block_count=" << thicknessBlocks
           << ", length_block_count=" << lengthBlocks
           << ", angular_increment_count=" << angularIncrements
           << ", thickness_increment_count=" << thicknessIncrements
           << ", length_increment_count=" << lengthIncrements
        << ")";

    return output.str();
}
