#include "Coil.h"
#include "LegendreMatrix.h"
#include "CUDAFunctions/ConstantsAndStructs/CoilDataStructs.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>

#include <sstream>


namespace
{
    const int g_defaultLegendreOrder = 12;
    const int g_defaultBlockCount = 1;
}

PrecisionArguments::PrecisionArguments() :
        PrecisionArguments(g_defaultBlockCount, g_defaultBlockCount, g_defaultBlockCount,
                           g_defaultLegendreOrder, g_defaultLegendreOrder, g_defaultLegendreOrder) {}

PrecisionArguments::PrecisionArguments(
        int angularBlocks, int thicknessBlocks, int lengthBlocks,
        int angularIncrements, int thicknessIncrements, int lengthIncrements) :
        angularBlockCount(angularBlocks), thicknessBlockCount(thicknessBlocks),
        lengthBlockCount(lengthBlocks), angularIncrementCount(angularIncrements),
        thicknessIncrementCount(thicknessIncrements), lengthIncrementCount(lengthIncrements)
{
    if (angularIncrements > Legendre::maxLegendreOrder || angularIncrements < 1)
        PrecisionArguments::angularIncrementCount = g_defaultLegendreOrder;

    if (thicknessIncrements > Legendre::maxLegendreOrder || thicknessIncrements < 1)
        PrecisionArguments::thicknessIncrementCount = g_defaultLegendreOrder;

    if (lengthIncrements > Legendre::maxLegendreOrder || lengthIncrements < 1)
        PrecisionArguments::lengthIncrementCount = g_defaultLegendreOrder;


    if (angularBlocks < 1)
        PrecisionArguments::angularBlockCount = g_defaultBlockCount;

    if (thicknessBlocks < 1)
        PrecisionArguments::thicknessBlockCount = g_defaultBlockCount;

    if (lengthIncrements < 1)
        PrecisionArguments::lengthBlockCount = g_defaultBlockCount;
}

PrecisionArguments PrecisionArguments::getCoilPrecisionArgumentsCPU(const Coil &coil, PrecisionFactor precisionFactor)
{
    int lengthArrayIndex, thicknessArrayIndex, angularArrayIndex;

    double angularRoot = std::sqrt(M_PI * (coil.getInnerRadius() + 0.5 * coil.getThickness()));
    double thicknessRoot = std::sqrt(coil.getThickness());
    double lengthRoot = std::sqrt(coil.getLength());

    int totalIncrements = g_baseLayerIncrementsCPU;
    int currentIncrements;

    switch (coil.getCoilType())
    {
        case CoilType::RECTANGULAR:
        {
            totalIncrements *= g_baseLayerIncrementsCPU;
            thicknessArrayIndex = g_minPrimThicknessIncrements - 1;
            lengthArrayIndex = g_minPrimLengthIncrements - 1;
            break;
        }
        case CoilType::THIN:
        {
            thicknessArrayIndex = 0;
            lengthArrayIndex = g_minPrimLengthIncrements - 1;
            break;
        }
        case CoilType::FLAT:
        {
            totalIncrements *= g_baseLayerIncrementsCPU;
            thicknessArrayIndex = g_minPrimThicknessIncrements - 1;
            lengthArrayIndex = 0;
            break;
        }
        case CoilType::FILAMENT:
        {
            thicknessArrayIndex = 0;
            lengthArrayIndex = 0;
            break;
        }
    }
    angularArrayIndex = g_minPrimAngularIncrements - 1;

    totalIncrements *= std::pow(2, precisionFactor.relativePrecision - 1.0);

    double angularStep, lengthStep, thicknessStep;

    do
    {
        angularStep = angularRoot /
                (blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex]);

        lengthStep = lengthRoot /
                (blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex]);

        thicknessStep = thicknessRoot /
                (blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex]);


        if (angularStep >= std::max(lengthStep, thicknessStep))
            angularArrayIndex++;
        else
        {
            if (thicknessStep >= lengthStep)
                thicknessArrayIndex++;
            else
                lengthArrayIndex++;
        }

        currentIncrements =
                blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex] *
                blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex];
    }
    while (currentIncrements <= totalIncrements);

    #if PRINT_ENABLED
        printf("%d : %d %d %d\n", currentIncrements,
               blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex],
               blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex],
               blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex]);
    #endif // PRINT_ENABLED

    return PrecisionArguments(blockPrecisionCPUArray[angularArrayIndex],
                              blockPrecisionCPUArray[thicknessArrayIndex],
                              blockPrecisionCPUArray[lengthArrayIndex],
                              incrementPrecisionCPUArray[angularArrayIndex],
                              incrementPrecisionCPUArray[thicknessArrayIndex],
                              incrementPrecisionCPUArray[lengthArrayIndex]);
}

PrecisionArguments PrecisionArguments::getCoilPrecisionArgumentsGPU(const Coil &coil, PrecisionFactor precisionFactor)
{
    int lengthIncrements, thicknessIncrements, angularIncrements;

    double angularRoot = std::sqrt(M_PI * (coil.getInnerRadius() + 0.5 * coil.getThickness()));
    double thicknessRoot = std::sqrt(coil.getThickness());
    double lengthRoot = std::sqrt(coil.getLength());

    int totalIncrements = g_baseLayerIncrementsGPU;
    int currentIncrements;

    auto coilType = coil.getCoilType();

    switch (coilType)
    {
        case CoilType::RECTANGULAR:
        {
            totalIncrements *= g_baseLayerIncrementsGPU;
            thicknessIncrements = g_minPrimThicknessIncrements;
            lengthIncrements = g_minPrimLengthIncrements;
            break;
        }
        case CoilType::THIN:
        {
            thicknessIncrements = 1;
            lengthIncrements = g_minPrimLengthIncrements;
            break;
        }
        case CoilType::FLAT:
        {
            totalIncrements *= g_baseLayerIncrementsGPU;
            thicknessIncrements = g_minPrimThicknessIncrements;
            lengthIncrements = 1;
            break;
        }
        case CoilType::FILAMENT:
        {
            thicknessIncrements = 1;
            lengthIncrements = 1;
            break;
        }
    }
    angularIncrements = g_minPrimAngularIncrements;

    totalIncrements *= std::pow(2, precisionFactor.relativePrecision - 1.0);

    double angularStep, lengthStep, thicknessStep;
    bool exitLoop = false;

    do
    {
        lengthStep = lengthRoot / lengthIncrements;

        if (angularIncrements >= GPU_INCREMENTS)
            angularStep = 0.0;
        else
            angularStep = angularRoot / angularIncrements;

        if (thicknessIncrements >= GPU_INCREMENTS)
            thicknessStep = 0.0;
        else
            thicknessStep = thicknessRoot / thicknessIncrements;


        if (angularStep >= std::max(lengthStep, thicknessStep))
            angularIncrements++;
        else
        {
            if (thicknessStep >= lengthStep)
                thicknessIncrements++;
            else
                lengthIncrements++;
        }

        currentIncrements = angularIncrements * thicknessIncrements;

        if (currentIncrements > totalIncrements)
            exitLoop = true;

        switch (coilType) {
            case CoilType::RECTANGULAR:
            {
                if (angularIncrements >= GPU_INCREMENTS && thicknessIncrements >= GPU_INCREMENTS)
                    exitLoop = true;
                break;
            }
            case CoilType::THIN:
            {
                if (angularIncrements >= GPU_INCREMENTS)
                    exitLoop = true;
                break;
            }
            case CoilType::FLAT:
            {
                if (angularIncrements >= GPU_INCREMENTS && thicknessIncrements >= GPU_INCREMENTS)
                    exitLoop = true;
                break;
            }
            case CoilType::FILAMENT:
            {
                if (angularIncrements >= GPU_INCREMENTS)
                    exitLoop = true;
                break;
            }
        }
    }
    while (!exitLoop);

    #if PRINT_ENABLED
        printf("%d : %d %d %d\n", currentIncrements, lengthIncrements, thicknessIncrements, angularIncrements);
    #endif // PRINT_ENABLED

    return PrecisionArguments(1,1,1, angularIncrements, thicknessIncrements, lengthIncrements);
}

PrecisionArguments::operator std::string() const
{
    std::stringstream output;

    output << "PrecisionArguments("
        << "angular_block_count=" << angularBlockCount
        << ", thickness_block_count=" << thicknessBlockCount
        << ", length_block_count=" << lengthBlockCount
        << ", angular_increment_count=" << angularIncrementCount
        << ", thickness_increment_count=" << thicknessIncrementCount
        << ", length_increment_count=" << lengthIncrementCount
        << ")";

    return output.str();
}
