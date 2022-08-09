#include "Coil.h"
#include "CUDAUtils/ConstantsAndStructs/CoilDataStructs.h"

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
    {
        if (!pureGPU)
            return CoilPairArguments::calculateCoilPairArgumentsGPU(primary, secondary, precisionFactor, zAxisCase);
        else
            return CoilPairArguments::calculateCoilPairArgumentsGPUPure(primary, secondary, precisionFactor, zAxisCase);
    }
    else
        return CoilPairArguments::calculateCoilPairArgumentsCPU(primary, secondary, precisionFactor, zAxisCase);
}


CoilPairArguments CoilPairArguments::calculateCoilPairArgumentsCPU(const Coil &primary, const Coil &secondary,
                                                                   PrecisionFactor precisionFactor, bool zAxisCase)
{
    int primLengthArrayIndex, primThicknessArrayIndex, primAngularArrayIndex;
    int secLengthArrayIndex, secThicknessArrayIndex, secAngularArrayIndex;

    int totalIncrements = 1;
    int currentIncrements;

    double primAngularRoot = std::sqrt(M_PI * (primary.getInnerRadius() + 0.5 * primary.getThickness()));
    double primThicknessRoot = std::sqrt(primary.getThickness());
    double primLengthRoot = std::sqrt(primary.getLength());

    double secAngularRoot = std::sqrt(2 * M_PI * (secondary.getInnerRadius() + 0.5 * secondary.getThickness()));
    double secThicknessRoot = std::sqrt(secondary.getThickness());
    double secLengthRoot = std::sqrt(secondary.getLength());

    // multiplying for primary and secondary angular increments in case they are present
    if (!zAxisCase)
        totalIncrements *= g_baseLayerIncrementsCPU * g_baseLayerIncrementsCPU;

    // if the slow methods are used in the z-axis case one extra layer is required
    else if (secondary.getCoilType() == CoilType::FILAMENT || secondary.getCoilType() == CoilType::FLAT)
        totalIncrements *= g_baseLayerIncrementsCPU;

    switch (primary.getCoilType())
    {
        case CoilType::RECTANGULAR:
        {
            totalIncrements *= g_baseLayerIncrementsCPU;
            primThicknessArrayIndex = g_minPrimThicknessIncrements - 1;
            primLengthArrayIndex = g_minPrimLengthIncrements - 1;
            break;
        }
        case CoilType::THIN:
        {
            primThicknessArrayIndex = 0;
            primLengthArrayIndex = g_minPrimThicknessIncrements - 1;
            break;
        }
        case CoilType::FLAT:
        {
            totalIncrements *= g_baseLayerIncrementsCPU;
            primThicknessArrayIndex = g_minPrimThicknessIncrements - 1;
            primLengthArrayIndex = 0;
            break;
        }
        case CoilType::FILAMENT:
        {
            primThicknessArrayIndex = 0;
            primLengthArrayIndex = 0;
            break;
        }
    }
    primAngularArrayIndex = g_minPrimAngularIncrements - 1;

    switch (secondary.getCoilType())
    {
        case CoilType::RECTANGULAR:
        {
            totalIncrements *= g_baseLayerIncrementsCPU * g_baseLayerIncrementsCPU;
            secLengthArrayIndex = g_minSecLengthIncrements - 1;
            secThicknessArrayIndex = g_minSecThicknessIncrements - 1;
            break;
        }
        case CoilType::THIN:
        {
            totalIncrements *= g_baseLayerIncrementsCPU;
            secThicknessArrayIndex = 0;
            secLengthArrayIndex = g_minSecLengthIncrements - 1;
            break;
        }
        case CoilType::FLAT:
        {
            totalIncrements *= g_baseLayerIncrementsCPU;
            secLengthArrayIndex = 0;
            secThicknessArrayIndex = g_minSecAngularIncrements - 1;
            break;
        }
        case CoilType::FILAMENT:
        {
            secLengthArrayIndex = 0;
            secThicknessArrayIndex = 0;
            break;
        }
    }
    secAngularArrayIndex = g_minSecAngularIncrements - 1;

    totalIncrements *= std::pow(2, precisionFactor.relativePrecision - 1.0);

    double primAngularStep, primLengthStep, primThicknessStep, primLinearStep;
    double secAngularStep, secLengthStep, secThicknessStep, secLinearStep;

    do
    {
        primAngularStep = g_primAngularWeightModifier * primAngularRoot /
                (blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex]);

        primLengthStep = g_primLinearWeightModifier * primLengthRoot /
                (blockPrecisionCPUArray[primLengthArrayIndex] * incrementPrecisionCPUArray[primLengthArrayIndex]);

        primThicknessStep = g_primLinearWeightModifier * primThicknessRoot /
                (blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex]);

        secAngularStep = g_secAngularWeightModifier * secAngularRoot /
                (blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);

        secLengthStep = g_secLinearWeightModifier * secLengthRoot /
                (blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex]);

        secThicknessStep = g_secLinearWeightModifier * secThicknessRoot /
                (blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

        primLinearStep = std::max(primLengthStep, primThicknessStep);
        secLinearStep = std::max(secLengthStep, secThicknessStep);

        if (std::max(primAngularStep, primLinearStep) >= std::max(secAngularStep, secLinearStep))
        {
            if (primAngularStep >= primLinearStep)
                primAngularArrayIndex++;
            else
            {
                if (primLengthStep >= primThicknessStep)
                    primLengthArrayIndex++;
                else
                    primThicknessArrayIndex++;
            }
        }
        else
        {
            if (secAngularStep >= secLinearStep)
                secAngularArrayIndex++;
            else
            {
                if (secLengthStep >= secThicknessStep)
                    secLengthArrayIndex++;
                else
                    secThicknessArrayIndex++;
            }
        }

        currentIncrements =
                blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex] *
                blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex] *
                blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex];

        if (!zAxisCase)
            currentIncrements *= blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex] *
                                 blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex];

    }
    while (currentIncrements < totalIncrements);

    auto primaryPrecision = PrecisionArguments(blockPrecisionCPUArray[primAngularArrayIndex],
                                               blockPrecisionCPUArray[primThicknessArrayIndex],
                                               blockPrecisionCPUArray[primLengthArrayIndex],
                                               incrementPrecisionCPUArray[primAngularArrayIndex],
                                               incrementPrecisionCPUArray[primThicknessArrayIndex],
                                               incrementPrecisionCPUArray[primLengthArrayIndex]);

    auto secondaryPrecision = PrecisionArguments(blockPrecisionCPUArray[secAngularArrayIndex],
                                                 blockPrecisionCPUArray[secThicknessArrayIndex],
                                                 blockPrecisionCPUArray[secLengthArrayIndex],
                                                 incrementPrecisionCPUArray[secAngularArrayIndex],
                                                 incrementPrecisionCPUArray[secThicknessArrayIndex],
                                                 incrementPrecisionCPUArray[secLengthArrayIndex]);

    #if PRINT_ENABLED
    printf("%d : %d %d %d | %d %d %d\n", currentIncrements,
           blockPrecisionCPUArray[primLengthArrayIndex] * incrementPrecisionCPUArray[primLengthArrayIndex],
           blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex],
           blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex],
           blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex],
           blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex],
           blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);

    printf("%.6g %.6g %.6g | %.6g %.6g %.6g\n",
           primLengthStep, primThicknessStep, primAngularStep, secLengthStep, secThicknessStep, secAngularStep);
    #endif // PRINT_ENABLED

    return CoilPairArguments(primaryPrecision, secondaryPrecision);
}

CoilPairArguments CoilPairArguments::calculateCoilPairArgumentsGPU(const Coil &primary, const Coil &secondary,
                                                                   PrecisionFactor precisionFactor, bool zAxisCase)
{
    int primLengthIncrements, primThicknessIncrements, primAngularIncrements;
    int secLengthArrayIndex, secThicknessArrayIndex, secAngularArrayIndex;

    int totalIncrements = 1;
    int currentIncrements;

    double primAngularRoot = std::sqrt(M_PI * (primary.getInnerRadius() + 0.5 * primary.getThickness()));
    double primThicknessRoot = std::sqrt(primary.getThickness());
    double primLengthRoot = std::sqrt(primary.getLength());

    double secAngularRoot = std::sqrt(2 * M_PI * (secondary.getInnerRadius() + 0.5 * secondary.getThickness()));
    double secThicknessRoot = std::sqrt(secondary.getThickness());
    double secLengthRoot = std::sqrt(secondary.getLength());

    // multiplying for primary and secondary angular increments in case they are present
    if (!zAxisCase)
        totalIncrements *= g_baseLayerIncrementsCPU * g_baseLayerIncrementsCPU;

        // if the slow methods are used in the z-axis case one extra layer is required
    else if (secondary.getCoilType() == CoilType::FILAMENT || secondary.getCoilType() == CoilType::FLAT)
        totalIncrements *= g_baseLayerIncrementsCPU;

    CoilType primType = primary.getCoilType();
    CoilType secType = secondary.getCoilType();

    switch (primType)
    {
        case CoilType::RECTANGULAR:
        {
            totalIncrements *= g_baseLayerIncrementsCPU;
            primThicknessIncrements = g_minPrimThicknessIncrements - 1;
            primLengthIncrements = g_minPrimLengthIncrements - 1;
            break;
        }
        case CoilType::THIN:
        {
            primThicknessIncrements = 0;
            primLengthIncrements = g_minPrimThicknessIncrements - 1;
            break;
        }
        case CoilType::FLAT:
        {
            totalIncrements *= g_baseLayerIncrementsCPU;
            primThicknessIncrements = g_minPrimThicknessIncrements - 1;
            primLengthIncrements = 0;
            break;
        }
        case CoilType::FILAMENT:
        {
            primThicknessIncrements = 0;
            primLengthIncrements = 0;
            break;
        }
    }
    primAngularIncrements = g_minPrimAngularIncrements - 1;

    switch (secType)
    {
        case CoilType::RECTANGULAR:
        {
            totalIncrements *= g_baseLayerIncrementsCPU * g_baseLayerIncrementsCPU;
            secLengthArrayIndex = g_minSecLengthIncrements - 1;
            secThicknessArrayIndex = g_minSecThicknessIncrements - 1;
            break;
        }
        case CoilType::THIN:
        {
            totalIncrements *= g_baseLayerIncrementsCPU;
            secThicknessArrayIndex = 0;
            secLengthArrayIndex = g_minSecLengthIncrements - 1;
            break;
        }
        case CoilType::FLAT:
        {
            totalIncrements *= g_baseLayerIncrementsCPU;
            secLengthArrayIndex = 0;
            secThicknessArrayIndex = g_minSecAngularIncrements - 1;
            break;
        }
        case CoilType::FILAMENT:
        {
            secLengthArrayIndex = 0;
            secThicknessArrayIndex = 0;
            break;
        }
    }
    secAngularArrayIndex = g_minSecAngularIncrements - 1;

    totalIncrements *= std::pow(2, precisionFactor.relativePrecision - 1.0);

    double primAngularStep, primLengthStep, primThicknessStep, primLinearStep;
    double secAngularStep, secLengthStep, secThicknessStep, secLinearStep;

    bool exitLoop = false;

    do
    {
        primLengthStep = g_primLinearWeightModifier * primLengthRoot / primLengthIncrements;

        if (primAngularIncrements >= GPU_INCREMENTS)
            primAngularStep = 0.0;
        else
            primAngularStep = g_primAngularWeightModifier * primAngularRoot / primAngularIncrements;

        if (primThicknessIncrements >= GPU_INCREMENTS)
            primThicknessStep = 0.0;
        else
            primThicknessStep = g_primLinearWeightModifier * primThicknessRoot / primThicknessIncrements;

        secAngularStep = g_secAngularWeightModifier * secAngularRoot /
                         (blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);

        secLengthStep = g_secLinearWeightModifier * secLengthRoot /
                        (blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex]);

        secThicknessStep = g_secLinearWeightModifier * secThicknessRoot /
                           (blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

        primLinearStep = std::max(primLengthStep, primThicknessStep);
        secLinearStep = std::max(secLengthStep, secThicknessStep);

        if (std::max(primAngularStep, primLinearStep) >= std::max(secAngularStep, secLinearStep))
        {
            if (primAngularStep >= primLinearStep)
                primAngularIncrements++;
            else
            {
                if (primLengthStep >= primThicknessStep)
                    primLengthIncrements++;
                else
                    primThicknessIncrements++;
            }
        }
        else
        {
            if (secAngularStep >= secLinearStep)
                secAngularArrayIndex++;
            else
            {
                if (secLengthStep >= secThicknessStep)
                    secLengthArrayIndex++;
                else
                    secThicknessArrayIndex++;
            }
        }

        currentIncrements =
                primThicknessIncrements *
                primAngularIncrements *
                blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex];

        if (!zAxisCase)
            currentIncrements *= blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex] *
                                 blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex];
        else
        {
            if ((secType == CoilType::THIN || secType == CoilType::FILAMENT))
            {
                switch (primType) {
                    case CoilType::RECTANGULAR:
                    {
                        if (primAngularIncrements >= GPU_INCREMENTS && primThicknessIncrements >= GPU_INCREMENTS)
                            exitLoop = true;
                        break;
                    }
                    case CoilType::THIN:
                    {
                        if (primAngularIncrements >= GPU_INCREMENTS)
                            exitLoop = true;
                        break;
                    }
                    case CoilType::FLAT:
                    {
                        if (primAngularIncrements >= GPU_INCREMENTS && primThicknessIncrements >= GPU_INCREMENTS)
                            exitLoop = true;
                        break;
                    }
                    case CoilType::FILAMENT:
                    {
                        if (primAngularIncrements >= GPU_INCREMENTS)
                            exitLoop = true;
                        break;
                    }
                }
            }
        }

        if (currentIncrements > totalIncrements)
            exitLoop = true;
    }
    while (!exitLoop);

    auto primaryPrecision = PrecisionArguments(1,1,1,
                                               primAngularIncrements,
                                               primThicknessIncrements,
                                               primLengthIncrements);

    auto secondaryPrecision = PrecisionArguments(blockPrecisionCPUArray[secAngularArrayIndex],
                                                 blockPrecisionCPUArray[secThicknessArrayIndex],
                                                 blockPrecisionCPUArray[secLengthArrayIndex],
                                                 incrementPrecisionCPUArray[secAngularArrayIndex],
                                                 incrementPrecisionCPUArray[secThicknessArrayIndex],
                                                 incrementPrecisionCPUArray[secLengthArrayIndex]);

    #if PRINT_ENABLED
        printf("%d : %d %d %d | %d %d %d\n", currentIncrements,
        primLengthIncrements,
        primThicknessIncrements,
        primAngularIncrements,
        blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex],
        blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex],
        blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);

        printf("%.6g %.6g %.6g | %.6g %.6g %.6g\n",
               primLengthStep, primThicknessStep, primAngularStep, secLengthStep, secThicknessStep, secAngularStep);
    #endif // PRINT_ENABLED

    return CoilPairArguments(primaryPrecision, secondaryPrecision);
}


CoilPairArguments CoilPairArguments::calculateCoilPairArgumentsGPUPure(const Coil &primary, const Coil &secondary,
                                                                       PrecisionFactor precisionFactor,
                                                                       bool zAxisCase)
{
    int primLengthIncrements, primThicknessIncrements, primAngularIncrements;
    int secLengthIncrements, secThicknessIncrements, secAngularIncrements;

    int totalIncrements = 1;
    int currentIncrements;

    double primAngularRoot = std::sqrt(M_PI * (primary.getInnerRadius() + 0.5 * primary.getThickness()));
    double primThicknessRoot = std::sqrt(primary.getThickness());
    double primLengthRoot = std::sqrt(primary.getLength());

    double secAngularRoot = std::sqrt(2 * M_PI * (secondary.getInnerRadius() + 0.5 * secondary.getThickness()));
    double secThicknessRoot = std::sqrt(secondary.getThickness());
    double secLengthRoot = std::sqrt(secondary.getLength());

    // multiplying for primary and secondary angular increments in case they are present
    if (!zAxisCase)
        totalIncrements *= g_baseLayerIncrementsGPU * g_baseLayerIncrementsGPU;

        // if the slow methods are used in the z-axis case one extra layer is required
    else if (secondary.getCoilType() == CoilType::FILAMENT || secondary.getCoilType() == CoilType::FLAT)
        totalIncrements *= g_baseLayerIncrementsGPU;

    CoilType primType = primary.getCoilType();
    CoilType secType = secondary.getCoilType();

    switch (primType)
    {
        case CoilType::RECTANGULAR:
        {
            totalIncrements *= g_baseLayerIncrementsGPU;
            primThicknessIncrements = g_minPrimThicknessIncrements - 1;
            primLengthIncrements = g_minPrimLengthIncrements - 1;
            break;
        }
        case CoilType::THIN:
        {
            primThicknessIncrements = 0;
            primLengthIncrements = g_minPrimThicknessIncrements - 1;
            break;
        }
        case CoilType::FLAT:
        {
            totalIncrements *= g_baseLayerIncrementsGPU;
            primThicknessIncrements = g_minPrimThicknessIncrements - 1;
            primLengthIncrements = 0;
            break;
        }
        case CoilType::FILAMENT:
        {
            primThicknessIncrements = 0;
            primLengthIncrements = 0;
            break;
        }
    }
    primAngularIncrements = g_minPrimAngularIncrements - 1;

    switch (secType)
    {
        case CoilType::RECTANGULAR:
        {
            totalIncrements *= g_baseLayerIncrementsGPU * g_baseLayerIncrementsGPU;
            secLengthIncrements = g_minSecLengthIncrements - 1;
            secThicknessIncrements = g_minSecThicknessIncrements - 1;
            break;
        }
        case CoilType::THIN:
        {
            totalIncrements *= g_baseLayerIncrementsGPU;
            secThicknessIncrements = 0;
            secLengthIncrements = g_minSecLengthIncrements - 1;
            break;
        }
        case CoilType::FLAT:
        {
            totalIncrements *= g_baseLayerIncrementsGPU;
            secLengthIncrements = 0;
            secThicknessIncrements = g_minSecAngularIncrements - 1;
            break;
        }
        case CoilType::FILAMENT:
        {
            secLengthIncrements = 0;
            secThicknessIncrements = 0;
            break;
        }
    }
    secAngularIncrements = g_minSecAngularIncrements - 1;

    totalIncrements *= std::pow(2, precisionFactor.relativePrecision - 1.0);

    double primAngularStep, primLengthStep, primThicknessStep, primLinearStep;
    double secAngularStep, secLengthStep, secThicknessStep, secLinearStep;

    bool exitLoop = false;

    do
    {
        if (primAngularIncrements >= GPU_INCREMENTS)
            primAngularStep = 0.0;
        else
            primAngularStep = g_primAngularWeightModifier * primAngularRoot / primAngularIncrements;

        if (primThicknessIncrements >= GPU_INCREMENTS)
            primThicknessStep = 0.0;
        else
            primThicknessStep = g_primLinearWeightModifier * primThicknessRoot / primThicknessIncrements;

        if (primLengthIncrements >= GPU_INCREMENTS)
            primLengthStep = 0.0;
        else
            primLengthStep = g_primAngularWeightModifier * primLengthRoot / primLengthIncrements;

        if (secAngularIncrements >= GPU_INCREMENTS)
            secAngularStep = 0.0;
        else
            secAngularStep = g_secAngularWeightModifier * secAngularRoot / secAngularIncrements;

        if (secThicknessIncrements >= GPU_INCREMENTS)
            secThicknessStep = 0.0;
        else
            secThicknessStep = g_secAngularWeightModifier * secThicknessRoot / secThicknessIncrements;

        if (secLengthIncrements >= GPU_INCREMENTS)
            secLengthStep = 0.0;
        else
            secLengthStep = g_secAngularWeightModifier * secLengthRoot / secLengthIncrements;

        primLinearStep = std::max(primLengthStep, primThicknessStep);
        secLinearStep = std::max(secLengthStep, secThicknessStep);

        if (std::max(primAngularStep, primLinearStep) >= std::max(secAngularStep, secLinearStep))
        {
            if (primAngularStep >= primLinearStep)
                primAngularIncrements++;
            else
            {
                if (primLengthStep >= primThicknessStep)
                    primLengthIncrements++;
                else
                    primThicknessIncrements++;
            }
        }
        else
        {
            if (secAngularStep >= secLinearStep)
                secAngularIncrements++;
            else
            {
                if (secLengthStep >= secThicknessStep)
                    secLengthIncrements++;
                else
                    secThicknessIncrements++;
            }
        }

        currentIncrements = primThicknessIncrements * primAngularIncrements * secThicknessIncrements;

        if (!zAxisCase)
            currentIncrements *= secAngularIncrements * secLengthIncrements;

        if (primAngularIncrements >= GPU_INCREMENTS &&
            primThicknessIncrements >= GPU_INCREMENTS &&
            primLengthIncrements >= GPU_INCREMENTS &&
            secAngularIncrements >= GPU_INCREMENTS &&
            secThicknessIncrements >= GPU_INCREMENTS &&
            secLengthIncrements >= GPU_INCREMENTS)
        {
            exitLoop = true;
        }

        if (currentIncrements > totalIncrements)
            exitLoop = true;
    }
    while (!exitLoop);

    auto primaryPrecision = PrecisionArguments(1,1,1,
                                               primAngularIncrements,
                                               primThicknessIncrements,
                                               primLengthIncrements);

    auto secondaryPrecision = PrecisionArguments(1,1,1,
                                                 secAngularIncrements,
                                                 secThicknessIncrements,
                                                 secLengthIncrements);

#if PRINT_ENABLED
    printf("%d : %d %d %d | %d %d %d\n", currentIncrements,
           primLengthIncrements, primThicknessIncrements, primAngularIncrements,
           secLengthIncrements, secThicknessIncrements, secAngularIncrements
    );

    printf("%.6g %.6g %.6g | %.6g %.6g %.6g\n",
           primLengthStep, primThicknessStep, primAngularStep, secLengthStep, secThicknessStep, secAngularStep);
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
