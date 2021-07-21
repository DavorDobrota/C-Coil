
#include <cmath>
#include <cstdio>
#include <vector>
#include <functional>

#include "Coil.h"
#include "ComputeMethod.h"
#include "hardware_acceleration.h"
#include "LegendreMatrix.h"

namespace
{
    const double PI = 3.14159265357989323;
    const double g_MiReduced = 0.0000001;

    const int g_defaultLegendreOrder = 12;
    const int g_defaultBlockCount = 1;

    const double g_defaultCurrent = 1.0;
    const double g_defaultResistivity = 1.63e-8;
    const double g_defaultSineFrequency = 50;
    const PrecisionArguments g_defaultPrecision = PrecisionArguments(2, 1, 1, 12, 12, 12);

    const double g_minPrecisionFactor = 1.0;
    const double g_maxPrecisionFactor = 8.0;
    const double g_defaultPrecisionFactor = 5.0;

    const int g_minPrimLengthIncrements = 6;
    const int g_minPrimThicknessIncrements = 6;
    const int g_minPrimAngularIncrements = 6;

    const int g_minSecLengthIncrements = 4;
    const int g_minSecThicknessIncrements = 4;
    const int g_minSecAngularIncrements = 4;

    const double g_thinCoilApproximationRatio = 1e-6;
}

namespace Precision
{
    const PrecisionArguments defaultPrecision_ULTRAFAST = PrecisionArguments(1, 1, 1, 12, 8, 8);
    const PrecisionArguments defaultPrecision_FAST = PrecisionArguments(1, 1, 1, 12, 12, 12);
    const PrecisionArguments defaultPrecision_NORMAL = PrecisionArguments(2, 1, 1, 12, 12, 12);
    const PrecisionArguments defaultPrecision_PRECISE = PrecisionArguments(1, 1, 1, 48, 16, 16);
}

const int blockPrecisionCPUArray[precisionArraySize] =
    {
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 9,
        9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14,
        15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20, 21, 21, 21, 22, 22, 22, 23,
        23, 23, 24, 24, 24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 32, 32, 33, 33, 34, 34, 35, 35,
        36, 36, 37, 37, 38, 38, 39, 39, 40, 40, 41, 41, 42, 42, 43, 43, 44, 44, 45, 45, 46, 46, 47, 47, 48, 48, 49,
        49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75,
        76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 102, 104,
        106, 108, 110, 112, 114, 116, 118, 120, 122, 124, 126, 128, 130, 132, 134, 136, 138, 140, 142, 144, 146, 148,
        150, 153, 156, 159, 162, 165, 168, 171, 174, 177, 180, 183, 186, 189, 192, 195, 198, 200, 204, 208, 212, 216,
        220, 224, 228, 232, 236, 240, 244, 248, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300, 306, 312, 318,
        324, 330, 336, 342, 348, 350, 357, 364, 371, 378, 385, 392, 399, 400, 408, 416, 424, 432, 440, 448, 450, 459,
        468, 477, 486, 495, 500, 510, 520, 530, 540, 550, 561, 572, 583, 594, 600, 612, 624, 636, 648, 650, 663, 676,
        689, 700, 714, 728, 742, 750, 765, 780, 795, 800
    };

const int incrementPrecisionCPUArray[precisionArraySize] =
    {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
        31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 26, 27, 28, 29, 30, 31, 32,
        33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 34, 35, 36, 37, 38, 39, 40, 41, 42,
        43, 44, 45, 46, 47, 48, 49, 50, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 41, 42, 43, 44, 45, 46,
        47, 48, 49, 50, 42, 43, 44, 45, 46, 47, 48, 49, 50, 43, 44, 45, 46, 47, 48, 49, 50, 44, 45, 46, 47, 48, 49,
        50, 45, 46, 47, 48, 49, 50, 46, 47, 48, 49, 50, 46, 47, 48, 49, 50, 46, 47, 48, 49, 50, 47, 48, 49, 50, 47,
        48, 49, 50, 47, 48, 49, 50, 47, 48, 49, 50, 48, 49, 50, 48, 49, 50, 48, 49, 50, 48, 49, 50, 48, 49, 50, 48,
        49, 50, 48, 49, 50, 48, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49,
        50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50, 49, 50,
        49, 50, 49, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
        50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
        50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
        50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
        50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
        50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
        50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50
    };


PrecisionFactor::PrecisionFactor() : PrecisionFactor(g_defaultPrecisionFactor){}

PrecisionFactor::PrecisionFactor(double relativePrecision)
{
    if (relativePrecision < g_minPrecisionFactor || relativePrecision > g_maxPrecisionFactor)
        PrecisionFactor::relativePrecision = g_defaultPrecisionFactor;
    else
        PrecisionFactor::relativePrecision = relativePrecision;
}


PrecisionArguments::PrecisionArguments() : PrecisionArguments(g_defaultPrecision) {}

PrecisionArguments::PrecisionArguments(
        int angularBlocks, int thicknessBlocks, int lengthBlocks,
        int angularIncrements, int thicknessIncrements, int lengthIncrements) :
        angularBlockCount(angularBlocks), thicknessBlockCount(thicknessBlocks),
        lengthBlockCount(lengthBlocks), angularIncrementCount(angularIncrements),
        thicknessIncrementCount(thicknessIncrements), lengthIncrementCount(lengthIncrements)
{
    //TODO - fix constructor calls from main
    if (angularIncrements > Legendre::maxLegendreOrder || angularIncrements < 1)
        PrecisionArguments::angularIncrementCount = g_defaultLegendreOrder;

    if (thicknessIncrements >= Legendre::maxLegendreOrder || thicknessIncrements < 1)
        PrecisionArguments::thicknessIncrementCount = g_defaultLegendreOrder;

    if (lengthIncrements >= Legendre::maxLegendreOrder || lengthIncrements < 1)
        PrecisionArguments::lengthIncrementCount = g_defaultLegendreOrder;


    if (angularBlocks < 1)
        PrecisionArguments::angularBlockCount = g_defaultBlockCount;

    if (thicknessBlocks < 1)
        PrecisionArguments::thicknessBlockCount = g_defaultBlockCount;

    if (lengthIncrements < 1)
        PrecisionArguments::lengthBlockCount = g_defaultBlockCount;
}

PrecisionArguments PrecisionArguments::getPrecisionArgumentsForCoilCPU(const Coil &coil, PrecisionFactor precisionFactor)
{
    int lengthIncrements = g_minPrimLengthIncrements;
    int thicknessIncrements = g_minPrimThicknessIncrements;
    int angularArrayIndex = g_minPrimAngularIncrements - 1;

    int totalIncrements = pow(2, 10 + precisionFactor.relativePrecision);
    int currentIncrements;

    double radius = coil.getInnerRadius();
    double thickness = coil.getThickness();
    double length = coil.getLength();

    do
    {
        double angularStep =
                PI * (radius + thickness * 0.5) /
                (blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex]);

        double lengthStep = sqrt(2) * length / lengthIncrements;
        double thicknessStep = sqrt(2) * thickness / thicknessIncrements;

        if (thickness / length < g_thinCoilApproximationRatio)
        {
            thicknessIncrements = 1;
            if (angularStep / lengthStep >= 1.0)
                angularArrayIndex++;
            else
                lengthIncrements++;
        }
        else if (length / thickness < g_thinCoilApproximationRatio)
        {
            lengthIncrements = 1;
            if (angularStep / thicknessStep >= 1.0)
                angularArrayIndex++;
            else
                thicknessIncrements++;
        }
        else
        {
            if (angularStep / sqrt(lengthStep * thicknessStep) >= 1.0)
                angularArrayIndex++;
            else
                lengthIncrements++; thicknessIncrements++;
        }

        currentIncrements = lengthIncrements * thicknessIncrements *
                            blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex];
    }
    while (currentIncrements < totalIncrements);

    return PrecisionArguments(blockPrecisionCPUArray[angularArrayIndex], 1, 1,
                              incrementPrecisionCPUArray[angularArrayIndex], thicknessIncrements, lengthIncrements);
}

PrecisionArguments PrecisionArguments::getPrecisionArgumentsForCoilGPU(const Coil &coil, PrecisionFactor precisionFactor)
{
    // TODO - implement GPU variant of the adaptive increment function
    return PrecisionArguments();
}


MInductanceArguments::MInductanceArguments() : MInductanceArguments(g_defaultPrecision, g_defaultPrecision) {}

MInductanceArguments::MInductanceArguments(const PrecisionArguments &primaryPrecision,
                                           const PrecisionArguments &secondaryPrecision)
{
    MInductanceArguments::primaryPrecision = primaryPrecision;
    MInductanceArguments::secondaryPrecision = secondaryPrecision;
}
MInductanceArguments MInductanceArguments::getMInductanceArgumentsZCPU(const Coil &primary, const Coil &secondary,
                                                                       PrecisionFactor precisionFactor)
{
    int primLengthIncrements = g_minPrimLengthIncrements;
    int primThicknessIncrements = g_minPrimThicknessIncrements;
    int secLengthIncrements = g_minSecLengthIncrements;
    int secThicknessIncrements = g_minSecThicknessIncrements;

    int primAngularArrayIndex = g_minPrimAngularIncrements - 1;

    int totalIncrements = 0;
    int currentIncrements;
    int caseIndex;

    getMInductanceCaseAndIncrements(primary, secondary, precisionFactor, caseIndex, totalIncrements);

    printf("case %d; ", caseIndex);
    do
    {
        double primAngularStep = PI * (primary.getInnerRadius() + primary.getThickness() * 0.5) /
                (blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex]);

        double primLengthStep = sqrt(2) * primary.getLength() / primLengthIncrements;
        double primThicknessStep = sqrt(2) * primary.getThickness() / primThicknessIncrements;

        double secLengthStep = secondary.getLength() / secLengthIncrements;
        double secThicknessStep = secondary.getThickness() / secThicknessIncrements;

        double primLinearStep = sqrt(primLengthStep * primThicknessStep);
        double secLinearStep = sqrt(secLengthStep * secThicknessStep);

        switch (caseIndex)
        {
            case (1):
                primThicknessIncrements = 1;
                secThicknessIncrements = 1; secLengthIncrements = 1;

                if (primAngularStep / primLengthStep >= 1.0)
                    primAngularArrayIndex++;
                else
                    primLengthIncrements++;
                break;
            case (2):
                primThicknessIncrements = 1;
                secThicknessIncrements = 1;

                if (primAngularStep / primLengthStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primLengthStep / secLengthStep >= 1.0)
                        primLengthIncrements++;
                    else
                        secLengthIncrements++;
                }
                break;
            case (3):
                primThicknessIncrements = 1;
                secLengthIncrements = 1;

                if (primAngularStep / primLengthStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primLengthStep / secThicknessStep >= 1.0)
                        primLengthIncrements++;
                    else
                        secThicknessIncrements++;
                }
                break;
            case (4):
                primThicknessIncrements = 1;

                if (primAngularStep / primLengthStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primLengthStep / secLinearStep >= 1.0)
                        primLengthIncrements++;
                    else
                        { secLengthIncrements++; secThicknessIncrements++; }
                }
                break;
            case (5):
                primLengthIncrements = 1;
                secThicknessIncrements = 1; secLengthIncrements = 1;

                if (primAngularStep / primThicknessStep >= 1.0)
                    primAngularArrayIndex++;
                else
                    primThicknessIncrements++;
                break;
            case (6):
                primLengthIncrements = 1;
                secThicknessIncrements = 1;

                if (primAngularStep / primThicknessStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primThicknessStep / secLengthStep >= 1.0)
                        primThicknessIncrements++;
                    else
                        secLengthIncrements++;
                }
                break;
            case (7):
                primLengthIncrements = 1;
                secLengthIncrements = 1;

                if (primAngularStep / primThicknessStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primThicknessStep / secThicknessStep >= 1.0)
                        primThicknessIncrements++;
                    else
                        secThicknessIncrements++;
                }
                break;
            case (8):
                primLengthIncrements = 1;

                if (primAngularStep / primThicknessStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primThicknessStep / secLinearStep >= 1.0)
                        primThicknessIncrements++;
                    else
                    { secLengthIncrements++; secThicknessIncrements++; }
                }
                break;
            case (9):
                primLengthIncrements = 1; primThicknessIncrements = 1;
                secLengthIncrements = 1; secThicknessIncrements = 1;

                primAngularArrayIndex++;
                break;
            case (10):
                primLengthIncrements = 1; primThicknessIncrements = 1;
                secThicknessIncrements = 1;

                if (primAngularStep / secLengthStep >= 1.0)
                    primAngularArrayIndex++;
                else
                    secLengthIncrements++;
                break;
            case (11):
                primLengthIncrements = 1; primThicknessIncrements = 1;
                secLengthIncrements = 1;

                if (primAngularStep / secThicknessStep >= 1.0)
                    primAngularArrayIndex++;
                else
                    secThicknessIncrements++;
                break;
            case (12):
                primLengthIncrements = 1; primThicknessIncrements = 1;

                if (primAngularStep / secLinearStep >= 1.9)
                    primAngularArrayIndex++;
                else
                    { secLengthIncrements++; secThicknessIncrements++; }
                break;
            case (13):
                secLengthIncrements = 1; secThicknessIncrements = 1;

                if (primAngularStep / primLinearStep >= 1.0)
                    primAngularArrayIndex++;
                else
                    { primLengthIncrements++; primThicknessIncrements++; }
                break;
            case (14):
                secThicknessIncrements = 1;

                if (primAngularStep / primLinearStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primLinearStep / secLengthStep >= 1.0)
                        { primLengthIncrements++; primThicknessIncrements++; }
                    else
                        secLengthIncrements++;
                }
                break;
            case (15):
                secLengthIncrements = 1;

                if (primAngularStep / primLinearStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primLinearStep / secThicknessStep >= 1.0)
                    { primLengthIncrements++; primThicknessIncrements++; }
                    else
                        secThicknessIncrements++;
                }
                break;
            default:
                if (primAngularStep / primLinearStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primLinearStep / secLinearStep >= 1.0)
                        { primLengthIncrements++; primThicknessIncrements++; }
                    else
                        { secLengthIncrements++; secThicknessIncrements++; }
                }
        }

        currentIncrements =
                primLengthIncrements * primThicknessIncrements * secLengthIncrements * secThicknessIncrements *
                blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex];
    }
    while (currentIncrements < totalIncrements);

    PrecisionArguments primaryPrecision = PrecisionArguments(blockPrecisionCPUArray[primAngularArrayIndex], 1, 1,
                                                             incrementPrecisionCPUArray[primAngularArrayIndex],
                                                             primThicknessIncrements, primLengthIncrements);

    PrecisionArguments secondaryPrecision = PrecisionArguments(0, 1, 1, 0,
                                                               secThicknessIncrements, secLengthIncrements);

    printf("%d : %d %d %d | %d %d\n", currentIncrements,
           primLengthIncrements, primThicknessIncrements,
           blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex],
           secLengthIncrements, secThicknessIncrements);

    return MInductanceArguments(primaryPrecision, secondaryPrecision);
}

MInductanceArguments MInductanceArguments::getMInductanceArgumentsGeneralCPU(const Coil &primary, const Coil &secondary,
                                                                             PrecisionFactor precisionFactor)
{
    int primLinearIncrements = g_minPrimLengthIncrements;
//    int primThicknessIncrements = g_minPrimThicknessIncrements;
    int secLinearIncrements = g_minSecLengthIncrements;

    int primAngularArrayIndex = g_minPrimAngularIncrements - 1;
    int secAngularArrayIndex = g_minSecAngularIncrements - 1;

    int totalIncrements = pow(2, 19 + precisionFactor.relativePrecision);
    int currentIncrements;

    do
    {
        double primAngularStep =
                PI * (primary.getInnerRadius() + primary.getThickness() * 0.5) /
                (blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex]);

        double primLinearStep = sqrt(2 * primary.getThickness() * primary.getLength()) / primLinearIncrements;

        double secAngularStep =
                PI * (secondary.getInnerRadius() + secondary.getThickness() * 0.5) /
                (blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);

        double secLinearStep = sqrt(secondary.getThickness() * secondary.getLength()) / secLinearIncrements;

        if ((primLinearStep * secLinearStep) / (primAngularStep * secAngularStep) <= 1.0)
        {
            if (secAngularStep / primAngularStep <= 1.0)
                primAngularArrayIndex++;
            else
                secAngularArrayIndex++;
        }
        else
        {
            if (secLinearStep / primLinearStep <= 1.0)
                primLinearIncrements++;
            else
                secLinearIncrements++;
        }
        currentIncrements =
                primLinearIncrements * primLinearIncrements * secLinearIncrements * secLinearIncrements *
                blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex] *
                blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex];

        printf("%d : %d %d %d %d\n", currentIncrements,
               primLinearIncrements,
               blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex],
               secLinearIncrements,
               blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);
    }
    while (currentIncrements < totalIncrements);

    PrecisionArguments primaryPrecision = PrecisionArguments(blockPrecisionCPUArray[primAngularArrayIndex], 1, 1,
                                                             incrementPrecisionCPUArray[primAngularArrayIndex],
                                                             primLinearIncrements, primLinearIncrements);

    PrecisionArguments secondaryPrecision = PrecisionArguments(blockPrecisionCPUArray[secAngularArrayIndex], 1, 1,
                                                               incrementPrecisionCPUArray[secAngularArrayIndex],
                                                               secLinearIncrements, secLinearIncrements);

    printf("%d : %d %d %d %d\n", currentIncrements,
           primLinearIncrements,
           blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex],
           secLinearIncrements,
           blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);

    return MInductanceArguments(primaryPrecision, secondaryPrecision);
}

void MInductanceArguments::getMInductanceCaseAndIncrements(const Coil &primary, const Coil &secondary,
                                                           PrecisionFactor precisionFactor,
                                                           int &caseIndex, int &totalIncrements)
{
    double primRadius = primary.getInnerRadius();
    double primThickness = primary.getThickness();
    double primLength = primary.getLength();

    double secRadius = secondary.getInnerRadius();
    double secThickness = secondary.getThickness();
    double secLength = secondary.getLength();

    if (primThickness / primLength < g_thinCoilApproximationRatio)
    {
        if (secThickness / secRadius < g_thinCoilApproximationRatio && secLength / secRadius < g_thinCoilApproximationRatio)
            { caseIndex = 1; totalIncrements = pow(2, 6 + precisionFactor.relativePrecision); }

        else if (secThickness / secLength < g_thinCoilApproximationRatio)
            { caseIndex = 2; totalIncrements = pow(2, 9 + precisionFactor.relativePrecision); }

        else if (secLength / secThickness < g_thinCoilApproximationRatio)
            { caseIndex = 3; totalIncrements = pow(2, 9 + precisionFactor.relativePrecision); }

        else
            { caseIndex = 4; totalIncrements = pow(2, 12 + precisionFactor.relativePrecision); }
    }
    else if (primLength / primThickness < g_thinCoilApproximationRatio)
    {
        if (secThickness / secRadius < g_thinCoilApproximationRatio && secLength / secRadius < g_thinCoilApproximationRatio)
            { caseIndex = 5; totalIncrements = pow(2, 6 + precisionFactor.relativePrecision); }

        else if (secThickness / secLength < g_thinCoilApproximationRatio)
            { caseIndex = 6; totalIncrements = pow(2, 9 + precisionFactor.relativePrecision); }

        else if (secLength / secThickness < g_thinCoilApproximationRatio)
            { caseIndex = 7; totalIncrements = pow(2, 9 + precisionFactor.relativePrecision); }

        else
            { caseIndex = 8; totalIncrements = pow(2, 12 + precisionFactor.relativePrecision); }
    }
    else if (primThickness / primRadius < g_thinCoilApproximationRatio && primLength / primRadius < g_thinCoilApproximationRatio)
    {
        if (secThickness / secRadius < g_thinCoilApproximationRatio && secLength / secRadius < g_thinCoilApproximationRatio)
            { caseIndex = 9; totalIncrements = pow(2, 3 + precisionFactor.relativePrecision); }

        else if (secThickness / secLength < g_thinCoilApproximationRatio)
            { caseIndex = 10; totalIncrements = pow(2, 6 + precisionFactor.relativePrecision); }

        else if (secLength / secThickness < g_thinCoilApproximationRatio)
            { caseIndex = 11; totalIncrements = pow(2, 6 + precisionFactor.relativePrecision); }

        else
            { caseIndex = 12; totalIncrements = pow(2, 9 + precisionFactor.relativePrecision); }
    }
    else
    {
        if (secThickness / secRadius < g_thinCoilApproximationRatio && secLength / secRadius < g_thinCoilApproximationRatio)
            { caseIndex = 13; totalIncrements = pow(2, 9 + precisionFactor.relativePrecision); }

        else if (secThickness / secLength < g_thinCoilApproximationRatio)
            { caseIndex = 14; totalIncrements = pow(2, 12 + precisionFactor.relativePrecision); }

        else if (secLength / secThickness < g_thinCoilApproximationRatio)
            { caseIndex = 15; totalIncrements = pow(2, 12 + precisionFactor.relativePrecision); }

        else
            { caseIndex = 16; totalIncrements = pow(2, 15 + precisionFactor.relativePrecision); }
    }
}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
           double wireResistivity, double sineFrequency, const PrecisionArguments &precisionSettings) :
           innerRadius(innerRadius), thickness(thickness), length(length), numOfTurns(numOfTurns),
           precisionSettings(precisionSettings)
{
    setCurrent(current);
    calculateAverageWireThickness();
    setWireResistivity(wireResistivity);
    setSineFrequency(sineFrequency);
    calculateSelfInductance();
}

Coil::Coil() : Coil(0.0, 0.0, 0.0, 3600, 0,
                    g_defaultResistivity, g_defaultSineFrequency, g_defaultPrecision){}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current, double sineFrequency) :
        Coil(innerRadius, thickness, length, numOfTurns, current,
             g_defaultResistivity, sineFrequency, g_defaultPrecision) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current, double sineFrequency,
           const PrecisionArguments &precisionSettings) :
           Coil(innerRadius, thickness, length, numOfTurns, current,
                g_defaultResistivity, sineFrequency, precisionSettings) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current) :
        Coil(innerRadius, thickness, length, numOfTurns, current,
             g_defaultResistivity, g_defaultSineFrequency, g_defaultPrecision) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
           const PrecisionArguments &precisionSettings) :
           Coil(innerRadius, thickness, length, numOfTurns, current,
                g_defaultResistivity, g_defaultSineFrequency, precisionSettings) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns) :
        Coil(innerRadius, thickness, length, numOfTurns,
             g_defaultCurrent, g_defaultResistivity, g_defaultSineFrequency, g_defaultPrecision) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns,
           const PrecisionArguments &precisionSettings) :
           Coil(innerRadius, thickness, length, numOfTurns, g_defaultCurrent,
                g_defaultResistivity, g_defaultSineFrequency, precisionSettings) {}


double Coil::getCurrentDensity() const { return currentDensity; }

double Coil::getCurrent() const { return current; }

int Coil::getNumOfTurns() const { return numOfTurns; }

double Coil::getInnerRadius() const { return innerRadius; }

double Coil::getThickness() const { return thickness; }

double Coil::getLength() const { return length; }

double Coil::getAverageWireThickness() const { return averageWireThickness; }

bool Coil::isSineDriven1() const { return isSineDriven; }

double Coil::getSineFrequency() const { return sineFrequency; }

double Coil::getSelfInductance() const { return selfInductance; }

double Coil::getMagneticMoment() const { return magneticMoment; }

double Coil::getWireResistivity() const { return wireResistivity; }

double Coil::getResistance() const { return resistance; }

double Coil::getReactance() const { return reactance; }

double Coil::getImpedance() const { return impedance; }

const PrecisionArguments &Coil::getPrecisionSettings() const { return precisionSettings; }


void Coil::setCurrentDensity(double currentDensity)
{
    Coil::currentDensity = currentDensity;
    current = currentDensity * length * thickness / numOfTurns;
    calculateMagneticMoment();
}

void Coil::setCurrent(double current)
{
    Coil::current = current;
    currentDensity = current * numOfTurns / (length * thickness);
    calculateMagneticMoment();
}

void Coil::setWireResistivity(double wireResistivity)
{
    Coil::wireResistivity = wireResistivity;
    calculateImpedance();
}

void Coil::setSineFrequency(double sineFrequency)
{
    if (sineFrequency > 0.0)
    {
        isSineDriven = true;
        Coil::sineFrequency = sineFrequency;
    }
    else
    {
        isSineDriven = false;
        Coil::sineFrequency = 0.0;
    }
    calculateImpedance();
}

void Coil::setPrecisionSettings(const PrecisionArguments &precisionSettings)
{
    Coil::precisionSettings = precisionSettings;
}

void Coil::calculateMagneticMoment()
{
    magneticMoment = PI * current * numOfTurns *
            (pow(innerRadius, 2) + pow(thickness, 2) + innerRadius * thickness / 3);
}

void Coil::calculateAverageWireThickness()
{
    averageWireThickness = sqrt(length * thickness / numOfTurns);
}

void Coil::calculateResistance()
{
    double wireRadius = averageWireThickness * 0.5;
    double ohmicResistance = wireResistivity * numOfTurns * 2*PI *
            (innerRadius + thickness * 0.5) / (wireRadius * wireRadius * PI);
    double skinDepth = sqrt(wireResistivity / (PI * sineFrequency * g_MiReduced));

    double ohmicSurface = PI * wireRadius * wireRadius;
    double sineSurface = 2*PI * (
            skinDepth * skinDepth * (exp(-wireRadius / skinDepth) - 1) +
            skinDepth * wireRadius);

    resistance = ohmicResistance * (ohmicSurface / sineSurface);
}

void Coil::calculateReactance()
{
    reactance = selfInductance * 2*PI * sineFrequency;
}

void Coil::calculateImpedance()
{
    calculateResistance();
    calculateReactance();
    impedance = sqrt(resistance * resistance + reactance * reactance);
}

void Coil::calculateSelfInductance()
{
    // TODO - firstGen solution applied: not very precise but error is less than 1%
//    selfInductance = computeMutualInductance(*this, *this, 0.0,
//                                             PrecisionFactor(g_maxPrecisionFactor));
    selfInductance = 0.0;
}

std::pair<double, double> Coil::calculateBField(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    double magneticFieldZ = 0.0;
    double magneticFieldH = 0.0;
    
    double lengthBlock = length / usedPrecision.lengthBlockCount;
    double thicknessBlock = thickness / usedPrecision.thicknessBlockCount;
    double angularBlock = PI / usedPrecision.angularBlockCount;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int lengthIncrements = usedPrecision.lengthIncrementCount - 1;
    int thicknessIncrements = usedPrecision.thicknessIncrementCount - 1;
    int angularIncrements = usedPrecision.angularIncrementCount - 1;

    // multiplication by 2 because cosine is an even function and by 0.125 for a triple change of interval (3 times 1/2)
    double constant = g_MiReduced * currentDensity * lengthBlock * thicknessBlock * angularBlock * 2 * 0.125;

    for (int indBlockL = 0; indBlockL < usedPrecision.lengthBlockCount; ++indBlockL)
    {
        for (int indBlockT = 0; indBlockT < usedPrecision.thicknessBlockCount; ++indBlockT)
        {
            double blockPositionL = (-1) * (length * 0.5) + lengthBlock * (indBlockL + 0.5);
            double blockPositionT = innerRadius + thicknessBlock * (indBlockT + 0.5);

            for (int incL = 0; incL <= lengthIncrements; ++incL)
            {
                for (int incT = 0; incT <= thicknessIncrements; ++incT)
                {
                    double incrementPositionL = blockPositionL +
                            (lengthBlock * 0.5) * Legendre::positionMatrix[lengthIncrements][incL];
                    double incrementPositionT = blockPositionT +
                            (thicknessBlock * 0.5) * Legendre::positionMatrix[thicknessIncrements][incT];

                    double incrementWeightS =
                            Legendre::weightsMatrix[lengthIncrements][incL] *
                            Legendre::weightsMatrix[thicknessIncrements][incT];

                    double tempConstA = incrementPositionT * incrementPositionT;
                    double tempConstB = incrementPositionT * (incrementPositionL + zAxis);
                    double tempConstC = incrementPositionT * rPolar;
                    double tempConstD = tempConstA + rPolar * rPolar +
                                        (incrementPositionL + zAxis) * (incrementPositionL + zAxis);
                    double tempConstE = constant * incrementWeightS;

                    for (int indBlockFi = 0; indBlockFi < usedPrecision.angularBlockCount; ++indBlockFi)
                    {
                        double blockPositionFi = angularBlock * (indBlockFi + 0.5);

                        for (int incFi = 0; incFi <= angularIncrements; ++incFi)
                        {
                            double incrementPositionFi = blockPositionFi +
                                    (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][incFi];

                            double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incFi];

                            double cosinePhi = cos(incrementPositionFi);
                            double tempConstF = 2 * tempConstC * cosinePhi;
                            double tempConstH = (tempConstD - tempConstF) * sqrt(tempConstD - tempConstF);
                            double tempConstG = tempConstE * incrementWeightFi / tempConstH;

                            magneticFieldZ += tempConstG * (tempConstA - tempConstC * cosinePhi);
                            magneticFieldH += tempConstG * (tempConstB * cosinePhi);
                        }
                    }
                }
            }
        }
    }
    std::pair<double, double> output;
    output.first = magneticFieldH;
    output.second = magneticFieldZ;

    return output;
}

double Coil::calculateAPotential(double zAxis, double rPolar, const PrecisionArguments &usedPrecision) const
{
    double magneticPotential = 0.0;

    double lengthBlock = length / usedPrecision.lengthBlockCount;
    double thicknessBlock = thickness / usedPrecision.thicknessBlockCount;
    double angularBlock = PI / usedPrecision.angularBlockCount;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int lengthIncrements = usedPrecision.lengthIncrementCount - 1;
    int thicknessIncrements = usedPrecision.thicknessIncrementCount - 1;
    int angularIncrements = usedPrecision.angularIncrementCount - 1;

    // multiplication by 2 because cosine is an even function and by 0.125 for a triple change of interval (3 times 1/2)
    double constant = g_MiReduced * currentDensity * lengthBlock * thicknessBlock * angularBlock * 2 * 0.125;

    for (int indBlockL = 0; indBlockL < usedPrecision.lengthBlockCount; ++indBlockL)
    {
        for (int indBlockT = 0; indBlockT < usedPrecision.thicknessBlockCount; ++indBlockT)
        {
            double blockPositionL = (-1) * (length * 0.5) + lengthBlock * (indBlockL + 0.5);
            double blockPositionT = innerRadius + thicknessBlock * (indBlockT + 0.5);

            for (int incL = 0; incL <= lengthIncrements; ++incL)
            {
                for (int incT = 0; incT <= thicknessIncrements; ++incT)
                {
                    double incrementPositionL = blockPositionL +
                                                (lengthBlock * 0.5) * Legendre::positionMatrix[lengthIncrements][incL];
                    double incrementPositionT = blockPositionT +
                                                (thicknessBlock * 0.5) * Legendre::positionMatrix[thicknessIncrements][incT];

                    double incrementWeightS =
                            Legendre::weightsMatrix[lengthIncrements][incL] *
                            Legendre::weightsMatrix[thicknessIncrements][incT];

                    double tempConstA = incrementPositionT;
                    double tempConstB = incrementPositionT * rPolar;
                    double tempConstC = incrementPositionT * incrementPositionT + rPolar * rPolar +
                                        (incrementPositionL + zAxis) * (incrementPositionL + zAxis);

                    for (int indBlockFi = 0; indBlockFi < usedPrecision.angularBlockCount; ++indBlockFi)
                    {
                        double blockPositionFi = angularBlock * (indBlockFi + 0.5);

                        for (int incFi = 0; incFi <= angularIncrements; ++incFi)
                        {
                            double incrementPositionFi = blockPositionFi +
                                                         (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][incFi];

                            double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incFi];

                            double cosinePhi = cos(incrementPositionFi);

                            magneticPotential += constant * incrementWeightS * incrementWeightFi *
                                    (tempConstA * cosinePhi) /sqrt(tempConstC - 2*tempConstB * cosinePhi);
                        }
                    }
                }
            }
        }
    }
    return magneticPotential;
}


double Coil::computeBFieldH(double cylindricalZ, double cylindricalR, const PrecisionArguments &usedPrecision) const
{
    std::pair<double, double> fields = calculateBField(cylindricalZ, cylindricalR, usedPrecision);
    return fields.first;
}

double Coil::computeBFieldH(double cylindricalZ, double cylindricalR) const
{
    return computeBFieldH(cylindricalZ, cylindricalR, precisionSettings);
}

double Coil::computeBFieldX(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                            const PrecisionArguments &usedPrecision) const
{
    return computeBFieldH(cylindricalZ, cylindricalR, usedPrecision) * cos(cylindricalPhi);
}

double Coil::computeBFieldX(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeBFieldX(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

double Coil::computeBFieldY(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                            const PrecisionArguments &usedPrecision) const
{
    return computeBFieldH(cylindricalZ, cylindricalR, usedPrecision) * sin(cylindricalPhi);
}

double Coil::computeBFieldY(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeBFieldY(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

double Coil::computeBFieldZ(double cylindricalZ, double cylindricalR, const PrecisionArguments &usedPrecision) const
{
    std::pair<double, double> fields = calculateBField(cylindricalZ, cylindricalR, usedPrecision);
    return fields.second;
}

double Coil::computeBFieldZ(double cylindricalZ, double cylindricalR) const
{
    return computeBFieldZ(cylindricalZ, cylindricalR, precisionSettings);
}

double Coil::computeBFieldAbs(double cylindricalZ, double cylindricalR, const PrecisionArguments &usedPrecision) const
{
    std::pair<double, double> fields = calculateBField(cylindricalZ, cylindricalR, usedPrecision);
    return sqrt(fields.first * fields.first + fields.second * fields.second);
}

double Coil::computeBFieldAbs(double cylindricalZ, double cylindricalR) const
{
    return computeBFieldAbs(cylindricalZ, cylindricalR, precisionSettings);
}

std::vector<double> Coil::computeBFieldVector(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                              const PrecisionArguments &usedPrecision) const
{
    std::vector<double> fieldVector;
    std::pair<double, double> fields = calculateBField(cylindricalZ, cylindricalR, usedPrecision);

    fieldVector.push_back(fields.first * cos(cylindricalPhi));
    fieldVector.push_back(fields.first * sin(cylindricalPhi));
    fieldVector.push_back(fields.second);

    return fieldVector;
}
std::vector<double> Coil::computeBFieldVector(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeBFieldVector(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}


double Coil::computeAPotentialX(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                const PrecisionArguments &usedPrecision) const
{
    return (-1) * calculateAPotential(cylindricalZ, cylindricalR, usedPrecision) * sin(cylindricalPhi);
}

double Coil::computeAPotentialX(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeAPotentialX(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

double Coil::computeAPotentialY(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                const PrecisionArguments &usedPrecision) const
{
    return calculateAPotential(cylindricalZ, cylindricalR, usedPrecision) * cos(cylindricalPhi);
}

double Coil::computeAPotentialY(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeAPotentialY(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

double Coil::computeAPotentialAbs(double cylindricalZ, double cylindricalR) const
{
    return calculateAPotential(cylindricalZ, cylindricalR, precisionSettings);
}

double Coil::computeAPotentialAbs(double cylindricalZ, double cylindricalR, const PrecisionArguments &usedPrecision) const
{
    return calculateAPotential(cylindricalZ, cylindricalR, usedPrecision);
}

std::vector<double> Coil::computeAPotentialVector(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                                  const PrecisionArguments &usedPrecision) const
{
    std::vector<double> potentialVector;
    double potential = calculateAPotential(cylindricalZ, cylindricalR, usedPrecision);

    potentialVector.push_back(potential * (-sin(cylindricalPhi)));
    potentialVector.push_back(potential * cos(cylindricalPhi));
    potentialVector.push_back(0.0);
    return potentialVector;
}

std::vector<double> Coil::computeAPotentialVector(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeAPotentialVector(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

double Coil::computeEFieldX(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                            const PrecisionArguments &usedPrecision) const
{
    return 2*PI * sineFrequency * computeAPotentialX(cylindricalZ, cylindricalR, cylindricalPhi, usedPrecision);
}

double Coil::computeEFieldX(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeEFieldX(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

double Coil::computeEFieldY(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                            const PrecisionArguments &usedPrecision) const
{
    return 2*PI * sineFrequency * computeAPotentialY(cylindricalZ, cylindricalR, cylindricalPhi, usedPrecision);
}

double Coil::computeEFieldY(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeEFieldY(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

double Coil::computeEFieldAbs(double cylindricalZ, double cylindricalR, const PrecisionArguments &usedPrecision) const
{
    return 2*PI * sineFrequency * computeAPotentialAbs(cylindricalZ, cylindricalR, usedPrecision);
}

double Coil::computeEFieldAbs(double cylindricalZ, double cylindricalR) const
{
    return computeEFieldAbs(cylindricalZ, cylindricalR, precisionSettings);
}

std::vector<double> Coil::computeEFieldVector(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                              const PrecisionArguments &usedPrecision) const
{
    std::vector<double> potentialVector;
    double potential = computeEFieldAbs(cylindricalZ, cylindricalR, usedPrecision);

    potentialVector.push_back(potential * (-sin(cylindricalPhi)));
    potentialVector.push_back(potential * cos(cylindricalPhi));
    potentialVector.push_back(0.0);
    return potentialVector;
}

std::vector<double> Coil::computeEFieldVector(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeEFieldVector(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}


void Coil::calculateAllBFieldST(const std::vector<double> &cylindricalZArr,
                                const std::vector<double> &cylindricalRArr,
                                std::vector<double> &computedFieldHArr,
                                std::vector<double> &computedFieldZArr,
                                const PrecisionArguments &usedPrecision) const
{
    computedFieldHArr.resize(0);
    computedFieldZArr.resize(0);

    for (int i = 0; i < cylindricalZArr.size(); ++i)
    {
        std::pair<double, double> values = calculateBField(cylindricalZArr[i], cylindricalRArr[i], usedPrecision);
        computedFieldHArr.push_back(values.first);
        computedFieldZArr.push_back(values.second);
    }
}

void Coil::calculateAllAPotentialST(const std::vector<double> &cylindricalZArr,
                                    const std::vector<double> &cylindricalRArr,
                                    std::vector<double> &computedPotentialArr,
                                    const PrecisionArguments &usedPrecision) const
{
    computedPotentialArr.resize(0);

    for (int i = 0; i < cylindricalZArr.size(); ++i)
    {
        double field = calculateAPotential(cylindricalZArr[i], cylindricalRArr[i], usedPrecision);
        computedPotentialArr.push_back(field);
    }
}

void Coil::calculateAllBFieldMT(const std::vector<double> &cylindricalZArr,
                                const std::vector<double> &cylindricalRArr,
                                std::vector<double> &computedFieldHArr,
                                std::vector<double> &computedFieldZArr,
                                const PrecisionArguments &usedPrecision) const
{
    // TODO - implement method using ThreadPool and calculateBField
}

void Coil::calculateAllAPotentialMT(const std::vector<double> &cylindricalZArr,
                                    const std::vector<double> &cylindricalRArr,
                                    std::vector<double> &computedPotentialArr,
                                    const PrecisionArguments &usedPrecision) const
{
    // TODO - implement method using ThreadPool and calculateAPotential
}

void Coil::calculateAllBFieldGPU(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 std::vector<double> &computedFieldHArr,
                                 std::vector<double> &computedFieldZArr,
                                 const PrecisionArguments &usedPrecision) const
{
    computedFieldHArr.resize(0);
    computedFieldZArr.resize(0);

    std::vector<float> polarR, polarTheta;

    for (int i = 0; i < cylindricalZArr.size(); ++i)
    {
        polarR.push_back(sqrt(cylindricalZArr[i] * cylindricalZArr[i] + cylindricalRArr[i] * cylindricalRArr[i]));
        polarTheta.push_back(atan2(cylindricalRArr[i], cylindricalZArr[i]));
    }

    std::vector<float> fieldHArr(polarR.size());
    std::vector<float> fieldZArr(polarR.size());

    Calculate_hardware_accelerated_b(polarR.size(), &polarTheta[0], &polarR[0],
                                     currentDensity, innerRadius, length, thickness,
                                     thickness/16, length/16, PI/48,
                                     &fieldHArr[0], &fieldZArr[0]);

    for (int i = 0; i < polarR.size(); ++i)
    {
        computedFieldHArr.push_back(fieldHArr[i]);
        computedFieldZArr.push_back(fieldZArr[i]);
    }
}

void Coil::calculateAllAPotentialGPU(const std::vector<double> &cylindricalZArr,
                                     const std::vector<double> &cylindricalRArr,
                                     std::vector<double> &computedPotentialArr,
                                     const PrecisionArguments &usedPrecision) const
{
    computedPotentialArr.resize(0);

    std::vector<float> polarR, polarTheta;

    for (int i = 0; i < cylindricalZArr.size(); ++i)
    {
        polarR.push_back(sqrt(cylindricalZArr[i] * cylindricalZArr[i] + cylindricalRArr[i] * cylindricalRArr[i]));
        polarTheta.push_back(atan2(cylindricalRArr[i], cylindricalZArr[i]));
    }
    std::vector<float> potentialArr(polarR.size());

    Calculate_hardware_accelerated_a(polarR.size(), &polarTheta[0], &polarR[0],
                                     currentDensity, innerRadius, length, thickness,
                                     thickness / 16, length / 16, PI / 48,
                                     nullptr, nullptr, &potentialArr[0]);

    // TODO - fix frequency in GPU potential calculation, current temporary fix
    for (int i = 0; i < polarR.size(); ++i)
        computedPotentialArr.push_back(potentialArr[i] / (2 * PI));
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


void Coil::computeAllBFieldX(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             const PrecisionArguments &usedPrecision,
                             ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedFieldArr.resize(0);
        std::vector<double> fieldH, fieldZ;

        calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision, method);

        for (int i = 0; i < fieldH.size(); ++i)
            computedFieldArr.push_back(fieldH[i] * cos(cylindricalPhiArr[i]));
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllBFieldX(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             ComputeMethod method) const
{
    computeAllBFieldX(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedFieldArr, precisionSettings, method);
}

void Coil::computeAllBFieldY(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             const PrecisionArguments &usedPrecision,
                             ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedFieldArr.resize(0);

        std::vector<double> fieldH, fieldZ;

        calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision, method);

        for (int i = 0; i < fieldH.size(); ++i)
            computedFieldArr.push_back(fieldH[i] * sin(cylindricalPhiArr[i]));
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllBFieldY(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             ComputeMethod method) const
{
    computeAllBFieldY(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedFieldArr, precisionSettings, method);
}

void Coil::computeAllBFieldH(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             const PrecisionArguments &usedPrecision,
                             ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        std::vector<double> fieldZ;

        calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, computedFieldArr, fieldZ, usedPrecision, method);
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllBFieldH(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             ComputeMethod method) const
{
    computeAllBFieldH(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedFieldArr, precisionSettings, method);
}

void Coil::computeAllBFieldZ(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             const PrecisionArguments &usedPrecision,
                             ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        std::vector<double> fieldH;

        calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr,  fieldH,computedFieldArr, usedPrecision, method);
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllBFieldZ(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             ComputeMethod method) const
{
    computeAllBFieldZ(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedFieldArr, precisionSettings, method);
}

void Coil::computeAllBFieldAbs(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedFieldArr,
                                 const PrecisionArguments &usedPrecision,
                                 ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedFieldArr.resize(0);

        std::vector<double> fieldH, fieldZ;

        calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision, method);

        for (int i = 0; i < fieldH.size(); ++i)
            computedFieldArr.push_back(std::sqrt(fieldH[i] * fieldH[i] + fieldZ[i] * fieldZ[i]));
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllBFieldAbs(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               ComputeMethod method) const
{
    computeAllBFieldAbs(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                        computedFieldArr, precisionSettings, method);
}

void Coil::computeAllBFieldComponents(const std::vector<double> &cylindricalZArr,
                                      const std::vector<double> &cylindricalRArr,
                                      const std::vector<double> &cylindricalPhiArr,
                                      std::vector<double> &computedFieldXArr,
                                      std::vector<double> &computedFieldYArr,
                                      std::vector<double> &computedFieldZArr,
                                      const PrecisionArguments &usedPrecision,
                                      ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedFieldXArr.resize(0);
        computedFieldYArr.resize(0);
        computedFieldZArr.resize(0);

        std::vector<double> fieldH, fieldZ;

        calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision, method);

        for (int i = 0; i < fieldH.size(); ++i)
        {
            computedFieldXArr.push_back(fieldH[i] * cos(cylindricalPhiArr[i]));
            computedFieldYArr.push_back(fieldH[i] * sin(cylindricalPhiArr[i]));
            computedFieldZArr.push_back(fieldZ[i]);
        }
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllBFieldComponents(const std::vector<double> &cylindricalZArr,
                                      const std::vector<double> &cylindricalRArr,
                                      const std::vector<double> &cylindricalPhiArr,
                                      std::vector<double> &computedFieldXArr,
                                      std::vector<double> &computedFieldYArr,
                                      std::vector<double> &computedFieldZArr,
                                      ComputeMethod method) const
{
    computeAllBFieldComponents(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
            computedFieldXArr, computedFieldYArr, computedFieldZArr,
            precisionSettings, method);
}

void Coil::computeAllAPotentialX(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedPotentialArr,
                                 const PrecisionArguments &usedPrecision,
                                 ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedPotentialArr.resize(0);

        std::vector<double> potentialA;

        calculateAllAPotentialSwitch(cylindricalZArr, cylindricalRArr, potentialA, usedPrecision, method);

        for (int i = 0; i < potentialA.size(); ++i)
            computedPotentialArr.push_back(potentialA[i] * (-1) * sin(cylindricalPhiArr[i]));

    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllAPotentialX(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedPotentialArr,
                                 ComputeMethod method) const
{
    computeAllAPotentialX(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedPotentialArr, precisionSettings, method);
}

void Coil::computeAllAPotentialY(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedPotentialArr,
                                 const PrecisionArguments &usedPrecision,
                                 ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedPotentialArr.resize(0);

        std::vector<double> potentialA;

        calculateAllAPotentialSwitch(cylindricalZArr, cylindricalRArr, potentialA, usedPrecision, method);

        for (int i = 0; i < potentialA.size(); ++i)
            computedPotentialArr.push_back(potentialA[i] * cos(cylindricalPhiArr[i]));
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllAPotentialY(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedPotentialArr,
                                 ComputeMethod method) const
{
    computeAllAPotentialY(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedPotentialArr, precisionSettings, method);
}

void
Coil::computeAllAPotentialAbs(const std::vector<double> &cylindricalZArr,
                              const std::vector<double> &cylindricalRArr,
                              std::vector<double> &computedPotentialArr,
                              const PrecisionArguments &usedPrecision,
                              ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size())
    {
        computedPotentialArr.resize(0);

        calculateAllAPotentialSwitch(cylindricalZArr, cylindricalRArr, computedPotentialArr, usedPrecision, method);
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void
Coil::computeAllAPotentialAbs(const std::vector<double> &cylindricalZArr,
                              const std::vector<double> &cylindricalRArr,
                              std::vector<double> &computedPotentialArr,
                              ComputeMethod method) const
{
    computeAllAPotentialAbs(
            cylindricalZArr, cylindricalRArr, computedPotentialArr, precisionSettings, method);
}

void Coil::computeAllAPotentialComponents(const std::vector<double> &cylindricalZArr,
                                          const std::vector<double> &cylindricalRArr,
                                          const std::vector<double> &cylindricalPhiArr,
                                          std::vector<double> &computedPotentialXArr,
                                          std::vector<double> &computedPotentialYArr,
                                          std::vector<double> &computedPotentialZArr,
                                          const PrecisionArguments &usedPrecision,
                                          ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedPotentialXArr.resize(0);
        computedPotentialYArr.resize(0);
        computedPotentialZArr.resize(0);

        std::vector<double> potentialA;

        calculateAllAPotentialSwitch(cylindricalZArr, cylindricalRArr, potentialA, usedPrecision, method);

        for (int i = 0; i < potentialA.size(); ++i)
        {
            computedPotentialXArr.push_back(potentialA[i] * (-1) * sin(cylindricalPhiArr[i]));
            computedPotentialYArr.push_back(potentialA[i] * cos(cylindricalPhiArr[i]));
            computedPotentialZArr.push_back(0.0);
        }
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllAPotentialComponents(const std::vector<double> &cylindricalZArr,
                                          const std::vector<double> &cylindricalRArr,
                                          const std::vector<double> &cylindricalPhiArr,
                                          std::vector<double> &computedPotentialXArr,
                                          std::vector<double> &computedPotentialYArr,
                                          std::vector<double> &computedPotentialZArr,
                                          ComputeMethod method) const
{
    computeAllAPotentialComponents(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
            computedPotentialXArr, computedPotentialYArr, computedPotentialZArr, precisionSettings, method);
}

void Coil::computeAllEFieldX(const std::vector<double> &cylindricalZArr, const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr, std::vector<double> &computedFieldArr,
                             const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    computeAllAPotentialX(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                          computedFieldArr, usedPrecision, method);

    for (double & i : computedFieldArr)
        i *= (2 * PI * sineFrequency);
}

void Coil::computeAllEFieldX(const std::vector<double> &cylindricalZArr, const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr, std::vector<double> &computedFieldArr,
                             ComputeMethod method) const
{
    computeAllEFieldX(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                      computedFieldArr, precisionSettings, method);
}

void Coil::computeAllEFieldY(const std::vector<double> &cylindricalZArr, const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr, std::vector<double> &computedFieldArr,
                             const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    computeAllAPotentialX(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                          computedFieldArr, usedPrecision, method);

    for (double & i : computedFieldArr)
        i *= (2 * PI * sineFrequency);
}

void Coil::computeAllEFieldY(const std::vector<double> &cylindricalZArr, const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr, std::vector<double> &computedFieldArr,
                             ComputeMethod method) const
{
    computeAllEFieldY(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                      computedFieldArr, precisionSettings, method);
}

void Coil::computeAllEFieldAbs(const std::vector<double> &cylindricalZArr, const std::vector<double> &cylindricalRArr,
                               std::vector<double> &computedFieldArr, const PrecisionArguments &usedPrecision,
                               ComputeMethod method) const
{
    computeAllAPotentialAbs(cylindricalZArr, cylindricalRArr, computedFieldArr, usedPrecision, method);

    for (double & i : computedFieldArr)
        i *= (2 * PI * sineFrequency);
}

void Coil::computeAllEFieldComponents(const std::vector<double> &cylindricalZArr,
                                      const std::vector<double> &cylindricalRArr,
                                      const std::vector<double> &cylindricalPhiArr,
                                      std::vector<double> &computedFieldXArr,
                                      std::vector<double> &computedFieldYArr,
                                      std::vector<double> &computedFieldZArr,
                                      const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    computeAllAPotentialComponents(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                                   computedFieldXArr, computedFieldYArr, computedFieldZArr,
                                   usedPrecision, method);

    for (int i = 0; i < computedFieldXArr.size(); ++i)
    {
        computedFieldXArr[i] *= (2 * PI * sineFrequency);
        computedFieldYArr[i] *= (2 * PI * sineFrequency);
        computedFieldZArr[i] *= (2 * PI * sineFrequency);
    }
}

void Coil::computeAllEFieldComponents(const std::vector<double> &cylindricalZArr,
                                      const std::vector<double> &cylindricalRArr,
                                      const std::vector<double> &cylindricalPhiArr,
                                      std::vector<double> &computedFieldXArr,
                                      std::vector<double> &computedFieldYArr,
                                      std::vector<double> &computedFieldZArr,
                                      ComputeMethod method) const
{
    computeAllEFieldComponents(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                               computedFieldXArr, computedFieldYArr, computedFieldZArr,
                               precisionSettings, method);
}

void Coil::computeAllEFieldAbs(const std::vector<double> &cylindricalZArr, const std::vector<double> &cylindricalRArr,
                               std::vector<double> &computedFieldArr, ComputeMethod method) const
{
    computeAllEFieldAbs(cylindricalZArr, cylindricalRArr, computedFieldArr, precisionSettings, method);
}


void Coil::calculateRingIncrementPosition(int angularBlocks, int angularIncrements,
                                          double alpha, double beta, double ringIntervalSize,
                                          std::vector<double> &ringXPosition,
                                          std::vector<double> &ringYPosition,
                                          std::vector<double> &ringZPosition,
                                          std::vector<double> &ringXTangent,
                                          std::vector<double> &ringYTangent,
                                          std::vector<double> &ringZTangent)
{
    ringXPosition.resize(0);
    ringYPosition.resize(0);
    ringZPosition.resize(0);

    ringXTangent.resize(0);
    ringYTangent.resize(0);
    ringZTangent.resize(0);

    double angularBlock = ringIntervalSize / angularBlocks;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    angularBlocks--;
    angularIncrements--;

    for (int phiBlock = 0; phiBlock <= angularBlocks; ++phiBlock)
    {
        double blockPositionPhi = angularBlock * (phiBlock + 0.5);

        for (int phiIndex = 0; phiIndex <= angularIncrements; ++phiIndex)
        {
            // PI/2 added to readjust to an even interval so a shortcut can be used
            double phi = PI/2 + blockPositionPhi +
                    (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][phiIndex];

            ringXPosition.push_back(cos(beta) * cos(phi) - sin(beta) * cos(alpha) * sin(phi));
            ringYPosition.push_back(sin(beta) * cos(phi) + cos(beta) * cos(alpha) * sin(phi));
            ringZPosition.push_back(sin(alpha) * sin(phi));

            ringXTangent.push_back((-1) * cos(beta) * sin(phi) - sin(beta) * cos(alpha) * cos(phi));
            ringYTangent.push_back((-1) * sin(beta) * sin(phi) + cos(beta) * cos(alpha) * cos(phi));
            ringZTangent.push_back(sin(alpha) * cos(phi));
        }
    }
}

double Coil::calculateMutualInductanceZAxis(const Coil &primary, const Coil &secondary, double zDisplacement,
                                            MInductanceArguments inductanceArguments, ComputeMethod method)
{
    PrecisionArguments primaryPrecisionArguments = inductanceArguments.primaryPrecision;
    int lengthIncrements = inductanceArguments.secondaryPrecision.lengthIncrementCount;
    int thicknessIncrements = inductanceArguments.secondaryPrecision.thicknessIncrementCount;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int maxLengthIndex = lengthIncrements - 1;
    int maxThicknessIndex = thicknessIncrements - 1;

    std::vector<double> zPositions;
    std::vector<double> rPositions;

    std::vector<double> weights;

    for (int zIndex = 0; zIndex < lengthIncrements; ++zIndex)
    {
        for (int rIndex = 0; rIndex < thicknessIncrements; ++rIndex)
        {
            zPositions.push_back(zDisplacement + (secondary.length * 0.5) *
                                                 Legendre::positionMatrix[maxLengthIndex][zIndex]);

            rPositions.push_back(secondary.innerRadius + secondary.thickness * 0.5 +
                                 (secondary.thickness * 0.5) * Legendre::positionMatrix[maxThicknessIndex][rIndex]);

            weights.push_back(
                    0.25 * Legendre::weightsMatrix[maxLengthIndex][zIndex] *
                    Legendre::weightsMatrix[maxThicknessIndex][rIndex]);
        }
    }

    std::vector<double> potentialA;
    double mutualInductance = 0.0;

    primary.computeAllAPotentialAbs(zPositions, rPositions, potentialA, primaryPrecisionArguments, method);

    for (int i = 0; i < potentialA.size(); ++i)
    {
        mutualInductance += 2*PI * rPositions[i] * potentialA[i] * weights[i];
    }

    return mutualInductance * secondary.numOfTurns / primary.current;
}

double Coil::calculateMutualInductanceGeneral(const Coil &primary, const Coil &secondary,
                                              double zDisplacement, double rDisplacement,
                                              double alphaAngle, double betaAngle,
                                              MInductanceArguments inductanceArguments, ComputeMethod method)
{
    if (rDisplacement == 0.0 && alphaAngle == 0.0)
    {
        return calculateMutualInductanceZAxis(primary, secondary, zDisplacement, inductanceArguments, method);
    }
    else {
        PrecisionArguments primaryPrecisionArguments = inductanceArguments.primaryPrecision;
        int lengthIncrements = inductanceArguments.secondaryPrecision.lengthIncrementCount;
        int thicknessIncrements = inductanceArguments.secondaryPrecision.thicknessIncrementCount;
        int angularBlocks = inductanceArguments.secondaryPrecision.angularBlockCount;
        int angularIncrements = inductanceArguments.secondaryPrecision.angularIncrementCount;

        int numElements = lengthIncrements * thicknessIncrements * angularBlocks * angularIncrements;

        // sometimes the function is even so a shortcut can be used to improve performance and efficiency
        double ringIntervalSize;

        if (rDisplacement == 0.0 || alphaAngle == 0.0 || betaAngle == 0.0)
            ringIntervalSize = PI;
        else
            ringIntervalSize = 2 * PI;

        std::vector<double> unitRingPointsX, unitRingPointsY, unitRingPointsZ;
        std::vector<double> unitRingTangentsX, unitRingTangentsY, unitRingTangentsZ;

        calculateRingIncrementPosition(angularBlocks, angularIncrements, alphaAngle, betaAngle, ringIntervalSize,
                                       unitRingPointsX, unitRingPointsY, unitRingPointsZ,
                                       unitRingTangentsX, unitRingTangentsY, unitRingTangentsZ);

        // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
        int maxLengthIndex = lengthIncrements - 1;
        int maxThicknessIndex = thicknessIncrements - 1;
        int maxAngularIncrementIndex = angularIncrements - 1;

        std::vector<double> zPositions;
        std::vector<double> rPositions;

        std::vector<double> weights;

        for (int zIndex = 0; zIndex < lengthIncrements; ++zIndex)
        {
            for (int rIndex = 0; rIndex < thicknessIncrements; ++rIndex)
            {
                double ringRadius = secondary.innerRadius + secondary.thickness * 0.5 +
                                    (secondary.thickness * 0.5) * Legendre::positionMatrix[maxThicknessIndex][rIndex];

                double lengthDisplacement = (secondary.length * 0.5) * Legendre::positionMatrix[maxLengthIndex][zIndex];

                for (int phiBlock = 0; phiBlock < angularBlocks; ++phiBlock)
                {
                    for (int phiIndex = 0; phiIndex < angularIncrements; ++phiIndex)
                    {
                        int phiPosition = phiBlock * angularIncrements + phiIndex;

                         double displacementX = lengthDisplacement * sin(alphaAngle) * sin(betaAngle) +
                                                ringRadius * unitRingPointsX[phiPosition];

                         double displacementY = rDisplacement - lengthDisplacement * sin(alphaAngle) * cos(betaAngle) +
                                                ringRadius * unitRingPointsY[phiPosition];

                         double displacementZ = zDisplacement + lengthDisplacement * cos(alphaAngle) +
                                                ringRadius * unitRingPointsZ[phiPosition];

                        zPositions.push_back(displacementZ);
                        rPositions.push_back(sqrt(displacementX * displacementX + displacementY * displacementY));

                        double rhoAngle = atan2(displacementY, displacementX);

                        double orientationFactor =
                                - sin(rhoAngle) * unitRingTangentsX[phiPosition] +
                                cos(rhoAngle) * unitRingTangentsY[phiPosition];

                        weights.push_back(
                                0.125 * orientationFactor * 2 * PI * ringRadius / angularBlocks *
                                Legendre::weightsMatrix[maxLengthIndex][zIndex] *
                                Legendre::weightsMatrix[maxThicknessIndex][rIndex] *
                                Legendre::weightsMatrix[maxAngularIncrementIndex][phiIndex]);
                    }
                }
            }
        }
        std::vector<double> potentialArray;
        double mutualInductance = 0.0;

        primary.computeAllAPotentialAbs(zPositions, rPositions, potentialArray, primaryPrecisionArguments, method);

        for (int i = 0; i < numElements; ++i)
            mutualInductance += potentialArray[i] * weights[i];

        return mutualInductance * secondary.numOfTurns / primary.current;
    }
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary, double zDisplacement,
                                     MInductanceArguments inductanceArguments, ComputeMethod method)
{
    return calculateMutualInductanceGeneral(primary, secondary, zDisplacement,
                                            0.0, 0.0, 0.0, inductanceArguments, method);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary, double zDisplacement,
                                     PrecisionFactor precisionFactor, ComputeMethod method)
{
    MInductanceArguments args = MInductanceArguments::getMInductanceArgumentsZCPU(primary, secondary, precisionFactor);
    return computeMutualInductance(primary, secondary, zDisplacement, args, method);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     double zDisplacement, double rDisplacement,
                                     MInductanceArguments inductanceArguments, ComputeMethod method)
{
    return calculateMutualInductanceGeneral(primary, secondary, zDisplacement, rDisplacement,
                                            0.0, 0.0, inductanceArguments, method);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     double zDisplacement, double rDisplacement,
                                     PrecisionFactor precisionFactor, ComputeMethod method)
{
    MInductanceArguments args = MInductanceArguments::getMInductanceArgumentsGeneralCPU(primary, secondary, precisionFactor);
    return computeMutualInductance(primary, secondary, zDisplacement, rDisplacement, args, method);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     double zDisplacement, double rDisplacement, double alphaAngle,
                                     MInductanceArguments inductanceArguments, ComputeMethod method)
{
    return calculateMutualInductanceGeneral(primary, secondary, zDisplacement, rDisplacement, alphaAngle,
                                            0.0, inductanceArguments, method);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     double zDisplacement, double rDisplacement, double alphaAngle,
                                     PrecisionFactor precisionFactor, ComputeMethod method)
{
    MInductanceArguments args = MInductanceArguments::getMInductanceArgumentsGeneralCPU(primary, secondary, precisionFactor);
    return computeMutualInductance(primary, secondary, zDisplacement, rDisplacement, alphaAngle, args, method);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     double zDisplacement, double rDisplacement, double alphaAngle, double betaAngle,
                                     MInductanceArguments inductanceArguments, ComputeMethod method)
{
    return calculateMutualInductanceGeneral(primary, secondary, zDisplacement, rDisplacement, alphaAngle, betaAngle,
                                            inductanceArguments, method);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     double zDisplacement, double rDisplacement, double alphaAngle, double betaAngle,
                                     PrecisionFactor precisionFactor, ComputeMethod method)
{
    MInductanceArguments args = MInductanceArguments::getMInductanceArgumentsGeneralCPU(primary, secondary, precisionFactor);
    return computeMutualInductance(primary, secondary, zDisplacement, rDisplacement, alphaAngle, betaAngle, args, method);
}

double Coil::computeInducedVoltageOn(const Coil &secondary, double zDisplacement,
                                     PrecisionFactor precisionFactor, ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, precisionFactor, method) *
           2 * PI * sineFrequency;
}

double Coil::computeInducedVoltageOn(const Coil &secondary, double zDisplacement,
                                     MInductanceArguments inductanceArguments, ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, inductanceArguments, method) *
           2 * PI * sineFrequency;;
}

double Coil::computeInducedVoltageOn(const Coil &secondary, double zDisplacement, double rDisplacement,
                                     PrecisionFactor precisionFactor, ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, rDisplacement, precisionFactor, method) *
           2 * PI * sineFrequency;
}

double Coil::computeInducedVoltageOn(const Coil &secondary, double zDisplacement, double rDisplacement,
                                     MInductanceArguments inductanceArguments, ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, rDisplacement, inductanceArguments, method) *
           2 * PI * sineFrequency;
}

double Coil::computeInducedVoltageOn(const Coil &secondary, double zDisplacement, double rDisplacement,
                                     double alphaAngle, PrecisionFactor precisionFactor, ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, rDisplacement, alphaAngle,
                                   precisionFactor, method) * 2 * PI * sineFrequency;
}

double Coil::computeInducedVoltageOn(const Coil &secondary, double zDisplacement, double rDisplacement,
                                     double alphaAngle, MInductanceArguments inductanceArguments,
                                     ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, rDisplacement, alphaAngle,
                                   inductanceArguments, method) * 2 * PI * sineFrequency;
}

double Coil::computeInducedVoltageOn(const Coil &secondary, double zDisplacement, double rDisplacement,
                                     double alphaAngle, double betaAngle, PrecisionFactor precisionFactor,
                                     ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, rDisplacement, alphaAngle, betaAngle,
                                   precisionFactor, method) * 2 * PI * sineFrequency;
}

double Coil::computeInducedVoltageOn(const Coil &secondary, double zDisplacement, double rDisplacement,
                                     double alphaAngle, double betaAngle, MInductanceArguments inductanceArguments,
                                     ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, rDisplacement, alphaAngle, betaAngle,
                                   inductanceArguments, method) * 2 * PI * sineFrequency;
}
