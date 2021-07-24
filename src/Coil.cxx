#include <cmath>
#include <cstdio>
#include <vector>
#include <functional>

#include "Coil.h"
#include "ComputeMethod.h"
#include "hardware_acceleration.h"
#include "Conversion_Util.h"
#include "LegendreMatrix.h"
#include "CoilData.h"

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

PrecisionArguments PrecisionArguments::getCoilPrecisionArgumentsCPU(const Coil &coil, PrecisionFactor precisionFactor)
{
    int lengthArrayIndex = g_minPrimLengthIncrements - 1;
    int thicknessArrayIndex = g_minPrimThicknessIncrements - 1;
    int angularArrayIndex = g_minPrimAngularIncrements - 1;

    int totalIncrements = 0;
    int currentIncrements;
    int caseIndex;

    if (coil.getThickness() / coil.getInnerRadius() < g_thinCoilApproximationRatio &&
        coil.getLength() / coil.getInnerRadius() < g_thinCoilApproximationRatio)
    {
        caseIndex = 1; totalIncrements = pow(2, 3 + precisionFactor.relativePrecision);
    }
    else if (coil.getThickness() / coil.getLength() < g_thinCoilApproximationRatio)
    {
        caseIndex = 2; totalIncrements = pow(2, 6 + precisionFactor.relativePrecision);
    }
    else if (coil.getLength() / coil.getThickness() < g_thinCoilApproximationRatio)
    {
        caseIndex = 3; totalIncrements = pow(2, 6 + precisionFactor.relativePrecision);
    }
    else
    {
        caseIndex = 4; totalIncrements = pow(2, 9 + precisionFactor.relativePrecision);
    }

    totalIncrements *= 2;

    do
    {
        double angularStep = PI * (coil.getInnerRadius() + coil.getThickness() * 0.5) /
                             (blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex]);

        double lengthStep = sqrt(2) * coil.getLength() /
                            (blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex]);

        double thicknessStep = sqrt(2) * coil.getThickness() /
                               (blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex]);

        double linearStep = sqrt(lengthStep * thicknessStep);

        switch (caseIndex)
        {
            case (1):
                thicknessArrayIndex = 0; lengthArrayIndex = 0;
                angularArrayIndex++;
                break;
            case (2):
                thicknessArrayIndex = 0;
                if (angularStep / lengthStep >= 1.0)
                    angularArrayIndex++;
                else
                    lengthArrayIndex++;
                break;
            case (3):
                lengthArrayIndex = 0;
                if (angularStep / thicknessStep >= 1.0)
                    angularArrayIndex++;
                else
                    thicknessArrayIndex++;
                break;
            default:
                if (angularStep / thicknessStep >= 1.0)
                    angularArrayIndex++;
                else
                { lengthArrayIndex++; thicknessArrayIndex++;}
        }

        currentIncrements =
                blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex] *
                blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex] *
                blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex];
    }
    while (currentIncrements < totalIncrements);

    printf("case %d - %d : %d %d %d\n", caseIndex, currentIncrements,
           blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex],
           blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex],
           blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex]);

    return PrecisionArguments(blockPrecisionCPUArray[angularArrayIndex],
                              blockPrecisionCPUArray[thicknessArrayIndex],
                              blockPrecisionCPUArray[lengthArrayIndex],
                              incrementPrecisionCPUArray[angularArrayIndex],
                              incrementPrecisionCPUArray[thicknessArrayIndex],
                              incrementPrecisionCPUArray[lengthArrayIndex]);
}

PrecisionArguments PrecisionArguments::getCoilPrecisionArgumentsGPU(const Coil &coil, PrecisionFactor precisionFactor)
{
    const int linearIncrements = arrSize;
    const int angularIncrements = arrSize;

    int linearBlocks = 1;
    int angularBlocks = 1;

    int totalIncrements = pow(2, 9 + precisionFactor.relativePrecision);
    int currentIncrements;

    do
    {
        double linearStep = sqrt(2 * coil.getLength() * coil.getThickness()) / (linearIncrements * linearBlocks);
        double angularStep = PI * (coil.getInnerRadius() + coil.getThickness() * 0.5) / (angularIncrements * angularBlocks);

        if (angularStep / linearStep >= 1.0)
            angularBlocks++;
        else
            linearBlocks++;

        currentIncrements =
                linearBlocks * linearIncrements * linearIncrements * linearBlocks * angularBlocks * angularIncrements;
    }
    while(currentIncrements < totalIncrements);

    return PrecisionArguments(angularBlocks, linearBlocks, linearBlocks, angularIncrements, linearIncrements, linearIncrements);
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
    int primLengthArrayIndex = g_minPrimLengthIncrements - 1;
    int primThicknessArrayIndex = g_minPrimThicknessIncrements -1;
    int primAngularArrayIndex = g_minPrimAngularIncrements - 1;

    int secLengthArrayIndex = g_minSecLengthIncrements - 1;
    int secThicknessArrayIndex = g_minSecThicknessIncrements - 1;

    int totalIncrements = 0;
    int currentIncrements;
    int caseIndex;

    getMInductanceCaseAndIncrements(primary, secondary, precisionFactor, caseIndex, totalIncrements);

    do
    {
        double primAngularStep = PI * (primary.getInnerRadius() + primary.getThickness() * 0.5) /
                (blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex]);

        double primLengthStep = sqrt(2) * primary.getLength() /
                (blockPrecisionCPUArray[primLengthArrayIndex] * incrementPrecisionCPUArray[primLengthArrayIndex]);

        double primThicknessStep = sqrt(2) * primary.getThickness() /
                (blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex]);

        double secLengthStep = secondary.getLength() /
                (blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex]);

        double secThicknessStep = secondary.getThickness() /
                (blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

        double primLinearStep = sqrt(primLengthStep * primThicknessStep);
        double secLinearStep = sqrt(secLengthStep * secThicknessStep);

        switch (caseIndex)
        {
            case (1):
                primThicknessArrayIndex = 0;
                secThicknessArrayIndex = 0; secLengthArrayIndex = 0;

                if (primAngularStep / primLengthStep >= 1.0)
                    primAngularArrayIndex++;
                else
                    primLengthArrayIndex++;
                break;
            case (2):
                primThicknessArrayIndex = 0;
                secThicknessArrayIndex = 0;

                if (primAngularStep / primLengthStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primLengthStep / secLengthStep >= 1.0)
                        primLengthArrayIndex++;
                    else
                        secLengthArrayIndex++;
                }
                break;
            case (3):
                primThicknessArrayIndex = 0;
                secLengthArrayIndex = 0;

                if (primAngularStep / primLengthStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primLengthStep / secThicknessStep >= 1.0)
                        primLengthArrayIndex++;
                    else
                        secThicknessArrayIndex++;
                }
                break;
            case (4):
                primThicknessArrayIndex = 0;

                if (primAngularStep / primLengthStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primLengthStep / secLinearStep >= 1.0)
                        primLengthArrayIndex++;
                    else
                        { secLengthArrayIndex++; secThicknessArrayIndex++; }
                }
                break;
            case (5):
                primLengthArrayIndex = 0;
                secThicknessArrayIndex = 0; secLengthArrayIndex = 0;

                if (primAngularStep / primThicknessStep >= 1.0)
                    primAngularArrayIndex++;
                else
                    primThicknessArrayIndex++;
                break;
            case (6):
                primLengthArrayIndex = 0;
                secThicknessArrayIndex = 0;

                if (primAngularStep / primThicknessStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primThicknessStep / secLengthStep >= 1.0)
                        primThicknessArrayIndex++;
                    else
                        secLengthArrayIndex++;
                }
                break;
            case (7):
                primLengthArrayIndex = 0;
                secLengthArrayIndex = 0;

                if (primAngularStep / primThicknessStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primThicknessStep / secThicknessStep >= 1.0)
                        primThicknessArrayIndex++;
                    else
                        secThicknessArrayIndex++;
                }
                break;
            case (8):
                primLengthArrayIndex = 0;

                if (primAngularStep / primThicknessStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primThicknessStep / secLinearStep >= 1.0)
                        primThicknessArrayIndex++;
                    else
                    { secLengthArrayIndex++; secThicknessArrayIndex++; }
                }
                break;
            case (9):
                primLengthArrayIndex = 0; primThicknessArrayIndex = 0;
                secLengthArrayIndex = 0; secThicknessArrayIndex = 0;

                primAngularArrayIndex++;
                break;
            case (10):
                primLengthArrayIndex = 0; primThicknessArrayIndex = 0;
                secThicknessArrayIndex = 0;

                if (primAngularStep / secLengthStep >= 1.0)
                    primAngularArrayIndex++;
                else
                    secLengthArrayIndex++;
                break;
            case (11):
                primLengthArrayIndex = 0; primThicknessArrayIndex = 0;
                secLengthArrayIndex = 0;

                if (primAngularStep / secThicknessStep >= 1.0)
                    primAngularArrayIndex++;
                else
                    secThicknessArrayIndex++;
                break;
            case (12):
                primLengthArrayIndex = 0; primThicknessArrayIndex = 0;

                if (primAngularStep / secLinearStep >= 1.0)
                    primAngularArrayIndex++;
                else
                    { secLengthArrayIndex++; secThicknessArrayIndex++; }
                break;
            case (13):
                secLengthArrayIndex = 0; secThicknessArrayIndex = 0;

                if (primAngularStep / primLinearStep >= 1.0)
                    primAngularArrayIndex++;
                else
                    { primLengthArrayIndex++; primThicknessArrayIndex++; }
                break;
            case (14):
                secThicknessArrayIndex = 0;

                if (primAngularStep / primLinearStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primLinearStep / secLengthStep >= 1.0)
                        { primLengthArrayIndex++; primThicknessArrayIndex++; }
                    else
                        secLengthArrayIndex++;
                }
                break;
            case (15):
                secLengthArrayIndex = 0;

                if (primAngularStep / primLinearStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primLinearStep / secThicknessStep >= 1.0)
                    { primLengthArrayIndex++; primThicknessArrayIndex++; }
                    else
                        secThicknessArrayIndex++;
                }
                break;
            default:
                if (primAngularStep / primLinearStep >= 1.0)
                    primAngularArrayIndex++;
                else
                {
                    if (primLinearStep / secLinearStep >= 1.0)
                        { primLengthArrayIndex++; primThicknessArrayIndex++; }
                    else
                        { secLengthArrayIndex++; secThicknessArrayIndex++; }
                }
        }

        currentIncrements =
                blockPrecisionCPUArray[primLengthArrayIndex] * incrementPrecisionCPUArray[primLengthArrayIndex] *
                blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex] *
                blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex] *
                blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex] *
                blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex];
    }
    while (currentIncrements < totalIncrements);

    auto primaryPrecision = PrecisionArguments(blockPrecisionCPUArray[primAngularArrayIndex],
                                               blockPrecisionCPUArray[primThicknessArrayIndex],
                                               blockPrecisionCPUArray[primLengthArrayIndex],
                                               incrementPrecisionCPUArray[primAngularArrayIndex],
                                               incrementPrecisionCPUArray[primThicknessArrayIndex],
                                               incrementPrecisionCPUArray[primLengthArrayIndex]);

    auto secondaryPrecision = PrecisionArguments(0,
                                                 blockPrecisionCPUArray[secThicknessArrayIndex],
                                                 blockPrecisionCPUArray[secLengthArrayIndex],
                                                 0,
                                                 incrementPrecisionCPUArray[secThicknessArrayIndex],
                                                 incrementPrecisionCPUArray[secLengthArrayIndex]);

//    printf("case %d - %d : %d %d %d | %d %d\n", caseIndex, currentIncrements,
//           blockPrecisionCPUArray[primLengthArrayIndex] * incrementPrecisionCPUArray[primLengthArrayIndex],
//           blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex],
//           blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex],
//           blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex],
//           blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

    return MInductanceArguments(primaryPrecision, secondaryPrecision);
}

MInductanceArguments MInductanceArguments::getMInductanceArgumentsGeneralCPU(const Coil &primary, const Coil &secondary,
                                                                             PrecisionFactor precisionFactor)
{
    int primLengthArrayIndex = g_minPrimLengthIncrements - 1;
    int primThicknessArrayIndex = g_minPrimThicknessIncrements - 1;
    int primAngularArrayIndex = g_minPrimAngularIncrements - 1;

    int secLengthArrayIndex = g_minSecLengthIncrements - 1;
    int secThicknessArrayIndex = g_minSecThicknessIncrements - 1;
    int secAngularArrayIndex = g_minSecAngularIncrements - 1;

    int totalIncrements = 0;
    int currentIncrements;
    int caseIndex;

    getMInductanceCaseAndIncrements(primary, secondary, precisionFactor, caseIndex, totalIncrements);
    // multiplying the number of increments by 16 because of one extra dimension of integration compared to z-axis
    totalIncrements *= 16;

    do
    {
        double primAngularStep = PI * (primary.getInnerRadius() + primary.getThickness() * 0.5) /
                (blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex]);

        double primLengthStep = sqrt(2) * primary.getLength() /
                (blockPrecisionCPUArray[primLengthArrayIndex] * incrementPrecisionCPUArray[primLengthArrayIndex]);

        double primThicknessStep = sqrt(2) * primary.getThickness() /
                (blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex]);

        double secLengthStep = secondary.getLength() /
                (blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex]);

        double secThicknessStep = secondary.getThickness() /
                (blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

        double secAngularStep = PI * (secondary.getInnerRadius() + secondary.getThickness() * 0.5) /
                (blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);

        double primLinearStep = sqrt(primLengthStep * primThicknessStep);
        double secLinearStep = sqrt(secLengthStep * secThicknessStep);

        switch (caseIndex)
        {
            case (1):
                primThicknessArrayIndex = 0;
                secThicknessArrayIndex = 0; secLengthArrayIndex = 0;

                if ((primLengthStep * primLengthStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                    primLengthArrayIndex++;
                break;
            case (2):
                primThicknessArrayIndex = 0;
                secThicknessArrayIndex = 0;

                if ((primLengthStep * secLengthStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secLengthStep / primLengthStep <= 1.0)
                        primLengthArrayIndex++;
                    else
                        secLengthArrayIndex++;
                }
                break;
            case (3):
                primThicknessArrayIndex = 0;
                secLengthArrayIndex = 0;

                if ((primLengthStep * secThicknessStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secThicknessStep / primLengthStep <= 1.0)
                        primLengthArrayIndex++;
                    else
                        secThicknessArrayIndex++;
                }
                break;
            case (4):
                primThicknessArrayIndex = 0;

                if ((primLengthStep * secLinearStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secLinearStep / primLengthStep <= 1.0)
                        primLengthArrayIndex++;
                    else
                        { secLengthArrayIndex++; secThicknessArrayIndex++; }
                }
                break;
            case (5):
                primLengthArrayIndex = 0;
                secThicknessArrayIndex = 0; secLengthArrayIndex = 0;

                if ((primThicknessStep * primThicknessStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                    primThicknessArrayIndex++;
                break;
            case (6):
                primLengthArrayIndex = 0;
                secThicknessArrayIndex = 0;

                if ((primThicknessStep * secLengthStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secLengthStep / primThicknessStep <= 1.0)
                        primThicknessArrayIndex++;
                    else
                        secLengthArrayIndex++;
                }
                break;
            case (7):
                primLengthArrayIndex = 0;
                secLengthArrayIndex = 0;

                if ((primThicknessStep * secThicknessStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secThicknessStep / primThicknessStep <= 1.0)
                        primThicknessArrayIndex++;
                    else
                        secThicknessArrayIndex++;
                }
                break;
            case (8):
                primLengthArrayIndex = 0;

                if ((primThicknessStep * secLinearStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secLinearStep / primThicknessStep <= 1.0)
                        primThicknessArrayIndex++;
                    else
                        { secLengthArrayIndex++; secThicknessArrayIndex++; }
                }
                break;
            case (9):
                primLengthArrayIndex = 0; primThicknessArrayIndex = 0;
                secLengthArrayIndex = 0; secThicknessArrayIndex = 0;

                if (secAngularStep / primAngularStep <= 1.0)
                    primAngularArrayIndex++;
                else
                    secAngularArrayIndex++;
                break;
            case (10):
                primLengthArrayIndex = 0; primThicknessArrayIndex = 0;
                secThicknessArrayIndex = 0;

                if ((secLengthStep * secLengthStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                    secLengthArrayIndex++;
                break;
            case (11):
                primLengthArrayIndex = 0; primThicknessArrayIndex = 0;
                secLengthArrayIndex = 0;

                if ((secThicknessStep * secThicknessStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                    secThicknessArrayIndex++;
                break;
            case (12):
                primLengthArrayIndex = 0; primThicknessArrayIndex = 0;

                if ((secLengthStep * secThicknessStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                    { secThicknessArrayIndex++; secThicknessArrayIndex++; }
                break;
            case (13):
                secLengthArrayIndex = 0; secThicknessArrayIndex = 0;

                if ((primLengthStep * primThicknessStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                    { primLengthArrayIndex++; primThicknessArrayIndex++; }
                break;
            case (14):
                secThicknessArrayIndex = 0;

                if ((primLinearStep * secLengthStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secLengthStep / primLinearStep <= 1.0)
                        { primLengthArrayIndex++; primThicknessArrayIndex++; }
                    else
                        secLengthArrayIndex++;
                }
                break;
            case (15):
                secLengthArrayIndex = 0;

                if ((primLinearStep * secThicknessStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularArrayIndex++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secThicknessStep / primLinearStep <= 1.0)
                    { primLengthArrayIndex++; primThicknessArrayIndex++; }
                    else
                        secThicknessArrayIndex++;
                }
                break;
            default:
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
                        { primLengthArrayIndex++; primThicknessArrayIndex++; }
                    else
                        { secLengthArrayIndex++; secThicknessArrayIndex++; }
                }
        }

        currentIncrements =
                blockPrecisionCPUArray[primLengthArrayIndex] * incrementPrecisionCPUArray[primLengthArrayIndex] *
                blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex] *
                blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex] *
                blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex] *
                blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex] *
                blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex];
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

//    printf("case %d - %d : %d %d %d | %d %d %d\n", caseIndex, currentIncrements,
//           blockPrecisionCPUArray[primLengthArrayIndex] * incrementPrecisionCPUArray[primLengthArrayIndex],
//           blockPrecisionCPUArray[primThicknessArrayIndex] * incrementPrecisionCPUArray[primThicknessArrayIndex],
//           blockPrecisionCPUArray[primAngularArrayIndex] * incrementPrecisionCPUArray[primAngularArrayIndex],
//           blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex],
//           blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex],
//           blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);

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

MInductanceArguments MInductanceArguments::getMInductanceArgumentsZGPU(const Coil &primary, const Coil &secondary,
                                                                       PrecisionFactor precisionFactor)
{
    const int primLinearIncrements = arrSize;
    const int primAngularIncrements = arrSize;

    int primLinearBlocks = 1;
    int primAngularBlocks = 1;

    int secLengthArrayIndex = g_minSecLengthIncrements - 1;
    int secThicknessArrayIndex = g_minSecThicknessIncrements - 1;

    double secRadius = secondary.getInnerRadius();
    double secThickness = secondary.getThickness();
    double secLength = secondary.getLength();

    int totalIncrements = 0;
    int currentIncrements;
    int caseIndex;

    if (secThickness / secRadius < g_thinCoilApproximationRatio && secLength / secRadius < g_thinCoilApproximationRatio)
        { caseIndex = 1; totalIncrements = pow(2, 10 + precisionFactor.relativePrecision); }

    else if (secThickness / secLength < g_thinCoilApproximationRatio)
        { caseIndex = 2; totalIncrements = pow(2, 13 + precisionFactor.relativePrecision); }

    else if (secLength / secThickness < g_thinCoilApproximationRatio)
        { caseIndex = 3; totalIncrements = pow(2, 13 + precisionFactor.relativePrecision); }

    else
        { caseIndex = 4; totalIncrements = pow(2, 16 + precisionFactor.relativePrecision); }

    do
    {
        double primAngularStep = PI * (primary.getInnerRadius() + primary.getThickness() * 0.5) /
                                 (primAngularBlocks * primAngularIncrements);

        double primLinearStep = sqrt(primary.getThickness() * primary.getLength()) /
                                (primLinearBlocks * primLinearIncrements);

        double secLengthStep = secondary.getLength() /
                (blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex]);

        double secThicknessStep = secondary.getThickness() /
                (blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

        double secLinearStep = sqrt(secLengthStep * secThicknessStep);

        switch (caseIndex)
        {
            case (1):
                secLengthArrayIndex = 0; secThicknessArrayIndex = 0;
                if (primAngularStep / primLinearStep <= 1.0)
                    primLinearBlocks++;
                else
                    primAngularBlocks++;
                break;
            case (2):
                secThicknessArrayIndex = 0;
                if (primAngularStep / sqrt(secLengthStep * primLinearStep) <= 1.0)
                {
                    if (secLengthStep / primLinearStep >= 1.0)
                        secLengthArrayIndex++;
                    else
                        primLinearBlocks++;
                }
                else
                {
                    primAngularBlocks++;
                }
                break;
            case (3):
                secLengthArrayIndex = 0;
                if (primAngularStep / sqrt(secThicknessStep * primLinearStep) <= 1.0)
                {
                    if (secThicknessStep / primLinearStep >= 1.0)
                        secThicknessArrayIndex++;
                    else
                        primLinearBlocks++;
                }
                else
                    primAngularBlocks++;
                break;
            default:
                if (primAngularStep / sqrt(secLinearStep * primLinearStep) <= 1.0)
                {
                    if (secLinearStep / primLinearStep >= 1.0)
                        { secLengthArrayIndex++; secThicknessArrayIndex++; }
                    else
                        primLinearBlocks++;
                }
                else
                    primAngularBlocks++;
        }
        currentIncrements = primAngularBlocks * primAngularIncrements *
                primLinearBlocks * primLinearIncrements * primLinearBlocks * primLinearIncrements *
                blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex] *
                blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex];
    }
    while (currentIncrements < totalIncrements);

    auto primaryPrecision = PrecisionArguments(primAngularBlocks, primLinearBlocks, primLinearBlocks,
                                               primAngularIncrements, primLinearIncrements, primLinearIncrements);

    auto secondaryPrecision = PrecisionArguments(0,
                                                 blockPrecisionCPUArray[secThicknessArrayIndex],
                                                 blockPrecisionCPUArray[secLengthArrayIndex],
                                                 0,
                                                 incrementPrecisionCPUArray[secThicknessArrayIndex],
                                                 incrementPrecisionCPUArray[secLengthArrayIndex]);

//    printf("case %d - %d : %d %d %d | %d %d\n", caseIndex, currentIncrements,
//           primLinearBlocks * primLinearIncrements,
//           primLinearBlocks * primLinearIncrements,
//           primAngularBlocks * primAngularIncrements,
//           blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex],
//           blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

    return MInductanceArguments(primaryPrecision, secondaryPrecision);
}

MInductanceArguments MInductanceArguments::getMInductanceArgumentsGeneralGPU(const Coil &primary, const Coil &secondary,
                                                                             PrecisionFactor precisionFactor)
{
    const int primLinearIncrements = arrSize;
    const int primAngularIncrements = arrSize;

    int primLinearBlocks = 1;
    int primAngularBlocks = 1;

    int secLengthArrayIndex = g_minSecLengthIncrements - 1;
    int secThicknessArrayIndex = g_minSecThicknessIncrements - 1;
    int secAngularArrayIndex = g_minSecAngularIncrements - 1;

    double secRadius = secondary.getInnerRadius();
    double secThickness = secondary.getThickness();
    double secLength = secondary.getLength();

    int totalIncrements = 0;
    int currentIncrements;
    int caseIndex;

    if (secThickness / secRadius < g_thinCoilApproximationRatio && secLength / secRadius < g_thinCoilApproximationRatio)
    { caseIndex = 1; totalIncrements = pow(2, 13     + precisionFactor.relativePrecision); }

    else if (secThickness / secLength < g_thinCoilApproximationRatio)
    { caseIndex = 2; totalIncrements = pow(2, 16 + precisionFactor.relativePrecision); }

    else if (secLength / secThickness < g_thinCoilApproximationRatio)
    { caseIndex = 3; totalIncrements = pow(2, 16 + precisionFactor.relativePrecision); }

    else
    { caseIndex = 4; totalIncrements = pow(2, 19 + precisionFactor.relativePrecision); }

    do
    {
        double primAngularStep = PI * (primary.getInnerRadius() + primary.getThickness() * 0.5) /
                                 (primAngularBlocks * primAngularIncrements);

        double primLinearStep = sqrt(primary.getThickness() * primary.getLength()) /
                                (primLinearBlocks * primLinearIncrements);

        double secLengthStep = secondary.getLength() /
                               (blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex]);

        double secThicknessStep = secondary.getThickness() /
                (blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex]);

        double secAngularStep = PI * (secondary.getInnerRadius() + secondary.getThickness() * 0.5) /
                (blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);

        double secLinearStep = sqrt(secLengthStep * secThicknessStep);

        switch (caseIndex)
        {
            case (1):
                secLengthArrayIndex = 0; secThicknessArrayIndex = 0;

                if ((primLinearStep * primLinearStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularBlocks++;
                    else
                        secAngularArrayIndex++;
                }
                else
                    primLinearBlocks++;
                break;
            case (2):
                secThicknessArrayIndex = 0;

                if ((primLinearStep * secLengthStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularBlocks++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secLengthStep / primLinearStep <= 1.0)
                        primLinearBlocks++;
                    else
                        secLengthArrayIndex++;
                }
                break;
            case (3):
                secLengthArrayIndex = 0;

                if ((primLinearStep * secThicknessStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularBlocks++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secThicknessStep / primLinearStep <= 1.0)
                        primLinearBlocks++;
                    else
                        secThicknessArrayIndex++;
                }
                break;
            default:
                if ((primLinearStep * secLinearStep) / (primAngularStep * secAngularStep) <= 1.0)
                {
                    if (secAngularStep / primAngularStep <= 1.0)
                        primAngularBlocks++;
                    else
                        secAngularArrayIndex++;
                }
                else
                {
                    if (secLinearStep / primLinearStep <= 1.0)
                        primLinearBlocks++;
                    else
                    { secLengthArrayIndex++; secThicknessArrayIndex++; }
                }
        }
        currentIncrements = primAngularBlocks * primAngularIncrements *
                            primLinearBlocks * primLinearIncrements * primLinearBlocks * primLinearIncrements *
                            blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex] *
                            blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex] *
                            blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex];
    }
    while (currentIncrements < totalIncrements);

    auto primaryPrecision = PrecisionArguments(primAngularBlocks, primLinearBlocks, primLinearBlocks,
                                               primAngularIncrements, primLinearIncrements, primLinearIncrements);

    auto secondaryPrecision = PrecisionArguments(blockPrecisionCPUArray[secAngularArrayIndex],
                                                 blockPrecisionCPUArray[secThicknessArrayIndex],
                                                 blockPrecisionCPUArray[secLengthArrayIndex],
                                                 incrementPrecisionCPUArray[secAngularArrayIndex],
                                                 incrementPrecisionCPUArray[secThicknessArrayIndex],
                                                 incrementPrecisionCPUArray[secLengthArrayIndex]);

//    printf("case %d - %d : %d %d %d | %d %d %d\n", caseIndex, currentIncrements,
//           primLinearBlocks * primLinearIncrements,
//           primLinearBlocks * primLinearIncrements,
//           primAngularBlocks * primAngularIncrements,
//           blockPrecisionCPUArray[secLengthArrayIndex] * incrementPrecisionCPUArray[secLengthArrayIndex],
//           blockPrecisionCPUArray[secThicknessArrayIndex] * incrementPrecisionCPUArray[secThicknessArrayIndex],
//           blockPrecisionCPUArray[secAngularArrayIndex] * incrementPrecisionCPUArray[secAngularArrayIndex]);

    return MInductanceArguments(primaryPrecision, secondaryPrecision);
}

MInductanceArguments MInductanceArguments::getSelfInductanceArguments(const Coil &coil, PrecisionFactor precisionFactor)
{
    int lengthArrayIndex = g_minPrimLengthIncrements - 1;
    int thicknessArrayIndex = g_minPrimThicknessIncrements -1;
    int angularArrayIndex = g_minPrimAngularIncrements - 1;

    int totalIncrements = 0;
    int currentIncrements;
    int caseIndex;

    if (coil.getThickness() / coil.getInnerRadius() < g_thinCoilApproximationRatio &&
        coil.getLength() / coil.getInnerRadius() < g_thinCoilApproximationRatio)
    {
        caseIndex = 1; totalIncrements = pow(2, 5 + precisionFactor.relativePrecision);
    }
    else if (coil.getThickness() / coil.getLength() < g_thinCoilApproximationRatio)
    {
        caseIndex = 2; totalIncrements = pow(2, 11 + precisionFactor.relativePrecision);
    }
    else if (coil.getLength() / coil.getThickness() < g_thinCoilApproximationRatio)
    {
        caseIndex = 3; totalIncrements = pow(2, 11 + precisionFactor.relativePrecision);
    }
    else
    {
        caseIndex = 4; totalIncrements = pow(2, 17 + precisionFactor.relativePrecision);
    }

    do
    {
        double angularStep = PI * (coil.getInnerRadius() + coil.getThickness() * 0.5) /
                (blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex]);

        double lengthStep = sqrt(2) * coil.getLength() /
                (blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex]);

        double thicknessStep = sqrt(2) * coil.getThickness() /
                (blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex]);

        double linearStep = sqrt(lengthStep * thicknessStep);

        switch (caseIndex)
        {
            case (1):
                thicknessArrayIndex = 0; lengthArrayIndex = 0;
                angularArrayIndex++;
                break;
            case (2):
                thicknessArrayIndex = 0;
                if (angularStep / lengthStep >= 1.0)
                    angularArrayIndex++;
                else
                    lengthArrayIndex++;
                break;
            case (3):
                lengthArrayIndex = 0;
                if (angularStep / thicknessStep >= 1.0)
                    angularArrayIndex++;
                else
                    thicknessArrayIndex++;
                break;
            default:
                if (angularStep / thicknessStep >= 1.0)
                    angularArrayIndex++;
                else
                    { lengthArrayIndex++; thicknessArrayIndex++;}
        }

        currentIncrements =
                blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex] *
                blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex] *
                blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex] *
                blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex] *
                blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex];
    }
    while (currentIncrements < totalIncrements);

    auto primaryPrecision = PrecisionArguments(blockPrecisionCPUArray[angularArrayIndex],
                                               blockPrecisionCPUArray[thicknessArrayIndex],
                                               blockPrecisionCPUArray[lengthArrayIndex],
                                               incrementPrecisionCPUArray[angularArrayIndex],
                                               incrementPrecisionCPUArray[thicknessArrayIndex],
                                               incrementPrecisionCPUArray[lengthArrayIndex]);

    auto secondaryPrecision = PrecisionArguments(0,
                                                 blockPrecisionCPUArray[thicknessArrayIndex],
                                                 blockPrecisionCPUArray[lengthArrayIndex],
                                                 0,
                                                 incrementPrecisionCPUArray[thicknessArrayIndex],
                                                 incrementPrecisionCPUArray[lengthArrayIndex]);

//    printf("case %d - %d : %d %d %d | %d %d\n", caseIndex, currentIncrements,
//           blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex],
//           blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex],
//           blockPrecisionCPUArray[angularArrayIndex] * incrementPrecisionCPUArray[angularArrayIndex],
//           blockPrecisionCPUArray[lengthArrayIndex] * incrementPrecisionCPUArray[lengthArrayIndex],
//           blockPrecisionCPUArray[thicknessArrayIndex] * incrementPrecisionCPUArray[thicknessArrayIndex]);

    return MInductanceArguments(primaryPrecision, secondaryPrecision);
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
    // TODO - firstGen solution applied: not very precise but error is significant
//    selfInductance = computeMutualInductance(*this, *this, 0.0, PrecisionFactor(8.0));
    selfInductance = 0.0;

    precisionSettings = PrecisionArguments::getCoilPrecisionArgumentsCPU(*this, PrecisionFactor(8.0));

    double lengthBlock = length / precisionSettings.lengthBlockCount;
    double thicknessBlock = thickness / precisionSettings.thicknessBlockCount;
    double angularBlock = PI / precisionSettings.angularBlockCount;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int lengthIncrements = precisionSettings.lengthIncrementCount - 1;
    int thicknessIncrements = precisionSettings.thicknessIncrementCount - 1;
    int angularIncrements = precisionSettings.angularIncrementCount - 1;

    // multiplication by 2 because cosine is an even function and by 0.125 for a triple change of interval (3 times 1/2)
    double constant = g_MiReduced * currentDensity * lengthBlock * thicknessBlock * angularBlock * 2 * 0.125;

    for (int zBlock = 0; zBlock < precisionSettings.lengthBlockCount; zBlock++)
    {
        for (int rBlock = 0; rBlock < precisionSettings.thicknessBlockCount; rBlock++)
        {
            double zBlockPosition = (-1) * (length * 0.5) + lengthBlock * (zBlock + 0.5);
            double rBlockPosition = innerRadius + thicknessBlock * (rBlock + 0.5);

            for (int zIndex = 0; zIndex <= lengthIncrements; ++zIndex)
            {
                for (int rIndex = 0; rIndex <= thicknessIncrements; ++rIndex)
                {
                    double incrementPositionZ = zBlockPosition +
                            (lengthBlock * 0.5) * Legendre::positionMatrix[lengthIncrements][zIndex];
                    double incrementPositionR = rBlockPosition +
                            (thicknessBlock * 0.5) * Legendre::positionMatrix[thicknessIncrements][rIndex];

                    double potential = 0.0;

                    for (int indBlockL = 0; indBlockL < precisionSettings.lengthBlockCount; ++indBlockL)
                    {
                        for (int indBlockT = 0; indBlockT < precisionSettings.thicknessBlockCount; ++indBlockT)
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

                                    if (!(std::fabs(incrementPositionL - incrementPositionZ) / length < 1e-14 &&
                                        std::fabs(incrementPositionT - incrementPositionR) / thickness < 1e-14))
                                    {
                                        double incrementWeightS =
                                                Legendre::weightsMatrix[lengthIncrements][incL] *
                                                Legendre::weightsMatrix[thicknessIncrements][incT];

                                        double tempConstA = incrementPositionT;
                                        double tempConstB = incrementPositionT * incrementPositionR;
                                        double tempConstC =
                                                incrementPositionT * incrementPositionT +
                                                incrementPositionR * incrementPositionR +
                                                (incrementPositionL + incrementPositionZ) *
                                                (incrementPositionL + incrementPositionZ);

                                        for (int indBlockFi = 0; indBlockFi < precisionSettings.angularBlockCount; ++indBlockFi)
                                        {
                                            double blockPositionFi = angularBlock * (indBlockFi + 0.5);

                                            for (int incFi = 0; incFi <= angularIncrements; ++incFi)
                                            {
                                                double incrementPositionFi = blockPositionFi +
                                                        (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][incFi];

                                                double incrementWeightFi = Legendre::weightsMatrix[angularIncrements][incFi];

                                                double cosinePhi = cos(incrementPositionFi);

                                                potential += constant * incrementWeightS * incrementWeightFi *
                                                        (tempConstA * cosinePhi) /sqrt(tempConstC - 2*tempConstB * cosinePhi);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    selfInductance += 2 * PI * rBlockPosition * potential * 0.25 *
                                      Legendre::weightsMatrix[lengthIncrements][zIndex] *
                                      Legendre::weightsMatrix[thicknessIncrements][rIndex];
                }
            }
        }
    }
    selfInductance *= (numOfTurns / current);
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

double Coil::computeAPotentialZ(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                const PrecisionArguments &usedPrecision) const
{
    // TODO - add functionality in the future
    return 0.0;
}

double Coil::computeAPotentialZ(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeAPotentialZ(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
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

double Coil::computeEFieldZ(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                            const PrecisionArguments &usedPrecision) const
{
    return 2*PI * sineFrequency * computeAPotentialZ(cylindricalZ, cylindricalR, cylindricalPhi, usedPrecision);
}

double Coil::computeEFieldZ(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeEFieldZ(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
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

void Coil::computeAllAPotentialZ(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedPotentialArr,
                                 const PrecisionArguments &usedPrecision,
                                 ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedPotentialArr.resize(0);

        for (int i = 0; i < cylindricalZArr.size(); ++i)
            computedPotentialArr.push_back(0.0);
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllAPotentialZ(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedPotentialArr,
                                 ComputeMethod method) const
{
    computeAllAPotentialZ(
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
    computeAllAPotentialY(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
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

void Coil::computeAllEFieldZ(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             const PrecisionArguments &usedPrecision,
                             ComputeMethod method) const
{
    computeAllAPotentialZ(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                          computedFieldArr, usedPrecision, method);

    for (double & i : computedFieldArr)
        i *= (2 * PI * sineFrequency);
}

void Coil::computeAllEFieldZ(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             ComputeMethod method) const
{
    computeAllEFieldZ(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedFieldArr, precisionSettings, method);
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

    int lengthBlocks = inductanceArguments.secondaryPrecision.lengthBlockCount;
    int lengthIncrements = inductanceArguments.secondaryPrecision.lengthIncrementCount;

    int thicknessBlocks = inductanceArguments.secondaryPrecision.thicknessBlockCount;
    int thicknessIncrements = inductanceArguments.secondaryPrecision.thicknessIncrementCount;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    int maxLengthIndex = lengthIncrements - 1;
    int maxThicknessIndex = thicknessIncrements - 1;

    double lengthBlockSize = secondary.length / lengthBlocks;
    double thicknessBlockSize = secondary.thickness / thicknessBlocks;

    std::vector<double> zPositions;
    std::vector<double> rPositions;

    std::vector<double> weights;

    for (int zBlock = 0; zBlock < lengthBlocks; ++zBlock)
    {
        for (int rBlock = 0; rBlock < thicknessBlocks; ++rBlock)
        {
            double zBlockPosition = (-1) * (secondary.length * 0.5) + lengthBlockSize * (zBlock + 0.5);
            double rBlockPosition = secondary.innerRadius + thicknessBlockSize * (rBlock + 0.5);

            for (int zIndex = 0; zIndex < lengthIncrements; ++zIndex)
            {
                for (int rIndex = 0; rIndex < thicknessIncrements; ++rIndex)
                {
                    double incrementPositionZ = zDisplacement + zBlockPosition +
                            (lengthBlockSize * 0.5) * Legendre::positionMatrix[maxLengthIndex][zIndex];
                    double incrementPositionR = rBlockPosition +
                            (thicknessBlockSize * 0.5) * Legendre::positionMatrix[maxThicknessIndex][rIndex];

                    zPositions.push_back(incrementPositionZ);
                    rPositions.push_back(incrementPositionR);

                    weights.push_back(
                            0.25 * Legendre::weightsMatrix[maxLengthIndex][zIndex] *
                            Legendre::weightsMatrix[maxThicknessIndex][rIndex]);
                }
            }
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

        int lengthBlocks = inductanceArguments.secondaryPrecision.lengthBlockCount;
        int lengthIncrements = inductanceArguments.secondaryPrecision.lengthIncrementCount;

        int thicknessBlocks = inductanceArguments.secondaryPrecision.thicknessBlockCount;
        int thicknessIncrements = inductanceArguments.secondaryPrecision.thicknessIncrementCount;

        int angularBlocks = inductanceArguments.secondaryPrecision.angularBlockCount;
        int angularIncrements = inductanceArguments.secondaryPrecision.angularIncrementCount;

        int numElements = lengthBlocks * lengthIncrements * thicknessBlocks * thicknessIncrements * angularBlocks * angularIncrements;

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

        double lengthBlockSize = secondary.length / lengthBlocks;
        double thicknessBlockSize = secondary.thickness / thicknessBlocks;

        std::vector<double> zPositions;
        std::vector<double> rPositions;

        std::vector<double> weights;

        for (int zBlock = 0; zBlock < lengthBlocks; ++zBlock)
        {
            for (int rBlock = 0; rBlock < thicknessBlocks; ++rBlock)
            {
                double zBlockPosition = (-1) * (secondary.length * 0.5) + lengthBlockSize * (zBlock + 0.5);
                double rBlockPosition = secondary.innerRadius + thicknessBlockSize * (rBlock + 0.5);

                for (int zIndex = 0; zIndex < lengthIncrements; ++zIndex)
                {
                    for (int rIndex = 0; rIndex < thicknessIncrements; ++rIndex)
                    {
                        double ringRadius = rBlockPosition +
                                (thicknessBlockSize * 0.5) * Legendre::positionMatrix[maxThicknessIndex][rIndex];

                        double lengthDisplacement = zBlockPosition +
                                (lengthBlockSize * 0.5) * Legendre::positionMatrix[maxLengthIndex][zIndex];

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

MInductanceArguments Coil::calculateAppropriateMInductanceArguments(const Coil &primary, const Coil &secondary,
                                                                    PrecisionFactor precisionFactor,
                                                                    ComputeMethod method, bool isGeneral)
{
    if (!isGeneral)
    {
        if (method == GPU)
            return MInductanceArguments::getMInductanceArgumentsZGPU(primary, secondary, precisionFactor);
        else
            return MInductanceArguments::getMInductanceArgumentsZCPU(primary, secondary, precisionFactor);
    }
    else
    {
        if (method == GPU)
            return MInductanceArguments::getMInductanceArgumentsGeneralGPU(primary, secondary, precisionFactor);
        else
            return MInductanceArguments::getMInductanceArgumentsGeneralCPU(primary, secondary, precisionFactor);
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
    auto args = calculateAppropriateMInductanceArguments(primary, secondary, precisionFactor, method, false);
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
    auto args = calculateAppropriateMInductanceArguments(primary, secondary, precisionFactor, method);
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
    auto args = calculateAppropriateMInductanceArguments(primary, secondary, precisionFactor, method);
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
    auto args = calculateAppropriateMInductanceArguments(primary, secondary, precisionFactor, method);
    return computeMutualInductance(primary, secondary, zDisplacement, rDisplacement, alphaAngle, betaAngle, args, method);
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement,
                                            PrecisionFactor precisionFactor, ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, precisionFactor, method) *
           2 * PI * sineFrequency;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement,
                                            MInductanceArguments inductanceArguments, ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, inductanceArguments, method) *
           2 * PI * sineFrequency;;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement, double rDisplacement,
                                            PrecisionFactor precisionFactor, ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, rDisplacement, precisionFactor, method) *
           2 * PI * sineFrequency;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement, double rDisplacement,
                                            MInductanceArguments inductanceArguments, ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, rDisplacement, inductanceArguments, method) *
           2 * PI * sineFrequency;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement, double rDisplacement,
                                            double alphaAngle, PrecisionFactor precisionFactor, ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, rDisplacement, alphaAngle,
                                   precisionFactor, method) * 2 * PI * sineFrequency;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement, double rDisplacement,
                                            double alphaAngle, MInductanceArguments inductanceArguments,
                                            ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, rDisplacement, alphaAngle,
                                   inductanceArguments, method) * 2 * PI * sineFrequency;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement, double rDisplacement,
                                            double alphaAngle, double betaAngle, PrecisionFactor precisionFactor,
                                            ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, rDisplacement, alphaAngle, betaAngle,
                                   precisionFactor, method) * 2 * PI * sineFrequency;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement, double rDisplacement,
                                            double alphaAngle, double betaAngle, MInductanceArguments inductanceArguments,
                                            ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, zDisplacement, rDisplacement, alphaAngle, betaAngle,
                                   inductanceArguments, method) * 2 * PI * sineFrequency;
}
