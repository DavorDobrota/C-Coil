#include "CoilGroup.h"
#include "CoilGroupAcceleration.h"

#define _USE_MATH_DEFINES
#include <math.h>

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"


std::vector<double> CoilGroup::calculateAllMutualInductanceArrangementsGPU(const Coil &secondary,
                                                                           const vec3::Vector3Array &secondaryPositions,
                                                                           const std::vector<double> &secondaryYAngles,
                                                                           const std::vector<double> &secondaryZAngles,
                                                                           PrecisionFactor precisionFactor) const
{
    size_t size = secondaryPositions.size();
    auto secPrecision = PrecisionArguments::getCoilPrecisionArgumentsGPU(secondary, precisionFactor);

    auto *configArr = static_cast<SecondaryCoilPositionData *>(calloc(size, sizeof(SecondaryCoilPositionData)));
    auto *resultArr = static_cast<TYPE *>(calloc(size, sizeof(TYPE)));

    if (!configArr || !resultArr)
        throw std::bad_alloc();

    for (long long i = 0; i < size; ++i)
    {
        vec3::Vector3 tempSecPos = secondaryPositions[i];

        configArr[i].positionVector[0] = tempSecPos.x;
        configArr[i].positionVector[1] = tempSecPos.y;
        configArr[i].positionVector[2] = tempSecPos.z;

        configArr[i].alphaAngle = secondaryYAngles[i];
        configArr[i].betaAngle = secondaryZAngles[i];
    }

    size_t coilArrSize = memberCoils.size();
    bool removeSpecificCoil = false;
    for (const auto& memberCoil : memberCoils)
        if (memberCoil->getId() == secondary.getId())
        {
            --coilArrSize; removeSpecificCoil = true;
            break;
        }

    auto *coilArr = static_cast<CoilData *>(calloc(coilArrSize, sizeof(CoilData)));
    generateCoilDataArray(coilArr, removeSpecificCoil, secondary.getId());

    SecondaryCoilData secondaryData;
    generateSecondaryData(secondary, secondaryData, secPrecision, false);

    long long pointCount = secPrecision.lengthIncrementCount *
                           secPrecision.thicknessIncrementCount *
                           secPrecision.angularIncrementCount;

    #if USE_GPU == 1
        Calculate_mutual_inductance_configurations_group(
            coilArrSize, size, pointCount,
            secondaryData, coilArr, configArr, resultArr
        );
    #else
        free(configArr);
        free(resultArr);
        free(coilArr);
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    free(configArr);
    free(coilArr);

    std::vector<double> outputArr;
    outputArr.reserve(size);

    for (long long i = 0; i < size; ++i)
        outputArr.emplace_back(resultArr[i]);

    free(resultArr);

    return outputArr;
}