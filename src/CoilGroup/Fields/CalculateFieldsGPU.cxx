#include "CoilGroup.h"
#include "CoilGroupAcceleration.h"

#define _USE_MATH_DEFINES
#include <math.h>

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"


vec3::Vector3Array CoilGroup::calculateAllAPotentialGPU(const vec3::Vector3Array &pointVectors) const
{
    long long size = pointVectors.size();
    long long coils = memberCoils.size();

    auto *coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *resultArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *coilDataArr = static_cast<CoilData *>(calloc(coils, sizeof(CoilData)));

    if(!coordinateArr || !resultArr || !coilDataArr)
        throw std::bad_alloc();

    for (long long i = 0; i < size; ++i)
    {
        vec3::Vector3 tempVec = pointVectors[i];

        coordinateArr[i].x = tempVec.x;
        coordinateArr[i].y = tempVec.y;
        coordinateArr[i].z = tempVec.z;
    }

    generateCoilDataArray(coilDataArr, this->defaultPrecisionFactor);

    #if USE_GPU == 1
        Calculate_hardware_accelerated_a_group(coils, size, coilDataArr, coordinateArr, resultArr);
    #else
        free(coordinateArr);
        free(resultArr);
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    free(coordinateArr);

    vec3::Vector3Array computedPotentialArr;
    computedPotentialArr.reserve(size);

    for (long long i = 0; i < pointVectors.size(); ++i)
        computedPotentialArr.append(resultArr[i].x, resultArr[i].y, resultArr[i].z);

    free(resultArr);

    return computedPotentialArr;
}


vec3::Vector3Array CoilGroup::calculateAllBFieldGPU(const vec3::Vector3Array &pointVectors) const
{
    long long size = pointVectors.size();
    long long coils = memberCoils.size();

    auto *coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *resultArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *coilDataArr = static_cast<CoilData *>(calloc(coils, sizeof(CoilData)));

    if(!coordinateArr || !resultArr || !coilDataArr)
        throw std::bad_alloc();

    for (long long i = 0; i < size; ++i)
    {
        vec3::Vector3 tempVec = pointVectors[i];

        coordinateArr[i].x = tempVec.x;
        coordinateArr[i].y = tempVec.y;
        coordinateArr[i].z = tempVec.z;
    }

    generateCoilDataArray(coilDataArr, this->defaultPrecisionFactor);

    #if USE_GPU == 1
        Calculate_hardware_accelerated_b_group(coils, size, coilDataArr, coordinateArr, resultArr);
    #else
        free(coordinateArr);
        free(resultArr);
        free(coilDataArr);
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    free(coordinateArr);

    vec3::Vector3Array computedFieldArr;
    computedFieldArr.reserve(size);

    for (long long i = 0; i < pointVectors.size(); ++i)
        computedFieldArr.append(resultArr[i].x, resultArr[i].y, resultArr[i].z);

    free(resultArr);

    return computedFieldArr;
}


vec3::Vector3Array CoilGroup::calculateAllEFieldGPU(const vec3::Vector3Array &pointVectors) const
{
    long long size = pointVectors.size();
    long long coils = memberCoils.size();

    auto *coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *resultArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *coilDataArr = static_cast<CoilData *>(calloc(coils, sizeof(CoilData)));

    if(!coordinateArr || !resultArr || !coilDataArr)
        throw std::bad_alloc();

    for (long long i = 0; i < size; ++i)
    {
        vec3::Vector3 tempVec = pointVectors[i];

        coordinateArr[i].x = tempVec.x;
        coordinateArr[i].y = tempVec.y;
        coordinateArr[i].z = tempVec.z;
    }

    generateCoilDataArray(coilDataArr, this->defaultPrecisionFactor);

    #if USE_GPU == 1
        Calculate_hardware_accelerated_e_group(coils, size, coilDataArr, coordinateArr, resultArr);
    #else
        free(coordinateArr);
        free(resultArr);
        free(coilDataArr);
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    free(coordinateArr);

    vec3::Vector3Array computedFieldArr;
    computedFieldArr.reserve(size);

    for (long long i = 0; i < pointVectors.size(); ++i)
        computedFieldArr.append(resultArr[i].x, resultArr[i].y, resultArr[i].z);

    free(resultArr);

    return computedFieldArr;
}


vec3::Matrix3Array CoilGroup::calculateAllBGradientGPU(const vec3::Vector3Array &pointVectors) const
{
    long long size = pointVectors.size();
    long long coils = memberCoils.size();

    auto *coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *resultArr = static_cast<DataMatrix *>(calloc(size, sizeof(DataMatrix)));
    auto *coilDataArr = static_cast<CoilData *>(calloc(coils, sizeof(CoilData)));

    if(!coordinateArr || !resultArr || !coilDataArr)
        throw std::bad_alloc();

    for (long long i = 0; i < size; ++i)
    {
        vec3::Vector3 tempVec = pointVectors[i];

        coordinateArr[i].x = tempVec.x;
        coordinateArr[i].y = tempVec.y;
        coordinateArr[i].z = tempVec.z;
    }

    generateCoilDataArray(coilDataArr, this->defaultPrecisionFactor);

    #if USE_GPU == 1
        Calculate_hardware_accelerated_g_group(coils, size, coilDataArr, coordinateArr, resultArr);
    #else
        free(coordinateArr);
        free(resultArr);
        free(coilDataArr);
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    free(coordinateArr);

    vec3::Matrix3Array computedGradientArr;
    computedGradientArr.reserve(size);

    for (long long i = 0; i < pointVectors.size(); ++i)
        computedGradientArr.append(resultArr[i].xx, resultArr[i].xy, resultArr[i].xz,
                                   resultArr[i].yx, resultArr[i].yy, resultArr[i].yz,
                                   resultArr[i].zx, resultArr[i].zy, resultArr[i].zz);

    free(resultArr);

    return computedGradientArr;

}

#pragma clang diagnostic pop
