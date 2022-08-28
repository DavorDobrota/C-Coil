#include "Coil.h"
#include "CoilAcceleration.h"


#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"
vec3::Vector3Array Coil::calculateAllAPotentialGPU(const vec3::Vector3Array &pointVectors,
                                                   const PrecisionArguments &usedPrecision) const
{
    long long size = pointVectors.size();

    auto *coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *resultArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));

    if(!coordinateArr || !resultArr)
        throw std::bad_alloc();

    for (long long i = 0; i < size; ++i)
    {
        vec3::Vector3 tempVec = pointVectors[i];

        coordinateArr[i].x = tempVec.x;
        coordinateArr[i].y = tempVec.y;
        coordinateArr[i].z = tempVec.z;
    }

    CoilData coilData;
    generateCoilData(coilData, usedPrecision);

    #if USE_GPU == 1
        Calculate_hardware_accelerated_a(size, coilData, coordinateArr, resultArr);
    #else
        free(coordinateArr);
        free(resultArr);
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    free(coordinateArr);

    vec3::Vector3Array computedPotentialArr;
    computedPotentialArr.reserve(size);
    std::vector<vec3::Vector3> &outputRef = computedPotentialArr.getItems();

    for (long long i = 0; i < pointVectors.size(); ++i)
        outputRef.emplace_back(resultArr[i].x, resultArr[i].y, resultArr[i].z);

    free(resultArr);

    return computedPotentialArr;
}

vec3::Vector3Array Coil::calculateAllBFieldGPU(const vec3::Vector3Array &pointVectors,
                                               const PrecisionArguments &usedPrecision) const
{
    long long size = pointVectors.size();

    auto *coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *resultArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));

    if(!coordinateArr || !resultArr)
        throw std::bad_alloc();

    for (long long i = 0; i < size; ++i)
    {
        vec3::Vector3 tempVec = pointVectors[i];

        coordinateArr[i].x = tempVec.x;
        coordinateArr[i].y = tempVec.y;
        coordinateArr[i].z = tempVec.z;
    }

    CoilData coilData;
    generateCoilData(coilData, usedPrecision);

    #if USE_GPU == 1
        Calculate_hardware_accelerated_b(size, coilData, coordinateArr, resultArr);
    #else
        free(coordinateArr);
        free(resultArr);
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    free(coordinateArr);

    vec3::Vector3Array computedFieldArr;
    computedFieldArr.reserve(pointVectors.size());
    std::vector<vec3::Vector3> &outputRef = computedFieldArr.getItems();

    for (long long i = 0; i < pointVectors.size(); ++i)
        computedFieldArr.append(resultArr[i].x, resultArr[i].y, resultArr[i].z);

    free(resultArr);

    return computedFieldArr;
}

vec3::Matrix3Array Coil::calculateAllBGradientGPU(const vec3::Vector3Array &pointVectors,
                                                  const PrecisionArguments &usedPrecision) const
{
    long long size = pointVectors.size();

    auto *coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *resultArr = static_cast<DataMatrix *>(calloc(size, sizeof(DataMatrix)));

    if(!coordinateArr || !resultArr)
        throw std::bad_alloc();

    for (long long i = 0; i < size; ++i)
    {
        vec3::Vector3 tempVec = pointVectors[i];

        coordinateArr[i].x = tempVec.x;
        coordinateArr[i].y = tempVec.y;
        coordinateArr[i].z = tempVec.z;
    }

    CoilData coilData;
    generateCoilData(coilData, usedPrecision);

    #if USE_GPU == 1
        Calculate_hardware_accelerated_g(size, coilData, coordinateArr, resultArr);
    #else
        free(coordinateArr);
        free(resultArr);
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    free(coordinateArr);

    vec3::Matrix3Array computedGradientArr;
    computedGradientArr.reserve(size);
    std::vector<vec3::Matrix3> &outputRef = computedGradientArr.getItems();

    for (long long i = 0; i < pointVectors.size(); ++i)
        outputRef.emplace_back(resultArr[i].xx, resultArr[i].xy, resultArr[i].xz,
                               resultArr[i].yx, resultArr[i].yy, resultArr[i].yz,
                               resultArr[i].zx, resultArr[i].zy, resultArr[i].zz);

    free(resultArr);

    return computedGradientArr;
}
#pragma clang diagnostic pop