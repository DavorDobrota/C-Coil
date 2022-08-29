#include "Coil.h"

#include "CoilAcceleration.h"
#include "ThreadPool.h"
#include "Timing.h"


namespace
{
    threadPool::ThreadPoolControl g_threadPool;
}


#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"
vec3::Vector3Array Coil::calculateAllAPotentialGPU(const vec3::Vector3Array &pointVectors,
                                                   const PrecisionArguments &usedPrecision) const
{
    #if PRINT_TIMINGS == 1
        double interval;

        recordStartPoint();
        recordStartPoint();
    #endif //PRINT_TIMINGS

    long long size = pointVectors.size();

    auto *coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *resultArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));

    if(!coordinateArr || !resultArr)
        throw std::bad_alloc();

    #if PRINT_TIMINGS == 1
        interval = getIntervalDuration();
        printf("\tInitialising input:      %.9g s | %.3g GB/s\n",
               interval, 1e-9 * double(size * (sizeof(DataVector) + sizeof(vec3::Vector3))) / interval
        );
        recordStartPoint();
    #endif //PRINT_TIMINGS

    for (long long i = 0; i < size; ++i)
    {
        vec3::Vector3 tempVec = pointVectors[i];

        coordinateArr[i].x = tempVec.x;
        coordinateArr[i].y = tempVec.y;
        coordinateArr[i].z = tempVec.z;
    }

    #if PRINT_TIMINGS == 1
        interval = getIntervalDuration();
        printf("\tPreparing input array:   %.9g s | %.3g GB/s\n",
               interval, 1e-9 * double(size * (sizeof(DataVector) + sizeof(vec3::Vector3))) / interval
        );
        recordStartPoint();
    #endif //PRINT_TIMINGS

    CoilData coilData;
    generateCoilData(coilData, usedPrecision);

    #if PRINT_TIMINGS == 1
        interval = getIntervalDuration();
        printf("\tPreparing coil data:     %.9g s\n", interval);
        recordStartPoint();
    #endif //PRINT_TIMINGS

    #if USE_GPU == 1
        Calculate_hardware_accelerated_a(size, coilData, coordinateArr, resultArr);
    #else
        free(coordinateArr);
        free(resultArr);
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    #if PRINT_TIMINGS == 1
        interval = getIntervalDuration();
        printf("\tCalculations:            %.9g s\n", interval);
        recordStartPoint();
    #endif //PRINT_TIMINGS

    free(coordinateArr);

    vec3::Vector3Array computedPotentialArr;
    computedPotentialArr.reserve(pointVectors.size());
    std::vector<vec3::Vector3> &outputRef = computedPotentialArr.getItems();

    #if PRINT_TIMINGS == 1
        interval = getIntervalDuration();
        printf("\tInitialising output:     %.9g s | %.3g GB/s\n",
               interval, 1e-9 * double(size * sizeof(vec3::Vector3)) / interval
        );
        recordStartPoint();
    #endif //PRINT_TIMINGS

    for (long long i = 0; i < pointVectors.size(); ++i)
        outputRef.emplace_back(resultArr[i].x, resultArr[i].y, resultArr[i].z);

    free(resultArr);

    #if PRINT_TIMINGS == 1
        interval = getIntervalDuration();
        printf("\tPreparing output array:  %.9g s | %.3g GB/s\n",
               interval, 1e-9 * double(size * (sizeof(DataVector) + sizeof(vec3::Vector3))) / interval
        );
        interval = getIntervalDuration();
        printf("\nTotal time:                  %.9g s\n", interval);
        printf("=======================================================\n");
    #endif //PRINT_TIMINGS

    return computedPotentialArr;
}

vec3::Vector3Array Coil::calculateAllBFieldGPU(const vec3::Vector3Array &pointVectors,
                                               const PrecisionArguments &usedPrecision) const
{
    #if PRINT_TIMINGS == 1
        double interval;

        recordStartPoint();
        recordStartPoint();
    #endif //PRINT_TIMINGS
    long long size = pointVectors.size();

    auto *coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *resultArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));

    if(!coordinateArr || !resultArr)
        throw std::bad_alloc();

    #if PRINT_TIMINGS == 1
        interval = getIntervalDuration();
        printf("\tInitialising input:      %.9g s | %.3g GB/s\n",
               interval, 1e-9 * double(size * (sizeof(DataVector) + sizeof(vec3::Vector3))) / interval
        );
        recordStartPoint();
    #endif //PRINT_TIMINGS

    for (long long i = 0; i < size; ++i)
    {
        vec3::Vector3 tempVec = pointVectors[i];

        coordinateArr[i].x = tempVec.x;
        coordinateArr[i].y = tempVec.y;
        coordinateArr[i].z = tempVec.z;
    }

    #if PRINT_TIMINGS == 1
        interval = getIntervalDuration();
        printf("\tPreparing input array:   %.9g s | %.3g GB/s\n",
               interval, 1e-9 * double(size * (sizeof(DataVector) + sizeof(vec3::Vector3))) / interval
        );
        recordStartPoint();
    #endif //PRINT_TIMINGS

    CoilData coilData;
    generateCoilData(coilData, usedPrecision);

    #if PRINT_TIMINGS == 1
        interval = getIntervalDuration();
        printf("\tPreparing coil data:     %.9g s\n", interval);
        recordStartPoint();
    #endif //PRINT_TIMINGS

    #if USE_GPU == 1
        Calculate_hardware_accelerated_b(size, coilData, coordinateArr, resultArr);
    #else
        free(coordinateArr);
        free(resultArr);
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    #if PRINT_TIMINGS == 1
        interval = getIntervalDuration();
        printf("\tCalculations:            %.9g s\n", interval);
        recordStartPoint();
    #endif //PRINT_TIMINGS

    free(coordinateArr);

    vec3::Vector3Array computedFieldArr;
    computedFieldArr.reserve(pointVectors.size());
    std::vector<vec3::Vector3> &outputRef = computedFieldArr.getItems();

    #if PRINT_TIMINGS == 1
        interval = getIntervalDuration();
        printf("\tInitialising output:     %.9g s | %.3g GB/s\n",
               interval, 1e-9 * double(size * sizeof(vec3::Vector3)) / interval
        );
        recordStartPoint();
    #endif //PRINT_TIMINGS

    for (long long i = 0; i < pointVectors.size(); ++i)
        outputRef.emplace_back(resultArr[i].x, resultArr[i].y, resultArr[i].z);

    free(resultArr);

    #if PRINT_TIMINGS == 1
        interval = getIntervalDuration();
        printf("\tPreparing output array:  %.9g s | %.3g GB/s\n",
               interval, 1e-9 * double(size * (sizeof(DataVector) + sizeof(vec3::Vector3))) / interval
        );
        interval = getIntervalDuration();
        printf("\nTotal time:                  %.9g s\n", interval);
        printf("=======================================================\n");
    #endif //PRINT_TIMINGS

    return computedFieldArr;
}

vec3::Matrix3Array Coil::calculateAllBGradientGPU(const vec3::Vector3Array &pointVectors,
                                                  const PrecisionArguments &usedPrecision) const
{
    #if PRINT_TIMINGS == 1
        double interval;

        recordStartPoint();
        recordStartPoint();
    #endif //PRINT_TIMINGS

    long long size = pointVectors.size();

    auto *coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *resultArr = static_cast<DataMatrix *>(calloc(size, sizeof(DataMatrix)));

    if(!coordinateArr || !resultArr)
        throw std::bad_alloc();

     #if PRINT_TIMINGS == 1
        interval = getIntervalDuration();
        printf("\tInitialising input:      %.9g s | %.3g GB/s\n",
               interval, 1e-9 * double(size * (sizeof(DataVector) + sizeof(vec3::Vector3))) / interval
        );
        recordStartPoint();
    #endif //PRINT_TIMINGS

    for (long long i = 0; i < size; ++i)
    {
        vec3::Vector3 tempVec = pointVectors[i];

        coordinateArr[i].x = tempVec.x;
        coordinateArr[i].y = tempVec.y;
        coordinateArr[i].z = tempVec.z;
    }

    #if PRINT_TIMINGS == 1
        interval = getIntervalDuration();
        printf("\tPreparing input array:   %.9g s | %.3g GB/s\n",
               interval, 1e-9 * double(size * (sizeof(DataVector) + sizeof(vec3::Vector3))) / interval
        );
        recordStartPoint();
    #endif //PRINT_TIMINGS

    CoilData coilData;
    generateCoilData(coilData, usedPrecision);

    #if PRINT_TIMINGS == 1
        interval = getIntervalDuration();
        printf("\tPreparing coil data:     %.9g s\n", interval);
        recordStartPoint();
    #endif //PRINT_TIMINGS

    #if USE_GPU == 1
        Calculate_hardware_accelerated_g(size, coilData, coordinateArr, resultArr);
    #else
        free(coordinateArr);
        free(resultArr);
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    #if PRINT_TIMINGS == 1
        interval = getIntervalDuration();
        printf("\tCalculations:            %.9g s\n", interval);
        recordStartPoint();
    #endif //PRINT_TIMINGS

    free(coordinateArr);

    vec3::Matrix3Array computedGradientArr;
    computedGradientArr.reserve(size);
    std::vector<vec3::Matrix3> &outputRef = computedGradientArr.getItems();

    #if PRINT_TIMINGS == 1
        interval = getIntervalDuration();
        printf("\tInitialising output:     %.9g s | %.3g GB/s\n",
               interval, 1e-9 * double(size * sizeof(vec3::Matrix3)) / interval
        );
        recordStartPoint();
    #endif //PRINT_TIMINGS

    for (long long i = 0; i < pointVectors.size(); ++i)
        outputRef.emplace_back(resultArr[i].xx, resultArr[i].xy, resultArr[i].xz,
                               resultArr[i].yx, resultArr[i].yy, resultArr[i].yz,
                               resultArr[i].zx, resultArr[i].zy, resultArr[i].zz);

    free(resultArr);

    #if PRINT_TIMINGS == 1
        interval = getIntervalDuration();
        printf("\tPreparing output array:  %.9g s | %.3g GB/s\n",
               interval, 1e-9 * double(size * (sizeof(DataMatrix) + sizeof(vec3::Matrix3))) / interval
        );
        interval = getIntervalDuration();
        printf("\nTotal time:                  %.9g s\n", interval);
        printf("=======================================================\n");
    #endif //PRINT_TIMINGS

    return computedGradientArr;
}
#pragma clang diagnostic pop