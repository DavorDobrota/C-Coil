#include "Coil.h"
#include "ThreadPool.h"
#include "CUDAFunctions/CoilKernels/CoilAcceleration.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <functional>

using namespace std::placeholders;


namespace
{
    threadPool::ThreadPoolControl g_threadPool;
}


std::vector<double> Coil::calculateAllMutualInductanceArrangementsMTD(Coil primary, Coil secondary,
                                                                      const vec3::Vector3Array &primaryPositions,
                                                                      const vec3::Vector3Array &secondaryPositions,
                                                                      const std::vector<double> &primaryYAngles,
                                                                      const std::vector<double> &primaryZAngles,
                                                                      const std::vector<double> &secondaryYAngles,
                                                                      const std::vector<double> &secondaryZAngles,
                                                                      PrecisionFactor precisionFactor)
{
    size_t arrangementCount = primaryPositions.size();

    std::vector<double> outputMInductances(arrangementCount);

    g_threadPool.setTaskCount(arrangementCount);
    g_threadPool.getCompletedTasks().store(0ull);

    auto calcThread = []
        (
            int idx,
            Coil primary,
            Coil secondary,
            vec3::Vector3 primaryPosition,
            vec3::Vector3 secondaryPosition,
            double primaryYAngle,
            double primaryZAngle,
            double secondaryYAngle,
            double secondaryZAngle,
            PrecisionFactor precisionFactor,
            double &mutualInductance
        )
    {
        primary.setPositionAndOrientation(primaryPosition, primaryYAngle, primaryZAngle);
        secondary.setPositionAndOrientation(secondaryPosition, secondaryYAngle, secondaryZAngle);

        mutualInductance = Coil::computeMutualInductance(primary, secondary, precisionFactor);

        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (int i = 0; i < arrangementCount; ++i)
    {
        g_threadPool.push(calcThread,
                          std::ref(primary),
                          std::ref(secondary),
                          primaryPositions[i],
                          secondaryPositions[i],
                          primaryYAngles[i],
                          primaryZAngles[i],
                          secondaryYAngles[i],
                          secondaryZAngles[i],
                          precisionFactor,
                          std::ref(outputMInductances[i]));
    }

    g_threadPool.synchronizeThreads();

    return outputMInductances;
}


#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"
std::vector<double> Coil::calculateAllMutualInductanceArrangementsGPU(Coil primary, Coil secondary,
                                                                      const vec3::Vector3Array &primaryPositions,
                                                                      const vec3::Vector3Array &secondaryPositions,
                                                                      const std::vector<double> &primaryYAngles,
                                                                      const std::vector<double> &primaryZAngles,
                                                                      const std::vector<double> &secondaryYAngles,
                                                                      const std::vector<double> &secondaryZAngles,
                                                                      PrecisionFactor precisionFactor)
{
    size_t size = primaryPositions.size();

    auto *configArr = static_cast<CoilPairPositionData *>(calloc(size, sizeof(CoilPairPositionData)));
    auto *resultArr = static_cast<TYPE *>(calloc(size, sizeof(TYPE)));

    if(!configArr || !resultArr)
        throw std::bad_alloc();

    for (long long i = 0; i < size; ++i)
    {
        vec3::Vector3 tempPrimPos = primaryPositions[i];
        vec3::Vector3 tempSecPos = secondaryPositions[i];

        configArr[i].primPositionVector[0] = tempPrimPos.x;
        configArr[i].primPositionVector[1] = tempPrimPos.y;
        configArr[i].primPositionVector[2] = tempPrimPos.z;

        configArr[i].secPositionVector[0] = tempSecPos.x;
        configArr[i].secPositionVector[1] = tempSecPos.y;
        configArr[i].secPositionVector[2] = tempSecPos.z;

        configArr[i].primAlphaAngle = primaryYAngles[i];
        configArr[i].primBetaAngle = primaryZAngles[i];

        configArr[i].secAlphaAngle = secondaryYAngles[i];
        configArr[i].secBetaAngle = secondaryZAngles[i];
    }

    CoilPairArguments inductanceArguments = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, GPU);
    CoilPairArgumentsData coilPairArgumentsData;

    generateCoilPairArgumentsData(primary, secondary, coilPairArgumentsData, inductanceArguments, false);
    long long pointCount = inductanceArguments.secondaryPrecision.lengthIncrementCount *
                          inductanceArguments.secondaryPrecision.thicknessIncrementCount *
                          inductanceArguments.secondaryPrecision.angularIncrementCount;

    #if USE_GPU == 1
        Calculate_mutual_inductance_configurations(size, pointCount, coilPairArgumentsData, configArr, resultArr);
    #else
        free(configArr);
        free(resultArr);
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    free(configArr);

    std::vector<double> outputArr;
    outputArr.reserve(size);

    for (long long i = 0; i < size; ++i)
        outputArr.emplace_back(resultArr[i]);

    free(resultArr);

    return outputArr;
}
#pragma clang diagnostic pop