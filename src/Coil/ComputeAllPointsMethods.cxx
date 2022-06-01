#include "Coil.h"
#include "ThreadPool.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <numeric>
#include <cstdio>

namespace
{
    const int pointMultiplier = 1024;
    threadPool::ThreadPoolControl g_threadPool;
}


void Coil::adaptInputVectorsForAllPoints(const std::vector<vec3::CoordVector3> &pointVectors,
                                         std::vector<double> &cylindricalZArr,
                                         std::vector<double> &cylindricalRArr,
                                         std::vector<double> &cylindricalPhiArr) const
{
    cylindricalZArr.resize(pointVectors.size());
    cylindricalRArr.resize(pointVectors.size());
    cylindricalPhiArr.resize(pointVectors.size());

    const int chunkSize = pointVectors.size() / (2 * threadCount);

    if (pointVectors.size() > pointMultiplier * threadCount)
    {
//        printf("Made it here\n");
        g_threadPool.setTaskCount(cylindricalZArr.size());
        g_threadPool.getCompletedTasks().store(0ull);

        auto calcThread = [](
                int idx,
                Coil coil,
                const std::vector<vec3::CoordVector3> &pointVectors,
                std::vector<double> &cylindricalZArr,
                std::vector<double> &cylindricalRArr,
                std::vector<double> &cylindricalPhiArr,
                size_t startIdx, size_t stopIdx
        ) -> void
        {
            for(size_t i = startIdx; i < stopIdx; i++)
            {
                auto result = coil.adaptInputVectorForPoint(pointVectors[i]);

                cylindricalZArr[i] = result.comp1;
                cylindricalRArr[i] = result.comp2;
                cylindricalPhiArr[i] = result.comp3;

                g_threadPool.getCompletedTasks().fetch_add(1ull);
            }
        };

        for(size_t i = 0; i < (size_t)std::ceil((double)cylindricalZArr.size() / (double)chunkSize); i++)
        {
            g_threadPool.push(
                    calcThread,
                    std::ref(*this),
                    std::ref(pointVectors),
                    std::ref(cylindricalZArr), std::ref(cylindricalRArr), std::ref(cylindricalPhiArr),
                    i * chunkSize, std::min((i + 1) * chunkSize, cylindricalZArr.size())
            );
        }

        g_threadPool.synchronizeThreads();
    }
    else {
        vec3::FieldVector3 positionVec = vec3::CoordVector3::convertToFieldVector(positionVector);

        for (int i = 0; i < pointVectors.size(); ++i) {
            vec3::FieldVector3 pointVec = vec3::CoordVector3::convertToFieldVector(pointVectors[i]);
            vec3::FieldVector3 transformedVec = inverseTransformationMatrix * (pointVec - positionVec);

            vec3::CoordVector3 finalVec = vec3::CoordVector3::convertToCoordVector(transformedVec);
            finalVec.convertToCylindrical();

            cylindricalZArr[i] = finalVec.comp1;
            cylindricalRArr[i] = finalVec.comp2;
            cylindricalPhiArr[i] = finalVec.comp3;
        }
    }
}

std::vector<vec3::FieldVector3> Coil::adaptOutputVectorsForAllPoints(const std::vector<vec3::FieldVector3> &computedVectorArr) const
{
    std::vector<vec3::FieldVector3> outputVectorArr(computedVectorArr.size());

    for (int i = 0; i < computedVectorArr.size(); ++i)
        outputVectorArr[i] = transformationMatrix * computedVectorArr[i];

    return outputVectorArr;
}


std::vector<vec3::FieldVector3> Coil::computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                 const PrecisionArguments &usedPrecision,
                                                                 ComputeMethod computeMethod) const
{
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;
    std::vector<double> fieldH, fieldZ;
    std::vector<vec3::FieldVector3> computedFieldArr(pointVectors.size());

    adaptInputVectorsForAllPoints(pointVectors, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);
    calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision, computeMethod);

    for (int i = 0; i < pointVectors.size(); ++i)
        computedFieldArr[i] = vec3::FieldVector3(fieldH[i] * std::cos(cylindricalPhiArr[i]),
                                                 fieldH[i] * std::sin(cylindricalPhiArr[i]),
                                                 fieldZ[i]);

    return adaptOutputVectorsForAllPoints(computedFieldArr);
}

std::vector<vec3::FieldVector3> Coil::computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                 ComputeMethod computeMethod) const
{
    return computeAllBFieldComponents(pointVectors, defaultPrecision, computeMethod);
}

std::vector<double>Coil::computeAllBFieldX(const std::vector<vec3::CoordVector3> &pointVectors,
                                           const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(pointVectors, usedPrecision, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].x;

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldX(const std::vector<vec3::CoordVector3> &pointVectors,
                                            ComputeMethod computeMethod) const
{
    return computeAllBFieldX(pointVectors, defaultPrecision, computeMethod);
}

std::vector<double>Coil::computeAllBFieldY(const std::vector<vec3::CoordVector3> &pointVectors,
                                           const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(pointVectors, usedPrecision, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].y;

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldY(const std::vector<vec3::CoordVector3> &pointVectors,
                                           ComputeMethod computeMethod) const
{
    return computeAllBFieldY(pointVectors, defaultPrecision, computeMethod);
}

std::vector<double>Coil::computeAllBFieldZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                           const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(pointVectors, usedPrecision, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].z;

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                            ComputeMethod computeMethod) const
{
    return computeAllBFieldZ(pointVectors, defaultPrecision, computeMethod);
}

std::vector<double>Coil::computeAllBFieldAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                             const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(pointVectors, usedPrecision, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = std::sqrt(computedFieldArr[i].x * computedFieldArr[i].x +
                                 computedFieldArr[i].y * computedFieldArr[i].y +
                                 computedFieldArr[i].z * computedFieldArr[i].z);

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                              ComputeMethod computeMethod) const
{
    return computeAllBFieldAbs(pointVectors, defaultPrecision, computeMethod);
}


std::vector<vec3::FieldVector3> Coil::computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                     const PrecisionArguments &usedPrecision,
                                                                     ComputeMethod computeMethod) const
{
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;
    std::vector<double> potentialArr;
    std::vector<vec3::FieldVector3> computedFieldArr(pointVectors.size());

    adaptInputVectorsForAllPoints(pointVectors, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);
    calculateAllAPotentialSwitch(cylindricalZArr, cylindricalRArr, potentialArr, usedPrecision, computeMethod);

    for (int i = 0; i < pointVectors.size(); ++i)
        computedFieldArr[i] = vec3::FieldVector3(potentialArr[i] * (-1) * std::sin(cylindricalPhiArr[i]),
                                                 potentialArr[i] * std::cos(cylindricalPhiArr[i]),
                                                 0.0);

    return adaptOutputVectorsForAllPoints(computedFieldArr);
}

std::vector<vec3::FieldVector3> Coil::computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                     ComputeMethod computeMethod) const
{
    return computeAllAPotentialComponents(pointVectors, defaultPrecision, computeMethod);
}

std::vector<double> Coil::computeAllAPotentialX(const std::vector<vec3::CoordVector3> &pointVectors,
                                                const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllAPotentialComponents(pointVectors, usedPrecision, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].x;

    return outputArr;
}

std::vector<double> Coil::computeAllAPotentialX(const std::vector<vec3::CoordVector3> &pointVectors,
                                                ComputeMethod computeMethod) const
{
    return computeAllAPotentialX(pointVectors, defaultPrecision, computeMethod);
}

std::vector<double> Coil::computeAllAPotentialY(const std::vector<vec3::CoordVector3> &pointVectors,
                                                const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllAPotentialComponents(pointVectors, usedPrecision, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].y;

    return outputArr;
}

std::vector<double> Coil::computeAllAPotentialY(const std::vector<vec3::CoordVector3> &pointVectors,
                                                ComputeMethod computeMethod) const
{
    return computeAllAPotentialY(pointVectors, defaultPrecision, computeMethod);
}

std::vector<double> Coil::computeAllAPotentialZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                                const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllAPotentialComponents(pointVectors, usedPrecision, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].z;

    return outputArr;
}

std::vector<double> Coil::computeAllAPotentialZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                                ComputeMethod computeMethod) const
{
    return computeAllAPotentialZ(pointVectors, defaultPrecision, computeMethod);
}

std::vector<double> Coil::computeAllAPotentialAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                                  const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllAPotentialComponents(pointVectors, usedPrecision, computeMethod);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = std::sqrt(computedFieldArr[i].x * computedFieldArr[i].x +
                                 computedFieldArr[i].y * computedFieldArr[i].y +
                                 computedFieldArr[i].z * computedFieldArr[i].z);

    return outputArr;
}

std::vector<double> Coil::computeAllAPotentialAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                                  ComputeMethod computeMethod) const
{
    return computeAllAPotentialAbs(pointVectors, defaultPrecision, computeMethod);
}


std::vector<double> Coil::computeAllEFieldX(const std::vector<vec3::CoordVector3> &pointVectors,
                                            const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<double> output = computeAllAPotentialX(pointVectors, usedPrecision, computeMethod);

    for (double & i : output)
        i *= (2 * M_PI * sineFrequency);

    return output;
}

std::vector<double> Coil::computeAllEFieldX(const std::vector<vec3::CoordVector3> &pointVectors,
                                            ComputeMethod computeMethod) const
{
    return computeAllEFieldX(pointVectors, defaultPrecision, computeMethod);
}

std::vector<double> Coil::computeAllEFieldY(const std::vector<vec3::CoordVector3> &pointVectors,
                                            const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<double> output = computeAllAPotentialY(pointVectors, usedPrecision, computeMethod);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (double & i : output)
        i *= frequencyFactor;

    return output;
}

std::vector<double> Coil::computeAllEFieldY(const std::vector<vec3::CoordVector3> &pointVectors,
                                            ComputeMethod computeMethod) const
{
    return computeAllEFieldY(pointVectors, defaultPrecision, computeMethod);
}

std::vector<double> Coil::computeAllEFieldZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                            const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<double> output = computeAllAPotentialZ(pointVectors, usedPrecision, computeMethod);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (double & i : output)
        i *= frequencyFactor;

    return output;
}

std::vector<double> Coil::computeAllEFieldZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                            ComputeMethod computeMethod) const
{
    return computeAllEFieldZ(pointVectors, defaultPrecision, computeMethod);
}

std::vector<double> Coil::computeAllEFieldAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                              const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<double> output = computeAllAPotentialAbs(pointVectors, usedPrecision, computeMethod);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (double & i : output)
        i *= frequencyFactor;

    return output;
}

std::vector<double> Coil::computeAllEFieldAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                              ComputeMethod computeMethod) const
{
    return computeAllEFieldAbs(pointVectors, defaultPrecision, computeMethod);
}

std::vector<vec3::FieldVector3> Coil::computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                 const PrecisionArguments &usedPrecision,
                                                                 ComputeMethod computeMethod) const
{
    std::vector<vec3::FieldVector3> output = computeAllAPotentialComponents(pointVectors, usedPrecision, computeMethod);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (auto & i : output)
        i *= frequencyFactor;

    return output;
}

std::vector<vec3::FieldVector3> Coil::computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                 ComputeMethod computeMethod) const
{
    return computeAllEFieldComponents(pointVectors, defaultPrecision, computeMethod);
}


std::vector<vec3::Matrix3> Coil::computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                            const PrecisionArguments &usedPrecision, ComputeMethod computeMethod) const
{
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;

    adaptInputVectorsForAllPoints(pointVectors, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);

    std::vector<double> gradientRPhi;
    std::vector<double> gradientRR;
    std::vector<double> gradientRZ;
    std::vector<double> gradientZZ;

    calculateAllBGradientSwitch(cylindricalZArr, cylindricalRArr,
                                gradientRPhi, gradientRR, gradientRZ, gradientZZ, usedPrecision, computeMethod);

    std::vector<vec3::Matrix3> computedGradientMatrix(gradientRPhi.size());

    for (int i = 0; i < gradientRPhi.size(); ++i)
    {
        if (cylindricalRArr[i] / innerRadius < g_zAxisApproximationRatio)
        {
            computedGradientMatrix[i] = vec3::Matrix3(gradientRR[i], 0.0, 0.0,
                                                      0.0, gradientRR[i], 0.0,
                                                      0.0, 0.0, gradientZZ[i]);
        }
        else
        {
            double cosPhi = cos(cylindricalPhiArr[i]);
            double sinPhi = sin(cylindricalPhiArr[i]);

            double xx = gradientRZ[i] * cosPhi * cosPhi + gradientRPhi[i] * sinPhi * sinPhi;
            double yy = gradientRZ[i] * sinPhi * sinPhi + gradientRPhi[i] * cosPhi * cosPhi;
            double zz = gradientZZ[i];

            double xy = 0.5 * sin(2 * cylindricalPhiArr[i]) * (gradientRR[i] - gradientRPhi[i]);
            double xz = gradientRZ[i] * cosPhi;
            double yz = gradientRZ[i] * sinPhi;
            // the matrix is symmetric
            vec3::Matrix3 baseMatrix = vec3::Matrix3(xx, xy, xz,
                                                     xy, yy, yz,
                                                     xz, yz, zz);
            computedGradientMatrix[i] = transformationMatrix * baseMatrix;
        }
    }
    return computedGradientMatrix;
}

std::vector<vec3::Matrix3> Coil::computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                            ComputeMethod computeMethod) const
{
    return computeAllBGradientTensors(pointVectors, defaultPrecision, computeMethod);
}