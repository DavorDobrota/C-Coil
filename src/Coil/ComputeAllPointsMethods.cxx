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
                                         std::vector<double> &cylindricalPhiArr,
                                         ComputeMethod computeMethod) const
{
    cylindricalZArr.resize(pointVectors.size());
    cylindricalRArr.resize(pointVectors.size());
    cylindricalPhiArr.resize(pointVectors.size());

    int numOps = pointVectors.size();
    int average = std::floor((double)numOps / (double)threadCount);
    std::vector<int> ends(threadCount + 1);
    int remaining = numOps;
    ends[0] = 0;

    for(int i = 0; i < threadCount; i++)
    {
        int temp = (remaining % (threadCount - i) == 0 ? average : average + 1);
        ends[i+1] = (numOps - remaining) + temp;
        remaining -= temp;
    }

    if (pointVectors.size() > pointMultiplier * threadCount && computeMethod == GPU)
    {

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

        for(size_t i = 0; i < threadCount; i++)
        {
            g_threadPool.push(
                    calcThread,
                    std::ref(*this),
                    std::ref(pointVectors),
                    std::ref(cylindricalZArr), std::ref(cylindricalRArr), std::ref(cylindricalPhiArr),
                    ends[i], ends[i + 1]
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
    if (computeMethod == GPU) {

        std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;
        std::vector<double> fieldH(pointVectors.size());
        std::vector<double> fieldZ(pointVectors.size());

        std::vector<vec3::FieldVector3> computedFieldArr(pointVectors.size());

        adaptInputVectorsForAllPoints(pointVectors, cylindricalZArr, cylindricalRArr, cylindricalPhiArr, GPU);

        calculateAllBFieldGPU(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision);

        for (int i = 0; i < pointVectors.size(); ++i)
            computedFieldArr[i] = vec3::FieldVector3(fieldH[i] * std::cos(cylindricalPhiArr[i]),
                                                     fieldH[i] * std::sin(cylindricalPhiArr[i]),
                                                     fieldZ[i]);

        return adaptOutputVectorsForAllPoints(computedFieldArr);
    }
    else if (computeMethod == CPU_MT) {

        return calculateAllBFieldMT(pointVectors, usedPrecision);
    }
    else
    {
        std::vector<vec3::FieldVector3> computedFieldArr(pointVectors.size());

        for (int i = 0; i < pointVectors.size(); ++i)
            computedFieldArr[i] = computeBFieldVector(pointVectors[i], usedPrecision);

        return computedFieldArr;
    }
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
    if (computeMethod == GPU)
    {
        std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;
        std::vector<double> potentialArr;
        std::vector<vec3::FieldVector3> computedPotentialArr(pointVectors.size());

        adaptInputVectorsForAllPoints(pointVectors, cylindricalZArr, cylindricalRArr, cylindricalPhiArr, GPU);
        calculateAllAPotentialGPU(cylindricalZArr, cylindricalRArr, potentialArr, usedPrecision);

        for (int i = 0; i < pointVectors.size(); ++i)
            computedPotentialArr[i] = vec3::FieldVector3(potentialArr[i] * (-1) * std::sin(cylindricalPhiArr[i]),
                                                 potentialArr[i] * std::cos(cylindricalPhiArr[i]),
                                                 0.0);

        return adaptOutputVectorsForAllPoints(computedPotentialArr);
    }
    else if (computeMethod == CPU_MT)
    {
        return calculateAllAPotentialMT(pointVectors, usedPrecision);
    }
    else
    {
        std::vector<vec3::FieldVector3> computedPotentialArr(pointVectors.size());

        for (int i = 0; i < pointVectors.size(); ++i)
            computedPotentialArr[i] = computeAPotentialVector(pointVectors[i], usedPrecision);

        return computedPotentialArr;
    }
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
                                                            const PrecisionArguments &usedPrecision,
                                                            ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
    {
        std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;

        adaptInputVectorsForAllPoints(pointVectors, cylindricalZArr, cylindricalRArr, cylindricalPhiArr, GPU);

        std::vector<double> gradientRPhi;
        std::vector<double> gradientRR;
        std::vector<double> gradientRZ;
        std::vector<double> gradientZZ;

        std::vector<vec3::Matrix3> computedGradientMatrix(cylindricalZArr.size());
        calculateAllBGradientGPU(cylindricalZArr, cylindricalRArr, gradientRPhi, gradientRR, gradientRZ, gradientZZ, usedPrecision);

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
    else if (computeMethod == CPU_MT)
    {
        return calculateAllBGradientMT(pointVectors, usedPrecision);
    }
    else
    {
        std::vector<vec3::Matrix3> computedGradientArr(pointVectors.size());

        for (int i = 0; i < pointVectors.size(); ++i)
            computedGradientArr[i] = computeBGradientTensor(pointVectors[i], usedPrecision);

        return computedGradientArr;
    }
}

std::vector<vec3::Matrix3> Coil::computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                            ComputeMethod computeMethod) const
{
    return computeAllBGradientTensors(pointVectors, defaultPrecision, computeMethod);
}