#include "Coil.h"
#include "PrecisionGlobalVars.h"

#include <cmath>
#include <cstdio>


void Coil::adaptInputVectorToCalculateMethods(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                              std::vector<double> &cylindricalZArr,
                                              std::vector<double> &cylindricalRArr,
                                              std::vector<double> &cylindricalPhiArr) const
{
    cylindricalZArr.resize(positionVectorArr.size());
    cylindricalRArr.resize(positionVectorArr.size());
    cylindricalPhiArr.resize(positionVectorArr.size());

    vec3::FieldVector3 positionVec = vec3::CoordVector3::convertToFieldVector(positionVector);
    positionVec.multiplyByConstant(-1.0);

    for (int i = 0; i < positionVectorArr.size(); ++i)
    {
        vec3::FieldVector3 pointVec = vec3::CoordVector3::convertToFieldVector(positionVectorArr[i]);
        vec3::FieldVector3 originVec = vec3::FieldVector3::addVectors(pointVec, positionVec);
        vec3::FieldVector3 transformedVec = vec3::Matrix3::matrixVectorMultiplication(inverseTransformationMatrix, originVec);
        vec3::CoordVector3 finalVec = vec3::CoordVector3::convertToCoordVector(transformedVec);

        finalVec.convertToCylindrical();

        cylindricalZArr[i] = finalVec.component1;
        cylindricalRArr[i] = finalVec.component2;
        cylindricalPhiArr[i] = finalVec.component3;
    }
}

std::vector<vec3::FieldVector3> Coil::adaptOutputVectorValues(const std::vector<vec3::FieldVector3> &computedVectorArr) const
{
    std::vector<vec3::FieldVector3> outputVectorArr(computedVectorArr.size());

    for (int i = 0; i < computedVectorArr.size(); ++i)
        outputVectorArr[i] = vec3::Matrix3::matrixVectorMultiplication(transformationMatrix, computedVectorArr[i]);

    return outputVectorArr;
}


std::vector<vec3::FieldVector3> Coil::computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                                 const PrecisionArguments &usedPrecision,
                                                                 ComputeMethod method) const
{
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;
    std::vector<double> fieldH, fieldZ;
    std::vector<vec3::FieldVector3> computedFieldArr(positionVectorArr.size());

    adaptInputVectorToCalculateMethods(positionVectorArr, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);
    calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision, method);

    for (int i = 0; i < positionVectorArr.size(); ++i)
        computedFieldArr[i] = vec3::FieldVector3(fieldH[i] * std::cos(cylindricalPhiArr[i]),
                                                 fieldH[i] * std::sin(cylindricalPhiArr[i]),
                                                 fieldZ[i]);

    return adaptOutputVectorValues(computedFieldArr);
}

std::vector<vec3::FieldVector3> Coil::computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                                 ComputeMethod method) const
                                                                 {
    return computeAllBFieldComponents(positionVectorArr, defaultPrecision, method);
                                                                 }

std::vector<double>Coil::computeAllBFieldX(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                           const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(positionVectorArr, usedPrecision, method);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].xComponent;

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldX(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                            ComputeMethod method) const
{
    return computeAllBFieldX(positionVectorArr, defaultPrecision, method);
}


std::vector<double>Coil::computeAllBFieldY(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                           const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(positionVectorArr, usedPrecision, method);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].yComponent;

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldY(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                           ComputeMethod method) const
{
    return computeAllBFieldY(positionVectorArr, defaultPrecision, method);
}

std::vector<double>Coil::computeAllBFieldH(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                           const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(positionVectorArr, usedPrecision, method);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = std::sqrt(computedFieldArr[i].xComponent * computedFieldArr[i].xComponent +
                                 computedFieldArr[i].yComponent * computedFieldArr[i].yComponent);

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldH(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                            ComputeMethod method) const
{
    return computeAllBFieldH(positionVectorArr, defaultPrecision, method);
}

std::vector<double>Coil::computeAllBFieldZ(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                           const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(positionVectorArr, usedPrecision, method);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].zComponent;

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldZ(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                            ComputeMethod method) const
{
    return computeAllBFieldZ(positionVectorArr, defaultPrecision, method);
}

std::vector<double>Coil::computeAllBFieldAbs(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                             const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(positionVectorArr, usedPrecision, method);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = std::sqrt(computedFieldArr[i].xComponent * computedFieldArr[i].xComponent +
                                 computedFieldArr[i].yComponent * computedFieldArr[i].yComponent +
                                 computedFieldArr[i].zComponent * computedFieldArr[i].zComponent);

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldAbs(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                              ComputeMethod method) const
{
    return computeAllBFieldAbs(positionVectorArr, defaultPrecision, method);
}


std::vector<vec3::FieldVector3> Coil::computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                                     const PrecisionArguments &usedPrecision,
                                                                     ComputeMethod method) const
{
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;
    std::vector<double> potentialArr;
    std::vector<vec3::FieldVector3> computedFieldArr(positionVectorArr.size());

    adaptInputVectorToCalculateMethods(positionVectorArr, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);
    calculateAllAPotentialSwitch(cylindricalZArr, cylindricalRArr, potentialArr, usedPrecision, method);

    for (int i = 0; i < positionVectorArr.size(); ++i)
        computedFieldArr[i] = vec3::FieldVector3(potentialArr[i] * (-1) * std::sin(cylindricalPhiArr[i]),
                                                 potentialArr[i] * std::cos(cylindricalPhiArr[i]),
                                                 0.0);

    return adaptOutputVectorValues(computedFieldArr);
}

std::vector<vec3::FieldVector3> Coil::computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                                     ComputeMethod method) const
{
    return computeAllAPotentialComponents(positionVectorArr, defaultPrecision, method);
}

std::vector<double> Coil::computeAllAPotentialX(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllAPotentialComponents(positionVectorArr, usedPrecision, method);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].xComponent;

    return outputArr;
}

std::vector<double> Coil::computeAllAPotentialX(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                ComputeMethod method) const
{
    return computeAllAPotentialX(positionVectorArr, defaultPrecision, method);
}

std::vector<double> Coil::computeAllAPotentialY(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllAPotentialComponents(positionVectorArr, usedPrecision, method);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].yComponent;

    return outputArr;
}

std::vector<double> Coil::computeAllAPotentialY(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                ComputeMethod method) const
{
    return computeAllAPotentialY(positionVectorArr, defaultPrecision, method);
}

std::vector<double> Coil::computeAllAPotentialZ(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllAPotentialComponents(positionVectorArr, usedPrecision, method);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].zComponent;

    return outputArr;
}

std::vector<double> Coil::computeAllAPotentialZ(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                ComputeMethod method) const
{
    return computeAllAPotentialZ(positionVectorArr, defaultPrecision, method);
}

std::vector<double> Coil::computeAllAPotentialAbs(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                  const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllAPotentialComponents(positionVectorArr, usedPrecision, method);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = std::sqrt(computedFieldArr[i].xComponent * computedFieldArr[i].xComponent +
                                 computedFieldArr[i].yComponent * computedFieldArr[i].yComponent +
                                 computedFieldArr[i].zComponent * computedFieldArr[i].zComponent);

    return outputArr;
}

std::vector<double> Coil::computeAllAPotentialAbs(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                  ComputeMethod method) const
{
    return computeAllAPotentialAbs(positionVectorArr, defaultPrecision, method);
}


std::vector<double> Coil::computeAllEFieldX(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                            const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<double> output = computeAllAPotentialX(positionVectorArr, usedPrecision, method);

    for (double & i : output)
        i *= (2 * M_PI * sineFrequency);

    return output;
}

std::vector<double> Coil::computeAllEFieldX(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                            ComputeMethod method) const
{
    return computeAllEFieldX(positionVectorArr, defaultPrecision, method);
}

std::vector<double> Coil::computeAllEFieldY(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                            const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<double> output = computeAllAPotentialY(positionVectorArr, usedPrecision, method);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (double & i : output)
        i *= frequencyFactor;

    return output;
}

std::vector<double> Coil::computeAllEFieldY(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                            ComputeMethod method) const
{
    return computeAllEFieldY(positionVectorArr, defaultPrecision, method);
}

std::vector<double> Coil::computeAllEFieldZ(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                            const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<double> output = computeAllAPotentialZ(positionVectorArr, usedPrecision, method);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (double & i : output)
        i *= frequencyFactor;

    return output;
}

std::vector<double> Coil::computeAllEFieldZ(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                            ComputeMethod method) const
{
    return computeAllEFieldZ(positionVectorArr, defaultPrecision, method);
}

std::vector<double> Coil::computeAllEFieldAbs(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                              const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<double> output = computeAllAPotentialAbs(positionVectorArr, usedPrecision, method);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (double & i : output)
        i *= frequencyFactor;

    return output;
}

std::vector<double> Coil::computeAllEFieldAbs(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                              ComputeMethod method) const
{
    return computeAllEFieldAbs(positionVectorArr, defaultPrecision, method);
}

std::vector<vec3::FieldVector3> Coil::computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                                 const PrecisionArguments &usedPrecision,
                                                                 ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> output = computeAllAPotentialComponents(positionVectorArr, usedPrecision, method);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (auto & i : output)
        i.multiplyByConstant(frequencyFactor);

    return output;
}

std::vector<vec3::FieldVector3> Coil::computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                                 ComputeMethod method) const
{
    return computeAllEFieldComponents(positionVectorArr, defaultPrecision, method);
}


std::vector<vec3::Matrix3> Coil::computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                            const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;

    adaptInputVectorToCalculateMethods(positionVectorArr, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);

    std::vector<double> gradientRPhi;
    std::vector<double> gradientRR;
    std::vector<double> gradientRZ;
    std::vector<double> gradientZZ;

    calculateAllBGradientSwitch(cylindricalZArr, cylindricalRArr,
                                gradientRPhi, gradientRR, gradientRZ, gradientZZ, usedPrecision, method);

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
            computedGradientMatrix[i] = vec3::Matrix3::matrixMultiplication(transformationMatrix, baseMatrix);
        }
    }
    return computedGradientMatrix;
}

std::vector<vec3::Matrix3> Coil::computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                            ComputeMethod method) const
{
    return computeAllBGradientTensors(positionVectorArr, defaultPrecision, method);
}