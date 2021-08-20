#include "Coil.h"
#include "PrecisionGlobalVars.h"

#include <cmath>
#include <cstdio>


void Coil::adaptInputVectorsForAllPoints(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                         std::vector<double> &cylindricalZArr,
                                         std::vector<double> &cylindricalRArr,
                                         std::vector<double> &cylindricalPhiArr) const
{
    cylindricalZArr.resize(pointVectorArr.size());
    cylindricalRArr.resize(pointVectorArr.size());
    cylindricalPhiArr.resize(pointVectorArr.size());

    vec3::FieldVector3 positionVec = vec3::CoordVector3::convertToFieldVector(positionVector);

    for (int i = 0; i < pointVectorArr.size(); ++i)
    {
        vec3::FieldVector3 pointVec = vec3::CoordVector3::convertToFieldVector(pointVectorArr[i]);
        vec3::FieldVector3 transformedVec = inverseTransformationMatrix * pointVec;
        vec3::FieldVector3 originVec = transformedVec - positionVec;
        //printf("%.15g %.15g %.15g\n", transformedVec.xComponent, transformedVec.yComponent, transformedVec.zComponent);
        vec3::CoordVector3 finalVec = vec3::CoordVector3::convertToCoordVector(originVec);

        finalVec.convertToCylindrical();

        cylindricalZArr[i] = finalVec.component1;
        cylindricalRArr[i] = finalVec.component2;
        cylindricalPhiArr[i] = finalVec.component3;
    }
}

std::vector<vec3::FieldVector3> Coil::adaptOutputVectorsForAllPoints(const std::vector<vec3::FieldVector3> &computedVectorArr) const
{
    std::vector<vec3::FieldVector3> outputVectorArr(computedVectorArr.size());

    for (int i = 0; i < computedVectorArr.size(); ++i)
        outputVectorArr[i] = transformationMatrix * computedVectorArr[i];

    return outputVectorArr;
}


std::vector<vec3::FieldVector3> Coil::computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                 const PrecisionArguments &usedPrecision,
                                                                 ComputeMethod method) const
{
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;
    std::vector<double> fieldH, fieldZ;
    std::vector<vec3::FieldVector3> computedFieldArr(pointVectorArr.size());

    adaptInputVectorsForAllPoints(pointVectorArr, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);
    calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision, method);

    for (int i = 0; i < pointVectorArr.size(); ++i)
        computedFieldArr[i] = vec3::FieldVector3(fieldH[i] * std::cos(cylindricalPhiArr[i]),
                                                 fieldH[i] * std::sin(cylindricalPhiArr[i]),
                                                 fieldZ[i]);

    return adaptOutputVectorsForAllPoints(computedFieldArr);
}

std::vector<vec3::FieldVector3> Coil::computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                 ComputeMethod method) const
{
    return computeAllBFieldComponents(pointVectorArr, defaultPrecision, method);
}

std::vector<double>Coil::computeAllBFieldX(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                           const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(pointVectorArr, usedPrecision, method);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].xComponent;

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldX(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                            ComputeMethod method) const
{
    return computeAllBFieldX(pointVectorArr, defaultPrecision, method);
}

std::vector<double>Coil::computeAllBFieldY(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                           const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(pointVectorArr, usedPrecision, method);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].yComponent;

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldY(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                           ComputeMethod method) const
{
    return computeAllBFieldY(pointVectorArr, defaultPrecision, method);
}

std::vector<double>Coil::computeAllBFieldH(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                           const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(pointVectorArr, usedPrecision, method);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = std::sqrt(computedFieldArr[i].xComponent * computedFieldArr[i].xComponent +
                                 computedFieldArr[i].yComponent * computedFieldArr[i].yComponent);

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldH(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                            ComputeMethod method) const
{
    return computeAllBFieldH(pointVectorArr, defaultPrecision, method);
}

std::vector<double>Coil::computeAllBFieldZ(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                           const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(pointVectorArr, usedPrecision, method);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].zComponent;

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldZ(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                            ComputeMethod method) const
{
    return computeAllBFieldZ(pointVectorArr, defaultPrecision, method);
}

std::vector<double>Coil::computeAllBFieldAbs(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                             const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllBFieldComponents(pointVectorArr, usedPrecision, method);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = std::sqrt(computedFieldArr[i].xComponent * computedFieldArr[i].xComponent +
                                 computedFieldArr[i].yComponent * computedFieldArr[i].yComponent +
                                 computedFieldArr[i].zComponent * computedFieldArr[i].zComponent);

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldAbs(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                              ComputeMethod method) const
{
    return computeAllBFieldAbs(pointVectorArr, defaultPrecision, method);
}


std::vector<vec3::FieldVector3> Coil::computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                     const PrecisionArguments &usedPrecision,
                                                                     ComputeMethod method) const
{
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;
    std::vector<double> potentialArr;
    std::vector<vec3::FieldVector3> computedFieldArr(pointVectorArr.size());

    adaptInputVectorsForAllPoints(pointVectorArr, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);
    calculateAllAPotentialSwitch(cylindricalZArr, cylindricalRArr, potentialArr, usedPrecision, method);

    for (int i = 0; i < pointVectorArr.size(); ++i)
        computedFieldArr[i] = vec3::FieldVector3(potentialArr[i] * (-1) * std::sin(cylindricalPhiArr[i]),
                                                 potentialArr[i] * std::cos(cylindricalPhiArr[i]),
                                                 0.0);

    return adaptOutputVectorsForAllPoints(computedFieldArr);
}

std::vector<vec3::FieldVector3> Coil::computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                     ComputeMethod method) const
{
    return computeAllAPotentialComponents(pointVectorArr, defaultPrecision, method);
}

std::vector<double> Coil::computeAllAPotentialX(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllAPotentialComponents(pointVectorArr, usedPrecision, method);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].xComponent;

    return outputArr;
}

std::vector<double> Coil::computeAllAPotentialX(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                ComputeMethod method) const
{
    return computeAllAPotentialX(pointVectorArr, defaultPrecision, method);
}

std::vector<double> Coil::computeAllAPotentialY(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllAPotentialComponents(pointVectorArr, usedPrecision, method);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].yComponent;

    return outputArr;
}

std::vector<double> Coil::computeAllAPotentialY(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                ComputeMethod method) const
{
    return computeAllAPotentialY(pointVectorArr, defaultPrecision, method);
}

std::vector<double> Coil::computeAllAPotentialZ(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllAPotentialComponents(pointVectorArr, usedPrecision, method);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = computedFieldArr[i].zComponent;

    return outputArr;
}

std::vector<double> Coil::computeAllAPotentialZ(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                ComputeMethod method) const
{
    return computeAllAPotentialZ(pointVectorArr, defaultPrecision, method);
}

std::vector<double> Coil::computeAllAPotentialAbs(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                  const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> computedFieldArr = computeAllAPotentialComponents(pointVectorArr, usedPrecision, method);
    std::vector<double> outputArr(computedFieldArr.size());

    for (int i = 0; i < computedFieldArr.size(); ++i)
        outputArr[i] = std::sqrt(computedFieldArr[i].xComponent * computedFieldArr[i].xComponent +
                                 computedFieldArr[i].yComponent * computedFieldArr[i].yComponent +
                                 computedFieldArr[i].zComponent * computedFieldArr[i].zComponent);

    return outputArr;
}

std::vector<double> Coil::computeAllAPotentialAbs(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                  ComputeMethod method) const
{
    return computeAllAPotentialAbs(pointVectorArr, defaultPrecision, method);
}


std::vector<double> Coil::computeAllEFieldX(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                            const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<double> output = computeAllAPotentialX(pointVectorArr, usedPrecision, method);

    for (double & i : output)
        i *= (2 * M_PI * sineFrequency);

    return output;
}

std::vector<double> Coil::computeAllEFieldX(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                            ComputeMethod method) const
{
    return computeAllEFieldX(pointVectorArr, defaultPrecision, method);
}

std::vector<double> Coil::computeAllEFieldY(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                            const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<double> output = computeAllAPotentialY(pointVectorArr, usedPrecision, method);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (double & i : output)
        i *= frequencyFactor;

    return output;
}

std::vector<double> Coil::computeAllEFieldY(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                            ComputeMethod method) const
{
    return computeAllEFieldY(pointVectorArr, defaultPrecision, method);
}

std::vector<double> Coil::computeAllEFieldZ(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                            const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<double> output = computeAllAPotentialZ(pointVectorArr, usedPrecision, method);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (double & i : output)
        i *= frequencyFactor;

    return output;
}

std::vector<double> Coil::computeAllEFieldZ(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                            ComputeMethod method) const
{
    return computeAllEFieldZ(pointVectorArr, defaultPrecision, method);
}

std::vector<double> Coil::computeAllEFieldAbs(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                              const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<double> output = computeAllAPotentialAbs(pointVectorArr, usedPrecision, method);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (double & i : output)
        i *= frequencyFactor;

    return output;
}

std::vector<double> Coil::computeAllEFieldAbs(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                              ComputeMethod method) const
{
    return computeAllEFieldAbs(pointVectorArr, defaultPrecision, method);
}

std::vector<vec3::FieldVector3> Coil::computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                 const PrecisionArguments &usedPrecision,
                                                                 ComputeMethod method) const
{
    std::vector<vec3::FieldVector3> output = computeAllAPotentialComponents(pointVectorArr, usedPrecision, method);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (auto & i : output)
        i *= frequencyFactor;

    return output;
}

std::vector<vec3::FieldVector3> Coil::computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                 ComputeMethod method) const
{
    return computeAllEFieldComponents(pointVectorArr, defaultPrecision, method);
}


std::vector<vec3::Matrix3> Coil::computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                            const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;

    adaptInputVectorsForAllPoints(pointVectorArr, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);

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
            computedGradientMatrix[i] = transformationMatrix * baseMatrix;
        }
    }
    return computedGradientMatrix;
}

std::vector<vec3::Matrix3> Coil::computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                            ComputeMethod method) const
{
    return computeAllBGradientTensors(pointVectorArr, defaultPrecision, method);
}