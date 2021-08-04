#include "Coil.h"

#include <cmath>


void Coil::adaptInputVectorToCalculateMethods(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                              std::vector<double> &cylindricalZArr,
                                              std::vector<double> &cylindricalRArr,
                                              std::vector<double> &cylindricalPhiArr)
{
    cylindricalZArr.resize(positionVectorArr.size());
    cylindricalRArr.resize(positionVectorArr.size());

    for (int i = 0; i < positionVectorArr.size(); ++i)
    {
        vec3::CoordVector3 vec = positionVectorArr[i];
        vec.convertToCylindrical();

        cylindricalZArr[i] = vec.component1;
        cylindricalRArr[i] = vec.component2;
        cylindricalPhiArr[i] = vec.component3;
    }
}


std::vector<double>Coil::computeAllBFieldX(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                           const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;
    std::vector<double> fieldH, fieldZ;
    std::vector<double> outputArr(positionVectorArr.size());

    adaptInputVectorToCalculateMethods(positionVectorArr, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);
    calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision, method);

    for (int i = 0; i < positionVectorArr.size(); ++i)
        outputArr[i] = fieldH[i] * std::cos(cylindricalPhiArr[i]);

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
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;
    std::vector<double> fieldH, fieldZ;
    std::vector<double> outputArr(positionVectorArr.size());

    adaptInputVectorToCalculateMethods(positionVectorArr, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);
    calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision, method);

    for (int i = 0; i < positionVectorArr.size(); ++i)
        outputArr[i] = fieldH[i] * std::sin(cylindricalPhiArr[i]);

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
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;
    std::vector<double> fieldH, fieldZ;

    adaptInputVectorToCalculateMethods(positionVectorArr, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);
    calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision, method);

    return fieldH;
}

std::vector<double> Coil::computeAllBFieldH(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                            ComputeMethod method) const
{
    return computeAllBFieldH(positionVectorArr, defaultPrecision, method);
}

std::vector<double>Coil::computeAllBFieldZ(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                           const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;
    std::vector<double> fieldH, fieldZ;

    adaptInputVectorToCalculateMethods(positionVectorArr, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);
    calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision, method);

    return fieldZ;
}

std::vector<double> Coil::computeAllBFieldZ(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                            ComputeMethod method) const
{
    return computeAllBFieldZ(positionVectorArr, defaultPrecision, method);
}

std::vector<double>Coil::computeAllBFieldAbs(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                             const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;
    std::vector<double> fieldH, fieldZ;
    std::vector<double> outputArr(positionVectorArr.size());

    adaptInputVectorToCalculateMethods(positionVectorArr, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);
    calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision, method);

    for (int i = 0; i < positionVectorArr.size(); ++i)
        outputArr[i] = std::sqrt(fieldH[i] * fieldH[i] + fieldZ[i] * fieldZ[i]);

    return outputArr;
}

std::vector<double> Coil::computeAllBFieldAbs(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                              ComputeMethod method) const
{
    return computeAllBFieldAbs(positionVectorArr, defaultPrecision, method);
}

std::vector<vec3::FieldVector3> Coil::computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                                 const PrecisionArguments &usedPrecision,
                                                                 ComputeMethod method) const
{
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;
    std::vector<double> fieldH, fieldZ;
    std::vector<vec3::FieldVector3> outputArr(positionVectorArr.size());

    adaptInputVectorToCalculateMethods(positionVectorArr, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);
    calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision, method);

    for (int i = 0; i < positionVectorArr.size(); ++i)
        outputArr[i] = vec3::FieldVector3(fieldH[i] * std::cos(cylindricalPhiArr[i]),
                                          fieldH[i] * std::sin(cylindricalPhiArr[i]),
                                          fieldZ[i]);

    return outputArr;
}

std::vector<vec3::FieldVector3> Coil::computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                                 ComputeMethod method) const
{
    return computeAllBFieldComponents(positionVectorArr, defaultPrecision, method);
}


std::vector<double> Coil::computeAllAPotentialX(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;
    std::vector<double> potentialArr;
    std::vector<double> outputArr(positionVectorArr.size());

    adaptInputVectorToCalculateMethods(positionVectorArr, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);
    calculateAllAPotentialSwitch(cylindricalZArr, cylindricalRArr, potentialArr, usedPrecision, method);

    for (int i = 0; i < positionVectorArr.size(); ++i)
        outputArr[i] = potentialArr[i] * (-1) * std::sin(cylindricalPhiArr[i]);

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
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;
    std::vector<double> potentialArr;
    std::vector<double> outputArr(positionVectorArr.size());

    adaptInputVectorToCalculateMethods(positionVectorArr, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);
    calculateAllAPotentialSwitch(cylindricalZArr, cylindricalRArr, potentialArr, usedPrecision, method);

    for (int i = 0; i < positionVectorArr.size(); ++i)
        outputArr[i] = potentialArr[i] * std::cos(cylindricalPhiArr[i]);

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
    std::vector<double> outputArr(positionVectorArr.size());
    std::fill(outputArr.begin(), outputArr.end(), 0.0);
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
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;
    std::vector<double> potentialArr;

    adaptInputVectorToCalculateMethods(positionVectorArr, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);
    calculateAllAPotentialSwitch(cylindricalZArr, cylindricalRArr, potentialArr, usedPrecision, method);

    return potentialArr;
}

std::vector<double> Coil::computeAllAPotentialAbs(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                  ComputeMethod method) const
{
    return computeAllAPotentialAbs(positionVectorArr, defaultPrecision, method);
}

std::vector<vec3::FieldVector3> Coil::computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                                     const PrecisionArguments &usedPrecision,
                                                                     ComputeMethod method) const
{
    std::vector<double> cylindricalZArr, cylindricalRArr, cylindricalPhiArr;
    std::vector<double> potentialArr;
    std::vector<vec3::FieldVector3> outputArr(positionVectorArr.size());

    adaptInputVectorToCalculateMethods(positionVectorArr, cylindricalZArr, cylindricalRArr, cylindricalPhiArr);
    calculateAllAPotentialSwitch(cylindricalZArr, cylindricalRArr, potentialArr, usedPrecision, method);

    for (int i = 0; i < positionVectorArr.size(); ++i)
        outputArr[i] = vec3::FieldVector3(potentialArr[i] * (-1) * std::sin(cylindricalPhiArr[i]),
                                          potentialArr[i] * std::cos(cylindricalPhiArr[i]),
                                          0.0);

    return outputArr;
}

std::vector<vec3::FieldVector3> Coil::computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                                     ComputeMethod method) const
{
    return computeAllAPotentialComponents(positionVectorArr, defaultPrecision, method);
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
    {
        i.xComponent *= frequencyFactor;
        i.yComponent *= frequencyFactor;
        i.zComponent *= frequencyFactor;
    }
}

std::vector<vec3::FieldVector3> Coil::computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                                                 ComputeMethod method) const
{
    return computeAllEFieldComponents(positionVectorArr, defaultPrecision, method);
}


void Coil::computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                      std::vector<double> &computedGradientXX,
                                      std::vector<double> &computedGradientXY,
                                      std::vector<double> &computedGradientXZ,
                                      std::vector<double> &computedGradientYX,
                                      std::vector<double> &computedGradientYY,
                                      std::vector<double> &computedGradientYZ,
                                      std::vector<double> &computedGradientZX,
                                      std::vector<double> &computedGradientZY,
                                      std::vector<double> &computedGradientZZ,
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

    computedGradientXX.resize(gradientRPhi.size());
    computedGradientXY.resize(gradientRPhi.size());
    computedGradientXZ.resize(gradientRPhi.size());

    computedGradientYX.resize(gradientRPhi.size());
    computedGradientYY.resize(gradientRPhi.size());
    computedGradientYZ.resize(gradientRPhi.size());

    computedGradientZX.resize(gradientRPhi.size());
    computedGradientZY.resize(gradientRPhi.size());
    computedGradientZZ.resize(gradientRPhi.size());

    for (int i = 0; i < gradientRPhi.size(); ++i)
    {
        if (cylindricalRArr[i] / innerRadius < 1e-14)
        {
            computedGradientXX[i] = gradientRR[i];
            computedGradientXY[i] = 0.0;
            computedGradientXZ[i] = 0.0;

            computedGradientYX[i] = 0.0;
            computedGradientYY[i] = gradientRR[i];
            computedGradientYZ[i] = 0.0;

            computedGradientZX[i] = 0.0;
            computedGradientZY[i] = 0.0;
            computedGradientZZ[i] = gradientZZ[i];
        }
        else
        {
            double cosPhi = cos(cylindricalPhiArr[i]);
            double sinPhi = sin(cylindricalPhiArr[i]);

            computedGradientXX[i] = gradientRZ[i] * cosPhi * cosPhi + gradientRPhi[i] * sinPhi * sinPhi;
            computedGradientXY[i] = 0.5 * sin(2 * cylindricalPhiArr[i]) * (gradientRR[i] - gradientRPhi[i]);
            computedGradientXZ[i] = gradientRZ[i] * cosPhi;

            computedGradientYX[i] = 0.5 * sin(2 * cylindricalPhiArr[i]) * (gradientRR[i] - gradientRPhi[i]);
            computedGradientYY[i] = gradientRZ[i] * sinPhi * sinPhi + gradientRPhi[i] * cosPhi * cosPhi;
            computedGradientYZ[i] = gradientRZ[i] * sinPhi;

            computedGradientZX[i] = gradientRZ[i] * cosPhi;
            computedGradientZY[i] = gradientRZ[i] * sinPhi;
            computedGradientZZ[i] = gradientZZ[i];
        }
    }
}

void Coil::computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &positionVectorArr,
                                      std::vector<double> &computedGradientXX,
                                      std::vector<double> &computedGradientXY,
                                      std::vector<double> &computedGradientXZ,
                                      std::vector<double> &computedGradientYX,
                                      std::vector<double> &computedGradientYY,
                                      std::vector<double> &computedGradientYZ,
                                      std::vector<double> &computedGradientZX,
                                      std::vector<double> &computedGradientZY,
                                      std::vector<double> &computedGradientZZ,
                                      ComputeMethod method) const
{
    computeAllBGradientTensors(positionVectorArr,
                               computedGradientXX, computedGradientXY, computedGradientXZ,
                               computedGradientYX, computedGradientYY, computedGradientYZ,
                               computedGradientZX, computedGradientZY, computedGradientZZ,
                               defaultPrecision, method);
}