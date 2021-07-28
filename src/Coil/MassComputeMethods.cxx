#include "Coil.h"

#include <cmath>

void Coil::computeAllBFieldX(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             const PrecisionArguments &usedPrecision,
                             ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedFieldArr.resize(0);
        std::vector<double> fieldH, fieldZ;

        calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision, method);

        for (int i = 0; i < fieldH.size(); ++i)
            computedFieldArr.push_back(fieldH[i] * cos(cylindricalPhiArr[i]));
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllBFieldX(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             ComputeMethod method) const
{
    computeAllBFieldX(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedFieldArr, precisionSettings, method);
}

void Coil::computeAllBFieldY(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             const PrecisionArguments &usedPrecision,
                             ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedFieldArr.resize(0);

        std::vector<double> fieldH, fieldZ;

        calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision, method);

        for (int i = 0; i < fieldH.size(); ++i)
            computedFieldArr.push_back(fieldH[i] * sin(cylindricalPhiArr[i]));
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllBFieldY(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             ComputeMethod method) const
{
    computeAllBFieldY(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedFieldArr, precisionSettings, method);
}

void Coil::computeAllBFieldH(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             std::vector<double> &computedFieldArr,
                             const PrecisionArguments &usedPrecision,
                             ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size())
    {
        std::vector<double> fieldZ;

        calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, computedFieldArr, fieldZ, usedPrecision, method);
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllBFieldH(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             std::vector<double> &computedFieldArr,
                             ComputeMethod method) const
{
    computeAllBFieldH(
            cylindricalZArr, cylindricalRArr, computedFieldArr, precisionSettings, method);
}

void Coil::computeAllBFieldZ(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             const PrecisionArguments &usedPrecision,
                             ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        std::vector<double> fieldH;

        calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr,  fieldH,computedFieldArr, usedPrecision, method);
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllBFieldZ(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             ComputeMethod method) const
{
    computeAllBFieldZ(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedFieldArr, precisionSettings, method);
}

void Coil::computeAllBFieldAbs(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               const PrecisionArguments &usedPrecision,
                               ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedFieldArr.resize(0);

        std::vector<double> fieldH, fieldZ;

        calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision, method);

        for (int i = 0; i < fieldH.size(); ++i)
            computedFieldArr.push_back(std::sqrt(fieldH[i] * fieldH[i] + fieldZ[i] * fieldZ[i]));
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllBFieldAbs(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               ComputeMethod method) const
{
    computeAllBFieldAbs(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                        computedFieldArr, precisionSettings, method);
}

void Coil::computeAllBFieldComponents(const std::vector<double> &cylindricalZArr,
                                      const std::vector<double> &cylindricalRArr,
                                      const std::vector<double> &cylindricalPhiArr,
                                      std::vector<double> &computedFieldXArr,
                                      std::vector<double> &computedFieldYArr,
                                      std::vector<double> &computedFieldZArr,
                                      const PrecisionArguments &usedPrecision,
                                      ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedFieldXArr.resize(0);
        computedFieldYArr.resize(0);
        computedFieldZArr.resize(0);

        std::vector<double> fieldH, fieldZ;

        calculateAllBFieldSwitch(cylindricalZArr, cylindricalRArr, fieldH, fieldZ, usedPrecision, method);

        for (int i = 0; i < fieldH.size(); ++i)
        {
            computedFieldXArr.push_back(fieldH[i] * cos(cylindricalPhiArr[i]));
            computedFieldYArr.push_back(fieldH[i] * sin(cylindricalPhiArr[i]));
            computedFieldZArr.push_back(fieldZ[i]);
        }
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllBFieldComponents(const std::vector<double> &cylindricalZArr,
                                      const std::vector<double> &cylindricalRArr,
                                      const std::vector<double> &cylindricalPhiArr,
                                      std::vector<double> &computedFieldXArr,
                                      std::vector<double> &computedFieldYArr,
                                      std::vector<double> &computedFieldZArr,
                                      ComputeMethod method) const
{
    computeAllBFieldComponents(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
            computedFieldXArr, computedFieldYArr, computedFieldZArr,
            precisionSettings, method);
}

void Coil::computeAllAPotentialX(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedPotentialArr,
                                 const PrecisionArguments &usedPrecision,
                                 ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedPotentialArr.resize(0);

        std::vector<double> potentialA;

        calculateAllAPotentialSwitch(cylindricalZArr, cylindricalRArr, potentialA, usedPrecision, method);

        for (int i = 0; i < potentialA.size(); ++i)
            computedPotentialArr.push_back(potentialA[i] * (-1) * sin(cylindricalPhiArr[i]));

    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllAPotentialX(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedPotentialArr,
                                 ComputeMethod method) const
{
    computeAllAPotentialX(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedPotentialArr, precisionSettings, method);
}

void Coil::computeAllAPotentialY(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedPotentialArr,
                                 const PrecisionArguments &usedPrecision,
                                 ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedPotentialArr.resize(0);

        std::vector<double> potentialA;

        calculateAllAPotentialSwitch(cylindricalZArr, cylindricalRArr, potentialA, usedPrecision, method);

        for (int i = 0; i < potentialA.size(); ++i)
            computedPotentialArr.push_back(potentialA[i] * cos(cylindricalPhiArr[i]));
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllAPotentialY(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedPotentialArr,
                                 ComputeMethod method) const
{
    computeAllAPotentialY(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedPotentialArr, precisionSettings, method);
}

void Coil::computeAllAPotentialZ(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedPotentialArr,
                                 const PrecisionArguments &usedPrecision,
                                 ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedPotentialArr.resize(0);

        for (int i = 0; i < cylindricalZArr.size(); ++i)
            computedPotentialArr.push_back(0.0);
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllAPotentialZ(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedPotentialArr,
                                 ComputeMethod method) const
{
    computeAllAPotentialZ(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedPotentialArr, precisionSettings, method);
}

void
Coil::computeAllAPotentialAbs(const std::vector<double> &cylindricalZArr,
                              const std::vector<double> &cylindricalRArr,
                              std::vector<double> &computedPotentialArr,
                              const PrecisionArguments &usedPrecision,
                              ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size())
    {
        computedPotentialArr.resize(0);

        calculateAllAPotentialSwitch(cylindricalZArr, cylindricalRArr, computedPotentialArr, usedPrecision, method);
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void
Coil::computeAllAPotentialAbs(const std::vector<double> &cylindricalZArr,
                              const std::vector<double> &cylindricalRArr,
                              std::vector<double> &computedPotentialArr,
                              ComputeMethod method) const
{
    computeAllAPotentialAbs(
            cylindricalZArr, cylindricalRArr, computedPotentialArr, precisionSettings, method);
}

void Coil::computeAllAPotentialComponents(const std::vector<double> &cylindricalZArr,
                                          const std::vector<double> &cylindricalRArr,
                                          const std::vector<double> &cylindricalPhiArr,
                                          std::vector<double> &computedPotentialXArr,
                                          std::vector<double> &computedPotentialYArr,
                                          std::vector<double> &computedPotentialZArr,
                                          const PrecisionArguments &usedPrecision,
                                          ComputeMethod method) const
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        computedPotentialXArr.resize(0);
        computedPotentialYArr.resize(0);
        computedPotentialZArr.resize(0);

        std::vector<double> potentialA;

        calculateAllAPotentialSwitch(cylindricalZArr, cylindricalRArr, potentialA, usedPrecision, method);

        for (int i = 0; i < potentialA.size(); ++i)
        {
            computedPotentialXArr.push_back(potentialA[i] * (-1) * sin(cylindricalPhiArr[i]));
            computedPotentialYArr.push_back(potentialA[i] * cos(cylindricalPhiArr[i]));
            computedPotentialZArr.push_back(0.0);
        }
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}

void Coil::computeAllAPotentialComponents(const std::vector<double> &cylindricalZArr,
                                          const std::vector<double> &cylindricalRArr,
                                          const std::vector<double> &cylindricalPhiArr,
                                          std::vector<double> &computedPotentialXArr,
                                          std::vector<double> &computedPotentialYArr,
                                          std::vector<double> &computedPotentialZArr,
                                          ComputeMethod method) const
{
    computeAllAPotentialComponents(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
            computedPotentialXArr, computedPotentialYArr, computedPotentialZArr, precisionSettings, method);
}

void Coil::computeAllEFieldX(const std::vector<double> &cylindricalZArr, const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr, std::vector<double> &computedFieldArr,
                             const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    computeAllAPotentialX(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                          computedFieldArr, usedPrecision, method);

    for (double & i : computedFieldArr)
        i *= (2 * M_PI * sineFrequency);
}

void Coil::computeAllEFieldX(const std::vector<double> &cylindricalZArr, const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr, std::vector<double> &computedFieldArr,
                             ComputeMethod method) const
{
    computeAllEFieldX(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                      computedFieldArr, precisionSettings, method);
}

void Coil::computeAllEFieldY(const std::vector<double> &cylindricalZArr, const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr, std::vector<double> &computedFieldArr,
                             const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    computeAllAPotentialY(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                          computedFieldArr, usedPrecision, method);

    for (double & i : computedFieldArr)
        i *= (2 * M_PI * sineFrequency);
}

void Coil::computeAllEFieldY(const std::vector<double> &cylindricalZArr, const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr, std::vector<double> &computedFieldArr,
                             ComputeMethod method) const
{
    computeAllEFieldY(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                      computedFieldArr, precisionSettings, method);
}

void Coil::computeAllEFieldZ(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             const PrecisionArguments &usedPrecision,
                             ComputeMethod method) const
{
    computeAllAPotentialZ(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                          computedFieldArr, usedPrecision, method);

    for (double & i : computedFieldArr)
        i *= (2 * M_PI * sineFrequency);
}

void Coil::computeAllEFieldZ(const std::vector<double> &cylindricalZArr,
                             const std::vector<double> &cylindricalRArr,
                             const std::vector<double> &cylindricalPhiArr,
                             std::vector<double> &computedFieldArr,
                             ComputeMethod method) const
{
    computeAllEFieldZ(
            cylindricalZArr, cylindricalRArr, cylindricalPhiArr, computedFieldArr, precisionSettings, method);
}

void Coil::computeAllEFieldAbs(const std::vector<double> &cylindricalZArr, const std::vector<double> &cylindricalRArr,
                               std::vector<double> &computedFieldArr, const PrecisionArguments &usedPrecision,
                               ComputeMethod method) const
{
    computeAllAPotentialAbs(cylindricalZArr, cylindricalRArr, computedFieldArr, usedPrecision, method);

    for (double & i : computedFieldArr)
        i *= (2 * M_PI * sineFrequency);
}

void Coil::computeAllEFieldComponents(const std::vector<double> &cylindricalZArr,
                                      const std::vector<double> &cylindricalRArr,
                                      const std::vector<double> &cylindricalPhiArr,
                                      std::vector<double> &computedFieldXArr,
                                      std::vector<double> &computedFieldYArr,
                                      std::vector<double> &computedFieldZArr,
                                      const PrecisionArguments &usedPrecision, ComputeMethod method) const
{
    computeAllAPotentialComponents(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                                   computedFieldXArr, computedFieldYArr, computedFieldZArr,
                                   usedPrecision, method);

    for (int i = 0; i < computedFieldXArr.size(); ++i)
    {
        computedFieldXArr[i] *= (2 * M_PI * sineFrequency);
        computedFieldYArr[i] *= (2 * M_PI * sineFrequency);
        computedFieldZArr[i] *= (2 * M_PI * sineFrequency);
    }
}

void Coil::computeAllEFieldComponents(const std::vector<double> &cylindricalZArr,
                                      const std::vector<double> &cylindricalRArr,
                                      const std::vector<double> &cylindricalPhiArr,
                                      std::vector<double> &computedFieldXArr,
                                      std::vector<double> &computedFieldYArr,
                                      std::vector<double> &computedFieldZArr,
                                      ComputeMethod method) const
{
    computeAllEFieldComponents(cylindricalZArr, cylindricalRArr, cylindricalPhiArr,
                               computedFieldXArr, computedFieldYArr, computedFieldZArr,
                               precisionSettings, method);
}

void Coil::computeAllEFieldAbs(const std::vector<double> &cylindricalZArr, const std::vector<double> &cylindricalRArr,
                               std::vector<double> &computedFieldArr, ComputeMethod method) const
{
    computeAllEFieldAbs(cylindricalZArr, cylindricalRArr, computedFieldArr, precisionSettings, method);
}

void Coil::computeAllBGradientTensors(const std::vector<double> &cylindricalZArr,
                                      const std::vector<double> &cylindricalRArr,
                                      const std::vector<double> &cylindricalPhiArr,
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
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        std::vector<double> gradientRPhi;
        std::vector<double> gradientRR;
        std::vector<double> gradientRZ;
        std::vector<double> gradientZZ;

        calculateAllBGradientSwitch(cylindricalZArr, cylindricalRArr,
                                    gradientRPhi, gradientRR, gradientRZ, gradientZZ, usedPrecision, method);

        for (int i = 0; i < gradientRPhi.size(); ++i)
        {
            if (cylindricalRArr[i] / innerRadius < 1e-14)
            {
                computedGradientXX.push_back(gradientRPhi[i]);
                computedGradientXY.push_back(0.0);
                computedGradientXZ.push_back(0.0);

                computedGradientYX.push_back(0.0);
                computedGradientYY.push_back(gradientRPhi[i]);
                computedGradientYZ.push_back(0.0);

                computedGradientZX.push_back(0.0);
                computedGradientZY.push_back(0.0);
                computedGradientZZ.push_back(gradientZZ[i]);
            }
            else
            {
                double cosPhi = cos(cylindricalPhiArr[i]);
                double sinPhi = sin(cylindricalPhiArr[i]);

                computedGradientXX.push_back(gradientRZ[i] * cosPhi * cosPhi + gradientRPhi[i] * sinPhi * sinPhi);
                computedGradientXY.push_back(0.5 * sin(2 * cylindricalPhiArr[i]) * (gradientRR[i] - gradientRPhi[i]));
                computedGradientXZ.push_back(gradientRZ[i] * cosPhi);

                computedGradientYX.push_back(0.5 * sin(2 * cylindricalPhiArr[i]) * (gradientRR[i] - gradientRPhi[i]));
                computedGradientYY.push_back(gradientRZ[i] * sinPhi * sinPhi + gradientRPhi[i] * cosPhi * cosPhi);
                computedGradientYZ.push_back(gradientRZ[i] * sinPhi);

                computedGradientZX.push_back(gradientRZ[i] * cosPhi);
                computedGradientZY.push_back(gradientRZ[i] * sinPhi);
                computedGradientZZ.push_back(gradientZZ[i]);
            }
        }
    }
    else
        throw "Number of elements in input data vectors is not the same!";
}