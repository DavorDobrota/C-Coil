#include "CoilGroup.h"
#include "CoilGroupAcceleration.h"
#include "LegendreMatrix.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace
{
    const double g_MiReduced = 0.0000001;
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"


void CoilGroup::generateCoilDataArray(CoilData *coilDataArr, bool removeSpecificCoil,
                                      unsigned long long specificCoilId) const
{
    size_t maxNum = memberCoils.size();

    for (int i = 0; i < maxNum; ++i)
    {
        if (removeSpecificCoil && specificCoilId == memberCoils[i].getId())
        {
            --i; --maxNum;
            continue;
        }

        if (memberCoils[i].isUsingFastMethod())
            coilDataArr[i].constFactor = g_MiReduced * memberCoils[i].getCurrentDensity() *
                                         memberCoils[i].getThickness() * M_PI * 0.5;
        else
            coilDataArr[i].constFactor = g_MiReduced * memberCoils[i].getCurrentDensity() *
                                         memberCoils[i].getThickness() * memberCoils[i].getLength() * M_PI * 0.5;

        coilDataArr[i].current = memberCoils[i].getCurrent();
        coilDataArr[i].frequency = memberCoils[i].getSineFrequency();
        coilDataArr[i].useFastMethod = memberCoils[i].isUsingFastMethod();

        coilDataArr[i].innerRadius = memberCoils[i].getInnerRadius();
        coilDataArr[i].thickness = memberCoils[i].getThickness();
        coilDataArr[i].length = memberCoils[i].getLength();

        PrecisionArguments coilPrecision = memberCoils[i].getPrecisionSettingsGPU();

        coilDataArr[i].lengthIncrements = coilPrecision.lengthIncrementCount;
        coilDataArr[i].thicknessIncrements = coilPrecision.thicknessIncrementCount;
        coilDataArr[i].angularIncrements = coilPrecision.angularIncrementCount;

        for (int j = 0; j < coilDataArr[i].angularIncrements; ++j)
        {
            double phiPosition = M_PI_2 * (1.0 + Legendre::positionMatrix[coilDataArr[i].angularIncrements - 1][j]);

            coilDataArr[i].cosPrecomputeArray[j] = std::cos(phiPosition);
            coilDataArr[i].angularWeightArray[j] = Legendre::weightsMatrix[coilDataArr[i].angularIncrements - 1][j];
        }

        for (int j = 0; j < coilDataArr[i].thicknessIncrements; ++j)
        {
            coilDataArr[i].thicknessPositionArray[j] = Legendre::positionMatrix[coilDataArr[i].thicknessIncrements - 1][j];
            coilDataArr[i].thicknessWeightArray[j] = Legendre::weightsMatrix[coilDataArr[i].thicknessIncrements - 1][j];
        }

        vec3::Vector3 tempVec = memberCoils[i].getPositionVector();

        coilDataArr[i].positionVector[0] = tempVec.x;
        coilDataArr[i].positionVector[1] = tempVec.y;
        coilDataArr[i].positionVector[2] = tempVec.z;

        vec3::Matrix3 transformMatrix = memberCoils[i].getTransformationMatrix();

        coilDataArr[i].transformArray[0] = transformMatrix.xx;
        coilDataArr[i].transformArray[1] = transformMatrix.xy;
        coilDataArr[i].transformArray[2] = transformMatrix.xz;
        coilDataArr[i].transformArray[3] = transformMatrix.yx;
        coilDataArr[i].transformArray[4] = transformMatrix.yy;
        coilDataArr[i].transformArray[5] = transformMatrix.yz;
        coilDataArr[i].transformArray[6] = transformMatrix.zx;
        coilDataArr[i].transformArray[7] = transformMatrix.zy;
        coilDataArr[i].transformArray[8] = transformMatrix.zz;

        vec3::Matrix3 invTransformMatrix = memberCoils[i].getInverseTransformationMatrix();

        coilDataArr[i].invTransformArray[0] = invTransformMatrix.xx;
        coilDataArr[i].invTransformArray[1] = invTransformMatrix.xy;
        coilDataArr[i].invTransformArray[2] = invTransformMatrix.xz;
        coilDataArr[i].invTransformArray[3] = invTransformMatrix.yx;
        coilDataArr[i].invTransformArray[4] = invTransformMatrix.yy;
        coilDataArr[i].invTransformArray[5] = invTransformMatrix.yz;
        coilDataArr[i].invTransformArray[6] = invTransformMatrix.zx;
        coilDataArr[i].invTransformArray[7] = invTransformMatrix.zy;
        coilDataArr[i].invTransformArray[8] = invTransformMatrix.zz;
    }
}

void CoilGroup::generateSecondaryData(const Coil &secondary, SecondaryCoilData &secondaryData, bool forceCalculation)
{
    secondaryData.innerRadius = secondary.getInnerRadius();
    secondaryData.thickness = secondary.getThickness();
    secondaryData.length = secondary.getLength();

    if (!forceCalculation)
        secondaryData.correctionFactor = 2*M_PI * secondary.getNumOfTurns();
    else
        secondaryData.correctionFactor = 2*M_PI * secondary.getNumOfTurns() * secondary.getCurrent();

    PrecisionArguments precisionArguments = secondary.getPrecisionSettingsGPU();

    secondaryData.lengthIncrements = precisionArguments.lengthIncrementCount;
    secondaryData.thicknessIncrements = precisionArguments.thicknessIncrementCount;
    secondaryData.angularIncrements = precisionArguments.angularIncrementCount;

    for (int i = 0; i < precisionArguments.lengthIncrementCount; ++i)
    {
        secondaryData.lengthPositionArray[i] = Legendre::positionMatrix[precisionArguments.lengthIncrementCount - 1][i];
        secondaryData.lengthWeightArray[i] = Legendre::weightsMatrix[precisionArguments.lengthIncrementCount - 1][i];
    }

    for (int i = 0; i < precisionArguments.thicknessIncrementCount; ++i)
    {
        secondaryData.thicknessPositionArray[i] = Legendre::positionMatrix[precisionArguments.thicknessIncrementCount - 1][i];
        secondaryData.thicknessWeightArray[i] = Legendre::weightsMatrix[precisionArguments.thicknessIncrementCount - 1][i];
    }

    for (int i = 0; i < precisionArguments.angularIncrementCount; ++i)
    {
        secondaryData.angularPositionArray[i] = Legendre::positionMatrix[precisionArguments.angularIncrementCount - 1][i];
        secondaryData.angularWeightArray[i] = Legendre::weightsMatrix[precisionArguments.angularIncrementCount - 1][i];
    }
}


#pragma clang diagnostic pop
