#include "CoilGroup.h"
#include "CoilData.h"
#include "hardware_acceleration.h"
#include "LegendreMatrix.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace
{
    const double g_MiReduced = 0.0000001;
}


void CoilGroup::generateCoilDataArray(CoilData *coilDataArr) const
{
    for (int i = 0; i < memberCoils.size(); ++i)
    {
        if (memberCoils[i].isUsingFastMethod())
            coilDataArr[i].constFactor = g_MiReduced * memberCoils[i].getCurrentDensity() *
                                         memberCoils[i].getThickness() * M_PI * 0.5;
        else
            coilDataArr[i].constFactor = g_MiReduced * memberCoils[i].getCurrentDensity() *
                                         memberCoils[i].getThickness() * memberCoils[i].getLength() * M_PI * 0.5;

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

        vec3::FieldVector3 tempVec = vec3::CoordVector3::convertToFieldVector(memberCoils[i].getPositionVector());

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

std::vector<vec3::FieldVector3>
CoilGroup::calculateAllAPotentialComponentsGPU(const std::vector<vec3::CoordVector3> &pointVectors) const
{
    long long size = pointVectors.size();
    long long coils = memberCoils.size();

    auto coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *resultArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *coilDataArr = static_cast<CoilData *>(calloc(coils, sizeof(CoilData)));

    if(!coordinateArr || !resultArr || !coilDataArr)
        throw std::bad_alloc();

    for (long long i = 0; i < size; ++i)
    {
        vec3::CoordVector3 vector = pointVectors[i];
        vector.convertToCartesian();

        coordinateArr[i].x = vector.comp1;
        coordinateArr[i].y = vector.comp2;
        coordinateArr[i].z = vector.comp3;
    }

    generateCoilDataArray(coilDataArr);

    #if USE_GPU == 1
        Calculate_hardware_accelerated_a_group(coils, size, coilDataArr, coordinateArr, resultArr);
    #else
        free(coordinateArr);
        free(resultArr);
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    free(coordinateArr);
    std::vector<vec3::FieldVector3> computedPotentialArr;
    computedPotentialArr.reserve(size);

    for (long long i = 0; i < pointVectors.size(); ++i)
        computedPotentialArr.emplace_back(resultArr[i].x, resultArr[i].y, resultArr[i].z);

    free(resultArr);

    return computedPotentialArr;
}

std::vector<vec3::FieldVector3>
CoilGroup::calculateAllBFieldComponentsGPU(const std::vector<vec3::CoordVector3> &pointVectors) const
{
    long long size = pointVectors.size();
    long long coils = memberCoils.size();

    auto coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *resultArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *coilDataArr = static_cast<CoilData *>(calloc(coils, sizeof(CoilData)));

    if(!coordinateArr || !resultArr || !coilDataArr)
        throw std::bad_alloc();

    for (long long i = 0; i < size; ++i)
    {
        vec3::CoordVector3 vector = pointVectors[i];
        vector.convertToCartesian();

        coordinateArr[i].x = vector.comp1;
        coordinateArr[i].y = vector.comp2;
        coordinateArr[i].z = vector.comp3;
    }

    generateCoilDataArray(coilDataArr);

    #if USE_GPU == 1
        Calculate_hardware_accelerated_b_group(coils, size, coilDataArr, coordinateArr, resultArr);
    #else
        free(coordinateArr);
        free(resultArr);
        free(coilDataArr);
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    free(coordinateArr);
    std::vector<vec3::FieldVector3> computedFieldArr;
    computedFieldArr.reserve(size);

    for (long long i = 0; i < pointVectors.size(); ++i)
        computedFieldArr.emplace_back(resultArr[i].x, resultArr[i].y, resultArr[i].z);

    free(resultArr);

    return computedFieldArr;
}

std::vector<vec3::Matrix3>
CoilGroup::calculateAllBGradientTensorsGPU(const std::vector<vec3::CoordVector3> &pointVectors) const
{
    long long size = pointVectors.size();
    long long coils = memberCoils.size();

    auto coordinateArr = static_cast<DataVector *>(calloc(size, sizeof(DataVector)));
    auto *resultArr = static_cast<DataMatrix *>(calloc(size, sizeof(DataMatrix)));
    auto *coilDataArr = static_cast<CoilData *>(calloc(coils, sizeof(CoilData)));

    if(!coordinateArr || !resultArr || !coilDataArr)
        throw std::bad_alloc();

    for (long long i = 0; i < size; ++i)
    {
        vec3::CoordVector3 vector = pointVectors[i];
        vector.convertToCartesian();

        coordinateArr[i].x = vector.comp1;
        coordinateArr[i].y = vector.comp2;
        coordinateArr[i].z = vector.comp3;
    }

    generateCoilDataArray(coilDataArr);

    #if USE_GPU == 1
        Calculate_hardware_accelerated_g_group(coils, size, coilDataArr, coordinateArr, resultArr);
    #else
        free(coordinateArr);
        free(resultArr);
        free(coilDataArr);
        throw std::logic_error("GPU functions are disabled. (rebuild the project with USE_GPU)");
    #endif // USE_GPU

    free(coordinateArr);
    std::vector<vec3::Matrix3> computedGradientArr;
    computedGradientArr.reserve(size);

    for (long long i = 0; i < pointVectors.size(); ++i)
        computedGradientArr.emplace_back(resultArr[i].xx, resultArr[i].xy, resultArr[i].xz,
                                         resultArr[i].yx, resultArr[i].yy, resultArr[i].yz,
                                         resultArr[i].zx, resultArr[i].zy, resultArr[i].zz);

    free(resultArr);

    return computedGradientArr;

}
