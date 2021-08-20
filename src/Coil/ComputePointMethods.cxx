#include "Coil.h"
#include "PrecisionGlobalVars.h"

#include <cmath>


vec3::CoordVector3 Coil::adaptInputVectorForPoint(const vec3::CoordVector3 &pointVector) const
{
    vec3::FieldVector3 positionVec = vec3::CoordVector3::convertToFieldVector(positionVector);
    vec3::FieldVector3 pointVec = vec3::CoordVector3::convertToFieldVector(pointVector);

    vec3::FieldVector3 transformedVec = inverseTransformationMatrix * pointVec;
    vec3::FieldVector3 originVec = transformedVec - positionVec;
    vec3::CoordVector3 finalVec = vec3::CoordVector3::convertToCoordVector(originVec);

    finalVec.convertToCylindrical();
    return finalVec;
}

vec3::FieldVector3 Coil::adaptOutputVectorForPoint(const vec3::FieldVector3 &computedVector) const
{
    return transformationMatrix * computedVector;
}


vec3::FieldVector3 Coil::computeBFieldVector(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    vec3::CoordVector3 transformedVector = adaptInputVectorForPoint(pointVector);
    std::pair fields = calculateBField(transformedVector.component1, transformedVector.component2, usedPrecision);

    vec3::FieldVector3 computedVector = vec3::FieldVector3(fields.second * cos(pointVector.component3),
                                                           fields.second * sin(pointVector.component3),
                                                           fields.first);
    return adaptOutputVectorForPoint(computedVector);
}

vec3::FieldVector3 Coil::computeBFieldVector(vec3::CoordVector3 pointVector) const
{
    return computeBFieldVector(pointVector, defaultPrecision);
}

double Coil::computeBFieldX(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return computeBFieldVector(pointVector, usedPrecision).xComponent;
}

double Coil::computeBFieldX(vec3::CoordVector3 pointVector) const
{
    return computeBFieldX(pointVector, defaultPrecision);
}

double Coil::computeBFieldY(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    pointVector.convertToCylindrical();
    return computeBFieldVector(pointVector, usedPrecision).yComponent;
}

double Coil::computeBFieldY(vec3::CoordVector3 pointVector) const
{
    return computeBFieldY(pointVector, defaultPrecision);
}

double Coil::computeBFieldH(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    vec3::FieldVector3 computedVector = computeBFieldVector(pointVector, usedPrecision);
    return std::sqrt(computedVector.xComponent * computedVector.xComponent + computedVector.yComponent * computedVector.yComponent);
}

double Coil::computeBFieldH(vec3::CoordVector3 pointVector) const
{
    return computeBFieldH(pointVector, defaultPrecision);
}

double Coil::computeBFieldZ(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return computeBFieldVector(pointVector, usedPrecision).zComponent;
}

double Coil::computeBFieldZ(vec3::CoordVector3 pointVector) const
{
    return computeBFieldZ(pointVector, defaultPrecision);
}

double Coil::computeBFieldAbs(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    vec3::FieldVector3 computedVector = computeBFieldVector(pointVector, usedPrecision);
    return std::sqrt(computedVector.xComponent * computedVector.xComponent +
                     computedVector.yComponent * computedVector.yComponent +
                     computedVector.zComponent * computedVector.zComponent);
}

double Coil::computeBFieldAbs(vec3::CoordVector3 pointVector) const
{
    return computeBFieldAbs(pointVector, defaultPrecision);
}


vec3::FieldVector3 Coil::computeAPotentialVector(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    vec3::CoordVector3 transformedVector = adaptInputVectorForPoint(pointVector);
    double potential = calculateAPotential(transformedVector.component1, transformedVector.component2, usedPrecision);

    vec3::FieldVector3 computedVector = vec3::FieldVector3(potential * (-1) * sin(pointVector.component3),
                                                           potential * cos(pointVector.component3),
                                                           0.0);
    return adaptOutputVectorForPoint(computedVector);
}

vec3::FieldVector3 Coil::computeAPotentialVector(vec3::CoordVector3 pointVector) const
{
    return computeAPotentialVector(pointVector, defaultPrecision);
}

double Coil::computeAPotentialX(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return computeAPotentialVector(pointVector, usedPrecision).xComponent;
}

double Coil::computeAPotentialX(vec3::CoordVector3 pointVector) const
{
    return computeAPotentialX(pointVector, defaultPrecision);
}

double Coil::computeAPotentialY(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return computeAPotentialVector(pointVector, usedPrecision).yComponent;
}

double Coil::computeAPotentialY(vec3::CoordVector3 pointVector) const
{
    return computeAPotentialY(pointVector, defaultPrecision);
}

double Coil::computeAPotentialZ(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return computeAPotentialVector(pointVector, usedPrecision).zComponent;
}

double Coil::computeAPotentialZ(vec3::CoordVector3 pointVector) const
{
    return computeAPotentialZ(pointVector, defaultPrecision);
}

double Coil::computeAPotentialAbs(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    vec3::FieldVector3 computedVector = computeAPotentialVector(pointVector, usedPrecision);
    return std::sqrt(computedVector.xComponent * computedVector.xComponent +
                     computedVector.yComponent * computedVector.yComponent +
                     computedVector.zComponent * computedVector.zComponent);
}

double Coil::computeAPotentialAbs(vec3::CoordVector3 pointVector) const
{
    return computeAPotentialAbs(pointVector, defaultPrecision);
}


double Coil::computeEFieldX(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialX(pointVector, usedPrecision);
}

double Coil::computeEFieldX(vec3::CoordVector3 pointVector) const
{
    return computeEFieldX(pointVector, defaultPrecision);
}

double Coil::computeEFieldY(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialY(pointVector, usedPrecision);
}

double Coil::computeEFieldY(vec3::CoordVector3 pointVector) const
{
    return computeEFieldY(pointVector, defaultPrecision);
}

double Coil::computeEFieldZ(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialZ(pointVector, usedPrecision);
}

double Coil::computeEFieldZ(vec3::CoordVector3 pointVector) const
{
    return computeEFieldZ(pointVector, defaultPrecision);
}

double Coil::computeEFieldAbs(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialAbs(pointVector, usedPrecision);
}

double Coil::computeEFieldAbs(vec3::CoordVector3 pointVector) const
{
    return computeEFieldAbs(pointVector, defaultPrecision);
}

vec3::FieldVector3 Coil::computeEFieldVector(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    vec3::FieldVector3 computedVector = computeAPotentialVector(pointVector, usedPrecision);
    computedVector *= 2*M_PI * sineFrequency;
    return computedVector;
}

vec3::FieldVector3 Coil::computeEFieldVector(vec3::CoordVector3 pointVector) const
{
    return computeEFieldVector(pointVector, defaultPrecision);
}


vec3::Matrix3 Coil::computeBGradientTensor(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    vec3::CoordVector3 transformedVector = adaptInputVectorForPoint(pointVector);
    double cylindricalZ = transformedVector.component1;
    double cylindricalR = transformedVector.component2;
    double cylindricalPhi = transformedVector.component3;

    std::vector<double> bufferValues = calculateBGradient(cylindricalZ, cylindricalR, usedPrecision);
    double bufferValueRPhi = bufferValues[0];
    double bufferValueRR = bufferValues[1];
    double bufferValueRZ = bufferValues[2];
    double bufferValueZZ = bufferValues[3];

    vec3::Matrix3 computedGradientMatrix;

    if (cylindricalR / innerRadius < g_zAxisApproximationRatio)
    {
        computedGradientMatrix = vec3::Matrix3(bufferValueRR, 0.0, 0.0,
                                               0.0, bufferValueRR, 0.0,
                                               0.0, 0.0, bufferValueZZ);
    }
    else
    {
        double cosPhi = cos(cylindricalPhi);
        double sinPhi = sin(cylindricalPhi);

        double xx = bufferValueRR * cosPhi * cosPhi + bufferValueRPhi * sinPhi * sinPhi;
        double yy = bufferValueRR * sinPhi * sinPhi + bufferValueRPhi * cosPhi * cosPhi;
        double zz = bufferValueZZ;

        double xy = 0.5 * sin(2 * cylindricalPhi) * (bufferValueRR - bufferValueRPhi);
        double xz = bufferValueRZ * cos(cylindricalPhi);
        double yz = bufferValueRZ * sin(cylindricalPhi);
        // the matrix is symmetric
        computedGradientMatrix = vec3::Matrix3(xx, xy, xz,
                                               xy, yy, yz,
                                               xz, yz, zz);
    }
    return transformationMatrix * computedGradientMatrix;
}

vec3::Matrix3 Coil::computeBGradientTensor(vec3::CoordVector3 pointVector) const
{
    return computeBGradientTensor(pointVector, defaultPrecision);
}