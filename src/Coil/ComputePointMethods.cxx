#include "Coil.h"

#define _USE_MATH_DEFINES
#include <math.h>


vec3::CoordVector3 Coil::adaptInputVectorForPoint(const vec3::CoordVector3 &pointVector) const
{
    vec3::Vector3 positionVec = vec3::CoordVector3::convertToFieldVector(positionVector);
    vec3::Vector3 pointVec = vec3::CoordVector3::convertToFieldVector(pointVector);

    vec3::Vector3 transformedVec = inverseTransformationMatrix * (pointVec - positionVec);
    vec3::CoordVector3 finalVec = vec3::CoordVector3::convertToCoordVector(transformedVec);

    finalVec.convertToCylindrical();
    return finalVec;
}

vec3::Vector3 Coil::adaptOutputVectorForPoint(const vec3::Vector3 &computedVector) const
{
    return transformationMatrix * computedVector;
}


vec3::Vector3 Coil::computeBFieldVector(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    vec3::CoordVector3 transformedVector = adaptInputVectorForPoint(pointVector);
    std::pair fields = calculateBField(transformedVector.comp1, transformedVector.comp2, usedPrecision);

    vec3::Vector3 computedVector = vec3::Vector3(fields.first * std::cos(transformedVector.comp3),
                                                           fields.first * std::sin(transformedVector.comp3),
                                                 fields.second);
    return adaptOutputVectorForPoint(computedVector);
}

vec3::Vector3 Coil::computeBFieldVector(vec3::CoordVector3 pointVector) const
{
    return computeBFieldVector(pointVector, defaultPrecisionCPU);
}

double Coil::computeBFieldX(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return computeBFieldVector(pointVector, usedPrecision).x;
}

double Coil::computeBFieldX(vec3::CoordVector3 pointVector) const
{
    return computeBFieldX(pointVector, defaultPrecisionCPU);
}

double Coil::computeBFieldY(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    pointVector.convertToCylindrical();
    return computeBFieldVector(pointVector, usedPrecision).y;
}

double Coil::computeBFieldY(vec3::CoordVector3 pointVector) const
{
    return computeBFieldY(pointVector, defaultPrecisionCPU);
}

double Coil::computeBFieldZ(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return computeBFieldVector(pointVector, usedPrecision).z;
}

double Coil::computeBFieldZ(vec3::CoordVector3 pointVector) const
{
    return computeBFieldZ(pointVector, defaultPrecisionCPU);
}

double Coil::computeBFieldAbs(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    vec3::Vector3 computedVector = computeBFieldVector(pointVector, usedPrecision);
    return std::sqrt(computedVector.x * computedVector.x +
                     computedVector.y * computedVector.y +
                     computedVector.z * computedVector.z);
}

double Coil::computeBFieldAbs(vec3::CoordVector3 pointVector) const
{
    return computeBFieldAbs(pointVector, defaultPrecisionCPU);
}


vec3::Vector3 Coil::computeAPotentialVector(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    vec3::CoordVector3 transformedVector = adaptInputVectorForPoint(pointVector);
    double potential = calculateAPotential(transformedVector.comp1, transformedVector.comp2, usedPrecision);

    vec3::Vector3 computedVector = vec3::Vector3(potential * (-1) * std::sin(transformedVector.comp3),
                                                           potential * std::cos(transformedVector.comp3),
                                                 0.0);
    return adaptOutputVectorForPoint(computedVector);
}

vec3::Vector3 Coil::computeAPotentialVector(vec3::CoordVector3 pointVector) const
{
    return computeAPotentialVector(pointVector, defaultPrecisionCPU);
}

double Coil::computeAPotentialX(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return computeAPotentialVector(pointVector, usedPrecision).x;
}

double Coil::computeAPotentialX(vec3::CoordVector3 pointVector) const
{
    return computeAPotentialX(pointVector, defaultPrecisionCPU);
}

double Coil::computeAPotentialY(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return computeAPotentialVector(pointVector, usedPrecision).y;
}

double Coil::computeAPotentialY(vec3::CoordVector3 pointVector) const
{
    return computeAPotentialY(pointVector, defaultPrecisionCPU);
}

double Coil::computeAPotentialZ(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return computeAPotentialVector(pointVector, usedPrecision).z;
}

double Coil::computeAPotentialZ(vec3::CoordVector3 pointVector) const
{
    return computeAPotentialZ(pointVector, defaultPrecisionCPU);
}

double Coil::computeAPotentialAbs(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    vec3::Vector3 computedVector = computeAPotentialVector(pointVector, usedPrecision);
    return std::sqrt(computedVector.x * computedVector.x +
                     computedVector.y * computedVector.y +
                     computedVector.z * computedVector.z);
}

double Coil::computeAPotentialAbs(vec3::CoordVector3 pointVector) const
{
    return computeAPotentialAbs(pointVector, defaultPrecisionCPU);
}


double Coil::computeEFieldX(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialX(pointVector, usedPrecision);
}

double Coil::computeEFieldX(vec3::CoordVector3 pointVector) const
{
    return computeEFieldX(pointVector, defaultPrecisionCPU);
}

double Coil::computeEFieldY(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialY(pointVector, usedPrecision);
}

double Coil::computeEFieldY(vec3::CoordVector3 pointVector) const
{
    return computeEFieldY(pointVector, defaultPrecisionCPU);
}

double Coil::computeEFieldZ(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialZ(pointVector, usedPrecision);
}

double Coil::computeEFieldZ(vec3::CoordVector3 pointVector) const
{
    return computeEFieldZ(pointVector, defaultPrecisionCPU);
}

double Coil::computeEFieldAbs(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialAbs(pointVector, usedPrecision);
}

double Coil::computeEFieldAbs(vec3::CoordVector3 pointVector) const
{
    return computeEFieldAbs(pointVector, defaultPrecisionCPU);
}

vec3::Vector3 Coil::computeEFieldVector(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    vec3::Vector3 computedVector = computeAPotentialVector(pointVector, usedPrecision);
    computedVector *= 2*M_PI * sineFrequency;
    return computedVector;
}

vec3::Vector3 Coil::computeEFieldVector(vec3::CoordVector3 pointVector) const
{
    return computeEFieldVector(pointVector, defaultPrecisionCPU);
}


vec3::Matrix3 Coil::computeBGradientTensor(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    vec3::CoordVector3 transformedVector = adaptInputVectorForPoint(pointVector);
    double cylindricalZ = transformedVector.comp1;
    double cylindricalR = transformedVector.comp2;
    double cylindricalPhi = transformedVector.comp3;

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
    return computeBGradientTensor(pointVector, defaultPrecisionCPU);
}