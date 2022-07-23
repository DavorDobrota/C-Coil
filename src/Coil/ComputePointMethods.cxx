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


vec3::Matrix3 Coil::computeBGradientMatrix(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const
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

vec3::Matrix3 Coil::computeBGradientMatrix(vec3::CoordVector3 pointVector) const
{
    return computeBGradientMatrix(pointVector, defaultPrecisionCPU);
}