#include "Coil.h"

#include <math.h>


vec3::Vector3 Coil::calculateAPotential(vec3::Vector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    vec3::Triplet transformedVector = (inverseTransformationMatrix * (pointVector - positionVector)).getAsCylindricalCoords();
    double potential;

    if (useFastMethod)
        potential = calculateAPotentialFast(transformedVector.first, transformedVector.second, usedPrecision);
    else
        potential = calculateAPotentialSlow(transformedVector.first, transformedVector.second, usedPrecision);

    vec3::Vector3 computedVector = vec3::Vector3(potential * (-1) * std::sin(transformedVector.third),
                                                 potential * std::cos(transformedVector.third),
                                                 0.0);
    return transformationMatrix * computedVector;

}

vec3::Vector3 Coil::calculateBField(vec3::Vector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    vec3::Triplet transformedVector = (inverseTransformationMatrix * (pointVector - positionVector)).getAsCylindricalCoords();
    std::pair<double, double> fields;

    if (useFastMethod)
        fields = calculateBFieldFast(transformedVector.first, transformedVector.second, usedPrecision);
    else
        fields = calculateBFieldSlow(transformedVector.first, transformedVector.second, usedPrecision);

    vec3::Vector3 computedVector = vec3::Vector3(fields.first * std::cos(transformedVector.third),
                                                 fields.first * std::sin(transformedVector.third),
                                                 fields.second);
    return transformationMatrix * computedVector;
}

vec3::Matrix3 Coil::calculateBGradient(vec3::Vector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    vec3::Triplet transformedVector = (inverseTransformationMatrix * (pointVector - positionVector)).getAsCylindricalCoords();
    double cylindricalZ = transformedVector.first;
    double cylindricalR = transformedVector.second;
    double cylindricalPhi = transformedVector.third;

    std::vector<double> bufferValues;

    if (useFastMethod)
        bufferValues = calculateBGradientFast(cylindricalZ, cylindricalR, usedPrecision);
    else
        bufferValues = calculateBGradientSlow(cylindricalZ, cylindricalR, usedPrecision);

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
