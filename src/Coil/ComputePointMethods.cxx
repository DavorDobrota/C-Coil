#include "Coil.h"
#include "PrecisionGlobalVars.h"

#include <cmath>


double Coil::computeBFieldH(vec3::CoordVector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    std::pair<double, double> fields = calculateBField(positionVector.component1, positionVector.component2, usedPrecision);
    return fields.first;
}

double Coil::computeBFieldH(vec3::CoordVector3 positionVector) const
{
    return computeBFieldH(positionVector, defaultPrecision);
}

double Coil::computeBFieldX(vec3::CoordVector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    return computeBFieldH(positionVector, usedPrecision) * std::cos(positionVector.component3);
}

double Coil::computeBFieldX(vec3::CoordVector3 positionVector) const
{
    return computeBFieldX(positionVector, defaultPrecision);
}

double Coil::computeBFieldY(vec3::CoordVector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    return computeBFieldH(positionVector, usedPrecision) * std::sin(positionVector.component3);
}

double Coil::computeBFieldY(vec3::CoordVector3 positionVector) const
{
    return computeBFieldY(positionVector, defaultPrecision);
}

double Coil::computeBFieldZ(vec3::CoordVector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    std::pair<double, double> fields = calculateBField(positionVector.component1, positionVector.component2, usedPrecision);
    return fields.second;
}

double Coil::computeBFieldZ(vec3::CoordVector3 positionVector) const
{
    return computeBFieldZ(positionVector, defaultPrecision);
}

double Coil::computeBFieldAbs(vec3::CoordVector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    std::pair<double, double> fields = calculateBField(positionVector.component1, positionVector.component2, usedPrecision);
    return std::sqrt(fields.first * fields.first + fields.second * fields.second);
}

double Coil::computeBFieldAbs(vec3::CoordVector3 positionVector) const
{
    return computeBFieldAbs(positionVector, defaultPrecision);
}

vec3::FieldVector3 Coil::computeBFieldVector(vec3::CoordVector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    std::pair<double, double> fields = calculateBField(positionVector.component1, positionVector.component2, usedPrecision);

    return vec3::FieldVector3(fields.second * cos(positionVector.component3),
                              fields.second * sin(positionVector.component3),
                              fields.first);
}
vec3::FieldVector3 Coil::computeBFieldVector(vec3::CoordVector3 positionVector) const
{
    return computeBFieldVector(positionVector, defaultPrecision);
}


double Coil::computeAPotentialX(vec3::CoordVector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    return (-1) * calculateAPotential(positionVector.component1, positionVector.component2, usedPrecision) * std::sin(positionVector.component3);
}

double Coil::computeAPotentialX(vec3::CoordVector3 positionVector) const
{
    return computeAPotentialX(positionVector, defaultPrecision);
}

double Coil::computeAPotentialY(vec3::CoordVector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    return calculateAPotential(positionVector.component1, positionVector.component2, usedPrecision) * std::cos(positionVector.component3);
}

double Coil::computeAPotentialY(vec3::CoordVector3 positionVector) const
{
    return computeAPotentialY(positionVector, defaultPrecision);
}

double Coil::computeAPotentialZ(vec3::CoordVector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    // TODO - add functionality in the future
    return 0.0;
}

double Coil::computeAPotentialZ(vec3::CoordVector3 positionVector) const
{
    return computeAPotentialZ(positionVector, defaultPrecision);
}

double Coil::computeAPotentialAbs(vec3::CoordVector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    return calculateAPotential(positionVector.component1, positionVector.component2, usedPrecision);
}

double Coil::computeAPotentialAbs(vec3::CoordVector3 positionVector) const
{
    return computeAPotentialAbs(positionVector, defaultPrecision);
}

vec3::FieldVector3 Coil::computeAPotentialVector(vec3::CoordVector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    double potential = calculateAPotential(positionVector.component1, positionVector.component2, usedPrecision);

    return vec3::FieldVector3(potential * (-1) * sin(positionVector.component3),
                              potential * cos(positionVector.component3),
                              0.0);
}

vec3::FieldVector3 Coil::computeAPotentialVector(vec3::CoordVector3 positionVector) const
{
    return computeAPotentialVector(positionVector, defaultPrecision);
}

double Coil::computeEFieldX(vec3::CoordVector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialX(positionVector, usedPrecision);
}

double Coil::computeEFieldX(vec3::CoordVector3 positionVector) const
{
    return computeEFieldX(positionVector, defaultPrecision);
}

double Coil::computeEFieldY(vec3::CoordVector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialY(positionVector, usedPrecision);
}

double Coil::computeEFieldY(vec3::CoordVector3 positionVector) const
{
    return computeEFieldY(positionVector, defaultPrecision);
}

double Coil::computeEFieldZ(vec3::CoordVector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialZ(positionVector, usedPrecision);
}

double Coil::computeEFieldZ(vec3::CoordVector3 positionVector) const
{
    return computeEFieldZ(positionVector, defaultPrecision);
}

double Coil::computeEFieldAbs(vec3::CoordVector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialAbs(positionVector, usedPrecision);
}

double Coil::computeEFieldAbs(vec3::CoordVector3 positionVector) const
{
    return computeEFieldAbs(positionVector, defaultPrecision);
}

vec3::FieldVector3 Coil::computeEFieldVector(vec3::CoordVector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    vec3::FieldVector3 output = computeAPotentialVector(positionVector, usedPrecision);

    double multiplier = 2*M_PI * sineFrequency;
    output.xComponent *= multiplier;
    output.yComponent *= multiplier;
    output.zComponent *= multiplier;

    return output;
}

vec3::FieldVector3 Coil::computeEFieldVector(vec3::CoordVector3 positionVector) const
{
    return computeEFieldVector(positionVector, defaultPrecision);
}


vec3::Matrix3 Coil::computeBGradientTensor(vec3::CoordVector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    double cylindricalZ = positionVector.component1;
    double cylindricalR = positionVector.component2;
    double cylindricalPhi = positionVector.component3;

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
    return computedGradientMatrix;
}

vec3::Matrix3 Coil::computeBGradientTensor(vec3::CoordVector3 positionVector) const

{
    return computeBGradientTensor(positionVector, defaultPrecision);
}