#include "Coil.h"

#include <cmath>


double Coil::computeBFieldH(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    std::pair<double, double> fields = calculateBField(positionVector.component1, positionVector.component2, usedPrecision);
    return fields.first;
}

double Coil::computeBFieldH(vec3::Vector3 positionVector) const
{
    return computeBFieldH(positionVector, defaultPrecision);
}

double Coil::computeBFieldX(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    return computeBFieldH(positionVector, usedPrecision) * std::cos(positionVector.component3);
}

double Coil::computeBFieldX(vec3::Vector3 positionVector) const
{
    return computeBFieldX(positionVector, defaultPrecision);
}

double Coil::computeBFieldY(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    return computeBFieldH(positionVector, usedPrecision) * std::sin(positionVector.component3);
}

double Coil::computeBFieldY(vec3::Vector3 positionVector) const
{
    return computeBFieldY(positionVector, defaultPrecision);
}

double Coil::computeBFieldZ(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    std::pair<double, double> fields = calculateBField(positionVector.component1, positionVector.component2, usedPrecision);
    return fields.second;
}

double Coil::computeBFieldZ(vec3::Vector3 positionVector) const
{
    return computeBFieldZ(positionVector, defaultPrecision);
}

double Coil::computeBFieldAbs(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    std::pair<double, double> fields = calculateBField(positionVector.component1, positionVector.component2, usedPrecision);
    return sqrt(fields.first * fields.first + fields.second * fields.second);
}

double Coil::computeBFieldAbs(vec3::Vector3 positionVector) const
{
    return computeBFieldAbs(positionVector, defaultPrecision);
}

vec3::Vector3 Coil::computeBFieldVector(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    std::pair<double, double> fields = calculateBField(positionVector.component1, positionVector.component2, usedPrecision);

    vec3::Vector3 output = vec3::Vector3::getVectorCylindrical(fields.second, fields.first, positionVector.component3);
    output.convertToCartesian();

    return output;
}
vec3::Vector3 Coil::computeBFieldVector(vec3::Vector3 positionVector) const
{
    return computeBFieldVector(positionVector, defaultPrecision);
}


double Coil::computeAPotentialX(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    return (-1) * calculateAPotential(positionVector.component1, positionVector.component2, usedPrecision) * std::sin(positionVector.component3);
}

double Coil::computeAPotentialX(vec3::Vector3 positionVector) const
{
    return computeAPotentialX(positionVector, defaultPrecision);
}

double Coil::computeAPotentialY(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    return calculateAPotential(positionVector.component1, positionVector.component2, usedPrecision) * std::cos(positionVector.component3);
}

double Coil::computeAPotentialY(vec3::Vector3 positionVector) const
{
    return computeAPotentialY(positionVector, defaultPrecision);
}

double Coil::computeAPotentialZ(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    // TODO - add functionality in the future
    return 0.0;
}

double Coil::computeAPotentialZ(vec3::Vector3 positionVector) const
{
    return computeAPotentialZ(positionVector, defaultPrecision);
}

double Coil::computeAPotentialAbs(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    return calculateAPotential(positionVector.component1, positionVector.component2, usedPrecision);
}

double Coil::computeAPotentialAbs(vec3::Vector3 positionVector) const
{
    return computeAPotentialAbs(positionVector, defaultPrecision);
}

vec3::Vector3 Coil::computeAPotentialVector(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    positionVector.convertToCylindrical();
    double potential = calculateAPotential(positionVector.component1, positionVector.component2, usedPrecision);

    return vec3::Vector3::getVectorCartesian(potential * (-1) * std::sin(positionVector.component3),
                                       potential * std::cos(positionVector.component3),
                                       0.0);
}

vec3::Vector3 Coil::computeAPotentialVector(vec3::Vector3 positionVector) const
{
    return computeAPotentialVector(positionVector, defaultPrecision);
}

double Coil::computeEFieldX(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialX(positionVector, usedPrecision);
}

double Coil::computeEFieldX(vec3::Vector3 positionVector) const
{
    return computeEFieldX(positionVector, defaultPrecision);
}

double Coil::computeEFieldY(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialY(positionVector, usedPrecision);
}

double Coil::computeEFieldY(vec3::Vector3 positionVector) const
{
    return computeEFieldY(positionVector, defaultPrecision);
}

double Coil::computeEFieldZ(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialZ(positionVector, usedPrecision);
}

double Coil::computeEFieldZ(vec3::Vector3 positionVector) const
{
    return computeEFieldZ(positionVector, defaultPrecision);
}

double Coil::computeEFieldAbs(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialAbs(positionVector, usedPrecision);
}

double Coil::computeEFieldAbs(Vector3 positionVector) const
{
    return computeEFieldAbs(positionVector, defaultPrecision);
}

Vector3 Coil::computeEFieldVector(Vector3 positionVector, const PrecisionArguments &usedPrecision) const
{
    Vector3 output = computeAPotentialVector(positionVector, usedPrecision);

    double multiplier = 2*M_PI * sineFrequency;
    output.component1 *= multiplier;
    output.component2 *= multiplier;
    output.component3 *= multiplier;

    return output;
}

Vector3 Coil::computeEFieldVector(Vector3 positionVector) const
{
    return computeEFieldVector(positionVector, defaultPrecision);
}

std::vector<double> Coil::computeBGradientTensor(Vector3 positionVector, const PrecisionArguments &usedPrecision) const
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

    std::vector<double> outputValues;

    if (cylindricalR / innerRadius < 1e-14)
    {
        outputValues.push_back(bufferValueRR);
        outputValues.push_back(0.0);
        outputValues.push_back(0.0);

        outputValues.push_back(0.0);
        outputValues.push_back(bufferValueRR);
        outputValues.push_back(0.0);

        outputValues.push_back(0.0);
        outputValues.push_back(0.0);
        outputValues.push_back(bufferValueZZ);
    }
    else
    {
        double cosPhi = cos(cylindricalPhi);
        double sinPhi = sin(cylindricalPhi);

        outputValues.push_back(bufferValueRR * cosPhi * cosPhi + bufferValueRPhi * sinPhi* sinPhi);
        outputValues.push_back(0.5 * sin(2 * cylindricalPhi) * (bufferValueRR - bufferValueRPhi));
        outputValues.push_back(bufferValueRZ * cos(cylindricalPhi));

        outputValues.push_back(0.5 * sin(2 * cylindricalPhi) * (bufferValueRR - bufferValueRPhi));
        outputValues.push_back(bufferValueRR * sinPhi * sinPhi + bufferValueRPhi * cosPhi * cosPhi);
        outputValues.push_back(bufferValueRZ * sin(cylindricalPhi));

        outputValues.push_back(bufferValueRZ * cosPhi);
        outputValues.push_back(bufferValueRZ * sinPhi);
        outputValues.push_back(bufferValueZZ);
    }
    return outputValues;
}

std::vector<double> Coil::computeBGradientTensor(Vector3 positionVector) const

{
    return computeBGradientTensor(positionVector, defaultPrecision);
}