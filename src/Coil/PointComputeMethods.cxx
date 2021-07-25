#include "Coil.h"

#include <cmath>


double Coil::computeBFieldH(double cylindricalZ, double cylindricalR, const PrecisionArguments &usedPrecision) const
{
    std::pair<double, double> fields = calculateBField(cylindricalZ, cylindricalR, usedPrecision);
    return fields.first;
}

double Coil::computeBFieldH(double cylindricalZ, double cylindricalR) const
{
    return computeBFieldH(cylindricalZ, cylindricalR, precisionSettings);
}

double Coil::computeBFieldX(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                            const PrecisionArguments &usedPrecision) const
{
    return computeBFieldH(cylindricalZ, cylindricalR, usedPrecision) * cos(cylindricalPhi);
}

double Coil::computeBFieldX(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeBFieldX(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

double Coil::computeBFieldY(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                            const PrecisionArguments &usedPrecision) const
{
    return computeBFieldH(cylindricalZ, cylindricalR, usedPrecision) * sin(cylindricalPhi);
}

double Coil::computeBFieldY(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeBFieldY(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

double Coil::computeBFieldZ(double cylindricalZ, double cylindricalR, const PrecisionArguments &usedPrecision) const
{
    std::pair<double, double> fields = calculateBField(cylindricalZ, cylindricalR, usedPrecision);
    return fields.second;
}

double Coil::computeBFieldZ(double cylindricalZ, double cylindricalR) const
{
    return computeBFieldZ(cylindricalZ, cylindricalR, precisionSettings);
}

double Coil::computeBFieldAbs(double cylindricalZ, double cylindricalR, const PrecisionArguments &usedPrecision) const
{
    std::pair<double, double> fields = calculateBField(cylindricalZ, cylindricalR, usedPrecision);
    return sqrt(fields.first * fields.first + fields.second * fields.second);
}

double Coil::computeBFieldAbs(double cylindricalZ, double cylindricalR) const
{
    return computeBFieldAbs(cylindricalZ, cylindricalR, precisionSettings);
}

std::vector<double> Coil::computeBFieldVector(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                              const PrecisionArguments &usedPrecision) const
{
    std::vector<double> fieldVector;
    std::pair<double, double> fields = calculateBField(cylindricalZ, cylindricalR, usedPrecision);

    fieldVector.push_back(fields.first * cos(cylindricalPhi));
    fieldVector.push_back(fields.first * sin(cylindricalPhi));
    fieldVector.push_back(fields.second);

    return fieldVector;
}
std::vector<double> Coil::computeBFieldVector(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeBFieldVector(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}


double Coil::computeAPotentialX(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                const PrecisionArguments &usedPrecision) const
{
    return (-1) * calculateAPotential(cylindricalZ, cylindricalR, usedPrecision) * sin(cylindricalPhi);
}

double Coil::computeAPotentialX(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeAPotentialX(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

double Coil::computeAPotentialY(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                const PrecisionArguments &usedPrecision) const
{
    return calculateAPotential(cylindricalZ, cylindricalR, usedPrecision) * cos(cylindricalPhi);
}

double Coil::computeAPotentialY(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeAPotentialY(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

double Coil::computeAPotentialZ(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                const PrecisionArguments &usedPrecision) const
{
    // TODO - add functionality in the future
    return 0.0;
}

double Coil::computeAPotentialZ(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeAPotentialZ(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

double Coil::computeAPotentialAbs(double cylindricalZ, double cylindricalR) const
{
    return calculateAPotential(cylindricalZ, cylindricalR, precisionSettings);
}

double Coil::computeAPotentialAbs(double cylindricalZ, double cylindricalR, const PrecisionArguments &usedPrecision) const
{
    return calculateAPotential(cylindricalZ, cylindricalR, usedPrecision);
}

std::vector<double> Coil::computeAPotentialVector(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                                  const PrecisionArguments &usedPrecision) const
{
    std::vector<double> potentialVector;
    double potential = calculateAPotential(cylindricalZ, cylindricalR, usedPrecision);

    potentialVector.push_back(potential * (-sin(cylindricalPhi)));
    potentialVector.push_back(potential * cos(cylindricalPhi));
    potentialVector.push_back(0.0);
    return potentialVector;
}

std::vector<double> Coil::computeAPotentialVector(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeAPotentialVector(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

double Coil::computeEFieldX(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                            const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialX(cylindricalZ, cylindricalR, cylindricalPhi, usedPrecision);
}

double Coil::computeEFieldX(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeEFieldX(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

double Coil::computeEFieldY(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                            const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialY(cylindricalZ, cylindricalR, cylindricalPhi, usedPrecision);
}

double Coil::computeEFieldY(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeEFieldY(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

double Coil::computeEFieldZ(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                            const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialZ(cylindricalZ, cylindricalR, cylindricalPhi, usedPrecision);
}

double Coil::computeEFieldZ(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeEFieldZ(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}

double Coil::computeEFieldAbs(double cylindricalZ, double cylindricalR, const PrecisionArguments &usedPrecision) const
{
    return 2*M_PI * sineFrequency * computeAPotentialAbs(cylindricalZ, cylindricalR, usedPrecision);
}

double Coil::computeEFieldAbs(double cylindricalZ, double cylindricalR) const
{
    return computeEFieldAbs(cylindricalZ, cylindricalR, precisionSettings);
}

std::vector<double> Coil::computeEFieldVector(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                              const PrecisionArguments &usedPrecision) const
{
    std::vector<double> potentialVector;
    double potential = computeEFieldAbs(cylindricalZ, cylindricalR, usedPrecision);

    potentialVector.push_back(potential * (-sin(cylindricalPhi)));
    potentialVector.push_back(potential * cos(cylindricalPhi));
    potentialVector.push_back(0.0);
    return potentialVector;
}

std::vector<double> Coil::computeEFieldVector(double cylindricalZ, double cylindricalR, double cylindricalPhi) const
{
    return computeEFieldVector(cylindricalZ, cylindricalR, cylindricalPhi, precisionSettings);
}