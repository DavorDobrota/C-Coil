#include "Coil.h"

double Coil::computeAmpereForceZAxis(const Coil &primary, const Coil &secondary, double zDisplacement,
                                     PrecisionFactor precisionFactor, ComputeMethod method)
{
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, method,
                                                                   false);
    return computeAmpereForceZAxis(primary, secondary, zDisplacement, args, method);
}

double Coil::computeAmpereForceZAxis(const Coil &primary, const Coil &secondary, double zDisplacement,
                                     CoilPairArguments forceArguments, ComputeMethod method)
{
    return calculateAmpereForceZAxis(primary, secondary, zDisplacement, forceArguments, method);
}

std::pair<vec3::FieldVector3, vec3::FieldVector3>
        Coil::computeAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                        double zDisplacement, double rDisplacement,
                                        PrecisionFactor precisionFactor, ComputeMethod method)
{
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, method);
    return computeAmpereForceGeneral(primary, secondary, zDisplacement, rDisplacement, args, method);
}

std::pair<vec3::FieldVector3, vec3::FieldVector3>
        Coil::computeAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                        double zDisplacement, double rDisplacement,
                                        CoilPairArguments forceArguments, ComputeMethod method)
{
    return calculateAmpereForceGeneral(primary, secondary, zDisplacement, rDisplacement,
                                       0.0, 0.0, forceArguments, method);
}

std::pair<vec3::FieldVector3, vec3::FieldVector3>
        Coil::computeAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                        double zDisplacement, double rDisplacement, double alphaAngle,
                                        PrecisionFactor precisionFactor, ComputeMethod method)
{
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, method);
    return computeAmpereForceGeneral(primary, secondary, zDisplacement, rDisplacement, alphaAngle, args, method);
}

std::pair<vec3::FieldVector3, vec3::FieldVector3>
        Coil::computeAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                        double zDisplacement, double rDisplacement, double alphaAngle,
                                        CoilPairArguments forceArguments, ComputeMethod method)
{
    return calculateAmpereForceGeneral(primary, secondary, zDisplacement, rDisplacement,
                                       alphaAngle, 0.0, forceArguments, method);
}

std::pair<vec3::FieldVector3, vec3::FieldVector3>
        Coil::computeAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                        double zDisplacement, double rDisplacement, double alphaAngle, double betaAngle,
                                        PrecisionFactor precisionFactor, ComputeMethod method)
{
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, method);
    return computeAmpereForceGeneral(primary, secondary, zDisplacement, rDisplacement, alphaAngle, betaAngle, args, method);
}

std::pair<vec3::FieldVector3, vec3::FieldVector3>
        Coil::computeAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                        double zDisplacement, double rDisplacement, double alphaAngle, double betaAngle,
                                        CoilPairArguments forceArguments, ComputeMethod method)
{
    return calculateAmpereForceGeneral(primary, secondary, zDisplacement, rDisplacement,
                                       alphaAngle, betaAngle, forceArguments, method);
}


std::pair<vec3::FieldVector3, vec3::FieldVector3>
Coil::computeForceOnDipoleMoment(vec3::CoordVector3 positionVector, vec3::FieldVector3 dipoleMoment,
                                 const PrecisionArguments &usedPrecision)
{
    vec3::FieldVector3 magneticField = computeBFieldVector(positionVector, usedPrecision);
    vec3::Matrix3 magneticGradient = computeBGradientTensor(positionVector, usedPrecision);

    vec3::FieldVector3 magneticTorque = vec3::FieldVector3::crossProduct(dipoleMoment, magneticField);
    vec3::FieldVector3 magneticForce = magneticGradient * dipoleMoment;

    return std::make_pair(magneticForce, magneticTorque);
}

std::pair<vec3::FieldVector3, vec3::FieldVector3>
Coil::computeForceOnDipoleMoment(vec3::CoordVector3 positionVector, vec3::FieldVector3 dipoleMoment)
{
    return computeForceOnDipoleMoment(positionVector, dipoleMoment, defaultPrecision);
}