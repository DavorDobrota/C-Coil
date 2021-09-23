#include "Coil.h"


std::pair<vec3::FieldVector3, vec3::FieldVector3>
Coil::computeAmpereForce(const Coil &primary, const Coil &secondary, CoilPairArguments forceArguments, ComputeMethod method)
{
    if (isZAxisCase(primary, secondary))
    {
        vec3::FieldVector3 secPositionVec = vec3::CoordVector3::convertToFieldVector(secondary.getPositionVector());
        double zForce = calculateAmpereForceZAxis(primary, secondary, secPositionVec.zComponent, forceArguments, method);

        return std::make_pair(vec3::FieldVector3(0.0, 0.0, zForce), vec3::FieldVector3());
    }
    else
        return calculateAmpereForceGeneral(primary, secondary, forceArguments, method);
}

std::pair<vec3::FieldVector3, vec3::FieldVector3>
Coil::computeAmpereForce(const Coil &primary, const Coil &secondary, PrecisionFactor precisionFactor, ComputeMethod method)
{
    bool zAxisCase = isZAxisCase(primary, secondary);
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, method, zAxisCase);

    return computeAmpereForce(primary, secondary, args, method);
}


std::pair<vec3::FieldVector3, vec3::FieldVector3>
Coil::computeForceOnDipoleMoment(vec3::CoordVector3 pointVector, vec3::FieldVector3 dipoleMoment,
                                 const PrecisionArguments &usedPrecision) const
{
    vec3::FieldVector3 magneticField = computeBFieldVector(pointVector, usedPrecision);
    vec3::Matrix3 magneticGradient = computeBGradientTensor(pointVector, usedPrecision);

    vec3::FieldVector3 magneticTorque = vec3::FieldVector3::crossProduct(dipoleMoment, magneticField);
    vec3::FieldVector3 magneticForce = magneticGradient * dipoleMoment;

    return std::make_pair(magneticForce, magneticTorque);
}

std::pair<vec3::FieldVector3, vec3::FieldVector3>
Coil::computeForceOnDipoleMoment(vec3::CoordVector3 pointVector, vec3::FieldVector3 dipoleMoment) const
{
    return computeForceOnDipoleMoment(pointVector, dipoleMoment, defaultPrecision);
}