#include "Coil.h"


std::pair<vec3::FieldVector3, vec3::FieldVector3>
Coil::computeAmpereForce(const Coil &primary, const Coil &secondary, CoilPairArguments forceArguments, ComputeMethod method)
{
    if (isZAxisCase(primary, secondary))
    {
        vec3::FieldVector3 secPositionVec = vec3::CoordVector3::convertToFieldVector(secondary.getPositionVector());
        double zForce = calculateAmpereForceZAxis(primary, secondary, secPositionVec.zComponent, forceArguments, method);

        return {vec3::FieldVector3(0.0, 0.0, zForce), vec3::FieldVector3()};
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

    return {magneticForce, magneticTorque};
}

std::pair<vec3::FieldVector3, vec3::FieldVector3>
Coil::computeForceOnDipoleMoment(vec3::CoordVector3 pointVector, vec3::FieldVector3 dipoleMoment) const
{
    return computeForceOnDipoleMoment(pointVector, dipoleMoment, defaultPrecision);
}

std::vector<std::pair<vec3::FieldVector3, vec3::FieldVector3>>
Coil::computeAllAmpereForceArrangements(Coil primary, Coil secondary,
                                        const std::vector<vec3::CoordVector3> &primaryPositions,
                                        const std::vector<vec3::CoordVector3> &secondaryPositions,
                                        const std::vector<double> &primaryYAngles, const std::vector<double> &primaryZAngles,
                                        const std::vector<double> &secondaryYAngles, const std::vector<double> &secondaryZAngles,
                                        PrecisionFactor precisionFactor, ComputeMethod method)
{
    std::vector<std::pair<vec3::FieldVector3, vec3::FieldVector3>> outputForcesAndTorques;

    if (primaryPositions.size() == secondaryPositions.size() ==
        primaryYAngles.size() == primaryZAngles.size() == secondaryYAngles.size() == secondaryZAngles.size())
    {
        unsigned long long size = primaryPositions.size();
        outputForcesAndTorques.resize(size);

        if (size < 4 * primary.getThreadCount() || method != CPU_MT)
        {
            for (int i = 0; i < size; ++i)
            {
                primary.setPositionAndOrientation(primaryPositions[i], primaryYAngles[i], primaryZAngles[i]);
                secondary.setPositionAndOrientation(secondaryPositions[i], secondaryYAngles[i], secondaryZAngles[i]);

                outputForcesAndTorques[i] = Coil::computeAmpereForce(primary, secondary, precisionFactor, method);
            }
        }
        else
        {
            //TODO - implement MTD code either here or in a dedicated method, I think it is alright here, same drill as in CoilGroup
        }
    }
    else
        throw "Array sized do not match";

    return outputForcesAndTorques;
}