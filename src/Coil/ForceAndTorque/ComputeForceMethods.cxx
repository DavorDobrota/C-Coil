#include "Coil.h"


std::pair<vec3::Vector3, vec3::Vector3>
Coil::computeAmpereForce(const Coil &primary, const Coil &secondary,
                         CoilPairArguments forceArguments, ComputeMethod computeMethod)
{
    if (isZAxisCase(primary, secondary))
    {
        vec3::Vector3 secPositionVec = secondary.getPositionVector();
        double zForce = 0.0;

        if ((primary.coilType == CoilType::THIN || primary.coilType == CoilType::RECTANGULAR) &&
            (secondary.coilType == CoilType::THIN || secondary.coilType == CoilType::RECTANGULAR))

            zForce = calculateAmpereForceZAxisFast(primary, secondary, secPositionVec.z, forceArguments, computeMethod);
        else
            zForce = calculateAmpereForceZAxisSlow(primary, secondary, secPositionVec.z, forceArguments, computeMethod);

        return {vec3::Vector3(0.0, 0.0, zForce), vec3::Vector3()};
    }
    else
        return calculateAmpereForceGeneral(primary, secondary, forceArguments, computeMethod);
}

std::pair<vec3::Vector3, vec3::Vector3>
Coil::computeAmpereForce(const Coil &primary, const Coil &secondary,
                         PrecisionFactor precisionFactor, ComputeMethod computeMethod)
{
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(
            primary, secondary, precisionFactor, computeMethod, isZAxisCase(primary, secondary)
    );

    return computeAmpereForce(primary, secondary, args, computeMethod);
}


std::pair<vec3::Vector3, vec3::Vector3>
Coil::computeForceOnDipoleMoment(vec3::Vector3 pointVector, vec3::Vector3 dipoleMoment,
                                 const PrecisionArguments &usedPrecision) const
{
    vec3::Vector3 magneticField = computeBFieldVector(pointVector, usedPrecision);
    vec3::Matrix3 magneticGradient = computeBGradientMatrix(pointVector, usedPrecision);

    vec3::Vector3 magneticTorque = vec3::Vector3::crossProduct(dipoleMoment, magneticField);
    vec3::Vector3 magneticForce = magneticGradient * dipoleMoment;

    return {magneticForce, magneticTorque};
}

std::pair<vec3::Vector3, vec3::Vector3>
Coil::computeForceOnDipoleMoment(vec3::Vector3 pointVector, vec3::Vector3 dipoleMoment) const
{
    return computeForceOnDipoleMoment(pointVector, dipoleMoment, defaultPrecisionCPU);
}


std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
Coil::computeAllAmpereForceArrangements(const Coil &primary, const Coil &secondary,
                                        const vec3::Vector3Array &primaryPositions,
                                        const vec3::Vector3Array &secondaryPositions,
                                        const std::vector<double> &primaryYAngles,
                                        const std::vector<double> &primaryZAngles,
                                        const std::vector<double> &secondaryYAngles,
                                        const std::vector<double> &secondaryZAngles,
                                        PrecisionFactor precisionFactor, ComputeMethod computeMethod)
{
    size_t arrangementCount = primaryPositions.size();

    if (arrangementCount == secondaryPositions.size() &&
        arrangementCount == primaryYAngles.size() &&
        arrangementCount == primaryZAngles.size() &&
        arrangementCount == secondaryYAngles.size() &&
        arrangementCount == secondaryZAngles.size())
    {
        if (computeMethod == GPU)
        {
            return calculateAllAmpereForceArrangementsGPU
            (
                primary, secondary,primaryPositions, secondaryPositions,
                primaryYAngles, primaryZAngles, secondaryYAngles, secondaryZAngles,
                precisionFactor
            );
        }
        else if (arrangementCount >= 2 * primary.getThreadCount() && computeMethod == CPU_MT)
        {
            return calculateAllAmpereForceArrangementsMTD
            (
                primary, secondary, primaryPositions, secondaryPositions,
                primaryYAngles, primaryZAngles, secondaryYAngles, secondaryZAngles,
                precisionFactor
            );
        }
        else
        {
            Coil prim = Coil(primary);
            Coil sec = Coil(secondary);

            std::vector<std::pair<vec3::Vector3, vec3::Vector3>> outputForcesAndTorques;

            for (int i = 0; i < arrangementCount; ++i)
            {
                prim.setPositionAndOrientation(
                    primaryPositions[i], primaryYAngles[i],primaryZAngles[i]
                );
                sec.setPositionAndOrientation(
                    secondaryPositions[i],secondaryYAngles[i],secondaryZAngles[i]
                );

                outputForcesAndTorques.emplace_back(computeAmpereForce(prim, sec, precisionFactor, computeMethod));
            }
            return outputForcesAndTorques;
        }
    }
    else
        throw std::logic_error("Array sizes do not match");
}