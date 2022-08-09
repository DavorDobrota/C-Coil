#include "CoilGroup.h"


std::pair<vec3::Vector3, vec3::Vector3>
CoilGroup::computeAmpereForce(const Coil &secondary, PrecisionFactor precisionFactor, ComputeMethod computeMethod) const
{
    if (memberCoils.size() >= 2 * threadCount && computeMethod == CPU_MT)
    {
        return calculateAmpereForceMTD(secondary, precisionFactor);
    }
    else
    {
        vec3::Vector3 totalForce;
        vec3::Vector3 totalTorque;
        std::pair<vec3::Vector3, vec3::Vector3> tempPair;

        for (const auto &memberCoil: memberCoils)
            if (memberCoil->getId() != secondary.getId())
            {
                tempPair = Coil::computeAmpereForce(*memberCoil, secondary, precisionFactor, computeMethod);
                totalForce += tempPair.first;
                totalTorque += tempPair.second;
            }

        return {totalForce, totalTorque};
    }
}


std::pair<vec3::Vector3, vec3::Vector3>
CoilGroup::computeForceOnDipoleMoment(vec3::Vector3 pointVector, vec3::Vector3 dipoleMoment) const
{
    vec3::Vector3 totalForce = vec3::Vector3();
    vec3::Vector3 totalTorque = vec3::Vector3();
    std::pair<vec3::Vector3, vec3::Vector3> tempPair;

    for (const auto& memberCoil : memberCoils)
    {
        tempPair = memberCoil->computeForceOnDipoleMoment(pointVector, dipoleMoment);
        totalForce += tempPair.first;
        totalTorque += tempPair.second;
    }
    return {totalForce, totalTorque};
}


std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
CoilGroup::computeAllAmpereForceArrangements(const Coil &secondary, const vec3::Vector3Array &secondaryPositions,
                                             const std::vector<double> &secondaryYAngles,
                                             const std::vector<double> &secondaryZAngles,
                                             PrecisionFactor precisionFactor, ComputeMethod computeMethod) const
{
    size_t arrangementCount = secondaryPositions.size();

    if (arrangementCount == secondaryPositions.size() &&
        arrangementCount == secondaryYAngles.size() &&
        arrangementCount == secondaryZAngles.size()) {
        if (computeMethod == GPU)
        {
            return calculateAllAmpereForceArrangementsGPU(
                secondary, secondaryPositions, secondaryYAngles, secondaryZAngles, precisionFactor
            );
        }
        else if (arrangementCount >= 2 * this->threadCount && computeMethod == CPU_MT)
        {
            return calculateAllAmpereForceArrangementsMTD(
                secondary, secondaryPositions, secondaryYAngles, secondaryZAngles, precisionFactor
            );
        }
        else
        {
            Coil sec = Coil(secondary);

            std::vector<std::pair<vec3::Vector3, vec3::Vector3>> outputMInductances;
            outputMInductances.reserve(arrangementCount);

            for (int i = 0; i < arrangementCount; ++i)
            {
                sec.setPositionAndOrientation(
                    secondaryPositions[i], secondaryYAngles[i], secondaryZAngles[i]
                );
                outputMInductances.emplace_back(computeAmpereForce(sec, precisionFactor, computeMethod));
            }
            return outputMInductances;
        }
    }
    else
        throw std::logic_error("Array sizes do not match!");
}
