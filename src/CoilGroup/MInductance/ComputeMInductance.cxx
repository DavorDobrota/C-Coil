#include "CoilGroup.h"


double CoilGroup::computeMutualInductance(const Coil &secondary, PrecisionFactor precisionFactor, ComputeMethod computeMethod) const
{
    if (memberCoils.size() >= 2 * threadCount && computeMethod == CPU_MT)
    {
        return calculateMutualInductanceMTD(secondary, precisionFactor);
    }
    else
    {
        double totalMutualInductance = 0.0;

        for (const auto &memberCoil: memberCoils)
            if (memberCoil->getId() != secondary.getId())
                totalMutualInductance += Coil::computeMutualInductance(*memberCoil, secondary, precisionFactor, computeMethod);

        return totalMutualInductance;
    }
}


std::vector<double> CoilGroup::computeAllMutualInductanceArrangements(const Coil &secondary,
                                                                      const vec3::Vector3Array &secondaryPositions,
                                                                      const std::vector<double> &secondaryYAngles,
                                                                      const std::vector<double> &secondaryZAngles,
                                                                      PrecisionFactor precisionFactor,
                                                                      ComputeMethod computeMethod) const
{
    size_t arrangementCount = secondaryPositions.size();

    if (arrangementCount == secondaryPositions.size() &&
        arrangementCount == secondaryYAngles.size() &&
        arrangementCount == secondaryZAngles.size()) {
        if (computeMethod == GPU)
        {
            return calculateAllMutualInductanceArrangementsGPU(
                    secondary, secondaryPositions, secondaryYAngles, secondaryZAngles,precisionFactor
            );
        }
        else if (arrangementCount >= 2 * this->threadCount && computeMethod == CPU_MT)
        {
            return calculateAllMutualInductanceArrangementsMTD(
                    secondary, secondaryPositions, secondaryYAngles, secondaryZAngles, precisionFactor
            );
        }
        else
        {
            Coil sec = Coil(secondary);

            std::vector<double> outputMInductances;
            outputMInductances.reserve(arrangementCount);

            for (int i = 0; i < arrangementCount; ++i)
            {
                sec.setPositionAndOrientation(secondaryPositions[i], secondaryYAngles[i], secondaryZAngles[i]);
                outputMInductances.emplace_back(computeMutualInductance(sec, precisionFactor, computeMethod));
            }
            return outputMInductances;
        }
    }
    else
        throw std::logic_error("Array sizes do not match!");
}