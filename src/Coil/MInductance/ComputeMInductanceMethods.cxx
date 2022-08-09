#include "Coil.h"

#define _USE_MATH_DEFINES
#include <math.h>


double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     const CoilPairArguments &inductanceArguments, ComputeMethod computeMethod)
{
    if (isZAxisCase(primary, secondary))
    {
        vec3::Vector3 secPositionVec = secondary.getPositionVector();

        if ((primary.coilType == CoilType::THIN || primary.coilType == CoilType::RECTANGULAR) &&
            (secondary.coilType == CoilType::THIN || secondary.coilType == CoilType::RECTANGULAR))
        {
            return calculateMutualInductanceZAxisFast(
                    primary, secondary, secPositionVec.z, inductanceArguments, computeMethod
            );
        }
        else
        {
            return calculateMutualInductanceZAxisSlow(
                    primary, secondary, secPositionVec.z, inductanceArguments, computeMethod
            );
        }
    }
    else
        return calculateMutualInductanceGeneral(primary, secondary, inductanceArguments, computeMethod);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     PrecisionFactor precisionFactor, ComputeMethod computeMethod)
{
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(
            primary, secondary, precisionFactor, computeMethod,isZAxisCase(primary, secondary)
    );

    return computeMutualInductance(primary, secondary, args, computeMethod);
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, const CoilPairArguments &inductanceArguments,
                                            ComputeMethod computeMethod) const
{
    return computeMutualInductance(*this, secondary, inductanceArguments, computeMethod) * 2*M_PI * sineFrequency;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, PrecisionFactor precisionFactor,
                                            ComputeMethod computeMethod) const
{
    return computeMutualInductance(*this, secondary, precisionFactor, computeMethod) * 2*M_PI * sineFrequency;
}


double Coil::computeAndSetSelfInductance(PrecisionFactor precisionFactor)
{
    if (coilType == CoilType::FILAMENT)
        throw std::logic_error("Coil loop self inductance is not defined, integral is divergent");

    // centering the coil at (0, 0, 0) and setting angles to 0 improves the accuracy by leveraging the z-axis formula
    vec3::Vector3 tempPosition = getPositionVector();
    std::pair tempAngles = getRotationAngles();
    setPositionAndOrientation();

    auto arguments = CoilPairArguments::getAppropriateCoilPairArguments(
            *this, *this, precisionFactor
    );
    double inductance;

    if (coilType == CoilType::FLAT)
        inductance = Coil::computeMutualInductance(*this, *this, arguments);
    else
        inductance = Coil::calculateSelfInductance(arguments);

    setSelfInductance(inductance);
    setPositionAndOrientation(tempPosition, tempAngles.first, tempAngles.second);

    return inductance;
}


std::vector<double> Coil::computeAllMutualInductanceArrangements(const Coil &primary, const Coil &secondary,
                                                                 const vec3::Vector3Array &primaryPositions,
                                                                 const vec3::Vector3Array &secondaryPositions,
                                                                 const std::vector<double> &primaryYAngles,
                                                                 const std::vector<double> &primaryZAngles,
                                                                 const std::vector<double> &secondaryYAngles,
                                                                 const std::vector<double> &secondaryZAngles,
                                                                 PrecisionFactor precisionFactor,
                                                                 ComputeMethod computeMethod)
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
            return calculateAllMutualInductanceArrangementsGPU
            (
                primary, secondary, primaryPositions, secondaryPositions,
                primaryYAngles, primaryZAngles, secondaryYAngles, secondaryZAngles,
                precisionFactor
            );
        }
        else if (arrangementCount >= 2 * primary.getThreadCount() && computeMethod == CPU_MT)
        {
            return calculateAllMutualInductanceArrangementsMTD
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

            std::vector<double> outputMInductances;
            outputMInductances.reserve(arrangementCount);

            for (int i = 0; i < arrangementCount; ++i)
            {
                prim.setPositionAndOrientation(
                    primaryPositions[i], primaryYAngles[i], primaryZAngles[i]
                );
                sec.setPositionAndOrientation(
                    secondaryPositions[i], secondaryYAngles[i], secondaryZAngles[i]
                );

                outputMInductances.emplace_back(Coil::computeMutualInductance(prim, sec, precisionFactor, computeMethod));
            }
            return outputMInductances;
        }
    }
    else
        throw std::logic_error("Array sizes do not match!");
}