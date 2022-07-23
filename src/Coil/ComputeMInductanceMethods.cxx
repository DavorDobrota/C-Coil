#include "Coil.h"
#include "ThreadPool.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <functional>

using namespace std::placeholders;


namespace
{
    threadPool::ThreadPoolControl g_threadPool;
}


double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     CoilPairArguments inductanceArguments, ComputeMethod computeMethod)
{
    if (isZAxisCase(primary, secondary))
    {
        vec3::Vector3 secPositionVec = vec3::CoordVector3::convertToFieldVector(secondary.getPositionVector());

        if ((primary.coilType == CoilType::THIN || primary.coilType == CoilType::RECTANGULAR) &&
            (secondary.coilType == CoilType::THIN || secondary.coilType == CoilType::RECTANGULAR))
        {
            return calculateMutualInductanceZAxisFast(primary, secondary, secPositionVec.z, inductanceArguments, computeMethod);
        } else
        {
            return calculateMutualInductanceZAxisSlow(primary, secondary, secPositionVec.z, inductanceArguments, computeMethod);
        }
    }
    else
        return calculateMutualInductanceGeneral(primary, secondary, inductanceArguments, computeMethod);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     PrecisionFactor precisionFactor, ComputeMethod computeMethod)
{
    bool zAxisCase = isZAxisCase(primary, secondary);
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, computeMethod, zAxisCase);

    return computeMutualInductance(primary, secondary, args, computeMethod);
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, CoilPairArguments inductanceArguments,
                                            ComputeMethod computeMethod) const
{
    return computeMutualInductance(*this, secondary, inductanceArguments, computeMethod) * 2*M_PI * sineFrequency;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, PrecisionFactor precisionFactor,
                                            ComputeMethod computeMethod) const
{
    return computeMutualInductance(*this, secondary, precisionFactor, computeMethod) * 2*M_PI * sineFrequency;
}

std::vector<double> Coil::computeAllMutualInductanceArrangements(Coil primary, Coil secondary,
                                                                 const std::vector<vec3::CoordVector3> &primaryPositions,
                                                                 const std::vector<vec3::CoordVector3> &secondaryPositions,
                                                                 const std::vector<double> &primaryYAngles,
                                                                 const std::vector<double> &primaryZAngles,
                                                                 const std::vector<double> &secondaryYAngles,
                                                                 const std::vector<double> &secondaryZAngles,
                                                                 PrecisionFactor precisionFactor,
                                                                 ComputeMethod computeMethod)
{
    std::vector<double> outputMInductances;

    if (primaryPositions.size() == secondaryPositions.size() &&
        primaryPositions.size() == primaryYAngles.size() &&
        primaryPositions.size() == primaryZAngles.size() &&
        primaryPositions.size() == secondaryYAngles.size() &&
        primaryPositions.size() == secondaryZAngles.size())
    {
        unsigned long long numArrangements = primaryPositions.size();
        outputMInductances.resize(numArrangements);

        if (numArrangements < 4 * primary.getThreadCount() || computeMethod != CPU_MT)
        {
            for (int i = 0; i < numArrangements; ++i)
            {
                primary.setPositionAndOrientation(primaryPositions[i], primaryYAngles[i], primaryZAngles[i]);
                secondary.setPositionAndOrientation(secondaryPositions[i], secondaryYAngles[i], secondaryZAngles[i]);

                outputMInductances[i] = Coil::computeMutualInductance(primary, secondary, precisionFactor, computeMethod);
            }
        }
        else
        {
            g_threadPool.setTaskCount(numArrangements);
            g_threadPool.getCompletedTasks().store(0ull);

            auto calcThread = []
                    (
                            int idx,
                            Coil primary,
                            Coil secondary,
                            vec3::CoordVector3 primaryPosition,
                            vec3::CoordVector3 secondaryPosition,
                            double primaryYAngle,
                            double primaryZAngle,
                            double secondaryYAngle,
                            double secondaryZAngle,
                            PrecisionFactor precisionFactor,
                            double &mutualInductance
                    )
            {
                primary.setPositionAndOrientation(primaryPosition, primaryYAngle, primaryZAngle);
                secondary.setPositionAndOrientation(secondaryPosition, secondaryYAngle, secondaryZAngle);

                mutualInductance = Coil::computeMutualInductance(primary, secondary, precisionFactor);

                g_threadPool.getCompletedTasks().fetch_add(1ull);
            };

            for (int i = 0; i < numArrangements; ++i)
            {
                g_threadPool.push(
                        calcThread, std::ref(primary), std::ref(secondary),
                        primaryPositions[i], secondaryPositions[i],
                        primaryYAngles[i], primaryZAngles[i], secondaryYAngles[i], secondaryZAngles[i],
                        precisionFactor, std::ref(outputMInductances[i]));
            }

            g_threadPool.synchronizeThreads();
        }
    }
    else
    {
        fprintf(stderr, "Array sizes do not match!\n");
        throw std::logic_error("Array sizes do not match");
    }

    return outputMInductances;
}
