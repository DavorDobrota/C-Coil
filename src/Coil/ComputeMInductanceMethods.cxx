#include "Coil.h"
#include "ThreadPool.h"

#include <cmath>
#include <functional>

using namespace std::placeholders;


namespace
{
    threadPool::ThreadPoolControl g_threadPool;
}


double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     CoilPairArguments inductanceArguments, ComputeMethod method)
{
    if (isZAxisCase(primary, secondary))
    {
        vec3::FieldVector3 secPositionVec = vec3::CoordVector3::convertToFieldVector(secondary.getPositionVector());

        if ((primary.coilType == CoilType::THIN || primary.coilType == CoilType::RECTANGULAR) &&
            (secondary.coilType == CoilType::THIN || secondary.coilType == CoilType::RECTANGULAR))
        {
            return calculateMutualInductanceZAxisFast(primary, secondary, secPositionVec.zComponent, inductanceArguments, method);
        } else
        {
            return calculateMutualInductanceZAxisSlow(primary, secondary, secPositionVec.zComponent, inductanceArguments, method);
        }
    }
    else
        return calculateMutualInductanceGeneral(primary, secondary, inductanceArguments, method);
}

double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     PrecisionFactor precisionFactor, ComputeMethod method)
{
    bool zAxisCase = isZAxisCase(primary, secondary);
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, method, zAxisCase);

    return computeMutualInductance(primary, secondary, args, method);
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, CoilPairArguments inductanceArguments,
                                            ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, inductanceArguments, method) * 2*M_PI * sineFrequency;
}

double Coil::computeSecondaryInducedVoltage(const Coil &secondary, PrecisionFactor precisionFactor,
                                            ComputeMethod method) const
{
    return computeMutualInductance(*this, secondary, precisionFactor, method) * 2*M_PI * sineFrequency;
}

std::vector<double> Coil::computeAllMutualInductanceArrangements(Coil primary, Coil secondary,
                                                                 const std::vector<vec3::CoordVector3> &primaryPositions,
                                                                 const std::vector<vec3::CoordVector3> &secondaryPositions,
                                                                 const std::vector<double> &primaryYAngles,
                                                                 const std::vector<double> &primaryZAngles,
                                                                 const std::vector<double> &secondaryYAngles,
                                                                 const std::vector<double> &secondaryZAngles,
                                                                 PrecisionFactor precisionFactor,
                                                                 ComputeMethod method)
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

        if (numArrangements < 4 * primary.getThreadCount() || method != CPU_MT)
        {
            for (int i = 0; i < numArrangements; ++i)
            {
                primary.setPositionAndOrientation(primaryPositions[i], primaryYAngles[i], primaryZAngles[i]);
                secondary.setPositionAndOrientation(secondaryPositions[i], secondaryYAngles[i], secondaryZAngles[i]);

                outputMInductances[i] = Coil::computeMutualInductance(primary, secondary, precisionFactor, method);
            }
        }
        else
        {
            //TODO - implement MTD code either here or in a dedicated method, I think it is alright here, same drill as in CoilGroup
            g_threadPool.setTaskCount(numArrangements);
            g_threadPool.getCompletedTasks().store(0ull);

            auto calcThread = []
                    (
                            int idx,
                            const Coil &coil,
                            const Coil &secondary,
                            PrecisionFactor precisionFactor,
                            double &mutualInductance
                    )
            {
                mutualInductance = Coil::computeMutualInductance(coil, secondary, precisionFactor);

                g_threadPool.getCompletedTasks().fetch_add(1ull);
            };

            for (int i = 0; i < numArrangements; ++i)
            {
                primary.setPositionAndOrientation(primaryPositions[i], primaryYAngles[i], primaryZAngles[i]);
                secondary.setPositionAndOrientation(secondaryPositions[i], secondaryYAngles[i], secondaryZAngles[i]);

                g_threadPool.push(
                        calcThread, std::ref(primary), std::ref(secondary), precisionFactor,
                        std::ref(outputMInductances[i]));
            }

            g_threadPool.synchronizeThreads();
        }
    }
    else
    {
        fprintf(stderr, "Array sizes do not match!\n");
        throw "Array sizes do not match";
    }

    return outputMInductances;
}