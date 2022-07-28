#include "Coil.h"
#include "ThreadPool.h"

#include <functional>

using namespace std::placeholders;


namespace
{
    threadPool::ThreadPoolControl g_threadPool;
}


std::pair<vec3::Vector3, vec3::Vector3>
Coil::computeAmpereForce(const Coil &primary, const Coil &secondary, CoilPairArguments forceArguments, ComputeMethod computeMethod)
{
    if (isZAxisCase(primary, secondary))
    {
        vec3::Vector3 secPositionVec = secondary.getPositionVector();
        double zForce = 0.0;

        if ((primary.coilType == CoilType::THIN || primary.coilType == CoilType::RECTANGULAR) &&
            (secondary.coilType == CoilType::THIN || secondary.coilType == CoilType::RECTANGULAR))
        {
            zForce = calculateAmpereForceZAxisFast(primary, secondary, secPositionVec.z, forceArguments, computeMethod);
        } else
        {
            zForce = calculateAmpereForceZAxisSlow(primary, secondary, secPositionVec.z, forceArguments, computeMethod);
        }

        return {vec3::Vector3(0.0, 0.0, zForce), vec3::Vector3()};
    }
    else
        return calculateAmpereForceGeneral(primary, secondary, forceArguments, computeMethod);
}

std::pair<vec3::Vector3, vec3::Vector3>
Coil::computeAmpereForce(const Coil &primary, const Coil &secondary, PrecisionFactor precisionFactor, ComputeMethod computeMethod)
{
    bool zAxisCase = isZAxisCase(primary, secondary);
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, computeMethod, zAxisCase);

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
Coil::computeAllAmpereForceArrangements(Coil primary, Coil secondary,
                                        const vec3::Vector3Array &primaryPositions,
                                        const vec3::Vector3Array &secondaryPositions,
                                        const std::vector<double> &primaryYAngles, const std::vector<double> &primaryZAngles,
                                        const std::vector<double> &secondaryYAngles, const std::vector<double> &secondaryZAngles,
                                        PrecisionFactor precisionFactor, ComputeMethod computeMethod)
{
    std::vector<std::pair<vec3::Vector3, vec3::Vector3>> outputForcesAndTorques;

    if (primaryPositions.size() == secondaryPositions.size() &&
        primaryPositions.size() == primaryYAngles.size() &&
        primaryPositions.size() == primaryZAngles.size() &&
        primaryPositions.size() == secondaryYAngles.size() &&
        primaryPositions.size() == secondaryZAngles.size())
    {
        unsigned long long numArrangements = primaryPositions.size();
        outputForcesAndTorques.resize(numArrangements);

        if (numArrangements < 4 * primary.getThreadCount() || computeMethod != CPU_MT)
        {
            for (int i = 0; i < numArrangements; ++i)
            {
                primary.setPositionAndOrientation(primaryPositions[i], primaryYAngles[i], primaryZAngles[i]);
                secondary.setPositionAndOrientation(secondaryPositions[i], secondaryYAngles[i], secondaryZAngles[i]);

                outputForcesAndTorques[i] = Coil::computeAmpereForce(primary, secondary, precisionFactor, computeMethod);
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
                            vec3::Vector3 primaryPosition,
                            vec3::Vector3 secondaryPosition,
                            double primaryYAngle,
                            double primaryZAngle,
                            double secondaryYAngle,
                            double secondaryZAngle,
                            PrecisionFactor precisionFactor,
                            std::pair<vec3::Vector3, vec3::Vector3> &ampereForce
                    )
            {
                primary.setPositionAndOrientation(primaryPosition, primaryYAngle, primaryZAngle);
                secondary.setPositionAndOrientation(secondaryPosition, secondaryYAngle, secondaryZAngle);

                ampereForce = Coil::computeAmpereForce(primary, secondary, precisionFactor);

                g_threadPool.getCompletedTasks().fetch_add(1ull);
            };

            for (int i = 0; i < numArrangements; ++i)
            {
                primary.setPositionAndOrientation(primaryPositions[i], primaryYAngles[i], primaryZAngles[i]);
                secondary.setPositionAndOrientation(secondaryPositions[i], secondaryYAngles[i], secondaryZAngles[i]);

                g_threadPool.push(
                        calcThread, std::ref(primary), std::ref(secondary),
                        primaryPositions[i], secondaryPositions[i],
                        primaryYAngles[i], primaryZAngles[i], secondaryYAngles[i], secondaryZAngles[i],
                        precisionFactor, std::ref(outputForcesAndTorques[i]));
            }

            g_threadPool.synchronizeThreads();
        }
    }
    else
    {
        throw std::logic_error("Array sizes do not match");
    }

    return outputForcesAndTorques;
}