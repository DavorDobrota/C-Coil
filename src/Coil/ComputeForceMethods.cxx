#include "Coil.h"
#include "ThreadPool.h"

#include <functional>

using namespace std::placeholders;


namespace
{
    threadPool::ThreadPoolControl g_threadPool;
}


std::pair<vec3::FieldVector3, vec3::FieldVector3>
Coil::computeAmpereForce(const Coil &primary, const Coil &secondary, CoilPairArguments forceArguments, ComputeMethod computeMethod)
{
    if (isZAxisCase(primary, secondary))
    {
        vec3::FieldVector3 secPositionVec = vec3::CoordVector3::convertToFieldVector(secondary.getPositionVector());
        double zForce = calculateAmpereForceZAxis(primary, secondary, secPositionVec.z, forceArguments, computeMethod);

        return {vec3::FieldVector3(0.0, 0.0, zForce), vec3::FieldVector3()};
    }
    else
        return calculateAmpereForceGeneral(primary, secondary, forceArguments, computeMethod);
}

std::pair<vec3::FieldVector3, vec3::FieldVector3>
Coil::computeAmpereForce(const Coil &primary, const Coil &secondary, PrecisionFactor precisionFactor, ComputeMethod computeMethod)
{
    bool zAxisCase = isZAxisCase(primary, secondary);
    auto args = CoilPairArguments::getAppropriateCoilPairArguments(primary, secondary, precisionFactor, computeMethod, zAxisCase);

    return computeAmpereForce(primary, secondary, args, computeMethod);
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

    if (primaryPositions.size() == secondaryPositions.size() &&
        primaryPositions.size() == primaryYAngles.size() &&
        primaryPositions.size() == primaryZAngles.size() &&
        primaryPositions.size() == secondaryYAngles.size() &&
        primaryPositions.size() == secondaryZAngles.size())
    {
        unsigned long long numArrangements = primaryPositions.size();
        outputForcesAndTorques.resize(numArrangements);

        if (numArrangements < 4 * primary.getThreadCount() || method != CPU_MT)
        {
            for (int i = 0; i < numArrangements; ++i)
            {
                primary.setPositionAndOrientation(primaryPositions[i], primaryYAngles[i], primaryZAngles[i]);
                secondary.setPositionAndOrientation(secondaryPositions[i], secondaryYAngles[i], secondaryZAngles[i]);

                outputForcesAndTorques[i] = Coil::computeAmpereForce(primary, secondary, precisionFactor, method);
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
                            std::pair<vec3::FieldVector3, vec3::FieldVector3> &ampereForce
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
        fprintf(stderr, "Array sizes do not match!\n");
        throw "Array sizes do not match";
    }

    return outputForcesAndTorques;
}