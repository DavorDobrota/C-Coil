#include "Coil.h"
#include "ThreadPool.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <sstream>


namespace
{
    const double g_MiReduced = 0.0000001;

    const double g_defaultCurrent = 1.0;
    const double g_defaultResistivity = 1.63e-8;
    const double g_defaultSineFrequency = 50;

    unsigned long long g_currentId = 0;

    threadPool::ThreadPoolControl g_threadPool;
}


Coil::Coil() : Coil(0.0, 0.0, 0.0, 3600, 0) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current, double wireResistivity,
           double sineFrequency, const PrecisionArguments &precisionSettingsCPU,
           const PrecisionArguments &precisionSettingsGPU, int threadCount, vec3::Vector3 coordinatePosition,
           double yAxisAngle, double zAxisAngle) :
        innerRadius(innerRadius),
        thickness(thickness / innerRadius < g_thinCoilApproximationRatio ? 1e-18 * innerRadius * g_thinCoilApproximationRatio : thickness),
        length(length / innerRadius < g_thinCoilApproximationRatio ? 1e-18 * innerRadius * g_thinCoilApproximationRatio : length),
        numOfTurns(numOfTurns),
        id(++g_currentId),
        defaultPrecisionCPU(precisionSettingsCPU),
        defaultPrecisionGPU(precisionSettingsGPU),
        useFastMethod(length / innerRadius >= g_thinCoilApproximationRatio),
        positionVector(coordinatePosition),
        yAxisAngle(yAxisAngle), zAxisAngle(zAxisAngle)
{
    calculateTransformationMatrices();
    calculateCoilType();
    setCurrent(current);
    calculateAverageWireThickness();
    setWireResistivity(wireResistivity);
    setSineFrequency(sineFrequency);
    setThreadCount(threadCount);
}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
           double wireResistivity, double sineFrequency, PrecisionFactor precisionFactor,
           int threadCount, vec3::Vector3 coordinatePosition, double yAxisAngle, double zAxisAngle) :
        innerRadius(innerRadius),
        thickness(thickness / innerRadius < g_thinCoilApproximationRatio ? 1e-18 * innerRadius * g_thinCoilApproximationRatio : thickness),
        length(length / innerRadius < g_thinCoilApproximationRatio ? 1e-18 * innerRadius * g_thinCoilApproximationRatio : length),
        numOfTurns(numOfTurns),
        id(++g_currentId),
        useFastMethod(length / innerRadius >= g_thinCoilApproximationRatio),
        positionVector(coordinatePosition),
        yAxisAngle(yAxisAngle), zAxisAngle(zAxisAngle)
{
    calculateTransformationMatrices();
    calculateCoilType();
    setCurrent(current);
    calculateAverageWireThickness();
    setWireResistivity(wireResistivity);
    setSineFrequency(sineFrequency);
    setDefaultPrecisionCPU(PrecisionArguments::getCoilPrecisionArgumentsCPU(*this, precisionFactor));
    setDefaultPrecisionGPU(PrecisionArguments::getCoilPrecisionArgumentsGPU(*this, precisionFactor));
    setThreadCount(threadCount);
}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current, double sineFrequency,
           PrecisionFactor precisionFactor, int threadCount, vec3::Vector3 coordinatePosition,
           double yAxisAngle, double zAxisAngle) :
           Coil(innerRadius, thickness, length, numOfTurns, current, g_defaultResistivity, sineFrequency,
                precisionFactor, threadCount, coordinatePosition, yAxisAngle, zAxisAngle) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current, double sineFrequency,
           const PrecisionArguments &precisionSettingsCPU, const PrecisionArguments &precisionSettingsGPU, int threadCount,
           vec3::Vector3 coordinatePosition, double yAxisAngle, double zAxisAngle) :
           Coil(innerRadius, thickness, length, numOfTurns, current,
                g_defaultResistivity, sineFrequency, precisionSettingsCPU, precisionSettingsGPU, threadCount,
                coordinatePosition, yAxisAngle, zAxisAngle) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
           PrecisionFactor precisionFactor, int threadCount, vec3::Vector3 coordinatePosition,
           double yAxisAngle, double zAxisAngle)  :
           Coil(innerRadius, thickness, length, numOfTurns, current, g_defaultResistivity,
                g_defaultSineFrequency, precisionFactor, threadCount) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
           const PrecisionArguments &precisionSettingsCPU, const PrecisionArguments &precisionSettingsGPU, int threadCount,
           vec3::Vector3 coordinatePosition, double yAxisAngle, double zAxisAngle) :
           Coil(innerRadius, thickness, length, numOfTurns, current, g_defaultResistivity,
                g_defaultSineFrequency, precisionSettingsCPU, precisionSettingsGPU, threadCount,
                coordinatePosition, yAxisAngle, zAxisAngle) {}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns, PrecisionFactor precisionFactor,
           int threadCount, vec3::Vector3 coordinatePosition, double yAxisAngle, double zAxisAngle) :
           Coil(innerRadius, thickness, length, numOfTurns, g_defaultCurrent, g_defaultResistivity,
                g_defaultSineFrequency, precisionFactor, threadCount, coordinatePosition, yAxisAngle, zAxisAngle){}

Coil::Coil(double innerRadius, double thickness, double length, int numOfTurns,
           const PrecisionArguments &precisionSettingsCPU, const PrecisionArguments &precisionSettingsGPU, int threadCount,
           vec3::Vector3 coordinatePosition, double yAxisAngle, double zAxisAngle) :
           Coil(innerRadius, thickness, length, numOfTurns, g_defaultCurrent, g_defaultResistivity,
                g_defaultSineFrequency, precisionSettingsCPU, precisionSettingsGPU, threadCount,
                coordinatePosition, yAxisAngle, zAxisAngle) {}


double Coil::getCurrentDensity() const { return currentDensity; }

double Coil::getCurrent() const { return current; }

int Coil::getNumOfTurns() const { return numOfTurns; }

unsigned long long Coil::getId() const { return id; }

double Coil::getInnerRadius() const { return innerRadius; }

double Coil::getThickness() const { return thickness; }

double Coil::getLength() const { return length; }

double Coil::getAverageWireThickness() const { return averageWireThickness; }

bool Coil::isSineDriven() const { return sineDriven; }

double Coil::getSineFrequency() const { return sineFrequency; }

vec3::Vector3 Coil::getMagneticMoment()
{
    calculateMagneticMoment();

    double sinA = std::sin(yAxisAngle); double cosA = std::cos(yAxisAngle);
    double sinB = std::sin(zAxisAngle); double cosB = std::cos(zAxisAngle);

    return vec3::Vector3(magneticMoment * sinA * cosB,magneticMoment * sinA * sinB,magneticMoment * cosA);
}

double Coil::getWireResistivity() const { return wireResistivity; }

double Coil::getSelfInductance() const { return selfInductance; }

double Coil::getResistance()
{
    calculateResistance();
    return resistance;
}

double Coil::getReactance()
{
    calculateReactance();
    return reactance;
}

double Coil::getImpedance()
{
    calculateImpedance();
    return impedance;
}

const PrecisionArguments &Coil::getPrecisionSettingsCPU() const { return defaultPrecisionCPU; }

const PrecisionArguments &Coil::getPrecisionSettingsGPU() const { return defaultPrecisionGPU; }

int Coil::getThreadCount() const { return threadCount; }

bool Coil::isUsingFastMethod() const { return useFastMethod; }

CoilType Coil::getCoilType() const { return coilType; }

vec3::Vector3 Coil::getPositionVector() const { return positionVector; }

std::pair<double, double> Coil::getRotationAngles() const { return {yAxisAngle, zAxisAngle}; }

vec3::Matrix3 Coil::getTransformationMatrix() const { return transformationMatrix; }

vec3::Matrix3 Coil::getInverseTransformationMatrix() const { return inverseTransformationMatrix; }


void Coil::setCurrentDensity(double currentDensity)
{
    this->currentDensity = currentDensity;
    current = currentDensity * length * thickness / numOfTurns;
}

void Coil::setCurrent(double current)
{
    this->current = current;
    currentDensity = current * numOfTurns / (length * thickness);
}

void Coil::setWireResistivity(double wireResistivity)
{
    this->wireResistivity = wireResistivity;
}

void Coil::setSineFrequency(double sineFrequency)
{
    if (sineFrequency > 0.0)
    {
        sineDriven = true;
        this->sineFrequency = sineFrequency;
    }
    else
    {
        sineDriven = false;
        this->sineFrequency = 0.0;
    }
    calculateImpedance();
}

void Coil::setDefaultPrecisionCPU(const PrecisionArguments &precisionSettings)
{
    this->defaultPrecisionCPU = precisionSettings;
}

void Coil::setDefaultPrecisionCPU(PrecisionFactor precisionFactor)
{
    this->defaultPrecisionCPU = PrecisionArguments::getCoilPrecisionArgumentsCPU(*this, precisionFactor);
}

void Coil::setDefaultPrecisionGPU(const PrecisionArguments &precisionSettings)
{
    this->defaultPrecisionGPU = precisionSettings;
}

void Coil::setDefaultPrecisionGPU(PrecisionFactor precisionFactor)
{
    this->defaultPrecisionGPU = PrecisionArguments::getCoilPrecisionArgumentsGPU(*this, precisionFactor);
}

void Coil::setDefaultPrecision(PrecisionFactor precisionFactor)
{
    this->defaultPrecisionCPU = PrecisionArguments::getCoilPrecisionArgumentsCPU(*this, precisionFactor);
    this->defaultPrecisionGPU = PrecisionArguments::getCoilPrecisionArgumentsGPU(*this, precisionFactor);
}


void Coil::setThreadCount(int threadCount)
{
    Coil::threadCount = threadCount;
    g_threadPool.setSize(threadCount);
}

void Coil::setPositionAndOrientation(vec3::Vector3 positionVector, double yAxisAngle, double zAxisAngle)
{
    this->positionVector = positionVector;
    this->yAxisAngle = yAxisAngle;
    this->zAxisAngle = zAxisAngle;

    calculateTransformationMatrices();
}

void Coil::setSelfInductance(double selfInductance)
{
    this->selfInductance = selfInductance;
}


vec3::Vector3 Coil::computeAPotentialVector(vec3::Vector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return calculateAPotential(pointVector, usedPrecision);
}

vec3::Vector3 Coil::computeAPotentialVector(vec3::Vector3 pointVector) const
{
    return computeAPotentialVector(pointVector, defaultPrecisionCPU);
}

vec3::Vector3 Coil::computeBFieldVector(vec3::Vector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return calculateBField(pointVector, usedPrecision);
}

vec3::Vector3 Coil::computeBFieldVector(vec3::Vector3 pointVector) const
{
    return computeBFieldVector(pointVector, defaultPrecisionCPU);
}

vec3::Vector3 Coil::computeEFieldVector(vec3::Vector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    vec3::Vector3 computedVector = computeAPotentialVector(pointVector, usedPrecision);
    computedVector *= 2*M_PI * sineFrequency;
    return computedVector;
}

vec3::Vector3 Coil::computeEFieldVector(vec3::Vector3 pointVector) const
{
    return computeEFieldVector(pointVector, defaultPrecisionCPU);
}

vec3::Matrix3 Coil::computeBGradientMatrix(vec3::Vector3 pointVector, const PrecisionArguments &usedPrecision) const
{
    return calculateBGradient(pointVector, usedPrecision);
}

vec3::Matrix3 Coil::computeBGradientMatrix(vec3::Vector3 pointVector) const
{
    return computeBGradientMatrix(pointVector, defaultPrecisionCPU);
}


vec3::Vector3Array Coil::computeAllAPotentialVectors(const vec3::Vector3Array &pointVectors,
                                                     const PrecisionArguments &usedPrecision,
                                                     ComputeMethod computeMethod) const
{
    if (computeMethod == CPU_MT)
        return calculateAllAPotentialMT(pointVectors, usedPrecision);

    else if (computeMethod == GPU)
        return calculateAllAPotentialGPU(pointVectors, usedPrecision);
    else
    {
        vec3::Vector3Array computedPotentialArr;
        computedPotentialArr.reserve(pointVectors.size());

        for (int i = 0; i < pointVectors.size(); ++i)
            computedPotentialArr += computeAPotentialVector(pointVectors[i], usedPrecision);

        return computedPotentialArr;
    }
}

vec3::Vector3Array Coil::computeAllAPotentialVectors(const vec3::Vector3Array &pointVectors,
                                                     ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllAPotentialVectors(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllAPotentialVectors(pointVectors, defaultPrecisionCPU, computeMethod);
}


vec3::Vector3Array Coil::computeAllBFieldVectors(const vec3::Vector3Array &pointVectors,
                                                 const PrecisionArguments &usedPrecision,
                                                 ComputeMethod computeMethod) const
{
    if (computeMethod == CPU_MT)
        return calculateAllBFieldMT(pointVectors, usedPrecision);

    else if (computeMethod == GPU)
        return calculateAllBFieldGPU(pointVectors, usedPrecision);
    else
    {
        vec3::Vector3Array computedFieldArr = vec3::Vector3Array();
        computedFieldArr.reserve(pointVectors.size());

        for (int i = 0; i < pointVectors.size(); ++i)
            computedFieldArr += computeBFieldVector(pointVectors[i], usedPrecision);

        return computedFieldArr;
    }
}

vec3::Vector3Array Coil::computeAllBFieldVectors(const vec3::Vector3Array &pointVectors,
                                                 ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllBFieldVectors(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllBFieldVectors(pointVectors, defaultPrecisionCPU, computeMethod);
}


vec3::Vector3Array Coil::computeAllEFieldVectors(const vec3::Vector3Array &pointVectors,
                                                 const PrecisionArguments &usedPrecision,
                                                 ComputeMethod computeMethod) const
{
    vec3::Vector3Array output = computeAllAPotentialVectors(pointVectors, usedPrecision, computeMethod);
    double frequencyFactor = 2 * M_PI * sineFrequency;

    for (int i = 0; i < output.size(); ++i)
        output[i] *= frequencyFactor;

    return output;
}

vec3::Vector3Array Coil::computeAllEFieldVectors(const vec3::Vector3Array &pointVectors,
                                                 ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllEFieldVectors(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllEFieldVectors(pointVectors, defaultPrecisionCPU, computeMethod);
}


vec3::Matrix3Array Coil::computeAllBGradientMatrices(const vec3::Vector3Array &pointVectors,
                                                     const PrecisionArguments &usedPrecision,
                                                     ComputeMethod computeMethod) const
{
    if (computeMethod == CPU_MT)
        return calculateAllBGradientMT(pointVectors, usedPrecision);

    else if (computeMethod == GPU)
        return calculateAllBGradientGPU(pointVectors, usedPrecision);
    else
    {
        vec3::Matrix3Array computedGradientArr;
        computedGradientArr.reserve(pointVectors.size());

        for (int i = 0; i < pointVectors.size(); ++i)
            computedGradientArr += computeBGradientMatrix(pointVectors[i], usedPrecision);

        return computedGradientArr;
    }
}

vec3::Matrix3Array Coil::computeAllBGradientMatrices(const vec3::Vector3Array &pointVectors,
                                                     ComputeMethod computeMethod) const
{
    if (computeMethod == GPU)
        return computeAllBGradientMatrices(pointVectors, defaultPrecisionGPU, computeMethod);
    else
        return computeAllBGradientMatrices(pointVectors, defaultPrecisionCPU, computeMethod);
}


double Coil::computeMutualInductance(const Coil &primary, const Coil &secondary,
                                     CoilPairArguments inductanceArguments, ComputeMethod computeMethod)
{
    if (isZAxisCase(primary, secondary))
    {
        vec3::Vector3 secPositionVec = secondary.getPositionVector();

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


double Coil::computeAndSetSelfInductance(PrecisionFactor precisionFactor)
{
    if (coilType == CoilType::FILAMENT)
        throw std::logic_error("Coil loop self inductance is not defined, integral is divergent");

    // centering the coil at (0, 0, 0) and setting angles to 0 improves the accuracy by leveraging the z-axis formula
    vec3::Vector3 tempPosition = getPositionVector();
    std::pair tempAngles = getRotationAngles();
    setPositionAndOrientation();

    auto arguments = CoilPairArguments::getAppropriateCoilPairArguments(*this, *this, precisionFactor);
    double inductance;

    if (coilType == CoilType::FLAT)
        inductance = Coil::computeMutualInductance(*this, *this, arguments);
    else
        inductance = Coil::calculateSelfInductance(arguments);

    setSelfInductance(inductance);
    setPositionAndOrientation(tempPosition, tempAngles.first, tempAngles.second);

    return inductance;
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

            zForce = calculateAmpereForceZAxisFast(primary, secondary, secPositionVec.z, forceArguments, computeMethod);
        else
            zForce = calculateAmpereForceZAxisSlow(primary, secondary, secPositionVec.z, forceArguments, computeMethod);

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


std::vector<double> Coil::computeAllMutualInductanceArrangements(Coil primary, Coil secondary,
                                                                 const vec3::Vector3Array &primaryPositions,
                                                                 const vec3::Vector3Array &secondaryPositions,
                                                                 const std::vector<double> &primaryYAngles,
                                                                 const std::vector<double> &primaryZAngles,
                                                                 const std::vector<double> &secondaryYAngles,
                                                                 const std::vector<double> &secondaryZAngles,
                                                                 PrecisionFactor precisionFactor,
                                                                 ComputeMethod computeMethod)
{
    size_t numArrangements = primaryPositions.size();

    if (numArrangements == secondaryPositions.size() &&
        numArrangements == primaryYAngles.size() &&
        numArrangements == primaryZAngles.size() &&
        numArrangements == secondaryYAngles.size() &&
        numArrangements == secondaryZAngles.size())
    {
        if (computeMethod == GPU) {
            return calculateAllMutualInductanceArrangementsGPU(primary, secondary, primaryPositions, secondaryPositions,
                                                               primaryYAngles, primaryZAngles, secondaryYAngles, secondaryZAngles,
                                                               precisionFactor);
        }
        else if (numArrangements >= 2 * primary.getThreadCount() && computeMethod == CPU_MT)
        {
            return calculateAllMutualInductanceArrangementsMTD(primary, secondary, primaryPositions, secondaryPositions,
                                                               primaryYAngles, primaryZAngles, secondaryYAngles, secondaryZAngles,
                                                               precisionFactor);
        } else
        {
            std::vector<double> outputMInductances;
            outputMInductances.reserve(numArrangements);

            for (int i = 0; i < numArrangements; ++i) {
                primary.setPositionAndOrientation(primaryPositions[i], primaryYAngles[i], primaryZAngles[i]);
                secondary.setPositionAndOrientation(secondaryPositions[i], secondaryYAngles[i], secondaryZAngles[i]);

                outputMInductances.emplace_back(Coil::computeMutualInductance(primary, secondary, precisionFactor, computeMethod));
            }
            return outputMInductances;
        }
    }
    else
        throw std::logic_error("Array sizes do not match!");
}

std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
Coil::computeAllAmpereForceArrangements(Coil primary, Coil secondary,
                                        const vec3::Vector3Array &primaryPositions,
                                        const vec3::Vector3Array &secondaryPositions,
                                        const std::vector<double> &primaryYAngles, const std::vector<double> &primaryZAngles,
                                        const std::vector<double> &secondaryYAngles, const std::vector<double> &secondaryZAngles,
                                        PrecisionFactor precisionFactor, ComputeMethod computeMethod)
{
    size_t numArrangements = primaryPositions.size();

    if (numArrangements == secondaryPositions.size() &&
        numArrangements == primaryYAngles.size() &&
        numArrangements == primaryZAngles.size() &&
        numArrangements == secondaryYAngles.size() &&
        numArrangements == secondaryZAngles.size())
    {
        if (computeMethod == GPU)
        {
            return calculateAllAmpereForceArrangementsGPU(primary, secondary, primaryPositions, secondaryPositions,
                                                          primaryYAngles, primaryZAngles, secondaryYAngles, secondaryZAngles,
                                                          precisionFactor);
        }
        else if (numArrangements >= 2 * primary.getThreadCount() && computeMethod == CPU_MT)
        {
            return calculateAllAmpereForceArrangementsMTD(primary, secondary, primaryPositions, secondaryPositions,
                                                          primaryYAngles, primaryZAngles, secondaryYAngles, secondaryZAngles,
                                                          precisionFactor);
        } else
        {
            std::vector<std::pair<vec3::Vector3, vec3::Vector3>> outputForcesAndTorques;

            for (int i = 0; i < numArrangements; ++i) {
                primary.setPositionAndOrientation(primaryPositions[i], primaryYAngles[i], primaryZAngles[i]);
                secondary.setPositionAndOrientation(secondaryPositions[i], secondaryYAngles[i], secondaryZAngles[i]);

                outputForcesAndTorques.emplace_back(Coil::computeAmpereForce(primary, secondary, precisionFactor, computeMethod));
            }
            return outputForcesAndTorques;
        }
    }
    else
        throw std::logic_error("Array sizes do not match");
}


Coil::operator std::string() const
{
    std::stringstream output;

    output << "Coil("
        << "id=" << id
        << ", inner_radius=" << innerRadius
        << ", thickness=" << thickness
        << ", length=" << length
        << ", num_of_turns=" << numOfTurns
        << ", current_density=" << currentDensity
        << ", current=" << current
        << ", wire_resistivity=" << wireResistivity
        << ", is_sine_driven=" << sineDriven
        << ", sine_frequency=" << sineFrequency
        << ", magnetic_moment=" << magneticMoment
        << ", average_wire_thickness=" << averageWireThickness
        << ", resistance=" << resistance
        << ", self_inductance=" << selfInductance
        << ", reactance=" << reactance
        << ", impedance=" << impedance
        << ", use_fast_method=" << useFastMethod
        << ", thread_count=" << threadCount
        << ", default_precision_CPU=" << std::string(defaultPrecisionCPU)
        << ", default_precision_GPU=" << std::string(defaultPrecisionGPU)
        << ", position_vector=" << std::string(positionVector)
        << ", y_axis_angle=" << yAxisAngle
        << ", z_axis_angle=" << zAxisAngle
        << ")";

    return output.str();
}
