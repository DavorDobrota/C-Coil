#include "Coil.h"
#include "ThreadPool.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <sstream>


namespace
{
    const double g_defaultCurrent = 1.0;
    const double g_defaultResistivity = 1.63e-8;
    const double g_defaultSineFrequency = 50;

    unsigned long long g_currentId = 0;

    threadPool::ThreadPoolControl g_threadPool;
}


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
                g_defaultSineFrequency, precisionFactor, threadCount,
                coordinatePosition, yAxisAngle, zAxisAngle) {}

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

bool Coil::isPointInside(vec3::Vector3 pointVector)
{
    vec3::Vector3 transformedVec = this->inverseTransformationMatrix * (pointVector - this->positionVector);
    vec3::Triplet cylindricalCoords = transformedVec.getAsCylindricalCoords();

    double zPos = cylindricalCoords.first;
    double rPos = cylindricalCoords.second;

    return (std::abs(zPos) < 0.5*length + g_thinCoilApproximationRatio * innerRadius) &&
           (std::abs(rPos - innerRadius - 0.5*thickness) < 0.5*thickness + g_thinCoilApproximationRatio * innerRadius);
}

void Coil::setSelfInductance(double selfInductance)
{
    this->selfInductance = selfInductance;
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
