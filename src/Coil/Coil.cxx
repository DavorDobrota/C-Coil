#include "Coil.h"
#include "LegendreMatrix.h"
#include "CoilData.h"

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

    return vec3::Vector3(magneticMoment * sinA * cosB,
                              magneticMoment * sinA * sinB,
                              magneticMoment * cosA);
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


void Coil::setSelfInductance(double selfInductance)
{
    this->selfInductance = selfInductance;
}


void Coil::setPositionAndOrientation(vec3::Vector3 positionVector, double yAxisAngle, double zAxisAngle)
{
    this->positionVector = positionVector;
    this->yAxisAngle = yAxisAngle;
    this->zAxisAngle = zAxisAngle;

    calculateTransformationMatrices();
}

void Coil::calculateMagneticMoment()
{
    magneticMoment = M_PI * current * numOfTurns *
            (innerRadius * innerRadius + innerRadius * thickness + thickness * thickness / 3);
}

void Coil::calculateAverageWireThickness()
{
    averageWireThickness = std::sqrt(length * thickness / numOfTurns);
}

void Coil::calculateResistance()
{
    double wireRadius = averageWireThickness * 0.5;
    double ohmicResistance = wireResistivity * numOfTurns * 2*M_PI *
            (innerRadius + thickness * 0.5) / (wireRadius * wireRadius * M_PI);
    double skinDepth = std::sqrt(wireResistivity / (M_PI * sineFrequency * g_MiReduced));

    double ohmicSurface = M_PI * wireRadius * wireRadius;
    double sineSurface = 2*M_PI * (
            skinDepth * skinDepth * (exp(-wireRadius / skinDepth) - 1) +
            skinDepth * wireRadius);

    resistance = ohmicResistance * (ohmicSurface / sineSurface);
}

void Coil::calculateReactance()
{
    reactance = selfInductance * 2*M_PI * sineFrequency;
}

void Coil::calculateImpedance()
{
    calculateResistance();
    calculateReactance();
    impedance = std::sqrt(resistance * resistance + reactance * reactance);
}

void Coil::calculateCoilType()
{
    if (thickness / innerRadius < g_thinCoilApproximationRatio && length / innerRadius < g_thinCoilApproximationRatio)
        coilType = CoilType::FILAMENT;

    else if (thickness / innerRadius < g_thinCoilApproximationRatio)
        coilType = CoilType::THIN;

    else if (length / innerRadius < g_thinCoilApproximationRatio)
        coilType = CoilType::FLAT;

    else
        coilType = CoilType::RECTANGULAR;
}

void Coil::calculateTransformationMatrices()
{

    double cosY = std::cos(yAxisAngle); double sinY = std::sin(yAxisAngle);
    double cosZ = std::cos(zAxisAngle); double sinZ = std::sin(zAxisAngle);

    transformationMatrix = vec3::Matrix3(cosZ * cosZ * cosY - sinZ * sinZ, -sinZ * cosZ * cosY - sinZ * cosZ, cosZ * sinY,
                                         sinZ * cosZ * cosY + sinZ * cosZ, cosZ * cosZ - sinZ * sinZ * cosY, sinZ * sinY,
                                         -sinY * cosZ, sinY * sinZ, cosY);

    inverseTransformationMatrix = vec3::Matrix3(cosZ * cosZ * cosY - sinZ * sinZ, sinZ * cosZ * cosY + sinZ * cosZ, -cosZ * sinY,
                                                -sinZ * cosZ * cosY - sinZ * cosZ, cosZ * cosZ - sinZ * sinZ * cosY, sinZ * sinY,
                                                sinY * cosZ, sinY * sinZ, cosY);
}

double Coil::computeAndSetSelfInductance(PrecisionFactor precisionFactor, ComputeMethod computeMethod)
{
    if (coilType == CoilType::FILAMENT)
    {
        fprintf(stderr, "ERROR: The integral of a filament is divergent, try a thin rectangular coil\n");
        throw std::logic_error("Coil loop calculation not supported");
    }
    // centering the coil at (0, 0, 0) and setting angles to 0 improves the accuracy by leveraging the z-axis formula
    vec3::Vector3 tempPosition = getPositionVector();
    std::pair tempAngles = getRotationAngles();
    setPositionAndOrientation();

    auto arguments = CoilPairArguments::getAppropriateCoilPairArguments(*this, *this, precisionFactor);
    double inductance;

    if (coilType == CoilType::FLAT)
        inductance = Coil::computeMutualInductance(*this, *this, arguments, computeMethod);
    else
        inductance = Coil::calculateSelfInductance(arguments, computeMethod);

    setSelfInductance(inductance);
    setPositionAndOrientation(tempPosition, tempAngles.first, tempAngles.second);

    return inductance;
}


std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
Coil::calculateRingIncrementPosition(int angularBlocks, int angularIncrements, double alpha, double beta)
{
    int numElements = angularBlocks * angularIncrements;

    std::vector<std::pair<vec3::Vector3, vec3::Vector3>> unitRingVector;
    unitRingVector.reserve(numElements);

    vec3::Vector3 ringPosition;
    vec3::Vector3 ringTangent;

    double angularBlock = 2*M_PI / angularBlocks;

    // subtracting 1 because n-th order Gauss quadrature has (n + 1) positions which here represent increments
    angularBlocks--;
    angularIncrements--;

    double sinA = std::sin(alpha); double cosA = std::cos(alpha);
    double sinB = std::sin(beta); double cosB = std::cos(beta);

    for (int phiBlock = 0; phiBlock <= angularBlocks; ++phiBlock)
    {
        double blockPositionPhi = angularBlock * (phiBlock + 0.5);

        for (int phiIndex = 0; phiIndex <= angularIncrements; ++phiIndex)
        {
            double phi = blockPositionPhi +
                         (angularBlock * 0.5) * Legendre::positionMatrix[angularIncrements][phiIndex];

            double sinPhi = std::sin(phi); double cosPhi = std::cos(phi);

            ringPosition = vec3::Vector3(cosB * cosA * cosPhi - sinB * sinPhi,
                                         sinB * cosA * cosPhi + cosB * sinPhi,
                                         (-1) * sinA * cosPhi);

            ringTangent = vec3::Vector3((-1) * cosB * cosA * sinPhi - sinB * cosPhi,
                                        (-1) * sinB * cosA * sinPhi + cosB * cosPhi,
                                        sinA * sinPhi);

            unitRingVector.emplace_back(ringPosition, ringTangent);
        }
    }
    return unitRingVector;
}


bool Coil::isZAxisCase(const Coil &primary, const Coil &secondary)
{
    vec3::Vector3 primPositionVec = primary.getPositionVector();
    vec3::Vector3 secPositionVec = secondary.getPositionVector();

    if (std::abs(primPositionVec.x / primary.innerRadius) < g_zAxisApproximationRatio &&
        std::abs(primPositionVec.y / primary.innerRadius) < g_zAxisApproximationRatio &&
        std::abs(secPositionVec.x / primary.innerRadius) < g_zAxisApproximationRatio &&
        std::abs(secPositionVec.y / primary.innerRadius) < g_zAxisApproximationRatio &&
        std::abs(primary.yAxisAngle / (2 * M_PI)) < g_zAxisApproximationRatio &&
        std::abs(secondary.yAxisAngle / (2 * M_PI)) < g_zAxisApproximationRatio)
    {
        return true;
    }
    return false;
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"
void Coil::generateCoilData(CoilData &coilData, const PrecisionArguments &usedPrecision) const
{
    if (useFastMethod)
        coilData.constFactor = g_MiReduced * currentDensity * thickness * M_PI * 0.5;
    else
        coilData.constFactor = g_MiReduced * currentDensity * thickness * length * M_PI * 0.5;

    coilData.useFastMethod = useFastMethod;

    coilData.innerRadius = innerRadius;
    coilData.thickness = thickness;
    coilData.length = length;

    coilData.lengthIncrements = usedPrecision.lengthIncrementCount;
    coilData.thicknessIncrements = usedPrecision.thicknessIncrementCount;
    coilData.angularIncrements = usedPrecision.angularIncrementCount;

    for (int i = 0; i < coilData.angularIncrements; ++i)
    {
        double phiPosition = M_PI_2 * (1.0 + Legendre::positionMatrix[coilData.angularIncrements - 1][i]);

        coilData.cosPrecomputeArray[i] = std::cos(phiPosition);
        coilData.angularWeightArray[i] = Legendre::weightsMatrix[coilData.angularIncrements - 1][i];
    }

    for (int i = 0; i < coilData.thicknessIncrements; ++i)
    {
        coilData.thicknessPositionArray[i] = Legendre::positionMatrix[coilData.thicknessIncrements - 1][i];
        coilData.thicknessWeightArray[i] = Legendre::weightsMatrix[coilData.thicknessIncrements - 1][i];
    }

    coilData.positionVector[0] = positionVector.x;
    coilData.positionVector[1] = positionVector.y;
    coilData.positionVector[2] = positionVector.z;

    coilData.transformArray[0] = transformationMatrix.xx;
    coilData.transformArray[1] = transformationMatrix.xy;
    coilData.transformArray[2] = transformationMatrix.xz;
    coilData.transformArray[3] = transformationMatrix.yx;
    coilData.transformArray[4] = transformationMatrix.yy;
    coilData.transformArray[5] = transformationMatrix.yz;
    coilData.transformArray[6] = transformationMatrix.zx;
    coilData.transformArray[7] = transformationMatrix.zy;
    coilData.transformArray[8] = transformationMatrix.zz;

    coilData.invTransformArray[0] = inverseTransformationMatrix.xx;
    coilData.invTransformArray[1] = inverseTransformationMatrix.xy;
    coilData.invTransformArray[2] = inverseTransformationMatrix.xz;
    coilData.invTransformArray[3] = inverseTransformationMatrix.yx;
    coilData.invTransformArray[4] = inverseTransformationMatrix.yy;
    coilData.invTransformArray[5] = inverseTransformationMatrix.yz;
    coilData.invTransformArray[6] = inverseTransformationMatrix.zx;
    coilData.invTransformArray[7] = inverseTransformationMatrix.zy;
    coilData.invTransformArray[8] = inverseTransformationMatrix.zz;
}
#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-narrowing-conversions"
void Coil::generateCoilPairArgumentsData(const Coil &primary, const Coil &secondary,
                                         CoilPairArgumentsData &coilPairArgumentsData,
                                         const CoilPairArguments &inductanceArguments, bool forceCalculation)
{
    if (primary.useFastMethod)
        coilPairArgumentsData.constFactor = g_MiReduced * primary.currentDensity * primary.thickness * M_PI * 0.5;
    else
        coilPairArgumentsData.constFactor = g_MiReduced * M_PI * 0.5 *
                                            primary.currentDensity * primary.thickness * primary.length;

    coilPairArgumentsData.useFastMethod = primary.useFastMethod;

    if (!forceCalculation)
        coilPairArgumentsData.correctionFactor = 2*M_PI * secondary.numOfTurns / primary.current;
    else
        coilPairArgumentsData.correctionFactor = 2*M_PI * secondary.numOfTurns * secondary.current;

    coilPairArgumentsData.primInnerRadius = primary.innerRadius;
    coilPairArgumentsData.primThickness = primary.thickness;
    coilPairArgumentsData.primLength = primary.length;

    coilPairArgumentsData.primLengthIncrements = inductanceArguments.primaryPrecision.lengthIncrementCount;
    coilPairArgumentsData.primThicknessIncrements = inductanceArguments.primaryPrecision.thicknessIncrementCount;
    coilPairArgumentsData.primAngularIncrements = inductanceArguments.primaryPrecision.angularIncrementCount;

    coilPairArgumentsData.secInnerRadius = secondary.innerRadius;
    coilPairArgumentsData.secThickness = secondary.thickness;
    coilPairArgumentsData.secLength = secondary.length;

    coilPairArgumentsData.secLengthIncrements = inductanceArguments.secondaryPrecision.lengthIncrementCount;
    coilPairArgumentsData.secThicknessIncrements = inductanceArguments.secondaryPrecision.thicknessIncrementCount;
    coilPairArgumentsData.secAngularIncrements = inductanceArguments.secondaryPrecision.angularIncrementCount;

    for (int i = 0; i < inductanceArguments.primaryPrecision.angularIncrementCount; ++i)
    {
        double phiPosition =
                M_PI_2 * (1.0 + Legendre::positionMatrix[inductanceArguments.primaryPrecision.angularIncrementCount - 1][i]);

        coilPairArgumentsData.primCosPrecomputeArray[i] = std::cos(phiPosition);
        coilPairArgumentsData.primAngularWeightArray[i] =
                Legendre::weightsMatrix[inductanceArguments.primaryPrecision.angularIncrementCount - 1][i];
    }
    for (int i = 0; i < inductanceArguments.primaryPrecision.thicknessIncrementCount; ++i)
    {
        coilPairArgumentsData.primThicknessPositionArray[i] =
                Legendre::positionMatrix[inductanceArguments.primaryPrecision.thicknessIncrementCount - 1][i];
        coilPairArgumentsData.primThicknessWeightArray[i] =
                Legendre::weightsMatrix[inductanceArguments.primaryPrecision.thicknessIncrementCount - 1][i];
    }

    for (int i = 0; i < inductanceArguments.secondaryPrecision.angularIncrementCount; ++i)
    {
        coilPairArgumentsData.secAngularPositionArray[i] =
                Legendre::positionMatrix[inductanceArguments.secondaryPrecision.angularIncrementCount - 1][i];
        coilPairArgumentsData.secAngularWeightArray[i] =
                Legendre::weightsMatrix[inductanceArguments.secondaryPrecision.angularIncrementCount - 1][i];
    }
    for (int i = 0; i < inductanceArguments.secondaryPrecision.thicknessIncrementCount; ++i)
    {
        coilPairArgumentsData.secThicknessPositionArray[i] =
                Legendre::positionMatrix[inductanceArguments.secondaryPrecision.thicknessIncrementCount - 1][i];
        coilPairArgumentsData.secThicknessWeightArray[i] =
                Legendre::weightsMatrix[inductanceArguments.secondaryPrecision.thicknessIncrementCount - 1][i];
    }
    for (int i = 0; i < inductanceArguments.secondaryPrecision.lengthIncrementCount; ++i)
    {
        coilPairArgumentsData.secLengthPositionArray[i] =
                Legendre::positionMatrix[inductanceArguments.secondaryPrecision.lengthIncrementCount - 1][i];
        coilPairArgumentsData.secLengthWeightArray[i] =
                Legendre::weightsMatrix[inductanceArguments.secondaryPrecision.lengthIncrementCount - 1][i];
    }

}
#pragma clang diagnostic pop

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
