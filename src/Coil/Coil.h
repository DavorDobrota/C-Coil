#ifndef GENERAL_COIL_PROGRAM_COIL_H
#define GENERAL_COIL_PROGRAM_COIL_H

#include "ComputeMethod.h"
#include "CoilType.h"
#include "Tensor.h"
#include "PrecisionGlobalVars.h"
#include "CoilData.h"

#include <vector>
#include <string>

#define PRINT_ENABLED 0


const int precisionArraySize = 864;
const int defaultThreadCount = 8;

const extern int blockPrecisionCPUArray[precisionArraySize];
const extern int incrementPrecisionCPUArray[precisionArraySize];

void generatePrecisionArrays();

class Coil;


struct PrecisionFactor
{
    PrecisionFactor();
    explicit PrecisionFactor(double relativePrecision);

    double relativePrecision;

    explicit operator std::string() const;
};

struct PrecisionArguments
{
    PrecisionArguments();
    explicit PrecisionArguments(int angularBlocks, int thicknessBlocks, int lengthBlocks,
                                int angularIncrements, int thicknessIncrements, int lengthIncrements);

    int angularBlockCount;
    int thicknessBlockCount;
    int lengthBlockCount;

    int angularIncrementCount;
    int thicknessIncrementCount;
    int lengthIncrementCount;

    static PrecisionArguments getCoilPrecisionArgumentsCPU(const Coil &coil, PrecisionFactor precisionFactor);

    static PrecisionArguments getCoilPrecisionArgumentsGPU(const Coil &coil, PrecisionFactor precisionFactor);

    explicit operator std::string() const;
};

struct CoilPairArguments
{
    CoilPairArguments();
    explicit CoilPairArguments(const PrecisionArguments &primaryPrecision,
                               const PrecisionArguments &secondaryPrecision);

    PrecisionArguments primaryPrecision;
    PrecisionArguments secondaryPrecision;

    static CoilPairArguments getAppropriateCoilPairArguments(const Coil &primary, const Coil &secondary,
                                                             PrecisionFactor precisionFactor,
                                                             ComputeMethod computeMethod = CPU_ST, bool zAxisCase = true);

    explicit operator std::string() const;

    private:
        static CoilPairArguments calculateCoilPairArgumentsCPU(const Coil &primary, const Coil &secondary,
                                                               PrecisionFactor precisionFactor, bool zAxisCase = false);

        static CoilPairArguments calculateCoilPairArgumentsGPU(const Coil &primary, const Coil &secondary,
                                                               PrecisionFactor precisionFactor, bool zAxisCase = false);
};

class Coil
{
    private:
        unsigned long long id;

        double innerRadius;
        double thickness;
        double length;
        int numOfTurns;

        double currentDensity{};
        double current{};

        double wireResistivity{};
        bool sineDriven{};
        double sineFrequency{};

        double magneticMoment{};
        double averageWireThickness{};

        double resistance{};
        double selfInductance{};
        double reactance{};
        double impedance{};

        CoilType coilType{};
        bool useFastMethod{};
        int threadCount{};

        PrecisionArguments defaultPrecisionCPU;
        PrecisionArguments defaultPrecisionGPU;

        vec3::CoordVector3 positionVector{};
        double yAxisAngle{};
        double zAxisAngle{};
        vec3::Matrix3 transformationMatrix{};
        vec3::Matrix3 inverseTransformationMatrix{};

    public:
        Coil();

        Coil(double innerRadius, double thickness, double length, int numOfTurns,
             double current, double wireResistivity, double sineFrequency,
             PrecisionFactor precisionFactor = PrecisionFactor(), int threadCount = defaultThreadCount,
             vec3::CoordVector3 coordinatePosition = vec3::CoordVector3(), double yAxisAngle = 0.0, double zAxisAngle = 0.0);
        Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
             double wireResistivity, double sineFrequency,
             const PrecisionArguments &precisionSettingsCPU, const PrecisionArguments &precisionSettingsGPU,
             int threadCount = defaultThreadCount, vec3::CoordVector3 coordinatePosition = vec3::CoordVector3(),
             double yAxisAngle = 0.0, double zAxisAngle = 0.0);
10
        Coil(double innerRadius, double thickness, double length, int numOfTurns, double current, double sineFrequency,
             PrecisionFactor precisionFactor = PrecisionFactor(), int threadCount = defaultThreadCount,
             vec3::CoordVector3 coordinatePosition = vec3::CoordVector3(), double yAxisAngle = 0.0, double zAxisAngle = 0.0);
        Coil(double innerRadius, double thickness, double length, int numOfTurns, double current, double sineFrequency,
             const PrecisionArguments &precisionSettingsCPU, const PrecisionArguments &precisionSettingsGPU,
             int threadCount = defaultThreadCount, vec3::CoordVector3 coordinatePosition = vec3::CoordVector3(),
             double yAxisAngle = 0.0, double zAxisAngle = 0.0);

        Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
             PrecisionFactor precisionFactor = PrecisionFactor(), int threadCount = defaultThreadCount,
             vec3::CoordVector3 coordinatePosition = vec3::CoordVector3(), double yAxisAngle = 0.0, double zAxisAngle = 0.0);
        Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
             const PrecisionArguments &precisionSettingsCPU, const PrecisionArguments &precisionSettingsGPU,
             int threadCount = defaultThreadCount, vec3::CoordVector3 coordinatePosition = vec3::CoordVector3(),
             double yAxisAngle = 0.0, double zAxisAngle = 0.0);

        Coil(double innerRadius, double thickness, double length, int numOfTurns,
             PrecisionFactor precisionFactor = PrecisionFactor(), int threadCount = defaultThreadCount,
             vec3::CoordVector3 coordinatePosition = vec3::CoordVector3(), double yAxisAngle = 0.0, double zAxisAngle = 0.0);
        Coil(double innerRadius, double thickness, double length, int numOfTurns,
             const PrecisionArguments &precisionSettingsCPU, const PrecisionArguments &precisionSettingsGPU,
             int threadCount = defaultThreadCount, vec3::CoordVector3 coordinatePosition = vec3::CoordVector3(),
             double yAxisAngle = 0.0, double zAxisAngle = 0.0);

        [[nodiscard]] unsigned long long getId() const;
        [[nodiscard]] double getInnerRadius() const;
        [[nodiscard]] double getThickness() const;
        [[nodiscard]] double getLength() const;
        [[nodiscard]] int getNumOfTurns() const;

        [[nodiscard]] double getCurrentDensity() const;
        [[nodiscard]] double getCurrent() const;

        [[nodiscard]] double getWireResistivity() const;
        [[nodiscard]] bool isSineDriven() const;
        [[nodiscard]] double getSineFrequency() const;

        [[nodiscard]] vec3::FieldVector3 getMagneticMoment();
        [[nodiscard]] double getAverageWireThickness() const;

        [[nodiscard]] double getSelfInductance() const;
        [[nodiscard]] double getResistance();
        [[nodiscard]] double getReactance();
        [[nodiscard]] double getImpedance();

        [[nodiscard]] const PrecisionArguments &getPrecisionSettingsCPU() const;
        [[nodiscard]] const PrecisionArguments &getPrecisionSettingsGPU() const;
        [[nodiscard]] int getThreadCount() const;
        [[nodiscard]] bool isUsingFastMethod() const;
        [[nodiscard]] CoilType getCoilType() const;

        [[nodiscard]] vec3::CoordVector3 getPositionVector() const;
        [[nodiscard]] std::pair<double, double> getRotationAngles() const;

        [[nodiscard]] vec3::Matrix3 getTransformationMatrix() const;
        [[nodiscard]] vec3::Matrix3 getInverseTransformationMatrix() const;

        void setCurrentDensity(double currentDensity);
        void setCurrent(double current);
        void setWireResistivity(double wireResistivity);
        void setSineFrequency(double sineFrequency);

        void setDefaultPrecisionCPU(const PrecisionArguments &precisionSettings);
        void setDefaultPrecisionCPU(PrecisionFactor precisionFactor = PrecisionFactor());
        void setDefaultPrecisionGPU(const PrecisionArguments &precisionSettings);
        void setDefaultPrecisionGPU(PrecisionFactor precisionFactor = PrecisionFactor());
        void setDefaultPrecision(PrecisionFactor precisionFactor = PrecisionFactor());

        void setThreadCount(int threadCount);
        void setPositionAndOrientation(vec3::CoordVector3 positionVector = vec3::CoordVector3(),
                                       double yAxisAngle = 0.0, double zAxisAngle = 0.0);

        void setSelfInductance(double selfInductance);


        [[nodiscard]] double computeBFieldX(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] double computeBFieldX(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeBFieldY(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] double computeBFieldY(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeBFieldZ(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] double computeBFieldZ(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeBFieldAbs(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] double computeBFieldAbs(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] vec3::FieldVector3 computeBFieldVector(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] vec3::FieldVector3 computeBFieldVector(vec3::CoordVector3 pointVector,
                                                             const PrecisionArguments &usedPrecision) const;


        [[nodiscard]] double computeAPotentialX(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] double computeAPotentialX(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeAPotentialY(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] double computeAPotentialY(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeAPotentialZ(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] double computeAPotentialZ(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeAPotentialAbs(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] double computeAPotentialAbs(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] vec3::FieldVector3 computeAPotentialVector(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] vec3::FieldVector3 computeAPotentialVector(vec3::CoordVector3 pointVector,
                                                                 const PrecisionArguments &usedPrecision) const;


        [[nodiscard]] double computeEFieldX(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] double computeEFieldX(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeEFieldY(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] double computeEFieldY(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeEFieldZ(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] double computeEFieldZ(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeEFieldAbs(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] double computeEFieldAbs(vec3::CoordVector3 pointVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] vec3::FieldVector3 computeEFieldVector(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] vec3::FieldVector3 computeEFieldVector(vec3::CoordVector3 pointVector,
                                                             const PrecisionArguments &usedPrecision) const;


        [[nodiscard]] vec3::Matrix3 computeBGradientTensor(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] vec3::Matrix3 computeBGradientTensor(vec3::CoordVector3 pointVector,
                                                           const PrecisionArguments &usedPrecision) const;


        [[nodiscard]] std::vector<double> computeAllBFieldX(const std::vector<vec3::CoordVector3> &pointVectors,
                                                            ComputeMethod computeMethod = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllBFieldX(const std::vector<vec3::CoordVector3> &pointVectors,
                                                            const PrecisionArguments &usedPrecision, ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllBFieldY(const std::vector<vec3::CoordVector3> &pointVectors,
                                                            ComputeMethod computeMethod = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllBFieldY(const std::vector<vec3::CoordVector3> &pointVectors,
                                                            const PrecisionArguments &usedPrecision, ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllBFieldZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                                            ComputeMethod computeMethod = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllBFieldZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                                            const PrecisionArguments &usedPrecision, ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllBFieldAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                                              ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllBFieldAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                                              const PrecisionArguments &usedPrecision, ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] std::vector<vec3::FieldVector3> computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                   ComputeMethod computeMethod = CPU_ST) const;
        [[nodiscard]] std::vector<vec3::FieldVector3> computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                   const PrecisionArguments &usedPrecision,
                                                                   ComputeMethod computeMethod = CPU_ST) const;


        [[nodiscard]] std::vector<double> computeAllAPotentialX(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                ComputeMethod computeMethod = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllAPotentialX(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                const PrecisionArguments &usedPrecision, ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllAPotentialY(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                ComputeMethod computeMethod = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllAPotentialY(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                const PrecisionArguments &usedPrecision, ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllAPotentialZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                ComputeMethod computeMethod = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllAPotentialZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                const PrecisionArguments &usedPrecision, ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllAPotentialAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                  ComputeMethod computeMethod = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllAPotentialAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                  const PrecisionArguments &usedPrecision, ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] std::vector<vec3::FieldVector3> computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                                     ComputeMethod computeMethod = CPU_ST) const;
        [[nodiscard]] std::vector<vec3::FieldVector3> computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                                     const PrecisionArguments &usedPrecision,
                                                                                     ComputeMethod computeMethod = CPU_ST) const;


        [[nodiscard]] std::vector<double> computeAllEFieldX(const std::vector<vec3::CoordVector3> &pointVectors,
                                                            ComputeMethod computeMethod = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllEFieldX(const std::vector<vec3::CoordVector3> &pointVectors,
                                                            const PrecisionArguments &usedPrecision, ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllEFieldY(const std::vector<vec3::CoordVector3> &pointVectors,
                                                            ComputeMethod computeMethod = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllEFieldY(const std::vector<vec3::CoordVector3> &pointVectors,
                                                            const PrecisionArguments &usedPrecision, ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllEFieldZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                                            ComputeMethod computeMethod = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllEFieldZ(const std::vector<vec3::CoordVector3> &pointVectors,
                                                            const PrecisionArguments &usedPrecision, ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllEFieldAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                                              ComputeMethod computeMethod = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllEFieldAbs(const std::vector<vec3::CoordVector3> &pointVectors,
                                                              const PrecisionArguments &usedPrecision, ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] std::vector<vec3::FieldVector3> computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                                 ComputeMethod computeMethod = CPU_ST) const;
        [[nodiscard]] std::vector<vec3::FieldVector3> computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                                 const PrecisionArguments &usedPrecision,
                                                                                 ComputeMethod computeMethod = CPU_ST) const;


        [[nodiscard]] std::vector<vec3::Matrix3> computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                            ComputeMethod computeMethod = CPU_ST) const;
        [[nodiscard]] std::vector<vec3::Matrix3> computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                            const PrecisionArguments &usedPrecision,
                                                                            ComputeMethod computeMethod = CPU_ST) const;


        static double computeMutualInductance(const Coil &primary, const Coil &secondary,
                                              PrecisionFactor precisionFactor = PrecisionFactor(),
                                              ComputeMethod computeMethod = CPU_ST);
        static double computeMutualInductance(const Coil &primary, const Coil &secondary,
                                              CoilPairArguments inductanceArguments, ComputeMethod computeMethod = CPU_ST);

        [[nodiscard]] double computeSecondaryInducedVoltage(const Coil &secondary, PrecisionFactor precisionFactor = PrecisionFactor(),
                                                            ComputeMethod computeMethod = CPU_ST) const;
        [[nodiscard]] double computeSecondaryInducedVoltage(const Coil &secondary, CoilPairArguments inductanceArguments,
                                                            ComputeMethod computeMethod = CPU_ST) const;

        double computeAndSetSelfInductance(PrecisionFactor precisionFactor = PrecisionFactor(), ComputeMethod computeMethod = CPU_ST);


        static std::pair<vec3::FieldVector3, vec3::FieldVector3>
        computeAmpereForce(const Coil &primary, const Coil &secondary,
                           PrecisionFactor precisionFactor = PrecisionFactor(), ComputeMethod computeMethod = CPU_ST);
        static std::pair<vec3::FieldVector3, vec3::FieldVector3>
        computeAmpereForce(const Coil &primary, const Coil &secondary,
                           CoilPairArguments forceArguments, ComputeMethod computeMethod = CPU_ST);

        [[nodiscard]] std::pair<vec3::FieldVector3, vec3::FieldVector3>
        computeForceOnDipoleMoment(vec3::CoordVector3 pointVector, vec3::FieldVector3 dipoleMoment) const;
        [[nodiscard]] std::pair<vec3::FieldVector3, vec3::FieldVector3>
        computeForceOnDipoleMoment(vec3::CoordVector3 pointVector, vec3::FieldVector3 dipoleMoment,
                                   const PrecisionArguments &usedPrecision) const;

        static std::vector<double>
        computeAllMutualInductanceArrangements(Coil primary, Coil secondary,
                                               const std::vector<vec3::CoordVector3> &primaryPositions,
                                               const std::vector<vec3::CoordVector3> &secondaryPositions,
                                               const std::vector<double> &primaryYAngles,
                                               const std::vector<double> &primaryZAngles,
                                               const std::vector<double> &secondaryYAngles,
                                               const std::vector<double> &secondaryZAngles,
                                               PrecisionFactor precisionFactor = PrecisionFactor(),
                                               ComputeMethod computeMethod = CPU_ST);

        static std::vector<std::pair<vec3::FieldVector3, vec3::FieldVector3>>
        computeAllAmpereForceArrangements(Coil primary, Coil secondary,
                                          const std::vector<vec3::CoordVector3> &primaryPositions,
                                          const std::vector<vec3::CoordVector3> &secondaryPositions,
                                          const std::vector<double> &primaryYAngles,
                                          const std::vector<double> &primaryZAngles,
                                          const std::vector<double> &secondaryYAngles,
                                          const std::vector<double> &secondaryZAngles,
                                          PrecisionFactor precisionFactor = PrecisionFactor(),
                                          ComputeMethod computeMethod = CPU_ST);

        explicit operator std::string() const;

    private:
        void calculateMagneticMoment();
        void calculateAverageWireThickness();
        void calculateResistance();
        void calculateReactance();
        void calculateImpedance();
        void calculateCoilType();
        void calculateTransformationMatrices();

        [[nodiscard]] std::pair<double, double> calculateBField(double zAxis, double rPolar,
                                                                const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double calculateAPotential(double zAxis, double rPolar,
                                                 const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::vector<double> calculateBGradient(double zAxis, double rPolar,
                                                             const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::pair<double, double> calculateBFieldSlow(double zCoord, double rCoord,
                                                                    const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double calculateAPotentialSlow(double zCoord, double rCoord,
                                                     const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::vector<double> calculateBGradientSlow(double zCoord, double rCoord,
                                                                 const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::pair<double, double> calculateBFieldFast(double zAxis, double rPolar,
                                                                    const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double calculateAPotentialFast(double zAxis, double rPolar,
                                                     const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::vector<double> calculateBGradientFast(double zAxis, double rPolar,
                                                                 const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] vec3::CoordVector3 adaptInputVectorForPoint(const vec3::CoordVector3 &pointVector) const;
        [[nodiscard]] vec3::FieldVector3 adaptOutputVectorForPoint(const vec3::FieldVector3 &computedVector) const;

        [[nodiscard]] std::vector<vec3::FieldVector3> calculateAllBFieldMT(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                           const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::vector<vec3::FieldVector3> calculateAllAPotentialMT(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                               const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::vector<vec3::Matrix3> calculateAllBGradientMT(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                         const PrecisionArguments &usedPrecision) const;

        void generateCoilData(CoilData &coilData, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::vector<vec3::FieldVector3> calculateAllAPotentialGPU(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                                const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::vector<vec3::FieldVector3> calculateAllBFieldGPU(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                            const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::vector<vec3::Matrix3> calculateAllBGradientGPU(const std::vector<vec3::CoordVector3> &pointVectors,
                                                                          const PrecisionArguments &usedPrecision) const;

        static std::vector<std::pair<vec3::FieldVector3, vec3::FieldVector3>>
        calculateRingIncrementPosition(int angularBlocks, int angularIncrements,
                                       double alpha, double beta, double ringIntervalSize);

        static bool isZAxisCase(const Coil &primary, const Coil &secondary);

        static double calculateMutualInductanceZAxisSlow(const Coil &primary, const Coil &secondary, double zDisplacement,
                                                         CoilPairArguments inductanceArguments,
                                                         ComputeMethod computeMethod = CPU_ST);

        static double calculateMutualInductanceZAxisFast(const Coil &primary, const Coil &secondary, double zDisplacement,
                                                         CoilPairArguments inductanceArguments,
                                                         ComputeMethod computeMethod = CPU_ST);

        static double calculateMutualInductanceGeneral(const Coil &primary, const Coil &secondary,
                                                       CoilPairArguments inductanceArguments, ComputeMethod computeMethod = CPU_ST);

        static double calculateAmpereForceZAxisSlow(const Coil &primary, const Coil &secondary, double zDisplacement,
                                                    CoilPairArguments forceArguments,
                                                    ComputeMethod computeMethod = CPU_ST);

        static double calculateAmpereForceZAxisFast(const Coil &primary, const Coil &secondary, double zDisplacement,
                                                    CoilPairArguments forceArguments,
                                                    ComputeMethod computeMethod = CPU_ST);

        static std::pair<vec3::FieldVector3, vec3::FieldVector3>
        calculateAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                    CoilPairArguments forceArguments, ComputeMethod computeMethod);

        [[nodiscard]] double calculateSelfInductance(CoilPairArguments inductanceArguments, ComputeMethod computeMethod) const;

        [[nodiscard]] std::vector<int> calculateChunkSize(int numOps) const;
};

#endif //GENERAL_COIL_PROGRAM_COIL_H
