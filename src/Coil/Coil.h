#ifndef GENERAL_COIL_PROGRAM_COIL_H
#define GENERAL_COIL_PROGRAM_COIL_H

#include "ComputeMethod.h"
#include "CoilType.h"
#include "Tensor/Tensor.h"
#include "PrecisionGlobalVars.h"

#include <vector>


#define PRINT_ENABLED 0

const int precisionArraySize = 500;
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
                                                             ComputeMethod method = CPU_ST, bool isGeneral = true);

    private:
        static void getGeometryCaseAndIncrementsSingleCoil(const Coil &coil, int &caseIndex, int &totalIncrements);

        static void getGeometryCaseAndIncrementsCoilPair(const Coil &primary, const Coil &secondary,
                                                         PrecisionFactor precisionFactor,
                                                         int &caseIndex, int &totalIncrements);

        static CoilPairArguments calculateCoilPairArgumentsZAxisCPU(const Coil &primary, const Coil &secondary,
                                                                    PrecisionFactor precisionFactor);

        static CoilPairArguments calculateCoilPairArgumentsGeneralCPU(const Coil &primary, const Coil &secondary,
                                                                      PrecisionFactor precisionFactor);

        static CoilPairArguments calculateCoilPairArgumentsZAxisGPU(const Coil &primary, const Coil &secondary,
                                                                    PrecisionFactor precisionFactor);

        static CoilPairArguments calculateCoilPairArgumentsGeneralGPU(const Coil &primary, const Coil &secondary,
                                                                      PrecisionFactor precisionFactor);
};

class Coil
{
    private:

        const double innerRadius;
        const double thickness;
        const double length;
        const int numOfTurns;

        double currentDensity{};
        double current{};

        double wireResistivity{};
        bool isSineDriven{};
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

        PrecisionArguments defaultPrecision;

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
        Coil(double innerRadius, double thickness, double length, int numOfTurns,
             double current, double wireResistivity, double sineFrequency,
             const PrecisionArguments &precisionSettings, int threadCount = defaultThreadCount,
             vec3::CoordVector3 coordinatePosition = vec3::CoordVector3(), double yAxisAngle = 0.0, double zAxisAngle = 0.0);

        Coil(double innerRadius, double thickness, double length, int numOfTurns, double current, double sineFrequency,
             PrecisionFactor precisionFactor = PrecisionFactor(), int threadCount = defaultThreadCount);
        Coil(double innerRadius, double thickness, double length, int numOfTurns, double current, double sineFrequency,
             const PrecisionArguments &precisionSettings, int threadCount = 1);

        Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
             PrecisionFactor precisionFactor = PrecisionFactor(), int threadCount = defaultThreadCount);
        Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
             const PrecisionArguments &precisionSettings, int threadCount = 1);

        Coil(double innerRadius, double thickness, double length, int numOfTurns,
             PrecisionFactor precisionFactor = PrecisionFactor(), int threadCount = defaultThreadCount);
        Coil(double innerRadius, double thickness, double length, int numOfTurns,
             const PrecisionArguments &precisionSettings, int threadCount = 1);


        [[nodiscard]] double getInnerRadius() const;
        [[nodiscard]] double getThickness() const;
        [[nodiscard]] double getLength() const;
        [[nodiscard]] int getNumOfTurns() const;

        [[nodiscard]] double getCurrentDensity() const;
        [[nodiscard]] double getCurrent() const;

        [[nodiscard]] double getWireResistivity() const;
        [[nodiscard]] bool isSineDriven1() const;
        [[nodiscard]] double getSineFrequency() const;

        [[nodiscard]] double getMagneticMoment();
        [[nodiscard]] double getAverageWireThickness() const;

        [[nodiscard]] double getSelfInductance() const;
        [[nodiscard]] double getResistance();
        [[nodiscard]] double getReactance();
        [[nodiscard]] double getImpedance();

        [[nodiscard]] const PrecisionArguments &getPrecisionSettings() const;
        [[nodiscard]] int getThreadCount() const;
        [[nodiscard]] bool isUsingFastMethod() const;
        [[nodiscard]] CoilType getCoilType() const;

        [[nodiscard]] vec3::CoordVector3 getPositionVector() const;
        [[nodiscard]] std::pair<double, double> getRotationAngles() const;

        void setCurrentDensity(double currentDensity);
        void setCurrent(double current);
        void setWireResistivity(double wireResistivity);
        void setSineFrequency(double sineFrequency);
        void setPrecisionSettings(const PrecisionArguments &precisionSettings);
        void setThreadCount(int threadCount);

        void setSelfInductance(double selfInductance);

        void setPositionAndOrientation(vec3::CoordVector3 positionVector = vec3::CoordVector3(),
                                       double yAxisAngle = 0.0, double zAxisAngle = 0.0);

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


        [[nodiscard]] std::vector<double> computeAllBFieldX(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                              ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllBFieldX(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                            const PrecisionArguments &usedPrecision, ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllBFieldY(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                            ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllBFieldY(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                            const PrecisionArguments &usedPrecision, ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllBFieldZ(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                              ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllBFieldZ(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                            const PrecisionArguments &usedPrecision, ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllBFieldAbs(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                              ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllBFieldAbs(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                              const PrecisionArguments &usedPrecision, ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<vec3::FieldVector3> computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                   ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<vec3::FieldVector3> computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                   const PrecisionArguments &usedPrecision,
                                                                   ComputeMethod method = CPU_ST) const;


        [[nodiscard]] std::vector<double> computeAllAPotentialX(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllAPotentialX(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                const PrecisionArguments &usedPrecision, ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllAPotentialY(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllAPotentialY(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                const PrecisionArguments &usedPrecision, ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllAPotentialZ(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllAPotentialZ(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                const PrecisionArguments &usedPrecision, ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllAPotentialAbs(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                  ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllAPotentialAbs(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                  const PrecisionArguments &usedPrecision, ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<vec3::FieldVector3> computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                                     ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<vec3::FieldVector3> computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                                     const PrecisionArguments &usedPrecision,
                                                                                     ComputeMethod method = CPU_ST) const;


        [[nodiscard]] std::vector<double> computeAllEFieldX(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                            ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllEFieldX(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                            const PrecisionArguments &usedPrecision, ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllEFieldY(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                            ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllEFieldY(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                            const PrecisionArguments &usedPrecision, ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllEFieldZ(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                            ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllEFieldZ(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                            const PrecisionArguments &usedPrecision, ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<double> computeAllEFieldAbs(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                              ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<double> computeAllEFieldAbs(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                              const PrecisionArguments &usedPrecision, ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<vec3::FieldVector3> computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                                 ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<vec3::FieldVector3> computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                                 const PrecisionArguments &usedPrecision,
                                                                                 ComputeMethod method = CPU_ST) const;


        [[nodiscard]] std::vector<vec3::Matrix3> computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                            ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<vec3::Matrix3> computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                                                            const PrecisionArguments &usedPrecision,
                                                                            ComputeMethod method = CPU_ST) const;


        static double computeMutualInductance(const Coil &primary, const Coil &secondary,
                                              PrecisionFactor precisionFactor = PrecisionFactor(),
                                              ComputeMethod method = CPU_ST);
        static double computeMutualInductance(const Coil &primary, const Coil &secondary,
                                              CoilPairArguments inductanceArguments, ComputeMethod method = CPU_ST);

        [[nodiscard]] double computeSecondaryInducedVoltage(const Coil &secondary, PrecisionFactor precisionFactor = PrecisionFactor(),
                                                            ComputeMethod method = CPU_ST) const;
        [[nodiscard]] double computeSecondaryInducedVoltage(const Coil &secondary, CoilPairArguments inductanceArguments,
                                                            ComputeMethod method = CPU_ST) const;

        double computeAndSetSelfInductance(PrecisionFactor precisionFactor, ComputeMethod method = CPU_ST);


        static std::pair<vec3::FieldVector3, vec3::FieldVector3>
        computeAmpereForce(const Coil &primary, const Coil &secondary,
                           PrecisionFactor precisionFactor = PrecisionFactor(), ComputeMethod method = CPU_ST);

        static std::pair<vec3::FieldVector3, vec3::FieldVector3>
        computeAmpereForce(const Coil &primary, const Coil &secondary,
                           CoilPairArguments forceArguments, ComputeMethod method = CPU_ST);

        std::pair<vec3::FieldVector3, vec3::FieldVector3>
        computeForceOnDipoleMoment(vec3::CoordVector3 pointVector, vec3::FieldVector3 dipoleMoment) const;

        std::pair<vec3::FieldVector3, vec3::FieldVector3>
        computeForceOnDipoleMoment(vec3::CoordVector3 pointVector, vec3::FieldVector3 dipoleMoment,
                                   const PrecisionArguments &usedPrecision) const;

    private:
        void calculateMagneticMoment();
        void calculateAverageWireThickness();
        void calculateResistance();
        void calculateReactance();
        void calculateImpedance();
        void calculateCoilType();
        void calculateTransformationMatrices();

        static void precomputeCosPhi(int numAngularBlocks, int numAngularIncrements, double *outputArray);

        [[nodiscard]] std::pair<double, double> calculateBField(double zAxis, double rPolar,
                                                                const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double calculateAPotential(double zAxis, double rPolar,
                                                 const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::vector<double> calculateBGradient(double zAxis, double rPolar,
                                                             const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::pair<double, double> calculateBFieldSlow(double zAxis, double rPolar,
                                                                    const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double calculateAPotentialSlow(double zAxis, double rPolar,
                                                     const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::vector<double> calculateBGradientSlow(double zAxis, double rPolar,
                                                                 const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::pair<double, double> calculateBFieldFast(double zAxis, double rPolar,
                                                                    const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double calculateAPotentialFast(double zAxis, double rPolar,
                                                     const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::vector<double> calculateBGradientFast(double zAxis, double rPolar,
                                                                 const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] vec3::CoordVector3 adaptInputVectorForPoint(const vec3::CoordVector3 &pointVector) const;
        [[nodiscard]] vec3::FieldVector3 adaptOutputVectorForPoint(const vec3::FieldVector3 &computedVector) const;

        void adaptInputVectorsForAllPoints(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                           std::vector<double> &cylindricalZArr,
                                           std::vector<double> &cylindricalRArr,
                                           std::vector<double> &cylindricalPhiArr) const;

        [[nodiscard]] std::vector<vec3::FieldVector3>
        adaptOutputVectorsForAllPoints(const std::vector<vec3::FieldVector3> &computedVectorArr) const;

        void calculateAllBFieldST(const std::vector<double> &cylindricalZArr,
                                  const std::vector<double> &cylindricalRArr,
                                  std::vector<double> &computedFieldHArr,
                                  std::vector<double> &computedFieldZArr,
                                  const PrecisionArguments &usedPrecision) const;

        void calculateAllAPotentialST(const std::vector<double> &cylindricalZArr,
                                      const std::vector<double> &cylindricalRArr,
                                      std::vector<double> &computedPotentialArr,
                                      const PrecisionArguments &usedPrecision) const;

        void calculateAllBGradientST(const std::vector<double> &cylindricalZArr,
                                     const std::vector<double> &cylindricalRArr,
                                     std::vector<double> &computedGradientRPhi,
                                     std::vector<double> &computedGradientRR,
                                     std::vector<double> &computedGradientRZ,
                                     std::vector<double> &computedGradientZZ,
                                     const PrecisionArguments &usedPrecision) const;

        void calculateAllBFieldMT(const std::vector<double> &cylindricalZArr,
                                  const std::vector<double> &cylindricalRArr,
                                  std::vector<double> &computedFieldHArr,
                                  std::vector<double> &computedFieldZArr,
                                  const PrecisionArguments &usedPrecision,
                                  int chunkSize = g_defaultChunkSize, bool async = false) const;

        void calculateAllAPotentialMT(const std::vector<double> &cylindricalZArr,
                                      const std::vector<double> &cylindricalRArr,
                                      std::vector<double> &computedPotentialArr,
                                      const PrecisionArguments &usedPrecision,
                                      int chunkSize = g_defaultChunkSize, bool async = false) const;

        void calculateAllBGradientMT(const std::vector<double> &cylindricalZArr,
                                     const std::vector<double> &cylindricalRArr,
                                     std::vector<double> &computedGradientRPhi,
                                     std::vector<double> &computedGradientRR,
                                     std::vector<double> &computedGradientRZ,
                                     std::vector<double> &computedGradientZZ,
                                     const PrecisionArguments &usedPrecision,
                                     int chunkSize = g_defaultChunkSize, bool async = false) const;

        void calculateAllBFieldGPU(const std::vector<double> &cylindricalZArr,
                                   const std::vector<double> &cylindricalRArr,
                                   std::vector<double> &computedFieldHArr,
                                   std::vector<double> &computedFieldZArr,
                                   const PrecisionArguments &usedPrecision) const;

        void calculateAllAPotentialGPU(const std::vector<double> &cylindricalZArr,
                                       const std::vector<double> &cylindricalRArr,
                                       std::vector<double> &computedPotentialArr,
                                       const PrecisionArguments &usedPrecision) const;

        void calculateAllBGradientGPU(const std::vector<double> &cylindricalZArr,
                                      const std::vector<double> &cylindricalRArr,
                                      std::vector<double> &computedGradientRPhi,
                                      std::vector<double> &computedGradientRR,
                                      std::vector<double> &computedGradientRZ,
                                      std::vector<double> &computedGradientZZ,
                                      const PrecisionArguments &usedPrecision) const;

        void calculateAllBFieldSwitch(const std::vector<double> &cylindricalZArr,
                                      const std::vector<double> &cylindricalRArr,
                                      std::vector<double> &computedFieldHArr,
                                      std::vector<double> &computedFieldZArr,
                                      const PrecisionArguments &usedPrecision,
                                      ComputeMethod method) const;

        void calculateAllAPotentialSwitch(const std::vector<double> &cylindricalZArr,
                                          const std::vector<double> &cylindricalRArr,
                                          std::vector<double> &computedPotentialArr,
                                          const PrecisionArguments &usedPrecision,
                                          ComputeMethod method) const;

        void calculateAllBGradientSwitch(const std::vector<double> &cylindricalZArr,
                                         const std::vector<double> &cylindricalRArr,
                                         std::vector<double> &computedGradientRPhi,
                                         std::vector<double> &computedGradientRR,
                                         std::vector<double> &computedGradientRZ,
                                         std::vector<double> &computedGradientZZ,
                                         const PrecisionArguments &usedPrecision,
                                         ComputeMethod method) const;

        static std::vector<std::pair<vec3::FieldVector3, vec3::FieldVector3>>
        calculateRingIncrementPosition(int angularBlocks, int angularIncrements,
                                       double alpha, double beta, double ringIntervalSize);

        static bool isZAxisCase(const Coil &primary, const Coil &secondary);

        static double calculateMutualInductanceZAxis(const Coil &primary, const Coil &secondary, double zDisplacement,
                                                     CoilPairArguments inductanceArguments,
                                                     ComputeMethod method = CPU_ST);

        static double calculateMutualInductanceGeneral(const Coil &primary, const Coil &secondary,
                                                       CoilPairArguments inductanceArguments, ComputeMethod method = CPU_ST);

        static double calculateAmpereForceZAxis(const Coil &primary, const Coil &secondary, double zDisplacement,
                                                CoilPairArguments forceArguments,
                                                ComputeMethod method = CPU_ST);

        static std::pair<vec3::FieldVector3, vec3::FieldVector3>
        calculateAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                    CoilPairArguments forceArguments, ComputeMethod method);

        void synchronizeThreads() const;
};

#endif //GENERAL_COIL_PROGRAM_COIL_H
