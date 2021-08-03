#ifndef GENERAL_COIL_PROGRAM_COIL_H
#define GENERAL_COIL_PROGRAM_COIL_H

#include <vector>

#include "ComputeMethod.h"
#include "Vector3/Vector3.h"

#define PRINT_ENABLED 0

const int precisionArraySize = 423;
const int defaultThreadCount = 4;

const extern int blockPrecisionCPUArray[precisionArraySize];
const extern int incrementPrecisionCPUArray[precisionArraySize];

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

    static CoilPairArguments getSelfInductanceArguments(const Coil &coil, PrecisionFactor precisionFactor);

    static CoilPairArguments getAppropriateCoilPairArguments(const Coil &primary, const Coil &secondary,
                                                             PrecisionFactor precisionFactor,
                                                             ComputeMethod method = CPU_ST, bool isGeneral = true);

    private:
        static void getGeometryCaseAndIncrementsSingleCoil(const Coil &coil, PrecisionFactor precisionFactor,
                                                           int &caseIndex, int &totalIncrements);

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

        int threadCount{};

        PrecisionArguments defaultPrecision;

    public:
        Coil();

        Coil(double innerRadius, double thickness, double length, int numOfTurns,
             double current, double wireResistivity, double sineFrequency,
             PrecisionFactor precisionFactor = PrecisionFactor(), int threadCount = defaultThreadCount);
        Coil(double innerRadius, double thickness, double length, int numOfTurns,
             double current, double wireResistivity, double sineFrequency,
             const PrecisionArguments &precisionSettings, int threadCount = defaultThreadCount);

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

        void setCurrentDensity(double currentDensity);
        void setCurrent(double current);
        void setWireResistivity(double wireResistivity);
        void setSineFrequency(double sineFrequency);
        void setPrecisionSettings(const PrecisionArguments &precisionSettings);

        void setThreadCount(int threadCount);

        void setSelfInductance(double selfInductance);

        [[nodiscard]] double computeBFieldX(vec3::Vector3 positionVector) const;
        [[nodiscard]] double computeBFieldX(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeBFieldY(vec3::Vector3 positionVector) const;
        [[nodiscard]] double computeBFieldY(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeBFieldH(vec3::Vector3 positionVector) const;
        [[nodiscard]] double computeBFieldH(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeBFieldZ(vec3::Vector3 positionVector) const;
        [[nodiscard]] double computeBFieldZ(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeBFieldAbs(vec3::Vector3 positionVector) const;
        [[nodiscard]] double computeBFieldAbs(vec3::Vector3 positionVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] vec3::Vector3 computeBFieldVector(vec3::Vector3 positionVector) const;
        [[nodiscard]] vec3::Vector3 computeBFieldVector(Vector3 positionVector, const PrecisionArguments &usedPrecision) const;


        [[nodiscard]] double computeAPotentialX(Vector3 positionVector) const;
        [[nodiscard]] double computeAPotentialX(Vector3 positionVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeAPotentialY(Vector3 positionVector) const;
        [[nodiscard]] double computeAPotentialY(Vector3 positionVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeAPotentialZ(Vector3 positionVector) const;
        [[nodiscard]] double computeAPotentialZ(Vector3 positionVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeAPotentialAbs(Vector3 positionVector) const;
        [[nodiscard]] double computeAPotentialAbs(Vector3 positionVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] Vector3 computeAPotentialVector(Vector3 positionVector) const;
        [[nodiscard]] Vector3 computeAPotentialVector(Vector3 positionVector, const PrecisionArguments &usedPrecision) const;


        [[nodiscard]] double computeEFieldX(Vector3 positionVector) const;
        [[nodiscard]] double computeEFieldX(Vector3 positionVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeEFieldY(Vector3 positionVector) const;
        [[nodiscard]] double computeEFieldY(Vector3 positionVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeEFieldZ(Vector3 positionVector) const;
        [[nodiscard]] double computeEFieldZ(Vector3 positionVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double computeEFieldAbs(Vector3 positionVector) const;
        [[nodiscard]] double computeEFieldAbs(Vector3 positionVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] Vector3 computeEFieldVector(Vector3 positionVector) const;
        [[nodiscard]] Vector3 computeEFieldVector(Vector3 positionVector, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::vector<double> computeBGradientTensor(Vector3 positionVector) const;
        [[nodiscard]] std::vector<double> computeBGradientTensor(Vector3 positionVector, const PrecisionArguments &usedPrecision) const;

        void computeAllBFieldX(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               ComputeMethod method = CPU_ST) const;
        void computeAllBFieldX(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               const PrecisionArguments &usedPrecision,
                               ComputeMethod method = CPU_ST) const;

        void computeAllBFieldY(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               ComputeMethod method = CPU_ST) const;
        void computeAllBFieldY(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               const PrecisionArguments &usedPrecision,
                               ComputeMethod method = CPU_ST) const;

        void computeAllBFieldH(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               std::vector<double> &computedFieldArr,
                               ComputeMethod method = CPU_ST) const;
        void computeAllBFieldH(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               std::vector<double> &computedFieldArr,
                               const PrecisionArguments &usedPrecision,
                               ComputeMethod method = CPU_ST) const;

        void computeAllBFieldZ(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               ComputeMethod method = CPU_ST) const;
        void computeAllBFieldZ(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               const PrecisionArguments &usedPrecision,
                               ComputeMethod method = CPU_ST) const;

        void computeAllBFieldAbs(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedFieldArr,
                                 ComputeMethod method = CPU_ST) const;

        void computeAllBFieldAbs(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 const std::vector<double> &cylindricalPhiArr,
                                 std::vector<double> &computedFieldArr,
                                 const PrecisionArguments &usedPrecision,
                                 ComputeMethod method = CPU_ST) const;

        void computeAllBFieldComponents(const std::vector<double> &cylindricalZArr,
                                        const std::vector<double> &cylindricalRArr,
                                        const std::vector<double> &cylindricalPhiArr,
                                        std::vector<double> &computedFieldXArr,
                                        std::vector<double> &computedFieldYArr,
                                        std::vector<double> &computedFieldZArr,
                                        ComputeMethod method = CPU_ST) const;
        void computeAllBFieldComponents(const std::vector<double> &cylindricalZArr,
                                        const std::vector<double> &cylindricalRArr,
                                        const std::vector<double> &cylindricalPhiArr,
                                        std::vector<double> &computedFieldXArr,
                                        std::vector<double> &computedFieldYArr,
                                        std::vector<double> &computedFieldZArr,
                                        const PrecisionArguments &usedPrecision,
                                        ComputeMethod method = CPU_ST) const;

        void computeAllAPotentialX(const std::vector<double> &cylindricalZArr,
                                   const std::vector<double> &cylindricalRArr,
                                   const std::vector<double> &cylindricalPhiArr,
                                   std::vector<double> &computedPotentialArr,
                                   ComputeMethod method = CPU_ST) const;
        void computeAllAPotentialX(const std::vector<double> &cylindricalZArr,
                                   const std::vector<double> &cylindricalRArr,
                                   const std::vector<double> &cylindricalPhiArr,
                                   std::vector<double> &computedPotentialArr,
                                   const PrecisionArguments &usedPrecision,
                                   ComputeMethod method = CPU_ST) const;

        void computeAllAPotentialY(const std::vector<double> &cylindricalZArr,
                                   const std::vector<double> &cylindricalRArr,
                                   const std::vector<double> &cylindricalPhiArr,
                                   std::vector<double> &computedPotentialArr,
                                   ComputeMethod method = CPU_ST) const;
        void computeAllAPotentialY(const std::vector<double> &cylindricalZArr,
                                   const std::vector<double> &cylindricalRArr,
                                   const std::vector<double> &cylindricalPhiArr,
                                   std::vector<double> &computedPotentialArr,
                                   const PrecisionArguments &usedPrecision,
                                   ComputeMethod method = CPU_ST) const;

        void computeAllAPotentialZ(const std::vector<double> &cylindricalZArr,
                                   const std::vector<double> &cylindricalRArr,
                                   const std::vector<double> &cylindricalPhiArr,
                                   std::vector<double> &computedPotentialArr,
                                   ComputeMethod method = CPU_ST) const;
        void computeAllAPotentialZ(const std::vector<double> &cylindricalZArr,
                                   const std::vector<double> &cylindricalRArr,
                                   const std::vector<double> &cylindricalPhiArr,
                                   std::vector<double> &computedPotentialArr,
                                   const PrecisionArguments &usedPrecision,
                                   ComputeMethod method = CPU_ST) const;

        void computeAllAPotentialAbs(const std::vector<double> &cylindricalZArr,
                                     const std::vector<double> &cylindricalRArr,
                                     std::vector<double> &computedPotentialArr,
                                     ComputeMethod method = CPU_ST) const;
        void computeAllAPotentialAbs(const std::vector<double> &cylindricalZArr,
                                     const std::vector<double> &cylindricalRArr,
                                     std::vector<double> &computedPotentialArr,
                                     const PrecisionArguments &usedPrecision,
                                     ComputeMethod method = CPU_ST) const;

        void computeAllAPotentialComponents(const std::vector<double> &cylindricalZArr,
                                            const std::vector<double> &cylindricalRArr,
                                            const std::vector<double> &cylindricalPhiArr,
                                            std::vector<double> &computedPotentialXArr,
                                            std::vector<double> &computedPotentialYArr,
                                            std::vector<double> &computedPotentialZArr,
                                            ComputeMethod method = CPU_ST) const;
        void computeAllAPotentialComponents(const std::vector<double> &cylindricalZArr,
                                            const std::vector<double> &cylindricalRArr,
                                            const std::vector<double> &cylindricalPhiArr,
                                            std::vector<double> &computedPotentialXArr,
                                            std::vector<double> &computedPotentialYArr,
                                            std::vector<double> &computedPotentialZArr,
                                            const PrecisionArguments &usedPrecision,
                                            ComputeMethod method = CPU_ST) const;

        void computeAllEFieldX(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               ComputeMethod method = CPU_ST) const;
        void computeAllEFieldX(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               const PrecisionArguments &usedPrecision,
                               ComputeMethod method = CPU_ST) const;

        void computeAllEFieldY(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               ComputeMethod method = CPU_ST) const;
        void computeAllEFieldY(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               const PrecisionArguments &usedPrecision,
                               ComputeMethod method = CPU_ST) const;

        void computeAllEFieldZ(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               ComputeMethod method = CPU_ST) const;
        void computeAllEFieldZ(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               const PrecisionArguments &usedPrecision,
                               ComputeMethod method = CPU_ST) const;


        void computeAllEFieldAbs(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 std::vector<double> &computedFieldArr,
                                 ComputeMethod method = CPU_ST) const;
        void computeAllEFieldAbs(const std::vector<double> &cylindricalZArr,
                                 const std::vector<double> &cylindricalRArr,
                                 std::vector<double> &computedFieldArr,
                                 const PrecisionArguments &usedPrecision,
                                 ComputeMethod method = CPU_ST) const;

        void computeAllEFieldComponents(const std::vector<double> &cylindricalZArr,
                                        const std::vector<double> &cylindricalRArr,
                                        const std::vector<double> &cylindricalPhiArr,
                                        std::vector<double> &computedFieldXArr,
                                        std::vector<double> &computedFieldYArr,
                                        std::vector<double> &computedFieldZArr,
                                        ComputeMethod method = CPU_ST) const;
        void computeAllEFieldComponents(const std::vector<double> &cylindricalZArr,
                                        const std::vector<double> &cylindricalRArr,
                                        const std::vector<double> &cylindricalPhiArr,
                                        std::vector<double> &computedFieldXArr,
                                        std::vector<double> &computedFieldYArr,
                                        std::vector<double> &computedFieldZArr,
                                        const PrecisionArguments &usedPrecision,
                                        ComputeMethod method = CPU_ST) const;

        void computeAllBGradientTensors(const std::vector<double> &cylindricalZArr,
                                        const std::vector<double> &cylindricalRArr,
                                        const std::vector<double> &cylindricalPhiArr,
                                        std::vector<double> &computedGradientXX,
                                        std::vector<double> &computedGradientXY,
                                        std::vector<double> &computedGradientXZ,
                                        std::vector<double> &computedGradientYX,
                                        std::vector<double> &computedGradientYY,
                                        std::vector<double> &computedGradientYZ,
                                        std::vector<double> &computedGradientZX,
                                        std::vector<double> &computedGradientZY,
                                        std::vector<double> &computedGradientZZ,
                                        ComputeMethod method = CPU_ST) const;
        void computeAllBGradientTensors(const std::vector<double> &cylindricalZArr,
                                        const std::vector<double> &cylindricalRArr,
                                        const std::vector<double> &cylindricalPhiArr,
                                        std::vector<double> &computedGradientXX,
                                        std::vector<double> &computedGradientXY,
                                        std::vector<double> &computedGradientXZ,
                                        std::vector<double> &computedGradientYX,
                                        std::vector<double> &computedGradientYY,
                                        std::vector<double> &computedGradientYZ,
                                        std::vector<double> &computedGradientZX,
                                        std::vector<double> &computedGradientZY,
                                        std::vector<double> &computedGradientZZ,
                                        const PrecisionArguments &usedPrecision,
                                        ComputeMethod method = CPU_ST) const;

        static double computeMutualInductance(const Coil &primary, const Coil &secondary, double zDisplacement,
                                              PrecisionFactor precisionFactor = PrecisionFactor(),
                                              ComputeMethod method = CPU_ST);
        static double computeMutualInductance(const Coil &primary, const Coil &secondary, double zDisplacement,
                                              CoilPairArguments inductanceArguments, ComputeMethod method = CPU_ST);

        static double computeMutualInductance(const Coil &primary, const Coil &secondary,
                                              double zDisplacement, double rDisplacement,
                                              PrecisionFactor precisionFactor = PrecisionFactor(),
                                              ComputeMethod method = CPU_ST);
        static double computeMutualInductance(const Coil &primary, const Coil &secondary,
                                              double zDisplacement, double rDisplacement,
                                              CoilPairArguments inductanceArguments, ComputeMethod method = CPU_ST);

        static double computeMutualInductance(const Coil &primary, const Coil &secondary,
                                              double zDisplacement, double rDisplacement, double alphaAngle,
                                              PrecisionFactor precisionFactor = PrecisionFactor(),
                                              ComputeMethod method = CPU_ST);
        static double computeMutualInductance(const Coil &primary, const Coil &secondary,
                                              double zDisplacement, double rDisplacement, double alphaAngle,
                                              CoilPairArguments inductanceArguments, ComputeMethod method = CPU_ST);

        static double computeMutualInductance(const Coil &primary, const Coil &secondary,
                                              double zDisplacement, double rDisplacement,
                                              double alphaAngle, double betaAngle,
                                              PrecisionFactor precisionFactor = PrecisionFactor(),
                                              ComputeMethod method = CPU_ST);
        static double computeMutualInductance(const Coil &primary, const Coil &secondary,
                                              double zDisplacement, double rDisplacement,
                                              double alphaAngle, double betaAngle,
                                              CoilPairArguments inductanceArguments, ComputeMethod method = CPU_ST);

        [[nodiscard]] double computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement,
                                                            PrecisionFactor precisionFactor = PrecisionFactor(),
                                                            ComputeMethod method = CPU_ST) const;
        [[nodiscard]] double computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement,
                                                            CoilPairArguments inductanceArguments,
                                                            ComputeMethod method = CPU_ST) const;

        [[nodiscard]] double computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement, double rDisplacement,
                                                            PrecisionFactor precisionFactor = PrecisionFactor(),
                                                            ComputeMethod method = CPU_ST) const;
        [[nodiscard]] double computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement, double rDisplacement,
                                                            CoilPairArguments inductanceArguments,
                                                            ComputeMethod method = CPU_ST) const;

        [[nodiscard]] double computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement, double rDisplacement,
                                                            double alphaAngle, PrecisionFactor precisionFactor = PrecisionFactor(),
                                                            ComputeMethod method = CPU_ST) const;
        [[nodiscard]] double computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement, double rDisplacement,
                                                            double alphaAngle, CoilPairArguments inductanceArguments,
                                                            ComputeMethod method = CPU_ST) const;

        [[nodiscard]] double computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement, double rDisplacement,
                                                            double alphaAngle, double betaAngle,
                                                            PrecisionFactor precisionFactor = PrecisionFactor(),
                                                            ComputeMethod method = CPU_ST) const;
        [[nodiscard]] double computeSecondaryInducedVoltage(const Coil &secondary, double zDisplacement, double rDisplacement,
                                                            double alphaAngle, double betaAngle,
                                                            CoilPairArguments inductanceArguments,
                                                            ComputeMethod method = CPU_ST) const;

        double computeAndSetSelfInductance(PrecisionFactor precisionFactor);

        double computeAndSetApproximateSelfInductance(PrecisionFactor precisionFactor, ComputeMethod method = CPU_ST);

        static double computeAmpereForceZAxis(const Coil &primary, const Coil &secondary, double zDisplacement,
                                              PrecisionFactor precisionFactor = PrecisionFactor(),
                                              ComputeMethod method = CPU_ST);
        static double computeAmpereForceZAxis(const Coil &primary, const Coil &secondary, double zDisplacement,
                                              CoilPairArguments inductanceArguments, ComputeMethod method = CPU_ST);

        static std::vector<double> computeAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                                             double zDisplacement, double rDisplacement,
                                                             PrecisionFactor precisionFactor = PrecisionFactor(),
                                                             ComputeMethod method = CPU_ST);
        static std::vector<double> computeAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                                             double zDisplacement, double rDisplacement,
                                                             CoilPairArguments forceArguments, ComputeMethod method = CPU_ST);

        static std::vector<double> computeAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                                             double zDisplacement, double rDisplacement, double alphaAngle,
                                                             PrecisionFactor precisionFactor = PrecisionFactor(),
                                                             ComputeMethod method = CPU_ST);
        static std::vector<double> computeAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                                             double zDisplacement, double rDisplacement, double alphaAngle,
                                                             CoilPairArguments forceArguments, ComputeMethod method = CPU_ST);

        static std::vector<double> computeAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                                             double zDisplacement, double rDisplacement,
                                                             double alphaAngle, double betaAngle,
                                                             PrecisionFactor precisionFactor = PrecisionFactor(),
                                                             ComputeMethod method = CPU_ST);
        static std::vector<double> computeAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                                             double zDisplacement, double rDisplacement,
                                                             double alphaAngle, double betaAngle,
                                                             CoilPairArguments forceArguments, ComputeMethod method = CPU_ST);

    private:
        void calculateMagneticMoment();
        void calculateAverageWireThickness();
        void calculateResistance();
        void calculateReactance();
        void calculateImpedance();

        [[nodiscard]] std::pair<double, double> calculateBField(double zAxis, double rPolar,
                                                                const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] double calculateAPotential(double zAxis, double rPolar,
                                                 const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::vector<double> calculateBGradient(double zAxis, double rPolar,
                                                             const PrecisionArguments &usedPrecision) const;

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
                                  const PrecisionArguments &usedPrecision) const;

        void calculateAllAPotentialMT(const std::vector<double> &cylindricalZArr,
                                      const std::vector<double> &cylindricalRArr,
                                      std::vector<double> &computedPotentialArr,
                                      const PrecisionArguments &usedPrecision) const;

        void calculateAllBGradientMT(const std::vector<double> &cylindricalZArr,
                                     const std::vector<double> &cylindricalRArr,
                                     std::vector<double> &computedGradientRPhi,
                                     std::vector<double> &computedGradientRR,
                                     std::vector<double> &computedGradientRZ,
                                     std::vector<double> &computedGradientZZ,
                                     const PrecisionArguments &usedPrecision) const;

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

        static void calculateRingIncrementPosition(int angularBlocks, int angularIncrements,
                                                   double alpha, double beta, double ringIntervalSize,
                                                   std::vector<double> &ringXPosition,
                                                   std::vector<double> &ringYPosition,
                                                   std::vector<double> &ringZPosition,
                                                   std::vector<double> &ringXTangent,
                                                   std::vector<double> &ringYTangent,
                                                   std::vector<double> &ringZTangent);

        static double calculateMutualInductanceZAxis(const Coil &primary, const Coil &secondary, double zDisplacement,
                                                     CoilPairArguments inductanceArguments,
                                                     ComputeMethod method = CPU_ST);

        static double calculateMutualInductanceGeneral(const Coil &primary, const Coil &secondary,
                                                       double zDisplacement, double rDisplacement,
                                                       double alphaAngle, double betaAngle,
                                                       CoilPairArguments inductanceArguments,
                                                       ComputeMethod method = CPU_ST);

        void calculateAndSetSelfInductance(PrecisionFactor precisionFactor);

        void calculateAndSetApproximateSelfInductance(PrecisionFactor precisionFactor, ComputeMethod method = CPU_ST);

        static double calculateAmpereForceZAxis(const Coil &primary, const Coil &secondary, double zDisplacement,
                                                CoilPairArguments forceArguments,
                                                ComputeMethod method = CPU_ST);

        static std::vector<double> calculateAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                                               double zDisplacement, double rDisplacement,
                                                               double alphaAngle, double betaAngle,
                                                               CoilPairArguments forceArguments, ComputeMethod method);
};

#endif //GENERAL_COIL_PROGRAM_COIL_H
