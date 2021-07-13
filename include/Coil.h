
#ifndef GENERAL_COIL_PROGRAM_COIL_H
#define GENERAL_COIL_PROGRAM_COIL_H

#include "ComputeMethod.h"

const int precisionArraySize = 300;

class Coil;

struct PrecisionArguments
{
    PrecisionArguments();
    explicit PrecisionArguments(int numOfAngularBlocks, int numOfThicknessBlocks, int numOfLengthBlocks,
                       int numOfAngularIncrements, int numOfThicknessIncrements, int numOfLengthIncrements);

    explicit PrecisionArguments(double precisionFactor);

    double precisionFactor;

    int angularBlockCount;
    int thicknessBlockCount;
    int lengthBlockCount;

    int angularIncrementCount;
    int thicknessIncrementCount;
    int lengthIncrementCount;

    private:
        static const int blockPrecisionCPUArray[precisionArraySize];
        static const int incrementPrecisionCPUArray[precisionArraySize];

    public:
        static void getMutualInductancePrecisionSettingsZCPU(const Coil &primCoil, const Coil &secCoil,
                                                             double precisionFactor,
                                                             PrecisionArguments &fieldPrecision,
                                                             int &linearIncrements);

        static void getMutualInductancePrecisionSettingsGeneralCPU(const Coil &primCoil, const Coil &secCoil,
                                                                   double precisionFactor,
                                                                   PrecisionArguments &fieldPrecision,
                                                                   int &linearIncrements, int &angularIncrements);

    private:
        void genParametersFromPrecision();

};

class Coil
{
    private:

        const double innerRadius;
        const double thickness;
        const double length;
        const int numOfTurns;

        double currentDensity;
        double current;

        double wireResistivity;
        bool isSineDriven;
        double sineFrequency;

        double magneticMoment;
        double averageWireThickness;

        double resistance;
        double selfInductance;
        double reactance;
        double impedance;

        PrecisionArguments precisionSettings;

    public:
        Coil();

        Coil(double innerRadius, double thickness, double length, int numOfTurns,
             double current, double wireResistivity, double sineFrequency, const PrecisionArguments &precisionSettings);

        Coil(double innerRadius, double thickness, double length, int numOfTurns, double current, double sineFrequency);
        Coil(double innerRadius, double thickness, double length, int numOfTurns, double current, double sineFrequency,
             const PrecisionArguments &precisionSettings);

        Coil(double innerRadius, double thickness, double length, int numOfTurns, double current);
        Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
             const PrecisionArguments &precisionSettings);

        Coil(double innerRadius, double thickness, double length, int numOfTurns);
        Coil(double innerRadius, double thickness, double length, int numOfTurns,
             const PrecisionArguments &precisionSettings);


        [[nodiscard]] double getInnerRadius() const;
        [[nodiscard]] double getThickness() const;
        [[nodiscard]] double getLength() const;
        [[nodiscard]] int getNumOfTurns() const;

        [[nodiscard]] double getCurrentDensity() const;
        [[nodiscard]] double getCurrent() const;

        [[nodiscard]] double getWireResistivity() const;
        [[nodiscard]] bool isSineDriven1() const;
        [[nodiscard]] double getSineFrequency() const;

        [[nodiscard]] double getMagneticMoment() const;
        [[nodiscard]] double getAverageWireThickness() const;

        [[nodiscard]] double getResistance() const;
        [[nodiscard]] double getSelfInductance() const;
        [[nodiscard]] double getReactance() const;
        [[nodiscard]] double getImpedance() const;

        [[nodiscard]] const PrecisionArguments &getPrecisionSettings() const;

        void setCurrentDensity(double currentDensity);
        void setCurrent(double current);
        void setWireResistivity(double wireResistivity);
        void setSineFrequency(double sineFrequency);
        void setPrecisionSettings(const PrecisionArguments &precisionSettings);

        double computeBFieldX(double cylindricalZ, double cylindricalR, double cylindricalPhi) const;
        double computeBFieldX(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                              const PrecisionArguments &usedPrecision) const;

        double computeBFieldY(double cylindricalZ, double cylindricalR, double cylindricalPhi) const;
        double computeBFieldY(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                              const PrecisionArguments &usedPrecision) const;

        double computeBFieldH(double cylindricalZ, double cylindricalR) const;
        double computeBFieldH(double cylindricalZ, double cylindricalR, const PrecisionArguments &usedPrecision) const;

        double computeBFieldZ(double cylindricalZ, double cylindricalR) const;
        double computeBFieldZ(double cylindricalZ, double cylindricalR, const PrecisionArguments &usedPrecision) const;

        std::vector<double> computeBFieldVector(double cylindricalZ, double cylindricalR, double cylindricalPhi) const;
        std::vector<double> computeBFieldVector(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                                const PrecisionArguments &usedPrecision) const;

        double computeAPotentialX(double cylindricalZ, double cylindricalR, double cylindricalPhi) const;
        double computeAPotentialX(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                  const PrecisionArguments &usedPrecision) const;

        double computeAPotentialY(double cylindricalZ, double cylindricalR, double cylindricalPhi) const;
        double computeAPotentialY(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                  const PrecisionArguments &usedPrecision) const;

        double computeAPotentialAbs(double cylindricalZ, double cylindricalR) const;
        double computeAPotentialAbs(double cylindricalZ, double cylindricalR, PrecisionArguments &usedPrecision) const;

        std::vector<double> computeAPotentialVector(double cylindricalZ, double cylindricalR, double cylindricalPhi) const;
        std::vector<double> computeAPotentialVector(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                                    const PrecisionArguments &usedPrecision) const;

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
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               ComputeMethod method = CPU_ST) const;
        void computeAllBFieldH(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
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

        static double computeMutualInductance(const Coil &primary, const Coil &secondary, double zDisplacement,
                                              double precisionFactor = 5.0, ComputeMethod method = CPU_ST);

        static double computeMutualInductance(const Coil &primary, const Coil &secondary,
                                             double zDisplacement, double rDisplacement,
                                             double precisionFactor = 5.0, ComputeMethod method = CPU_ST);

        static double computeMutualInductance(const Coil &primary, const Coil &secondary,
                                              double zDisplacement, double rDisplacement, double primaryRotationAngle,
                                              double precisionFactor = 5.0, ComputeMethod method = CPU_ST);

        static double computeMutualInductance(const Coil &primary, const Coil &secondary,
                                              double zDisplacement, double rDisplacement,
                                              double primaryRotationAngle, double secondaryRotationAngle,
                                              double precisionFactor = 5.0, ComputeMethod method = CPU_ST);



    private:
        void calculateMagneticMoment();
        void calculateAverageWireThickness();
        void calculateResistance();
        void calculateReactance();
        void calculateImpedance();
        void calculateSelfInductance();

        std::pair<double, double> calculateBField(double zAxis, double rPolar,
                                                  const PrecisionArguments &precisionSettings) const;
        double calculateBFieldVertical(double zAxis, double rPolar,
                                       const PrecisionArguments &precisionSettings) const;
        double calculateBFieldHorizontal(double zAxis, double rPolar,
                                         const PrecisionArguments &precisionSettings) const;
        double calculateAPotential(double zAxis, double rPolar,
                                   const PrecisionArguments &precisionSettings) const;

        void calculateAllBFieldSINGLE(const std::vector<double> &cylindricalZArr,
                                      const std::vector<double> &cylindricalRArr,
                                      std::vector<double> &computedFieldHArr,
                                      std::vector<double> &computedFieldZArr,
                                      const PrecisionArguments &usedPrecision) const;
        void calculateAllBFieldVerticalSINGLE(const std::vector<double> &cylindricalZArr,
                                              const std::vector<double> &cylindricalRArr,
                                              std::vector<double> &computedFieldZArr,
                                              const PrecisionArguments &usedPrecision) const;
        void calculateAllBFieldHorizontalSINGLE(const std::vector<double> &cylindricalZArr,
                                                const std::vector<double> &cylindricalRArr,
                                                std::vector<double> &computedFieldHArr,
                                                const PrecisionArguments &usedPrecision) const;
        void calculateAllAPotentialSINGLE(const std::vector<double> &cylindricalZArr,
                                          const std::vector<double> &cylindricalRArr,
                                          std::vector<double> &computedPotentialArr,
                                          const PrecisionArguments &usedPrecision) const;

        void calculateAllBFieldACCELERATED(const std::vector<double> &cylindricalZArr,
                                           const std::vector<double> &cylindricalRArr,
                                           std::vector<float> &computedFieldHArr,
                                           std::vector<float> &computedFieldZArr,
                                           const PrecisionArguments &usedPrecision) const;

        void calculateAllAPotentialACCELERATED(const std::vector<double> &cylindricalZArr,
                                               const std::vector<double> &cylindricalRArr,
                                               std::vector<float> &computedPotentialArr,
                                               const PrecisionArguments &usedPrecision) const;
};

#endif //GENERAL_COIL_PROGRAM_COIL_H
