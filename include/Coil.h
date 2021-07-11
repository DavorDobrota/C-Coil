
#ifndef GENERAL_COIL_PROGRAM_COIL_H
#define GENERAL_COIL_PROGRAM_COIL_H

#include "ComputeMethod.h"

struct PrecisionArguments
{
    PrecisionArguments(int numOfAngularBlocks, int numOfThicknessBlocks, int numOfLengthBlocks,
                       int numOfAngularIncrements, int numOfThicknessIncrements, int numOfLengthIncrements);

    explicit PrecisionArguments(double precisionFactor);

    double precisionFactor;

    int angularBlockCount;
    int thicknessBlockCount;
    int lengthBlockCount;

    int angularIncrementCount;
    int thicknessIncrementCount;
    int lengthIncrementCount;

    static int blockPrecisionArray[];
    static int incrementPrecisionArray[];

    private:
        void genParametersFromPrecision();

        void getMutualInductancePrecisionSettings(PrecisionArguments &fieldPrecision,
                                                  int &zIncrements, int &rIncrements);

        void getMutualInductancePrecisionSettings(PrecisionArguments &fieldPrecision,
                                                  int &zIncrements, int &rIncrements, int &phiIncrements);
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

        double computeBFieldX(double cylindricalZ, double cylindricalR, double cylindricalPhi);
        double computeBFieldX(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                              const PrecisionArguments &usedPrecision);

        double computeBFieldY(double cylindricalZ, double cylindricalR, double cylindricalPhi);
        double computeBFieldY(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                              const PrecisionArguments &usedPrecision);

        double computeBFieldH(double cylindricalZ, double cylindricalR);
        double computeBFieldH(double cylindricalZ, double cylindricalR, const PrecisionArguments &usedPrecision);

        double computeBFieldZ(double cylindricalZ, double cylindricalR);
        double computeBFieldZ(double cylindricalZ, double cylindricalR, const PrecisionArguments &usedPrecision);

        std::vector<double> computeBFieldVector(double cylindricalZ, double cylindricalR, double cylindricalPhi);
        std::vector<double> computeBFieldVector(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                                const PrecisionArguments &usedPrecision);

        double computeAPotentialX(double cylindricalZ, double cylindricalR, double cylindricalPhi);
        double computeAPotentialX(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                  const PrecisionArguments &usedPrecision);

        double computeAPotentialY(double cylindricalZ, double cylindricalR, double cylindricalPhi);
        double computeAPotentialY(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                  const PrecisionArguments &usedPrecision);

        double computeAPotentialAbs(double cylindricalZ, double cylindricalR);
        double computeAPotentialAbs(double cylindricalZ, double cylindricalR, PrecisionArguments &usedPrecision);

        std::vector<double> computeAPotentialVector(double cylindricalZ, double cylindricalR, double cylindricalPhi);
        std::vector<double> computeAPotentialVector(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                                    const PrecisionArguments &usedPrecision);

        void computeAllBFieldX(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               ComputeMethod method = CPU_ST);
        void computeAllBFieldX(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               const PrecisionArguments &usedPrecision,
                               ComputeMethod method = CPU_ST);

        void computeAllBFieldY(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               ComputeMethod method = CPU_ST);
        void computeAllBFieldY(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               const PrecisionArguments &usedPrecision,
                               ComputeMethod method = CPU_ST);

        void computeAllBFieldH(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               ComputeMethod method = CPU_ST);
        void computeAllBFieldH(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               const PrecisionArguments &usedPrecision,
                               ComputeMethod method = CPU_ST);

        void computeAllBFieldZ(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               ComputeMethod method = CPU_ST);
        void computeAllBFieldZ(const std::vector<double> &cylindricalZArr,
                               const std::vector<double> &cylindricalRArr,
                               const std::vector<double> &cylindricalPhiArr,
                               std::vector<double> &computedFieldArr,
                               const PrecisionArguments &usedPrecision,
                               ComputeMethod method = CPU_ST);

        void computeAllBFieldComponents(const std::vector<double> &cylindricalZArr,
                                        const std::vector<double> &cylindricalRArr,
                                        const std::vector<double> &cylindricalPhiArr,
                                        std::vector<double> &computedFieldXArr,
                                        std::vector<double> &computedFieldYArr,
                                        std::vector<double> &computedFieldZArr,
                                        ComputeMethod method = CPU_ST);
        void computeAllBFieldComponents(const std::vector<double> &cylindricalZArr,
                                        const std::vector<double> &cylindricalRArr,
                                        const std::vector<double> &cylindricalPhiArr,
                                        std::vector<double> &computedFieldXArr,
                                        std::vector<double> &computedFieldYArr,
                                        std::vector<double> &computedFieldZArr,
                                        const PrecisionArguments &usedPrecision,
                                        ComputeMethod method = CPU_ST);

        void computeAllAPotentialX(const std::vector<double> &cylindricalZArr,
                                   const std::vector<double> &cylindricalRArr,
                                   const std::vector<double> &cylindricalPhiArr,
                                   std::vector<double> &computedPotentialArr,
                                   ComputeMethod method = CPU_ST);
        void computeAllAPotentialX(const std::vector<double> &cylindricalZArr,
                                   const std::vector<double> &cylindricalRArr,
                                   const std::vector<double> &cylindricalPhiArr,
                                   std::vector<double> &computedPotentialArr,
                                   const PrecisionArguments &usedPrecision,
                                   ComputeMethod method = CPU_ST);

        void computeAllAPotentialY(const std::vector<double> &cylindricalZArr,
                                   const std::vector<double> &cylindricalRArr,
                                   const std::vector<double> &cylindricalPhiArr,
                                   std::vector<double> &computedPotentialArr,
                                   ComputeMethod method = CPU_ST);
        void computeAllAPotentialY(const std::vector<double> &cylindricalZArr,
                                   const std::vector<double> &cylindricalRArr,
                                   const std::vector<double> &cylindricalPhiArr,
                                   std::vector<double> &computedPotentialArr,
                                   const PrecisionArguments &usedPrecision,
                                   ComputeMethod method = CPU_ST);

        void computeAllAPotentialAbs(const std::vector<double> &cylindricalZArr,
                                     const std::vector<double> &cylindricalRArr,
                                     std::vector<double> &computedPotentialArr,
                                     ComputeMethod method = CPU_ST);
        void computeAllAPotentialAbs(const std::vector<double> &cylindricalZArr,
                                     const std::vector<double> &cylindricalRArr,
                                     std::vector<double> &computedPotentialArr,
                                     const PrecisionArguments &usedPrecision,
                                     ComputeMethod method = CPU_ST);

        void computeAllAPotentialComponents(const std::vector<double> &cylindricalZArr,
                                            const std::vector<double> &cylindricalRArr,
                                            const std::vector<double> &cylindricalPhiArr,
                                            std::vector<double> &computedPotentialXArr,
                                            std::vector<double> &computedPotentialYArr,
                                            std::vector<double> &computedPotentialZArr,
                                            ComputeMethod method = CPU_ST);
        void computeAllAPotentialComponents(const std::vector<double> &cylindricalZArr,
                                            const std::vector<double> &cylindricalRArr,
                                            const std::vector<double> &cylindricalPhiArr,
                                            std::vector<double> &computedPotentialXArr,
                                            std::vector<double> &computedPotentialYArr,
                                            std::vector<double> &computedPotentialZArr,
                                            const PrecisionArguments &usedPrecision,
                                            ComputeMethod method = CPU_ST);

        double computeMutualInductance(double zDisplacement, Coil secondary, ComputeMethod method = CPU_ST);

    private:
        void calculateMagneticMoment();
        void calculateAverageWireThickness();
        void calculateResistance();
        void calculateReactance();
        void calculateImpedance();
        void calculateSelfInductance();

        std::pair<double, double> calculateBField(double zAxis, double rPolar,
                                                  const PrecisionArguments &precisionSettings);
        double calculateBFieldVertical(double zAxis, double rPolar,
                                       const PrecisionArguments &precisionSettings);
        double calculateBFieldHorizontal(double zAxis, double rPolar,
                                         const PrecisionArguments &precisionSettings);
        double calculateAPotential(double zAxis, double rPolar,
                                   const PrecisionArguments &precisionSettings);

        void calculateAllBFieldSINGLE(const std::vector<double> &cylindricalZArr,
                                      const std::vector<double> &cylindricalRArr,
                                      std::vector<double> &computedFieldHArr,
                                      std::vector<double> &computedFieldZArr,
                                      const PrecisionArguments &usedPrecision);
        void calculateAllBFieldVerticalSINGLE(const std::vector<double> &cylindricalZArr,
                                              const std::vector<double> &cylindricalRArr,
                                              std::vector<double> &computedFieldZArr,
                                              const PrecisionArguments &usedPrecision);
        void calculateAllBFieldHorizontalSINGLE(const std::vector<double> &cylindricalZArr,
                                                const std::vector<double> &cylindricalRArr,
                                                std::vector<double> &computedFieldHArr,
                                                const PrecisionArguments &usedPrecision);
        void calculateAllAPotentialSINGLE(const std::vector<double> &cylindricalZArr,
                                          const std::vector<double> &cylindricalRArr,
                                          std::vector<double> &computedPotentialArr,
                                          const PrecisionArguments &usedPrecision);

        void convertPolarToCylindrical(double polarR, double polarTheta, double polarPhi,
                                       double &cylindricalZ, double &cylindricalR, double &cylindricalPhi);

        void convertAllPolarToCylindrical(const std::vector<double> &polarRArr,
                                          const std::vector<double> &polarThetaArr,
                                          const std::vector<double> &polarPhiArr,
                                          std::vector<double> &cylindricalZArr,
                                          std::vector<double> &cylindricalRArr,
                                          std::vector<double> &cylindricalPhiArr);

        void convertCylindricalToPolar(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                       double &polarR, double &polarTheta, double &polarPhi);

        void convertAllCylindricalToPolar(const std::vector<double> &cylindricalZArr,
                                          const std::vector<double> &cylindricalRArr,
                                          const std::vector<double> &cylindricalPhiArr,
                                          std::vector<double> &polarRArr,
                                          std::vector<double> &polarThetaArr,
                                          std::vector<double> &polarPhiArr);

        void calculateAllBFieldACCELERATED(const std::vector<double> &cylindricalZArr,
                                           const std::vector<double> &cylindricalRArr,
                                           std::vector<float> &computedFieldHArr,
                                           std::vector<float> &computedFieldZArr,
                                           const PrecisionArguments &usedPrecision);

        void calculateAllAPotentialACCELERATED(const std::vector<double> &cylindricalZArr,
                                               const std::vector<double> &cylindricalRArr,
                                               std::vector<float> &computedPotentialArr,
                                               const PrecisionArguments &usedPrecision);
};

#endif //GENERAL_COIL_PROGRAM_COIL_H
