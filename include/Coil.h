
#ifndef GENERAL_COIL_PROGRAM_COIL_H
#define GENERAL_COIL_PROGRAM_COIL_H

struct PrecisionArguments
{
    PrecisionArguments(int numOfAngularBlocks, int numOfThicknessBlocks, int numOfLengthBlocks,
                       int numOfAngularIncrements, int numOfThicknessIncrements, int numOfLengthIncrements);

    explicit PrecisionArguments(double precisionFactor);

    double precisionFactor;

    int numOfAngularBlocks;
    int numOfThicknessBlocks;
    int numOfLengthBlocks;

    int numOfAngularIncrements;
    int numOfThicknessIncrements;
    int numOfLengthIncrements;

    std::vector<double> angularIncrementPositions;
    std::vector<double> angularIncrementWeights;

    std::vector<double> thicknessIncrementPositions;
    std::vector<double> thicknessIncrementWeights;

    std::vector<double> lengthIncrementPositions;
    std::vector<double> lengthIncrementWeights;

    private:
        void genParametersFromPrecision();
        void genPrecisionVectors();
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


        void calculateMagneticMoment();
        void calculateAverageWireThickness();
        void calculateResistance();
        void calculateReactance();
        void calculateImpedance();
        void calculateSelfInductance();

        std::pair<double, double> calculateBField(double zAxis, double rPolar);
        double calculateBFieldVertical(double zAxis, double rPolar);
        double calculateBFieldHorizontal(double zAxis, double rPolar);

        double calculateAPotential(double zAxis, double rPolar);

    public:
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
        double computeBFieldY(double cylindricalZ, double cylindricalR, double cylindricalPhi);
        double computeBFieldH(double cylindricalZ, double cylindricalR);
        double computeBFieldZ(double cylindricalZ, double cylindricalR);
        std::vector<double> computeBFieldVector(double cylindricalZ, double cylindricalR, double cylindricalPhi);

        double computeAPotentialX(double cylindricalZ, double cylindricalR, double cylindricalPhi);
        double computeAPotentialY(double cylindricalZ, double cylindricalR, double cylindricalPhi);
        double computeAPotentialAbs(double cylindricalZ, double cylindricalR);
        std::vector<double> computeAPotentialVector(double cylindricalZ, double cylindricalR, double cylindricalPhi);
};


#endif //GENERAL_COIL_PROGRAM_COIL_H
