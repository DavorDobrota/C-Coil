
#ifndef GENERAL_COIL_PROGRAM_COIL_H
#define GENERAL_COIL_PROGRAM_COIL_H

#define calc_type double

class Coil
{
    private:
        double currentDensity;
        double current;
        int numOfTurns;

        double internalRadius;
        double thickness;
        double length;
        double averageWireThickness;

        bool isSineDriven;
        double sineFrequency;
        double selfInductance;
        double magneticMoment;

        double wireResistivity;
        double resistance;
        double reactance;
        double impedance;

        double precisionFactor;

        int numOfAngularBlocks;
        int numOfThicknessBlocks;
        int numOfLengthBlocks;

        int numOfAngularIncrements;
        int numOfThicknessIncrements;
        int numOfLengthIncrements;

        std::vector<calc_type> angularIncrementPositions;
        std::vector<calc_type> angularIncrementWeights;

        std::vector<calc_type> thicknessIncrementPositions;
        std::vector<calc_type> thicknessIncrementWeights;

        std::vector<calc_type> lengthIncrementPositions;
        std::vector<calc_type> lengthIncrementWeights;

    public:

};


#endif //GENERAL_COIL_PROGRAM_COIL_H
