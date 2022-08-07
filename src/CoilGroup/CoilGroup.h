#ifndef GENERAL_COIL_PROGRAM_COILGROUP_H
#define GENERAL_COIL_PROGRAM_COILGROUP_H

#include "Coil.h"

#include <string>


class CoilGroup
{
    private:
        std::vector<Coil> memberCoils;
        PrecisionFactor defaultPrecisionFactor;
        int threadCount{};
        
    public:
        explicit CoilGroup(std::vector<Coil> memberCoils = std::vector<Coil>(),
                           PrecisionFactor precisionFactor = PrecisionFactor(),
                           int threadCount = defaultThreadCount);
        
        [[nodiscard]] PrecisionFactor getDefaultPrecisionFactor() const;
        [[nodiscard]] int getThreadCount() const;
        [[nodiscard]] const std::vector<Coil> &getMemberCoils() const;

        void setDefaultPrecisionFactor(PrecisionFactor precisionFactor = PrecisionFactor());
        void setThreadCount(int threadCount);
        void addCoil(Coil coil);

        [[nodiscard]] vec3::Vector3 computeBFieldVector(vec3::Vector3 pointVector) const;
        [[nodiscard]] vec3::Vector3 computeAPotentialVector(vec3::Vector3 pointVector) const;
        [[nodiscard]] vec3::Vector3 computeEFieldVector(vec3::Vector3 pointVector) const;
        [[nodiscard]] vec3::Matrix3 computeBGradientMatrix(vec3::Vector3 pointVector) const;

        [[nodiscard]] vec3::Vector3Array computeAllAPotentialVectors(const vec3::Vector3Array &pointVectors,
                                                                     ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] vec3::Vector3Array computeAllBFieldVectors(const vec3::Vector3Array &pointVectors,
                                                                 ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] vec3::Vector3Array computeAllEFieldVectors(const vec3::Vector3Array &pointVectors,
                                                                 ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] vec3::Matrix3Array computeAllBGradientMatrices(const vec3::Vector3Array &pointVectors,
                                                                     ComputeMethod computeMethod = CPU_ST) const;


        [[nodiscard]] double computeMutualInductance(const Coil &secondary,
                                                     PrecisionFactor precisionFactor = PrecisionFactor(),
                                                     ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] std::pair<vec3::Vector3, vec3::Vector3>
        computeAmpereForce(const Coil &secondary, PrecisionFactor precisionFactor = PrecisionFactor(),
                           ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] std::pair<vec3::Vector3, vec3::Vector3>
        computeForceOnDipoleMoment(vec3::Vector3 pointVector, vec3::Vector3 dipoleMoment) const;

        [[nodiscard]] std::vector<double>
        computeAllMutualInductanceArrangements(Coil secondary,
                                               const vec3::Vector3Array &secondaryPositions,
                                               const std::vector<double> &secondaryYAngles,
                                               const std::vector<double> &secondaryZAngles,
                                               PrecisionFactor precisionFactor = PrecisionFactor(),
                                               ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
        computeAllAmpereForceArrangements(Coil secondary,
                                          const vec3::Vector3Array &secondaryPositions,
                                          const std::vector<double> &secondaryYAngles,
                                          const std::vector<double> &secondaryZAngles,
                                          PrecisionFactor precisionFactor = PrecisionFactor(),
                                          ComputeMethod computeMethod = CPU_ST) const;

        explicit operator std::string() const;

    private:
        // MTD stands for Multithreading Distributed - useful when there are many coils, each is given its own thread
        [[nodiscard]] vec3::Vector3Array calculateAllAPotentialMTD(const vec3::Vector3Array &pointVectors) const;
        [[nodiscard]] vec3::Vector3Array calculateAllBFieldMTD(const vec3::Vector3Array &pointVectors) const;
        [[nodiscard]] vec3::Vector3Array calculateAllEFieldMTD(const vec3::Vector3Array &pointVectors) const;
        [[nodiscard]] vec3::Matrix3Array calculateAllBGradientMTD(const vec3::Vector3Array &pointVectors) const;

        void generateCoilDataArray(CoilData *coilDataArr,
                                   bool removeSpecificCoil = false,
                                   unsigned long long specificCoilId = 0) const;
        static void generateSecondaryData(const Coil &secondary, SecondaryCoilData &secondaryData,
                                          bool forceCalculation = false);

        [[nodiscard]] vec3::Vector3Array calculateAllAPotentialGPU(const vec3::Vector3Array &pointVectors) const;
        [[nodiscard]] vec3::Vector3Array calculateAllBFieldGPU(const vec3::Vector3Array &pointVectors) const;
        [[nodiscard]] vec3::Vector3Array calculateAllEFieldGPU(const vec3::Vector3Array &pointVectors) const;
        [[nodiscard]] vec3::Matrix3Array calculateAllBGradientGPU(const vec3::Vector3Array &pointVectors) const;


        [[nodiscard]] double calculateMutualInductanceMTD(const Coil &secondary,
                                                          PrecisionFactor precisionFactor = PrecisionFactor()) const;

        [[nodiscard]] std::pair<vec3::Vector3, vec3::Vector3>
        calculateAmpereForceMTD(const Coil &secondary, PrecisionFactor precisionFactor = PrecisionFactor()) const;

        [[nodiscard]] std::vector<double>
        calculateAllMutualInductanceArrangementsMTD(Coil secondary,
                                                    const vec3::Vector3Array &secondaryPositions,
                                                    const std::vector<double> &secondaryYAngles,
                                                    const std::vector<double> &secondaryZAngles,
                                                    PrecisionFactor precisionFactor = PrecisionFactor()) const;
        [[nodiscard]] std::vector<double>
        calculateAllMutualInductanceArrangementsGPU(Coil secondary,
                                                    const vec3::Vector3Array &secondaryPositions,
                                                    const std::vector<double> &secondaryYAngles,
                                                    const std::vector<double> &secondaryZAngles,
                                                    PrecisionFactor precisionFactor = PrecisionFactor()) const;

        [[nodiscard]] std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
        calculateAllAmpereForceArrangementsMTD(Coil secondary,
                                               const vec3::Vector3Array &secondaryPositions,
                                               const std::vector<double> &secondaryYAngles,
                                               const std::vector<double> &secondaryZAngles,
                                               PrecisionFactor precisionFactor = PrecisionFactor()) const;

        [[nodiscard]] std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
        calculateAllAmpereForceArrangementsGPU(Coil secondary,
                                               const vec3::Vector3Array &secondaryPositions,
                                               const std::vector<double> &secondaryYAngles,
                                               const std::vector<double> &secondaryZAngles,
                                               PrecisionFactor precisionFactor = PrecisionFactor()) const;
};

#endif //GENERAL_COIL_PROGRAM_COILGROUP_H
