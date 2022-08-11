#ifndef GENERAL_COIL_PROGRAM_COILGROUP_H
#define GENERAL_COIL_PROGRAM_COILGROUP_H

#include "Coil.h"

#include <string>
#include <memory>

class CoilGroup
{
    private:
        std::vector<std::shared_ptr<Coil>> memberCoils;
        PrecisionFactor defaultPrecisionFactor;
        int threadCount{};
        
    public:
        explicit CoilGroup(std::vector<std::shared_ptr<Coil>> memberCoils = std::vector<std::shared_ptr<Coil>>(),
                           PrecisionFactor precisionFactor = PrecisionFactor(),
                           int threadCount = g_defaultThreadCount);
        
        [[nodiscard]] PrecisionFactor getDefaultPrecisionFactor() const;
        [[nodiscard]] int getThreadCount() const;
        [[nodiscard]] const std::vector<std::shared_ptr<Coil>> &getMemberCoils() const;

        void setDefaultPrecisionFactor(PrecisionFactor precisionFactor = PrecisionFactor());
        void setThreadCount(int threadCount);
        void addCoil(double innerRadius, double thickness, double length, int numOfTurns, double current = 1.0,
                     PrecisionFactor precisionFactor = PrecisionFactor(), int coilThreads = g_defaultThreadCount,
                     vec3::Vector3 coordinatePosition = vec3::Vector3(), double yAxisAngle = 0.0, double zAxisAngle = 0.0);

        Coil& operator[](size_t index) const;

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
        computeAllMutualInductanceArrangements(const Coil &secondary,
                                               const vec3::Vector3Array &secondaryPositions,
                                               const std::vector<double> &secondaryYAngles,
                                               const std::vector<double> &secondaryZAngles,
                                               PrecisionFactor precisionFactor = PrecisionFactor(),
                                               ComputeMethod computeMethod = CPU_ST) const;

        [[nodiscard]] std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
        computeAllAmpereForceArrangements(const Coil &secondary,
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

        void
        generateCoilDataArray(CoilData *coilDataArr, PrecisionFactor precisionFactor, bool removeSpecificCoil = false,
                              unsigned long long int specificCoilId = 0) const;
        static void generateSecondaryData(const Coil &secondary, SecondaryCoilData &secondaryData,
                                          const PrecisionArguments &precision, bool forceCalculation);

        [[nodiscard]] vec3::Vector3Array calculateAllAPotentialGPU(const vec3::Vector3Array &pointVectors) const;
        [[nodiscard]] vec3::Vector3Array calculateAllBFieldGPU(const vec3::Vector3Array &pointVectors) const;
        [[nodiscard]] vec3::Vector3Array calculateAllEFieldGPU(const vec3::Vector3Array &pointVectors) const;
        [[nodiscard]] vec3::Matrix3Array calculateAllBGradientGPU(const vec3::Vector3Array &pointVectors) const;


        [[nodiscard]] double calculateMutualInductanceMTD(const Coil &secondary,
                                                          PrecisionFactor precisionFactor = PrecisionFactor()) const;

        [[nodiscard]] std::pair<vec3::Vector3, vec3::Vector3>
        calculateAmpereForceMTD(const Coil &secondary, PrecisionFactor precisionFactor = PrecisionFactor()) const;

        [[nodiscard]] std::vector<double>
        calculateAllMutualInductanceArrangementsMTD(const Coil &secondary,
                                                    const vec3::Vector3Array &secondaryPositions,
                                                    const std::vector<double> &secondaryYAngles,
                                                    const std::vector<double> &secondaryZAngles,
                                                    PrecisionFactor precisionFactor = PrecisionFactor()) const;
        [[nodiscard]] std::vector<double>
        calculateAllMutualInductanceArrangementsGPU(const Coil &secondary,
                                                    const vec3::Vector3Array &secondaryPositions,
                                                    const std::vector<double> &secondaryYAngles,
                                                    const std::vector<double> &secondaryZAngles,
                                                    PrecisionFactor precisionFactor = PrecisionFactor()) const;

        [[nodiscard]] std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
        calculateAllAmpereForceArrangementsMTD(const Coil &secondary,
                                               const vec3::Vector3Array &secondaryPositions,
                                               const std::vector<double> &secondaryYAngles,
                                               const std::vector<double> &secondaryZAngles,
                                               PrecisionFactor precisionFactor = PrecisionFactor()) const;

        [[nodiscard]] std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
        calculateAllAmpereForceArrangementsGPU(const Coil &secondary,
                                               const vec3::Vector3Array &secondaryPositions,
                                               const std::vector<double> &secondaryYAngles,
                                               const std::vector<double> &secondaryZAngles,
                                               PrecisionFactor precisionFactor = PrecisionFactor()) const;

};

#endif //GENERAL_COIL_PROGRAM_COILGROUP_H
