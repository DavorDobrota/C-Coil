#ifndef GENERAL_COIL_PROGRAM_COILGROUP_H
#define GENERAL_COIL_PROGRAM_COILGROUP_H

#include "Coil.h"

class CoilGroup
{
    private:
        std::vector<Coil> memberCoils;
        PrecisionFactor defaultPrecisionFactor;
        int threadCount{};
        
    public:
        explicit CoilGroup(std::vector<Coil> memberCoils = std::vector<Coil>(), PrecisionFactor precisionFactor = PrecisionFactor(), 
                           int threadCount = defaultThreadCount);
        
        [[nodiscard]] PrecisionFactor getDefaultPrecisionFactor() const;
        [[nodiscard]] int getThreadCount() const;
        [[nodiscard]] const std::vector<Coil> &getMemberCoils() const;

        void setDefaultPrecisionFactor(PrecisionFactor precisionFactor = PrecisionFactor(), ComputeMethod method = CPU_ST);
        void setThreadCount(int threadCount);
        void addCoil(Coil coil);

        [[nodiscard]] vec3::FieldVector3 computeBFieldVector(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] vec3::FieldVector3 computeAPotentialVector(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] vec3::FieldVector3 computeEFieldVector(vec3::CoordVector3 pointVector) const;
        [[nodiscard]] vec3::Matrix3 computeBGradientTensor(vec3::CoordVector3 pointVector) const;

        [[nodiscard]] std::vector<vec3::FieldVector3>
        computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                   ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<vec3::FieldVector3>
        computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                       ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<vec3::FieldVector3>
        computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                   ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<vec3::Matrix3>
        computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                   ComputeMethod method = CPU_ST) const;

        [[nodiscard]] double computeMutualInductance(const Coil &secondary,
                                                     PrecisionFactor precisionFactor = PrecisionFactor(),
                                                     ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::pair<vec3::FieldVector3, vec3::FieldVector3>
        computeAmpereForce(const Coil &secondary, PrecisionFactor precisionFactor = PrecisionFactor(),
                           ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::pair<vec3::FieldVector3, vec3::FieldVector3>
        computeForceOnDipoleMoment(vec3::CoordVector3 pointVector, vec3::FieldVector3 dipoleMoment) const;

    private:
        // MTD stands for Multithreading Distributed - useful when there are many coils, each is given its own thread
        [[nodiscard]] std::vector<vec3::FieldVector3>
        calculateAllBFieldComponentsMTD(const std::vector<vec3::CoordVector3> &pointVectorArr) const;

        [[nodiscard]] std::vector<vec3::FieldVector3>
        calculateAllAPotentialComponentsMTD(const std::vector<vec3::CoordVector3> &pointVectorArr) const;

        [[nodiscard]] std::vector<vec3::FieldVector3>
        calculateAllEFieldComponentsMTD(const std::vector<vec3::CoordVector3> &pointVectorArr) const;

        [[nodiscard]] std::vector<vec3::Matrix3>
        calculateAllBGradientTensorsMTD(const std::vector<vec3::CoordVector3> &pointVectorArr) const;

        [[nodiscard]] double computeMutualInductanceMTD(const Coil &secondary,
                                                        PrecisionFactor precisionFactor = PrecisionFactor()) const;

        [[nodiscard]] std::pair<vec3::FieldVector3, vec3::FieldVector3>
        computeAmpereForceMTD(const Coil &secondary, PrecisionFactor precisionFactor = PrecisionFactor()) const;
};

#endif //GENERAL_COIL_PROGRAM_COILGROUP_H
