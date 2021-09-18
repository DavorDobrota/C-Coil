#ifndef GENERAL_COIL_PROGRAM_COILGROUP_H
#define GENERAL_COIL_PROGRAM_COILGROUP_H

#include "Coil.h"

class CoilGroup
{
    public:

        std::vector<Coil> memberCoils;
        PrecisionArguments defaultPrecision;
        int threadCount{};
        
    public:
        explicit CoilGroup(std::vector<Coil> memberCoils = std::vector<Coil>(), PrecisionFactor precisionFactor = PrecisionFactor(), 
                           int threadCount = defaultThreadCount);
        explicit CoilGroup(std::vector<Coil> memberCoils, PrecisionArguments defaultPrecision,int threadCount = defaultThreadCount);
        
        [[nodiscard]] PrecisionArguments getDefaultPrecision() const;
        [[nodiscard]] int getThreadCount() const;
        [[nodiscard]] const std::vector<Coil> &getMemberCoils() const;
        
        void setDefaultPrecision(const PrecisionArguments &defaultPrecision);
        void setDefaultPrecision(PrecisionFactor precisionFactor = PrecisionFactor(), ComputeMethod method = CPU_ST);
        void setThreadCount(int threadCount);
        void addCoil(Coil coil);

        [[nodiscard]] std::vector<vec3::FieldVector3>
        computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                   ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<vec3::FieldVector3>
        computeAllBFieldComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                   const PrecisionArguments &usedPrecision, ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<vec3::FieldVector3>
        computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                       ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<vec3::FieldVector3>
        computeAllAPotentialComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                       const PrecisionArguments &usedPrecision, ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<vec3::FieldVector3>
        computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                   ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<vec3::FieldVector3>
        computeAllEFieldComponents(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                   const PrecisionArguments &usedPrecision, ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::vector<vec3::Matrix3>
        computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                   ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::vector<vec3::Matrix3>
        computeAllBGradientTensors(const std::vector<vec3::CoordVector3> &pointVectorArr,
                                   const PrecisionArguments &usedPrecision, ComputeMethod method = CPU_ST) const;

        [[nodiscard]] double computeMutualInductance(const Coil &secondary,
                                                     PrecisionFactor precisionFactor = PrecisionFactor(),
                                                     ComputeMethod method = CPU_ST) const;
        [[nodiscard]] double computeMutualInductance(const Coil &secondary,
                                                     CoilPairArguments inductanceArguments,
                                                     ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::pair<vec3::FieldVector3, vec3::FieldVector3>
        computeAmpereForce(const Coil &secondary,
                           PrecisionFactor precisionFactor = PrecisionFactor(), ComputeMethod method = CPU_ST) const;
        [[nodiscard]] std::pair<vec3::FieldVector3, vec3::FieldVector3>
        computeAmpereForce(const Coil &secondary, CoilPairArguments forceArguments, ComputeMethod method = CPU_ST) const;

        [[nodiscard]] std::pair<vec3::FieldVector3, vec3::FieldVector3>
        computeForceOnDipoleMoment(vec3::CoordVector3 pointVector, vec3::FieldVector3 dipoleMoment) const;
        [[nodiscard]] std::pair<vec3::FieldVector3, vec3::FieldVector3>
        computeForceOnDipoleMoment(vec3::CoordVector3 pointVector, vec3::FieldVector3 dipoleMoment,
                                   const PrecisionArguments &usedPrecision) const;

    private:
        // MTD stands for Multithreading Distributed - useful when there are many coils, each is given its own thread
        [[nodiscard]] std::vector<vec3::FieldVector3>
        calculateAllBFieldComponentsMTD(const std::vector<vec3::CoordVector3> &pointVectorArr, 
                                        const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::vector<vec3::FieldVector3>
        calculateAllAPotentialComponentsMTD(const std::vector<vec3::CoordVector3> &pointVectorArr, 
                                            const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] std::vector<vec3::FieldVector3>
        calculateAllBGradientTensorsMTD(const std::vector<vec3::CoordVector3> &pointVectorArr, 
                                        const PrecisionArguments &usedPrecision) const;
};

#endif //GENERAL_COIL_PROGRAM_COILGROUP_H
