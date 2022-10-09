#ifndef GENERAL_COIL_PROGRAM_COILGROUP_H
#define GENERAL_COIL_PROGRAM_COILGROUP_H

#include "Coil.h"

#include <string>
#include <memory>

/**
 * @brief Represents a collection of unique Coil instances and is useful for representing multi-coil systems.
 * Enables faster field and interaction calculations, especially when using the GPU.
 * @details As all coils are unique, shared pointers are used to ensure there is only once instance of each Coil.
 * The coils are accessed individually as if this class were a list. The default PrecisionFactor is defined as well as
 * the number of threads used for calculations. When there are many coils, calculations are accelerated with MTD
 * (Multi-Threading Distributed, or more commonly Coarse-grained parallelism) methods for 10-100% more performance,
 * especially when using high core counts. Most methods are similar to Coil methods. When calculating MInductance or
 * Force, a Coil from the group, can be passed as an argument and its contribution is then ignored.
 */
class CoilGroup
{
    private:
        std::vector<std::shared_ptr<Coil>> memberCoils;
        PrecisionFactor defaultPrecisionFactor;
        int threadCount{};
        
    public:
        /// @brief All arguments have defaults so it is a also a default constructor, best used that way.
        explicit CoilGroup(std::vector<std::shared_ptr<Coil>> memberCoils = std::vector<std::shared_ptr<Coil>>(),
                           PrecisionFactor precisionFactor = PrecisionFactor(),
                           int threadCount = g_defaultThreadCount);

        /// @brief Returns the default PrecisionFactor according to which default arguments for all members are generated.
        [[nodiscard]] PrecisionFactor getDefaultPrecisionFactor() const;
        /// @brief Returns the number of threads used in CPU_MT calculations.
        [[nodiscard]] int getThreadCount() const;
        /// @brief Returns a constant reference to internal std::vector
        [[nodiscard]] const std::vector<std::shared_ptr<Coil>> &getMemberCoils() const;

        /// @brief Setts the default PrecisionFactor which is immediately applied to all members.
        void setDefaultPrecisionFactor(PrecisionFactor precisionFactor = PrecisionFactor());
        /// @brief Setts the default number of threads which is immediately applied to all members.
        void setThreadCount(int threadCount);

        /// @brief Adds a new member Coil to the back, most common Coil constructor is imitated for simplicity.
        void addCoil(double innerRadius, double thickness, double length, int numOfTurns, double current = 1.0,
                     PrecisionFactor precisionFactor = PrecisionFactor(), int coilThreads = g_defaultThreadCount,
                     vec3::Vector3 coordinatePosition = vec3::Vector3(), double yAxisAngle = 0.0, double zAxisAngle = 0.0);
         ///@brief Removes the Coil at the selected index from the CoilGroup.
        void removeCoil(size_t index);

        Coil& operator[](size_t index) const;

        /// @brief Checks if a given point lies inside the circular coil if coil type is RECTANGULAR,
        /// or if it is exceedingly close to a FILAMENT, THIN, or FLAT coil.
        [[nodiscard]] bool isPointInside(vec3::Vector3 pointVector);

        /**
         * @brief Calculates vector potential A of the magnetic field at the specified point from all member Coils.
         * Uses internally defined defaultPrecisionFactor.
         * @param pointVector Radius vector from the origin to the point where the field is calculated.
         * @return 3D Cartesian vector which represents vector potential A.
         */
        [[nodiscard]] vec3::Vector3 computeAPotentialVector(vec3::Vector3 pointVector) const;
        /**
         * @brief Calculates magnetic flux density B (magnetic field) at the specified point from all member Coils.
         * Uses internally defined defaultPrecisionFactor.
         * @param pointVector Radius vector from the origin to the point where the field is calculated.
         * @return 3D Cartesian vector which represents vector potential A.
         */
        [[nodiscard]] vec3::Vector3 computeBFieldVector(vec3::Vector3 pointVector) const;
        /**
         * @brief Calculates the amplitude vector of electric field E in sinusoidal steady-state at the specified point
         * from all member Coils. Uses internally defined defaultPrecisionFactor.
         * @param pointVector Radius vector from the origin to the point where the field is calculated.
         * @return 3D Cartesian vector which represents vector potential A.
         */
        [[nodiscard]] vec3::Vector3 computeEFieldVector(vec3::Vector3 pointVector) const;
        /**
         * @brief Calculates the gradient G of the magnetic field (total derivative of B) at the specified point
         * from all member Coils. Uses internally defined defaultPrecisionFactor.
         * @param pointVector Radius vector from the origin to the point where the field is calculated.
         * @return 3D Cartesian vector which represents vector potential A.
         */
        [[nodiscard]] vec3::Matrix3 computeBGradientMatrix(vec3::Vector3 pointVector) const;

        /**
         * @brief Calculates vector potential A of the magnetic field from all member Coils for a number of
         * specified points. Uses internally defined defaultPrecisionFactor.
         * @details There are multiple compute methods, GPU acceleration is best suited for over 100 points.
         * MTD is used if the number of coils is two or more times the number of threads and CPU_MT is selected.
         * @param pointVectors An array of radius vectors wrapped in class Vector3Array.
         * @param usedPrecision Custom precision settings used for this particular calculation.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Array of Cartesian vectors which represent vector potential A at specified points.
         */
        [[nodiscard]] vec3::Vector3Array computeAllAPotentialVectors(const vec3::Vector3Array &pointVectors,
                                                                     ComputeMethod computeMethod = CPU_ST) const;
        /**
         * @brief Calculates magnetic flux density B (magnetic field) from all member Coils for a number of
         * specified points. Uses internally defined defaultPrecisionFactor.
         * @details There are multiple compute methods, GPU acceleration is best suited for over 100 points.
         * MTD is used if the number of coils is two or more times the number of threads and CPU_MT is selected.
         * @param pointVectors An array of radius vectors wrapped in class Vector3Array.
         * @param usedPrecision Custom precision settings used for this particular calculation.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Array of Cartesian vectors which represent magnetic flux B at specified points.
         */
        [[nodiscard]] vec3::Vector3Array computeAllBFieldVectors(const vec3::Vector3Array &pointVectors,
                                                                 ComputeMethod computeMethod = CPU_ST) const;
        /**
         * @brief Calculates the amplitude vector of electric field E in sinusoidal steady-state from all member Coils
         * for a number of specified points. Uses internally defined defaultPrecisionFactor.
         * @details There are multiple compute methods, GPU acceleration is best suited for over 100 points.
         * MTD is used if the number of coils is two or more times the number of threads and CPU_MT is selected.
         * @param pointVectors An array of radius vectors wrapped in class Vector3Array.
         * @param usedPrecision Custom precision settings used for this particular calculation.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Array of Cartesian vectors which represent the amplitude of electric field E at specified points.
         */
        [[nodiscard]] vec3::Vector3Array computeAllEFieldVectors(const vec3::Vector3Array &pointVectors,
                                                                 ComputeMethod computeMethod = CPU_ST) const;
        /**
         * @brief Calculates the gradient G of the magnetic field (total derivative of B)
         * for a number of specified points. Uses internally defined defaultPrecisionFactor.
         * @details There are multiple compute methods, GPU acceleration is best suited for over 100 points.
         * MTD is used if the number of coils is two or more times the number of threads and CPU_MT is selected.
         * @param pointVectors An array of radius vectors wrapped in class Vector3Array.
         * @param usedPrecision Custom precision settings used for this particular calculation.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Array of 3x3 matrices which represents the magnetic gradient matrix G at specified points.
         */
        [[nodiscard]] vec3::Matrix3Array computeAllBGradientMatrices(const vec3::Vector3Array &pointVectors,
                                                                     ComputeMethod computeMethod = CPU_ST) const;

        /**
         * @brief Calculates the mutual inductance M between a provided coil and the rest of the member coils.
         * Uses separate CoilPairArguments for every pair of coils, generated with the given PrecisionFactor.
         * @details Using CPU_MT compute method is highly advisable, especially for higher precision factors.
         * GPU is a good option when an error of 1e-5 is good enough for the application (usual error of order 1e-6).
         * MTD is used if the number of coils is two or more times the number of threads and CPU_MT is selected.
         * @param secondary The coil represented with a number of points for which the potential is calculated.
         * Can be a member of CoilGroup or another Coil.
         * @param precisionFactor Determines the precision of given calculation.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Mutual inductance between the system and secondary coil.
         */
        [[nodiscard]] double computeMutualInductance(const Coil &secondary,
                                                     PrecisionFactor precisionFactor = PrecisionFactor(),
                                                     ComputeMethod computeMethod = CPU_ST) const;
        /**
         * @brief Calculates the force F and torque T on a provided coil from the rest of the member coils.
         * Uses separate CoilPairArguments for every pair of coils, generated with the given PrecisionFactor.
         * @details Using CPU_MT compute method is highly advisable, especially for higher precision factors.
         * GPU is a good option when an error of 1e-5 is good enough for the application (usual error of order 1e-6).
         * MTD is used if the number of coils is two or more times the number of threads and CPU_MT is selected.
         * @param secondary The coil represented with a number of points for which the magnetic field is calculated.
         * Can be a member of CoilGroup or another Coil.
         * @param precisionFactor Determines the precision of given calculation.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Pair of Cartesian vectors which represent force (first) and torque (second).
         */
        [[nodiscard]] std::pair<vec3::Vector3, vec3::Vector3>
        computeForceTorque(const Coil &secondary, PrecisionFactor precisionFactor = PrecisionFactor(),
                           ComputeMethod computeMethod = CPU_ST) const;
        /**
         * @brief Calculates force F and torque T between a coil and magnetostatic object with a dipole moment.
         * Uses precision internally defined by defaultPrecisionCPU.
         * @details This method can prove particularly useful for approximating the force and and torque
         * of CoilGroup on a coil which is sufficiently far apart, or when the coil is very small.
         * It can also be useful in particle simulations where the magnetic dipole moment is not negligible
         * @param pointVector Radius vector from the origin to the point where the magnetic dipole is located.
         * @param dipoleMoment Magnetic dipole moment vector of a secondary coil or another object (magnet, particle).
         * @return Pair of Cartesian vectors which represent force (first) and torque (second).
         */
        [[nodiscard]] std::pair<vec3::Vector3, vec3::Vector3>
        computeForceOnDipoleMoment(vec3::Vector3 pointVector, vec3::Vector3 dipoleMoment) const;

        /**
         * @brief Calculates the mutual inductance M between a provided coil and the rest of the member coils for
         * multiple positions and orientations of the secondary coil.
         * @details Uses separate CoilPairArguments for every pair of coils, generated with the given PrecisionFactor,
         * when using the CPU. When using the GPU every Coil is assigned a calculated precision factor which makes
         * the precision equivalent, in terms of total increments, to given PrecisionFactor.
         * This method is exceptionally powerful because it can more efficiently utilise the CPU with
         * MTD, and especially the GPU with a special pure GPU implementation of mutual inductance calculation.
         * @param secondary The coil that is represented with a number of points for which the potential is calculated.
         * @param secondaryPositions Positions of the center of the secondary coil.
         * @param secondaryYAngles Secondary coil rotation angles along the y-axis.
         * @param secondaryZAngles Secondary coil rotation angles along the z-axis.
         * @param precisionFactor Determines the precision of given calculation.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Array of mutual inductance values, one for each appropriate configuration.
         */
        [[nodiscard]] std::vector<double>
        computeAllMutualInductanceArrangements(const Coil &secondary,
                                               const vec3::Vector3Array &secondaryPositions,
                                               const std::vector<double> &secondaryYAngles,
                                               const std::vector<double> &secondaryZAngles,
                                               PrecisionFactor precisionFactor = PrecisionFactor(),
                                               ComputeMethod computeMethod = CPU_ST) const;
        /**
         * @brief Calculates force F and torque T on a provided coil from the rest of the member coils for
         * multiple positions and orientations of the secondary coil.
         * @details Uses separate CoilPairArguments for every pair of coils, generated with the given PrecisionFactor,
         * when using the CPU. When using the GPU every Coil is assigned a calculated precision factor which makes
         * the precision equivalent, in terms of total increments, to given PrecisionFactor.
         * This method is exceptionally powerful because it can more efficiently utilise the CPU with
         * MTD, and especially the GPU with a special pure GPU implementation of mutual inductance calculation.
         * @param secondary The coil that is represented with a number of points for which the magnetic field is calculated.
         * @param secondaryPositions Positions of the center of the secondary coil.
         * @param secondaryYAngles Secondary coil rotation angles along the y-axis.
         * @param secondaryZAngles Secondary coil rotation angles along the z-axis.
         * @param precisionFactor Determines the precision of given calculation.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * Array of pairs of force (first) and torque (second) vectors, one for each appropriate configuration.
         */
        [[nodiscard]] std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
        computeAllForceTorqueArrangements(const Coil &secondary,
                                          const vec3::Vector3Array &secondaryPositions,
                                          const std::vector<double> &secondaryYAngles,
                                          const std::vector<double> &secondaryZAngles,
                                          PrecisionFactor precisionFactor = PrecisionFactor(),
                                          ComputeMethod computeMethod = CPU_ST) const;

        /// @brief Generates a string object with all properties of the CoilGroup instance (and appropriate member Coils).
        explicit operator std::string() const;

    private:
        [[nodiscard]] vec3::Vector3Array calculateAllAPotentialMT(const vec3::Vector3Array &pointVectors) const;
        [[nodiscard]] vec3::Vector3Array calculateAllBFieldMT(const vec3::Vector3Array &pointVectors) const;
        [[nodiscard]] vec3::Vector3Array calculateAllEFieldMT(const vec3::Vector3Array &pointVectors) const;
        [[nodiscard]] vec3::Matrix3Array calculateAllBGradientMT(const vec3::Vector3Array &pointVectors) const;

        [[nodiscard]] vec3::Vector3Array calculateAllAPotentialMTD(const vec3::Vector3Array &pointVectors) const;
        [[nodiscard]] vec3::Vector3Array calculateAllBFieldMTD(const vec3::Vector3Array &pointVectors) const;
        [[nodiscard]] vec3::Vector3Array calculateAllEFieldMTD(const vec3::Vector3Array &pointVectors) const;
        [[nodiscard]] vec3::Matrix3Array calculateAllBGradientMTD(const vec3::Vector3Array &pointVectors) const;

        void generateCoilDataArray(CoilData *coilDataArr, PrecisionFactor precisionFactor,
                                   bool removeSpecificCoil = false, unsigned long long int specificCoilId = 0) const;
        static void generateSecondaryData(const Coil &secondary, SecondaryCoilData &secondaryData,
                                          const PrecisionArguments &precision, bool forceCalculation);

        [[nodiscard]] vec3::Vector3Array calculateAllAPotentialGPU(const vec3::Vector3Array &pointVectors) const;
        [[nodiscard]] vec3::Vector3Array calculateAllBFieldGPU(const vec3::Vector3Array &pointVectors) const;
        [[nodiscard]] vec3::Vector3Array calculateAllEFieldGPU(const vec3::Vector3Array &pointVectors) const;
        [[nodiscard]] vec3::Matrix3Array calculateAllBGradientGPU(const vec3::Vector3Array &pointVectors) const;


        [[nodiscard]] double calculateMutualInductanceMTD(const Coil &secondary,
                                                          PrecisionFactor precisionFactor = PrecisionFactor()) const;
        [[nodiscard]] std::pair<vec3::Vector3, vec3::Vector3>
        calculateForceTorqueMTD(const Coil &secondary, PrecisionFactor precisionFactor = PrecisionFactor()) const;

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
        calculateAllForceTorqueArrangementsMTD(const Coil &secondary,
                                               const vec3::Vector3Array &secondaryPositions,
                                               const std::vector<double> &secondaryYAngles,
                                               const std::vector<double> &secondaryZAngles,
                                               PrecisionFactor precisionFactor = PrecisionFactor()) const;

        [[nodiscard]] std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
        calculateAllForceTorqueArrangementsGPU(const Coil &secondary,
                                               const vec3::Vector3Array &secondaryPositions,
                                               const std::vector<double> &secondaryYAngles,
                                               const std::vector<double> &secondaryZAngles,
                                               PrecisionFactor precisionFactor = PrecisionFactor()) const;

};

#endif //GENERAL_COIL_PROGRAM_COILGROUP_H
