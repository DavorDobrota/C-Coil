#ifndef GENERAL_COIL_PROGRAM_COIL_H
#define GENERAL_COIL_PROGRAM_COIL_H

#include "Coil/EnumsAndConstants/ComputeMethod.h"
#include "Coil/EnumsAndConstants/CoilType.h"
#include "Coil/EnumsAndConstants/PrecisionGlobalVars.h"

#include "Tensor.h"
#include "CUDAUtils/ConstantsAndStructs/CoilDataStructs.h"

#include <vector>
#include <string>

#define PRINT_PRECISION 0
#define PRINT_TIMINGS 0


const int g_defaultThreadCount = 8;


class Coil;


struct PrecisionFactor
{
    PrecisionFactor();
    explicit PrecisionFactor(double relativePrecision);

    double relativePrecision;

    explicit operator std::string() const;
};

struct PrecisionArguments
{
    PrecisionArguments();
    explicit PrecisionArguments(int angularBlocks, int thicknessBlocks, int lengthBlocks,
                                int angularIncrements, int thicknessIncrements, int lengthIncrements);

    int angularBlocks;
    int thicknessBlocks;
    int lengthBlocks;

    int angularIncrements;
    int thicknessIncrements;
    int lengthIncrements;

    static PrecisionArguments getCoilPrecisionArgumentsCPU(const Coil &coil, PrecisionFactor precisionFactor);
    static PrecisionArguments getCoilPrecisionArgumentsGPU(const Coil &coil, PrecisionFactor precisionFactor);

    static PrecisionArguments getSecondaryCoilPrecisionArgumentsGPU(const Coil &coil, PrecisionFactor precisionFactor);

    explicit operator std::string() const;

    private:

        static PrecisionArguments calculatePrecisionArguments(const Coil &coil, PrecisionFactor precisionFactor,
                                                              bool useGPU = false);
};

struct CoilPairArguments
{
    friend PrecisionArguments;

    CoilPairArguments();
    explicit CoilPairArguments(const PrecisionArguments &primaryPrecision,
                               const PrecisionArguments &secondaryPrecision);

    PrecisionArguments primaryPrecision;
    PrecisionArguments secondaryPrecision;

    static CoilPairArguments
    getAppropriateCoilPairArguments(const Coil &primary, const Coil &secondary, PrecisionFactor precisionFactor,
                                    ComputeMethod computeMethod = CPU_ST, bool zAxisCase = false, bool pureGPU = false);

    explicit operator std::string() const;

    private:

        static std::vector<std::pair<int, int>> balanceIncrements(int totalIncrements,
                                                                  const std::vector<std::pair<int, double>> &components);
};


class CoilGroup;

/**
 * @brief Primary class in this project.
 * Models a circular coil with a rectangular cross section and uniform current density.
 * @details The model works best for a solid block of material and when the effects of windings are negligible.
 * Primary attributes are length (b), thickness (a), the inner radius (R) of the internal cylindrical hole,
 * and the number of turns of wire (windings).
 */
class Coil
{
    friend CoilGroup;

    private:
        unsigned long long id;

        double innerRadius;
        double thickness;
        double length;
        int numOfTurns;

        double currentDensity{};
        double current{};

        double wireResistivity{};
        bool sineDriven{};
        double sineFrequency{};

        double magneticMoment{};
        double averageWireThickness{};

        double resistance{};
        double selfInductance{};
        double reactance{};
        double impedance{};

        CoilType coilType{};
        bool useFastMethod{};
        int threadCount{};

        PrecisionArguments defaultPrecisionCPU;
        PrecisionArguments defaultPrecisionGPU;

        vec3::Vector3 positionVector{};
        double yAxisAngle{};
        double zAxisAngle{};
        vec3::Matrix3 transformationMatrix{};
        vec3::Matrix3 inverseTransformationMatrix{};

    public:

        /// @name CoilConstructors
        ///@{
        /**
         *
         */
        Coil();

        Coil(double innerRadius, double thickness, double length, int numOfTurns,
             double current, double wireResistivity, double sineFrequency,
             PrecisionFactor precisionFactor = PrecisionFactor(), int threadCount = g_defaultThreadCount,
             vec3::Vector3 coordinatePosition = vec3::Vector3(), double yAxisAngle = 0.0, double zAxisAngle = 0.0);
        Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
             double wireResistivity, double sineFrequency,
             const PrecisionArguments &precisionSettingsCPU, const PrecisionArguments &precisionSettingsGPU,
             int threadCount = g_defaultThreadCount, vec3::Vector3 coordinatePosition = vec3::Vector3(),
             double yAxisAngle = 0.0, double zAxisAngle = 0.0);

        Coil(double innerRadius, double thickness, double length, int numOfTurns, double current, double sineFrequency,
             PrecisionFactor precisionFactor = PrecisionFactor(), int threadCount = g_defaultThreadCount,
             vec3::Vector3 coordinatePosition = vec3::Vector3(), double yAxisAngle = 0.0, double zAxisAngle = 0.0);
        Coil(double innerRadius, double thickness, double length, int numOfTurns, double current, double sineFrequency,
             const PrecisionArguments &precisionSettingsCPU, const PrecisionArguments &precisionSettingsGPU,
             int threadCount = g_defaultThreadCount, vec3::Vector3 coordinatePosition = vec3::Vector3(),
             double yAxisAngle = 0.0, double zAxisAngle = 0.0);

        Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
             PrecisionFactor precisionFactor = PrecisionFactor(), int threadCount = g_defaultThreadCount,
             vec3::Vector3 coordinatePosition = vec3::Vector3(), double yAxisAngle = 0.0, double zAxisAngle = 0.0);
        Coil(double innerRadius, double thickness, double length, int numOfTurns, double current,
             const PrecisionArguments &precisionSettingsCPU, const PrecisionArguments &precisionSettingsGPU,
             int threadCount = g_defaultThreadCount, vec3::Vector3 coordinatePosition = vec3::Vector3(),
             double yAxisAngle = 0.0, double zAxisAngle = 0.0);

        Coil(double innerRadius, double thickness, double length, int numOfTurns,
             PrecisionFactor precisionFactor = PrecisionFactor(), int threadCount = g_defaultThreadCount,
             vec3::Vector3 coordinatePosition = vec3::Vector3(), double yAxisAngle = 0.0, double zAxisAngle = 0.0);
        Coil(double innerRadius, double thickness, double length, int numOfTurns,
             const PrecisionArguments &precisionSettingsCPU, const PrecisionArguments &precisionSettingsGPU,
             int threadCount = g_defaultThreadCount, vec3::Vector3 coordinatePosition = vec3::Vector3(),
             double yAxisAngle = 0.0, double zAxisAngle = 0.0);
        /// @}

        [[nodiscard]] unsigned long long getId() const;
        [[nodiscard]] double getInnerRadius() const;
        [[nodiscard]] double getThickness() const;
        [[nodiscard]] double getLength() const;
        [[nodiscard]] int getNumOfTurns() const;

        [[nodiscard]] double getCurrentDensity() const;
        [[nodiscard]] double getCurrent() const;

        [[nodiscard]] double getWireResistivity() const;
        [[nodiscard]] bool isSineDriven() const;
        [[nodiscard]] double getSineFrequency() const;

        [[nodiscard]] vec3::Vector3 getMagneticMoment();
        [[nodiscard]] double getAverageWireThickness() const;

        [[nodiscard]] double getSelfInductance() const;
        [[nodiscard]] double getResistance();
        [[nodiscard]] double getReactance();
        [[nodiscard]] double getImpedance();

        [[nodiscard]] const PrecisionArguments &getPrecisionSettingsCPU() const;
        [[nodiscard]] const PrecisionArguments &getPrecisionSettingsGPU() const;
        [[nodiscard]] int getThreadCount() const;
        [[nodiscard]] bool isUsingFastMethod() const;
        [[nodiscard]] CoilType getCoilType() const;

        [[nodiscard]] vec3::Vector3 getPositionVector() const;
        [[nodiscard]] std::pair<double, double> getRotationAngles() const;

        [[nodiscard]] vec3::Matrix3 getTransformationMatrix() const;
        [[nodiscard]] vec3::Matrix3 getInverseTransformationMatrix() const;

        void setCurrentDensity(double currentDensity);
        void setCurrent(double current);
        void setWireResistivity(double wireResistivity);
        void setSineFrequency(double sineFrequency);

        void setDefaultPrecisionCPU(const PrecisionArguments &precisionSettings);
        void setDefaultPrecisionCPU(PrecisionFactor precisionFactor = PrecisionFactor());
        void setDefaultPrecisionGPU(const PrecisionArguments &precisionSettings);
        void setDefaultPrecisionGPU(PrecisionFactor precisionFactor = PrecisionFactor());
        void setDefaultPrecision(PrecisionFactor precisionFactor = PrecisionFactor());

        void setThreadCount(int threadCount);
        void setPositionAndOrientation(vec3::Vector3 positionVector = vec3::Vector3(),
                                       double yAxisAngle = 0.0, double zAxisAngle = 0.0);

        void setSelfInductance(double selfInductance);


        /**
         * @brief Calculates vector potential A of the magnetic field at the specified point.
         * Uses precision internally defined by defaultPrecisionCPU.
         * @param pointVector Radius vector from the origin to the point where the field is calculated.
         * @return 3D Cartesian vector which represents vector potential A.
         */
        [[nodiscard]] vec3::Vector3 computeAPotentialVector(vec3::Vector3 pointVector) const;
        /**
         * @brief Calculates vector potential A of the magnetic field at the specified point.
         * Uses provided PrecisionArguments for precision settings.
         * @param pointVector Radius vector from the origin to the point where the field is calculated.
         * @param usedPrecision Custom precision settings used for this particular calculation.
         * @return 3D Cartesian vector which represents vector potential A.
         */

        [[nodiscard]] vec3::Vector3 computeAPotentialVector(vec3::Vector3 pointVector,
                                                            const PrecisionArguments &usedPrecision) const;
        /**
         * @brief Calculates magnetic flux density B (magnetic field) at the specified point.
         * Uses precision internally defined by defaultPrecisionCPU.
         * @param pointVector Radius vector from the origin to the point where the field is calculated.
         * @return 3D Cartesian vector which represents magnetic flux density B.
         */
        [[nodiscard]] vec3::Vector3 computeBFieldVector(vec3::Vector3 pointVector) const;
        /**
         * @brief Calculates magnetic flux density B (magnetic field) at the specified point.
         * Uses provided PrecisionArguments for precision settings.
         * @param pointVector Radius vector from the origin to the point where the field is calculated.
         * @param usedPrecision Custom precision settings used for this particular calculation.
         * @return 3D Cartesian vector which represents magnetic flux density B.
         */

        [[nodiscard]] vec3::Vector3 computeBFieldVector(vec3::Vector3 pointVector,
                                                        const PrecisionArguments &usedPrecision) const;
        /**
         * @brief Calculates the amplitude vector of electric field E in sinusoidal steady-state at the specified point.
         * Uses precision internally defined by defaultPrecisionCPU.
         * @param pointVector Radius vector from the origin to the point where the field is calculated.
         * @return 3D Cartesian vector which represents the amplitude of electric field E.
         */
        [[nodiscard]] vec3::Vector3 computeEFieldVector(vec3::Vector3 pointVector) const;

        /**
         * @brief Calculates the amplitude vector of electric field E in sinusoidal steady-state at the specified point.
         * Uses provided PrecisionArguments for precision settings.
         * @param pointVector Radius vector from the origin to the point where the field is calculated.
         * @param usedPrecision Custom precision settings used for this particular calculation.
         * @return 3D Cartesian vector which represents the amplitude of electric field E.
         */
        [[nodiscard]] vec3::Vector3 computeEFieldVector(vec3::Vector3 pointVector,
                                                        const PrecisionArguments &usedPrecision) const;

        /**
         * @brief Calculates the gradient G of the magnetic field (total derivative of B) at the specified point.
         * Uses precision internally defined by defaultPrecisionCPU.
         * @param pointVector Radius vector from the origin to the point where the field is calculated.
         * @return 3x3 Matrix which represents the magnetic gradient matrix G.
         */
        [[nodiscard]] vec3::Matrix3 computeBGradientMatrix(vec3::Vector3 pointVector) const;
        /**
         * @brief Calculates the gradient G of the magnetic field (total derivative of B) at the specified point.
         * Uses provided PrecisionArguments for precision settings.
         * @param pointVector Radius vector from the origin to the point where the field is calculated.
         * @param usedPrecision Custom precision settings used for this particular calculation.
         * @return 3x3 matrix which represents the magnetic gradient matrix G.
         */
        [[nodiscard]] vec3::Matrix3 computeBGradientMatrix(vec3::Vector3 pointVector,
                                                           const PrecisionArguments &usedPrecision) const;

        /**
         * @brief Calculates vector potential A of the magnetic field for a number of specified points.
         * Uses precision internally defined by defaultPrecisionCPU.
         * @details There are multiple compute methods, GPU acceleration is best suited for a large number of points.
         * @param pointVectors An array of radius vectors wrapped in class Vector3Array
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Array of Cartesian vectors which represent vector potential A at specified points.
         */
        [[nodiscard]] vec3::Vector3Array computeAllAPotentialVectors(const vec3::Vector3Array &pointVectors,
                                                                     ComputeMethod computeMethod = CPU_ST) const;
        /**
         * @brief Calculates vector potential A of the magnetic field for a number of specified points.
         * Uses provided PrecisionArguments for precision settings.
         * @details There are multiple compute methods, GPU acceleration is best suited for a large number of points.
         * @param pointVectors An array of radius vectors wrapped in class Vector3Array.
         * @param usedPrecision Custom precision settings used for this particular calculation.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Array of Cartesian vectors which represent vector potential A at specified points.
         */
        [[nodiscard]] vec3::Vector3Array computeAllAPotentialVectors(const vec3::Vector3Array &pointVectors,
                                                                     const PrecisionArguments &usedPrecision,
                                                                     ComputeMethod computeMethod = CPU_ST) const;

        /**
         * @brief Calculates magnetic flux density B (magnetic field) for a number of specified points.
         * Uses precision internally defined by defaultPrecisionCPU.
         * @details There are multiple compute methods, GPU acceleration is best suited for a large number of points.
         * @param pointVectors An array of radius vectors wrapped in class Vector3Array.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Array of Cartesian vectors which represent magnetic flux B at specified points.
         */
        [[nodiscard]] vec3::Vector3Array computeAllBFieldVectors(const vec3::Vector3Array &pointVectors,
                                                                 ComputeMethod computeMethod = CPU_ST) const;
        /**
         * @brief Calculates magnetic flux density B (magnetic field) for a number of specified points.
         * Uses provided PrecisionArguments for precision settings.
         * @details There are multiple compute methods, GPU acceleration is best suited for a large number of points.
         * @param pointVectors An array of radius vectors wrapped in class Vector3Array.
         * @param usedPrecision Custom precision settings used for this particular calculation.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Array of Cartesian vectors which represent magnetic flux B at specified points.
         */
        [[nodiscard]] vec3::Vector3Array computeAllBFieldVectors(const vec3::Vector3Array &pointVectors,
                                                                 const PrecisionArguments &usedPrecision,
                                                                 ComputeMethod computeMethod = CPU_ST) const;

        /**
         * @brief Calculates the amplitude vector of electric field E in sinusoidal steady-state for a number of specified points.
         * Uses precision internally defined by defaultPrecisionCPU.
         * @details There are multiple compute methods, GPU acceleration is best suited for a large number of points.
         * @param pointVectors An array of radius vectors wrapped in class Vector3Array.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Array of Cartesian vectors which represent the amplitude of electric field E at specified points.
         */
        [[nodiscard]] vec3::Vector3Array computeAllEFieldVectors(const vec3::Vector3Array &pointVectors,
                                                                 ComputeMethod computeMethod = CPU_ST) const;
        /**
         * @brief Calculates the amplitude vector of electric field E in sinusoidal steady-state for a number of specified points.
         * Uses provided PrecisionArguments for precision settings.
         * @details There are multiple compute methods, GPU acceleration is best suited for a large number of points.
         * @param pointVectors An array of radius vectors wrapped in class Vector3Array.
         * @param usedPrecision Custom precision settings used for this particular calculation.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Array of Cartesian vectors which represent the amplitude of electric field E at specified points.
         */
        [[nodiscard]] vec3::Vector3Array computeAllEFieldVectors(const vec3::Vector3Array &pointVectors,
                                                                 const PrecisionArguments &usedPrecision,
                                                                 ComputeMethod computeMethod = CPU_ST) const;

        /**
         * @brief Calculates the gradient G of the magnetic field (total derivative of B) for a number of specified points.
         * Uses precision internally defined by defaultPrecisionCPU.
         * @details There are multiple compute methods, GPU acceleration is best suited for a large number of points.
         * @param pointVectors An array of radius vectors wrapped in class Vector3Array.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Array of 3x3 matrices which represents the magnetic gradient matrix G at specified points.
         */
        [[nodiscard]] vec3::Matrix3Array computeAllBGradientMatrices(const vec3::Vector3Array &pointVectors,
                                                                     ComputeMethod computeMethod = CPU_ST) const;
        /**
         * @brief Calculates the gradient G of the magnetic field (total derivative of B) for a number of specified points.
         * Uses provided PrecisionArguments for precision settings.
         * @details There are multiple compute methods, GPU acceleration is best suited for a large number of points.
         * @param pointVectors An array of radius vectors wrapped in class Vector3Array.
         * @param usedPrecision Custom precision settings used for this particular calculation.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Array of 3x3 matrices which represents the magnetic gradient matrix G at specified points.
         */
        [[nodiscard]] vec3::Matrix3Array computeAllBGradientMatrices(const vec3::Vector3Array &pointVectors,
                                                                     const PrecisionArguments &usedPrecision,
                                                                     ComputeMethod computeMethod = CPU_ST) const;

        /**
         * @brief Calculates the Mutual induction M between two given coils.
         * Generates CoilPairArguments from given precisionFactor.
         * @details For better precision, the primary coil should be the bigger one, length is the most important parameter.
         * There are more performant implementations if both coils lie on the z-axis and have rotation angles set to 0.
         * Using CPU_MT compute method is highly advisable, especially for higher precision factors.
         * GPU is a good option when an error of 1e-5 is good enough for the application (usual error of order 1e-6).
         * @param primary The coil that generates the vector potential.
         * @param secondary The coil that is represented with a number of points for which the potential is calculated.
         * @param precisionFactor Determines the precision of given calculation, increments = baseIncrements * 2^precision.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Mutual inductance of the system of two coils.
         */
        static double computeMutualInductance(const Coil &primary, const Coil &secondary,
                                              PrecisionFactor precisionFactor = PrecisionFactor(),
                                              ComputeMethod computeMethod = CPU_ST);
        /**
         * @brief Calculates the Mutual induction M between two given coils.
         * Uses provided CoilPairArguments for precision settings.
         * @details For better precision, the primary coil should be the bigger one, length is the most important parameter.
         * There are more performant implementations if both coils lie on the z-axis and have rotation angles set to 0.
         * Using CPU_MT compute method is highly advisable, especially for higher precision factors.
         * GPU is a good option when an error of 1e-5 is good enough for the application (usual error of order 1e-6).
         * @param primary The coil that generates the vector potential.
         * @param secondary The coil that is represented with a number of points for which the potential is calculated.
         * @param inductanceArguments Custom precision settings used for this particular calculation.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Mutual inductance of the system of two coils.
         */
        static double computeMutualInductance(const Coil &primary, const Coil &secondary,
                                              const CoilPairArguments &inductanceArguments,
                                              ComputeMethod computeMethod = CPU_ST);

        /**
         * @brief Calculates the magnitude of sinusoidal steady-state voltage induced in the secondary coil.
         * Generates CoilPairArguments from given precisionFactor. Similar to computeMutualInductance.
         * @param secondary The coil that is represented with a number of points for which the E field is calculated.
         * @param precisionFactor Determines the precision of given calculation, increments = baseIncrements * 2^precision.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Magnitude of the voltage induced on the secondary coil.
         */
        [[nodiscard]] double computeSecondaryInducedVoltage(const Coil &secondary,
                                                            PrecisionFactor precisionFactor = PrecisionFactor(),
                                                            ComputeMethod computeMethod = CPU_ST) const;
        /**
         * @brief Calculates the magnitude of sinusoidal steady-state voltage induced in the secondary coil.
         * Uses provided CoilPairArguments for precision settings. Similar to computeMutualInductance.
         * @param secondary The coil that is represented with a number of points for which the E field is calculated.
         * @param inductanceArguments Custom precision settings used for this particular calculation.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Magnitude of the voltage induced on the secondary coil.
         */
        [[nodiscard]] double computeSecondaryInducedVoltage(const Coil &secondary,
                                                            const CoilPairArguments &inductanceArguments,
                                                            ComputeMethod computeMethod = CPU_ST) const;
        /**
         * @brief Special method which returns self inductance L of the given coil and sets it internally.
         * @details This method is exclusively single threaded and represents a shortcoming of this approach.
         * Low precision factors, below 5.0, are not advisable and good precision (error of order 1e-6)
         * can be achieved with precision factor 10.0. It works well for thick and thin coils,
         * but poorly for flat coils, and does not work for filaments (loops) as the integral is inherently divergent.
         * @param precisionFactor
         * @return Magnitude of the voltage induced on the secondary coil.
         */
        double computeAndSetSelfInductance(PrecisionFactor precisionFactor);

        /**
        * @brief Calculates the Ampere force F and torque T between two coils.
        * Generates CoilPairArguments from given precisionFactor.
        * @details For better precision, the primary coil should be the bigger one, length is the most important parameter.
        * There are more performant implementations if both coils lie on the z-axis and have rotation angles set to 0.
        * Using CPU_MT compute method is highly advisable, especially for higher precision factors and
        * GPU is a good option when an error of 1e-5 is good enough for the application (usual error of order 1e-6).
        * @param primary The coil that generates the magnetic field.
        * @param secondary The coil that is represented with a number of points for which the field is calculated.
        * @param precisionFactor Determines the precision of given calculation, increments = baseIncrements * 2^precision.
        * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
        * @return Pair of Cartesian vectors which represent force (first) and torque (second).
        */
        static std::pair<vec3::Vector3, vec3::Vector3>
        computeAmpereForce(const Coil &primary, const Coil &secondary,
                           PrecisionFactor precisionFactor = PrecisionFactor(), ComputeMethod computeMethod = CPU_ST);
        /**
        * @brief Calculates the Ampere force F and torque T between two coils.
        * Generates CoilPairArguments from given precisionFactor
        * @details For better precision, the primary coil should be the bigger one, length is the most important parameter.
        * There are more performant implementations if both coils lie on the z-axis and have rotation angles set to 0.
        * Using CPU_MT compute method is highly advisable, especially for higher precision factors and
        * GPU is a good option when an error of 1e-5 is good enough for the application (usual error of order 1e-6).
        * @param primary The coil that generates the magnetic field.
        * @param secondary The coil that is represented with a number of points for which the field is calculated.
        * @param inductanceArguments Custom precision settings used for this particular calculation.
        * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
        * @return Pair of Cartesian vectors which represent force (first) and torque (second).
        */
        static std::pair<vec3::Vector3, vec3::Vector3>
        computeAmpereForce(const Coil &primary, const Coil &secondary,
                           const CoilPairArguments &forceArguments, ComputeMethod computeMethod = CPU_ST);

        /**
         * @brief Calculates force F and torque T between a coil and magnetostatic object with a dipole moment.
         * Uses precision internally defined by defaultPrecisionCPU.
         * @details This method can prove particularly useful for approximating the force and and torque
         * between two coils which are sufficiently far apart, or when the secondary coil is very small.
         * It can also be useful in particle simulations where the magnetic dipole moment is not negligible
         * @param pointVector Radius vector from the origin to the point where the magnetic dipole is located.
         * @param dipoleMoment Magnetic dipole moment vector of a secondary coil or another object (magnet, particle).
         * @return Pair of Cartesian vectors which represent force (first) and torque (second).
         */
        [[nodiscard]] std::pair<vec3::Vector3, vec3::Vector3>
        computeForceOnDipoleMoment(vec3::Vector3 pointVector, vec3::Vector3 dipoleMoment) const;
        /**
         * @brief Calculates force F and torque T between a coil and magnetostatic object with a dipole moment.
         * Uses provided PrecisionArguments for precision settings.
         * @details This method can prove particularly useful for approximating the force and and torque
         * between two coils which are sufficiently far apart, or when the secondary coil is very small.
         * It can also be useful in particle simulations where the magnetic dipole moment is not negligible
         * @param pointVector Radius vector from the origin to the point where the magnetic dipole is located.
         * @param dipoleMoment Magnetic dipole moment vector of a secondary coil or another object (magnet, particle).
         * @param usedPrecision Custom precision settings used for this particular calculation.
         * @return Pair of Cartesian vectors which represent force (first) and torque (second).
         */
        [[nodiscard]] std::pair<vec3::Vector3, vec3::Vector3>
        computeForceOnDipoleMoment(vec3::Vector3 pointVector, vec3::Vector3 dipoleMoment,
                                   const PrecisionArguments &usedPrecision) const;

        /**
         * @brief Calculates mutual inductance M between two coils for different coil configurations.
         * All positional arguments can be changed. Generates CoilPairArguments from given precisionFactor.
         * @details This method is exceptionally powerful because it can more efficiently utilise the CPU with
         * distributed (coarse-grained) multithreading, and especially the GPU with a special pure GPU implementation
         * of mutual inductance calculation. Computation times can be as low as several microseconds per configuration.
         * @param primary The coil that generates the vector potential
         * @param secondary The coil that is represented with a number of points for which the potential is calculated.
         * @param primaryPositions Positions of the center of the primary coil.
         * @param secondaryPositions Positions of the center of the secondary coil.
         * @param primaryYAngles Primary coil rotation angles along the y-axis.
         * @param primaryZAngles Primary coil rotation angles along the z-axis.
         * @param secondaryYAngles Secondary coil rotation angles along the y-axis.
         * @param secondaryZAngles Secondary coil rotation angles along the z-axis.
         * @param precisionFactor Determines the precision of given calculation, increments = baseIncrements * 2^precision.
         * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
         * @return Array of mutual inductance values, one for each appropriate configuration.
         */
        static std::vector<double>
        computeAllMutualInductanceArrangements(const Coil &primary, const Coil &secondary,
                                               const vec3::Vector3Array &primaryPositions,
                                               const vec3::Vector3Array &secondaryPositions,
                                               const std::vector<double> &primaryYAngles,
                                               const std::vector<double> &primaryZAngles,
                                               const std::vector<double> &secondaryYAngles,
                                               const std::vector<double> &secondaryZAngles,
                                               PrecisionFactor precisionFactor = PrecisionFactor(),
                                               ComputeMethod computeMethod = CPU_ST);
        /**
        * @brief Calculates force F and torque T between two coils for different coil configurations.
        * All positional arguments can be changed. Generates CoilPairArguments from given precisionFactor.
        * @details This method is exceptionally powerful because it can more efficiently utilise the CPU with
        * distributed (coarse-grained) multithreading, and especially the GPU with a special pure GPU implementation
        * of force and torque calculation. Computation times can be as low as several microseconds per configuration.
        * @param primary The coil that generates the magnetic field
        * @param secondary The coil that is represented with a number of points for which the field is calculated.
        * @param primaryPositions Positions of the center of the primary coil.
        * @param secondaryPositions Positions of the center of the secondary coil.
        * @param primaryYAngles Primary coil rotation angles along the y-axis.
        * @param primaryZAngles Primary coil rotation angles along the z-axis.
        * @param secondaryYAngles Secondary coil rotation angles along the y-axis.
        * @param secondaryZAngles Secondary coil rotation angles along the z-axis.
        * @param precisionFactor Determines the precision of given calculation, increments = baseIncrements * 2^precision.
        * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
        * @return Array of pairs of force (first) and torque (second) vectors, one for each appropriate configuration.
        */
        static std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
        computeAllAmpereForceArrangements(const Coil &primary, const Coil &secondary,
                                          const vec3::Vector3Array &primaryPositions,
                                          const vec3::Vector3Array &secondaryPositions,
                                          const std::vector<double> &primaryYAngles,
                                          const std::vector<double> &primaryZAngles,
                                          const std::vector<double> &secondaryYAngles,
                                          const std::vector<double> &secondaryZAngles,
                                          PrecisionFactor precisionFactor = PrecisionFactor(),
                                          ComputeMethod computeMethod = CPU_ST);
        /**
         * @brief Generates a string object with all properties of the coil.
         * @return Coil string
         */
        explicit operator std::string() const;

    private:
        Coil (const Coil &original);

        void calculateMagneticMoment();
        void calculateAverageWireThickness();
        void calculateResistance();
        void calculateReactance();
        void calculateImpedance();
        void calculateCoilType();
        void calculateTransformationMatrices();


        [[nodiscard]] vec3::Vector3 calculateAPotential(vec3::Vector3 pointVector,
                                                        const PrecisionArguments &usedPrecision) const;
        [[nodiscard]] double calculateAPotentialSlow(double zCoord, double rCoord,
                                                     const PrecisionArguments &usedPrecision) const;
        [[nodiscard]] double calculateAPotentialFast(double zCoord, double rCoord,
                                                     const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] vec3::Vector3 calculateBField(vec3::Vector3 pointVector,
                                                    const PrecisionArguments &usedPrecision) const;
        [[nodiscard]] std::pair<double, double> calculateBFieldSlow(double zCoord, double rCoord,
                                                                    const PrecisionArguments &usedPrecision) const;
        [[nodiscard]] std::pair<double, double> calculateBFieldFast(double zCoord, double rCoord,
                                                                    const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] vec3::Matrix3 calculateBGradient(vec3::Vector3 pointVector,
                                                             const PrecisionArguments &usedPrecision) const;
        [[nodiscard]] std::vector<double> calculateBGradientSlow(double zCoord, double rCoord,
                                                                 const PrecisionArguments &usedPrecision) const;
        [[nodiscard]] std::vector<double> calculateBGradientFast(double zCoord, double rCoord,
                                                                 const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] static std::vector<size_t> calculateChunkSize(size_t opCount, int threads);

        [[nodiscard]] vec3::Vector3Array calculateAllBFieldMT(const vec3::Vector3Array &pointVectors,
                                                              const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] vec3::Vector3Array calculateAllAPotentialMT(const vec3::Vector3Array &pointVectors,
                                                                  const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] vec3::Matrix3Array calculateAllBGradientMT(const vec3::Vector3Array &pointVectors,
                                                                 const PrecisionArguments &usedPrecision) const;

        void generateCoilData(CoilData &coilData, const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] vec3::Vector3Array calculateAllAPotentialGPU(const vec3::Vector3Array &pointVectors,
                                                                   const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] vec3::Vector3Array calculateAllBFieldGPU(const vec3::Vector3Array &pointVectors,
                                                               const PrecisionArguments &usedPrecision) const;

        [[nodiscard]] vec3::Matrix3Array calculateAllBGradientGPU(const vec3::Vector3Array &pointVectors,
                                                                  const PrecisionArguments &usedPrecision) const;

        static bool isZAxisCase(const Coil &primary, const Coil &secondary);
        static std::pair<bool, double> improvedPrecisionCase(const Coil &primary, const Coil &secondary);

        static std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
        calculateRingIncrementPosition(int angularBlocks, int angularIncrements, const Coil &sec,
                                       bool improvedPrecision,
                                       double offset);

        static void generateCoilPairArgumentsData(const Coil &primary, const Coil &secondary,
                                                  CoilPairArgumentsData &coilPairArgumentsData,
                                                  const CoilPairArguments &inductanceArguments,
                                                  bool forceCalculation);

        static double calculateMutualInductanceZAxisSlow(const Coil &primary, const Coil &secondary, double zDisplacement,
                                                         const CoilPairArguments &inductanceArguments,
                                                         ComputeMethod computeMethod = CPU_ST);

        static double calculateMutualInductanceZAxisFast(const Coil &primary, const Coil &secondary, double zDisplacement,
                                                         const CoilPairArguments &inductanceArguments,
                                                         ComputeMethod computeMethod = CPU_ST);

        static double calculateMutualInductanceGeneral(const Coil &primary, const Coil &secondary,
                                                       const CoilPairArguments &inductanceArguments,
                                                       ComputeMethod computeMethod = CPU_ST);

        static double calculateAmpereForceZAxisSlow(const Coil &primary, const Coil &secondary, double zDisplacement,
                                                    const CoilPairArguments &forceArguments,
                                                    ComputeMethod computeMethod = CPU_ST);

        static double calculateAmpereForceZAxisFast(const Coil &primary, const Coil &secondary, double zDisplacement,
                                                    const CoilPairArguments &forceArguments,
                                                    ComputeMethod computeMethod = CPU_ST);

        static std::pair<vec3::Vector3, vec3::Vector3>
        calculateAmpereForceGeneral(const Coil &primary, const Coil &secondary,
                                    const CoilPairArguments &forceArguments, ComputeMethod computeMethod);

        static std::vector<double>
        calculateAllMutualInductanceArrangementsMTD(const Coil &primary, const Coil &secondary,
                                                    const vec3::Vector3Array &primaryPositions,
                                                    const vec3::Vector3Array &secondaryPositions,
                                                    const std::vector<double> &primaryYAngles,
                                                    const std::vector<double> &primaryZAngles,
                                                    const std::vector<double> &secondaryYAngles,
                                                    const std::vector<double> &secondaryZAngles,
                                                    PrecisionFactor precisionFactor = PrecisionFactor());
        static std::vector<double>
        calculateAllMutualInductanceArrangementsGPU(const Coil &primary, const Coil &secondary,
                                                    const vec3::Vector3Array &primaryPositions,
                                                    const vec3::Vector3Array &secondaryPositions,
                                                    const std::vector<double> &primaryYAngles,
                                                    const std::vector<double> &primaryZAngles,
                                                    const std::vector<double> &secondaryYAngles,
                                                    const std::vector<double> &secondaryZAngles,
                                                    PrecisionFactor precisionFactor = PrecisionFactor());

        static std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
        calculateAllAmpereForceArrangementsMTD(const Coil &primary, const Coil &secondary,
                                               const vec3::Vector3Array &primaryPositions,
                                               const vec3::Vector3Array &secondaryPositions,
                                               const std::vector<double> &primaryYAngles,
                                               const std::vector<double> &primaryZAngles,
                                               const std::vector<double> &secondaryYAngles,
                                               const std::vector<double> &secondaryZAngles,
                                               PrecisionFactor precisionFactor = PrecisionFactor());
        static std::vector<std::pair<vec3::Vector3, vec3::Vector3>>
        calculateAllAmpereForceArrangementsGPU(const Coil &primary, const Coil &secondary,
                                               const vec3::Vector3Array &primaryPositions,
                                               const vec3::Vector3Array &secondaryPositions,
                                               const std::vector<double> &primaryYAngles,
                                               const std::vector<double> &primaryZAngles,
                                               const std::vector<double> &secondaryYAngles,
                                               const std::vector<double> &secondaryZAngles,
                                               PrecisionFactor precisionFactor = PrecisionFactor());


        [[nodiscard]] double calculateSelfInductance(CoilPairArguments inductanceArguments) const;
};

#endif //GENERAL_COIL_PROGRAM_COIL_H
