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

/**
 * @brief Structure used to represent universal calculation precision.
 * A custom precision measure from interval [1.0, 15.0]. Increasing the factor by 1.0 doubles the performance.
 * @details Choosing a forbidden value results in the default value of 5.0, the same one used in the default constructor.
 * The total number of increments, where k is the number of integration dimensions, m the base number of increments,
 * and p the relativePrecision, is given as m^k * 2^(p-1). Currently, m = 10 for both CPU and GPU.
 */
struct PrecisionFactor
{
    ///@brief Default constructor, sets relativePrecision to 5.0.
    PrecisionFactor();
    ///@brief Sets relativePrecision to the given double value
    explicit PrecisionFactor(double relativePrecision);

    double relativePrecision;

    explicit operator std::string() const;
};

/**
 * @brief Structure used to store precision data (block and increment count) for an individual Coil.
 * @details There are 3 integration layers and thus 6 values arranged in 3 block-increment pairs.
 * The number of blocks is the number of sub-intervals in integration and increments determine quadrature order
 */
struct PrecisionArguments
{
    ///@brief Default constructor, sets all blocks to 1, and all increments to default quadrature value, currently 20.
    PrecisionArguments();
    ///@brief Takes 6 integer values specifying the number of blocks and increments assigned to each layer.
    ///@details Angular increments are along phi, thickness along a, and length along b (length is not utilised)
    explicit PrecisionArguments(int angularBlocks, int thicknessBlocks, int lengthBlocks,
                                int angularIncrements, int thicknessIncrements, int lengthIncrements);

    int angularBlocks;
    int thicknessBlocks;
    int lengthBlocks;

    int angularIncrements;
    int thicknessIncrements;
    int lengthIncrements;

    /**
     * @brief Returns results of ordinary (CPU) increment balancing algorithm for one coil.
     * @details The CPU can divide the interval of integration into multiple blocks with
     * @param coil Reference to the coil for which PrecisionArguments are generated
     * @param precisionFactor Relative precision, determines the total number of increments
     */
    static PrecisionArguments getCoilPrecisionArgumentsCPU(const Coil &coil, PrecisionFactor precisionFactor);
    /**
     * @brief Returns results of specialised GPU increment balancing algorithm for one coil.
     * @details The GPU uses only 1 block and a maximum of GPU_INCREMENTS increments per layer.
     * Precision factor definition may therefore be invalid as less increments are assigned than specified.
     * @param coil Reference to the coil for which PrecisionArguments are generated
     * @param precisionFactor Relative precision, determines the total number of increments
     */
    static PrecisionArguments getCoilPrecisionArgumentsGPU(const Coil &coil, PrecisionFactor precisionFactor);
    /**
     * @brief Usually only used for CoilGroup computeAll MInductance and Force.
     * Returns results of specialised GPU increment balancing algorithm for secondary coil.
     * @details The GPU uses only 1 block and a maximum of GPU_INCREMENTS increments per layer.
     * Precision factor definition may therefore be invalid as less increments are assigned than specified.
     * @param coil Secondary coil which is represented as a number of points
     * @param precisionFactor Relative precision, determines the total number of increments
     */
    static PrecisionArguments getSecondaryCoilPrecisionArgumentsGPU(const Coil &coil, PrecisionFactor precisionFactor);

    explicit operator std::string() const;

    private:

        static PrecisionArguments calculatePrecisionArguments(const Coil &coil, PrecisionFactor precisionFactor,
                                                              bool useGPU = false);
};

/**
 * @brief Structure used to store precision data (block and increment count) for a system of two coils.
 * Used for custom precision arguments in interaction calculations.
 * @details There are 3 integration layers and thus 6 values arranged in 3 block-increment pairs, for each coil.
 * Effectively stores two PrecisionArguments, one for the  primary and one for the secondary coil
 */
struct CoilPairArguments
{
    friend PrecisionArguments;
    ///@brief Default constructor, initialises primary and secondary coil PrecisionArguments to default values
    CoilPairArguments();
    ///@brief Takes two PrecisionArguments arguments and creates new CoilPairArguments
    explicit CoilPairArguments(const PrecisionArguments &primaryPrecision,
                               const PrecisionArguments &secondaryPrecision);

    PrecisionArguments primaryPrecision;
    PrecisionArguments secondaryPrecision;

    /**
     * @brief Applies increment balancing for a pair of coils and returns appropriate precision arguments.
     * Takes coil type, z-axis case, and GPU implementation limitations into consideration.
     * @param primary Reference to the primary coil, usually a larger coil is selected
     * @param secondary Reference to a secondary coil which is represented as a number of points
     * @param precisionFactor Relative precision, determines the total number of increments
     * @param computeMethod Three calculation options: CPU_ST, CPU_MT, and GPU. CPU_ST is default.
     * @param zAxisCase Are the arguments generated for two coils that have a common axis (x = 0, y = 0 for both).
     * @param pureGPU Is the method purely GPU based, such as in computeAllMInductance with GPU computeMethod.
     */
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
 * @brief Primary class in this project. Has a unique identifier.
 * Models a circular coil with a rectangular cross section and uniform current density.
 * @details The model works best for a solid block of material and when the effects of windings are negligible.
 * Primary attributes are length (b), thickness (a), the inner radius (R) of the internal cylindrical hole,
 * and the number of turns of wire (windings). Length and thickness can be set to 0.0 and that is interpreted as
 * having a thin coil (b = 0.0), flat coil (a = 0.0), or filament (a = 0.0, b = 0.0).
 * The coil is oriented like a spherical vector, when angles are (0.0, 0.0) the coil axis is along the z-axis.
 * The first angle is the rotation around the y-axis, and the second around the z-axis.
 * Default precision settings are stored and used if custom ones are not provided.
 * Calculating fields inside the coil is not recommended.
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

        ///@brief Returns the unique coil identifier.
        [[nodiscard]] unsigned long long getId() const;
        ///@brief Returns the radius of inner cylinder of the circular coil.
        [[nodiscard]] double getInnerRadius() const;
        ///@brief Returns the thickness of windings of the circular coil.
        [[nodiscard]] double getThickness() const;
        ///@brief Returns the length of the circular coil.
        [[nodiscard]] double getLength() const;
        ///@brief Returns the number of windings (turns) of the circular coil.
        [[nodiscard]] int getNumOfTurns() const;

        ///@brief Returns the current density of a circular coil. Ill defined for thin coils, pancakes and filaments.
        [[nodiscard]] double getCurrentDensity() const;
        ///@brief Returns a current passing through each winding of a circular coil.
        [[nodiscard]] double getCurrent() const;

        ///@brief Returns current wire resistivity, determined by the used material. By default, copper is used.
        [[nodiscard]] double getWireResistivity() const;
        ///@brief Returns if the coil is sine wave (AC) driven or DC driven.
        [[nodiscard]] bool isSineDriven() const;
        ///@brief Returns the frequency of the AC sine wave driving the coil. 0.0 if the it is DC driven.
        [[nodiscard]] double getSineFrequency() const;

        ///@brief Uses lazy loading and returns an equivalent magnetic dipole moment (approximation at large distance)
        [[nodiscard]] vec3::Vector3 getMagneticMoment();
        ///@brief Calculates average wire thickness supposing the winding is orthogonal.
        [[nodiscard]] double getAverageWireThickness() const;

        ///@brief Returns last set or calculated value of self inductance
        [[nodiscard]] double getSelfInductance() const;
        ///@brief Uses lazy loading and returns coil resistance of the coil with skin effect compensation.
        [[nodiscard]] double getResistance();
        ///@brief Uses lazy loading and returns inductive reactance of the coil, capacitance not included.
        [[nodiscard]] double getReactance();
        ///@brief Uses lazy loading and returns the magnitude of the coil impedance.
        [[nodiscard]] double getImpedance();

        ///@brief Returns the default PrecisionArguments used for CPU calculations.
        [[nodiscard]] const PrecisionArguments &getPrecisionSettingsCPU() const;
        ///@brief Returns the default PrecisionArguments used for GPU calculations.
        [[nodiscard]] const PrecisionArguments &getPrecisionSettingsGPU() const;
        ///@brief Returns the default number of threads.
        [[nodiscard]] int getThreadCount() const;
        ///@brief Returns the type of methods the coil is using. Thin and rectangular coils use fast methods.
        [[nodiscard]] bool isUsingFastMethod() const;
        ///@brief Returns the type of circular coil with rectangular cross section.
        [[nodiscard]] CoilType getCoilType() const;

        ///@brief Returns position of the coil in external Cartesian coordinate system
        [[nodiscard]] vec3::Vector3 getPositionVector() const;
        ///@brief Returns a pair of angles <yAxisAngle, zAxisAngle> which represent the coil orientation
        [[nodiscard]] std::pair<double, double> getRotationAngles() const;

        ///@brief Returns the inverse transformation matrix used to simplify field calculation.
        [[nodiscard]] vec3::Matrix3 getTransformationMatrix() const;
        ///@brief Returns the transformation matrix used to adapt field tensors to the external coordinate system.
        [[nodiscard]] vec3::Matrix3 getInverseTransformationMatrix() const;

        ///@brief Sets current density and calculates appropriate current. Not suitable for thin coils, pancakes and filaments.
        void setCurrentDensity(double currentDensity);
        ///@brief Sets the current through windings and calculates the appropriate current density.
        void setCurrent(double current);
        ///@brief Sets wire resistivity, necessary for resistance calculation.
        void setWireResistivity(double wireResistivity);
        ///@brief Sets AC sine driven frequency. 0.0 is used for DC.
        void setSineFrequency(double sineFrequency);

        ///@brief Sets the given, custom PrecisionArguments as default for CPU calculations.
        void setDefaultPrecisionCPU(const PrecisionArguments &precisionSettings);
        ///@brief Calculates the PrecisionArguments for the appropriate PrecisionFactor
        /// and sets them as default for CPU calculations
        void setDefaultPrecisionCPU(PrecisionFactor precisionFactor = PrecisionFactor());
        ///@brief Sets the given, custom PrecisionArguments as default for GPU calculations.
        void setDefaultPrecisionGPU(const PrecisionArguments &precisionSettings);
        ///@brief Calculates the PrecisionArguments for the appropriate PrecisionFactor
        /// and sets them as default for GPU calculations
        void setDefaultPrecisionGPU(PrecisionFactor precisionFactor = PrecisionFactor());
        ///@brief Calculates the PrecisionArguments for the appropriate PrecisionFactor
        /// and sets them as default for both CPU and GPU
        void setDefaultPrecision(PrecisionFactor precisionFactor = PrecisionFactor());

        ///@brief Sets the number of threads used in field and interaction calculations.
        void setThreadCount(int threadCount);
        /**
         * @brief Repositions and reorients the coil in the external coordinate system.
         * @param positionVector Position of the coil center.
         * @param yAxisAngle Rotation angle around the y-axis from [0, PI].
         * @param zAxisAngle Rotation angle around the z-axis from [0, 2PI].
         */
        void setPositionAndOrientation(vec3::Vector3 positionVector = vec3::Vector3(),
                                       double yAxisAngle = 0.0, double zAxisAngle = 0.0);
        ///@brief Sets self inductance of the coil to the provided value, overriding all previous calculations.
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
         * @param precisionFactor Determines the precision of given calculation.
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
         * @param precisionFactor Determines the precision of given calculation.
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
        * @param precisionFactor Determines the precision of given calculation.
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
        * Using CPU_MT compute method is highly advisable, especially for higher precision factors, and the
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
         * @param precisionFactor Determines the precision of given calculation.
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
        * @param precisionFactor Determines the precision of given calculation.
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
