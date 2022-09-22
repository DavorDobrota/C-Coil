#ifndef GENERAL_COIL_PROGRAM_COMPARE_H
#define GENERAL_COIL_PROGRAM_COMPARE_H

#include "Coil.h"
#include "CoilGroup.h"
#include "Benchmark.h"

/// @brief Contains functions that are used to compare precision of Coil and CoilGroup methods with relevant literature,
/// as well as some miscellaneous CPU, GPU, and MTD value generation.
namespace Compare
{
    /// @brief Compares vector potential A and magnetic field B computed by the CPU and GPU. Relative errors are printed.
    void fieldsPrecisionCPUvsGPU();
    /// @brief Compares the precision of CPU and GPU methods for computing mutual inductance, force, and torque.
    void mutualInductanceAndForceTorquePrecisionCPUvsGPU();
    /// @brief Prints mutual inductance parallel case values from paper K. Song, J. Feng, R. Zhao, X. Wu,
    /// ``A general mutual inductance formula for parallel non-coaxial circular coils,'' in ACES Journal, vol. 34,
    /// no. 9, pp. 1385-1390, September 2019.
    void mutualInductanceSpecialCase();

    /// @brief Prints force between a pair of filaments placed on the z-axis which is used to asses the precision of
    /// approach in paper Z. J. Wang and Y. Ren, ``Magnetic Force and Torque Calculation Between Circular
    /// Coils With Nonparallel Axes,'' in IEEE Trans. Appl. Supercond., vol. 24, no. 4, pp. 1-5, Aug. 2014,
    /// Art no. 4901505.
    void forceTorqueFilamentsZAxis();
    /// @brief Prints values from paper Z. J. Wang and Y. Ren, ``Magnetic Force and Torque Calculation Between Circular
    /// Coils With Nonparallel Axes,'' in IEEE Trans. Appl. Supercond., vol. 24, no. 4, pp. 1-5, Aug. 2014,
    /// Art no. 4901505.
    void forceTorqueThickCoilsGeneral();
    /// @brief Prints values from paper S. I. Babic and C. Akyel, ``Magnetic Force Calculation Between Thin Coaxial
    /// Circular Coils in Air,'' IEEE Trans. Magn., vol. 44, no. 4, pp. 445-452, April 2008.
    void forceTorqueThinCoilsZAxis();
    /// @brief Prints values from S. Babic and C. Akyel, ``Magnetic Force Between Inclined Circular Filaments Placed
    /// in Any Desired Position,'' IEEE Trans. Magn., vol. 48, no. 1, pp. 69-80, Jan. 2012.
    void forceTorqueFilamentsGeneral();
    /// @brief Prints values for z-axis force on a system of custom coils used in our prior research.
    void forceTorqueZAxis();
    /// @brief Evaluates the approximation of a Coil with a dipole moment for appropriate cases: coils far apart and
    /// one coil much smaller than the other.
    void forceOnDipoleVsForceTorque();

    /// @brief Prints mutual inductance general case values from paper S. Babic, C. Akyel and S. J. Salon,
    /// ``New procedures for calculating the mutual inductance of the system: filamentary circular coil-massive
    /// circular solenoid,'' in IEEE Trans. Magn., vol. 39, no. 3, pp. 1131-1134, May 2003.
    void mutualInductanceMisalignedCoils();
    /// @brief Outputs and prints mutual inductance parallel case values from paper J. T. Conway, ``Inductance
    /// Calculations for Circular Coils of Rectangular Cross Section and Parallel Axes Using Bessel and Struve Functions,''
    /// IEEE Trans. Magn., vol. 46, no. 1, pp. 75-81, Jan. 2010.
    void mutualInductanceParallelAxesGraphs();
    /// @brief Outputs and prints mutual inductance parallel case values from paper Y. Luo, X. Wang and X. Zhou,
    /// ``Inductance Calculations for Circular Coils With Rectangular Cross Section and Parallel Axes Using Inverse
    /// Mellin Transform and Generalized Hypergeometric Functions,'' IEEE Trans. Power Electron., vol. 32, no. 2,
    /// pp. 1367-1374, Feb. 2017.
    void mutualInductanceParallelAxes();
    /// @brief Outputs and prints mutual inductance general case values from paper J. T. Conway, ``Mutual inductance
    /// of thick coils for arbitrary relative orientation and position,'' 2017 Progress in Electromagnetics Research
    /// Symposium - Fall (PIERS - FALL), 2017, pp. 1388-1395.
    void mutualInductanceGeneralCase();
    /// @brief Outputs mutual inductance general case values from paper Y. Wang, X. Xie, Y. Zhou and W. Huan,
    /// ``Calculation and Modeling Analysis of Mutual Inductance Between Coreless Circular Coils With Rectangular
    /// Cross Section in Arbitrary Spatial Position,'' 2020 IEEE 5th Information Technology and Mechatronics
    /// Engineering Conference (ITOEC), 2020, pp. 1258-1267
    void mutualInductanceGeneralGraphs();
    /// @brief Prints mutual inductance general case values which we found to show reduced precision (edge cases).
    void mutualInductanceGeneralEdgeCases();

    /// @brief Prints mutual inductance z-axis case values from paper T. Župan, Ž. Štih and B. Trkulja, ``Fast and
    /// Precise Method for Inductance Calculation of Coaxial Circular Coils With Rectangular Cross Section Using the
    /// One-Dimensional Integration of Elementary Functions Applicable to Superconducting Magnets,'' IEEE Trans.
    /// Appl. Supercond., vol. 24, no. 2, pp. 81-89, April 2014, Art no. 4901309.
    void mutualInductanceZAxis();
    /// @brief Outputs and prints self inductance values from paper J. T. Conway, ``Inductance
    /// Calculations for Circular Coils of Rectangular Cross Section and Parallel Axes Using Bessel and Struve
    /// Functions,'' IEEE Trans. Magn., vol. 46, no. 1, pp. 75-81, Jan. 2010.
    void selfInductance();

    ///@brief Calculates and can print values of magnetic field for a number of coils in toroidal arrangement.
    void fieldsCoilGroupMTD(int coilCount = 100, int pointCount = 10'000,
                            int threadCount = g_defaultThreadCount, bool print = true);

}

#endif //GENERAL_COIL_PROGRAM_COMPARE_H
