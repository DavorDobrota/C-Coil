#ifndef GENERAL_COIL_PROGRAM_TEST_H
#define GENERAL_COIL_PROGRAM_TEST_H

#include "Coil.h"


/// @brief Contains functions that are used to test whether compute methods are implemented correctly. This module is
/// accessible only from C++ and is not included in Python.
namespace Test
{
    /// @brief Tests basic constructors and lazy loading of different parameters.
    void testNewCoilParameters();
    /// @brief Tests if the rotation matrix properly transforms fields (coordinate conversion).
    void testCoilPositionAndRotation();
    /// @brief Tests if coarse-grained multithreading (MTD) is implemented correctly (CPU_ST is reference).
    void testCoilGroupComputeAllMTD();

    /// @brief Tests if general mutual inductance method returns values similar to z-axis mutual inductance method.
    void testMutualInductanceGeneralForZAxis(ComputeMethod computeMethod);
    /// @brief Tests the increment balancing algorithm for z-axis CoilPairArguments generation.
    void testMInductanceZAxisArgumentGeneration();
    /// @brief Tests 16 possible configurations of coils for z-axis case. There are two circular coils,
    /// each either a filament, flat coil, thin coil, or a rectangular coil.
    void testMInductanceZAxisDifferentGeometries();

    /// @brief Tests the increment balancing algorithm for general CoilPairArguments generation.
    void testMInductanceGeneralArgumentGeneration();
    /// @brief Tests 16 possible configurations of coils for general case. There are two circular coils,
    /// each either a filament, flat coil, thin coil, or a rectangular coil.
    void testMInductanceGeneralDifferentGeometries();

    /// @brief Tests if general force and torque method returns values similar to z-axis force method.
    void testAmpereForceGeneralForZAxis();
    /// @brief Tests if the magnetic gradient tensor behaves properly, especially for z-axis positions (singular case).
    void testGradientTensor();

    /// @brief Tests if CoilGroup coarse-grained multithreading (MTD) is implemented properly.
    void testCoilGroupFieldsMTD();

    /// @brief Tests if Coil::computeAllMutualInductanceArrangements is implemented correctly (CPU_ST is reference)
    void testCoilMInductanceArrangements();
    /// @brief Tests if Coil::computeAllForceTorqueArrangements is implemented correctly (CPU_ST is reference)
    void testCoilForceArrangements();
    /// @brief Tests if CoilGroup::computeAllMutualInductanceArrangements is implemented correctly (CPU_ST is reference)
    void testGroupMInductanceArrangements();
    /// @brief Tests if CoilGroup::computeAllForceTorqueArrangements is implemented correctly (CPU_ST is reference)
    void testGroupForceArrangements();

}


#endif //GENERAL_COIL_PROGRAM_TEST_H
