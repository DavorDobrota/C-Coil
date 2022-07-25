#include "Test.h"
#include "Coil.h"

#include <cstdio>


void testMInductanceZAxisDifferentGeometries()
{
    vec3::Vector3 defaultVec = vec3::Vector3(0.0, 0.0, 0.2);

    printf("Testing mutual inductance z-axis in 16 different combinations\n\n");

    Coil prim1 = Coil(0.03, 0.0, 0.0, 1);
    Coil sec1 = Coil(0.03, 0.0, 0.0, 1);
    sec1.setPositionAndOrientation(defaultVec);
    printf(" 1:  %.15g\n", Coil::computeMutualInductance(prim1, sec1));

    Coil prim2 = Coil(0.03, 0.03, 0.0, 30);
    Coil sec2 = Coil(0.03, 0.0, 0.0, 1);
    sec2.setPositionAndOrientation(defaultVec);
    printf(" 2:  %.15g\n", Coil::computeMutualInductance(prim2, sec2));

    Coil prim3 = Coil(0.03, 0.0, 0.0, 1);
    Coil sec3 = Coil(0.03, 0.03, 0.0, 30);
    sec3.setPositionAndOrientation(defaultVec);
    printf(" 3:  %.15g\n", Coil::computeMutualInductance(prim3, sec3));

    Coil prim4 = Coil(0.03, 0.03, 0.0, 30);
    Coil sec4 = Coil(0.03, 0.03, 0.0, 30);
    sec4.setPositionAndOrientation(defaultVec);
    printf(" 4:  %.15g\n", Coil::computeMutualInductance(prim4, sec4));

    Coil prim5 = Coil(0.03, 0.0, 0.12, 120);
    Coil sec5 = Coil(0.03,  0.0, 0.0, 1);
    sec5.setPositionAndOrientation(defaultVec);
    printf(" 5:  %.15g\n", Coil::computeMutualInductance(prim5, sec5));

    Coil prim6 = Coil(0.03, 0.0, 0.0, 1);
    Coil sec6 = Coil(0.03,  0.0, 0.12, 120);
    sec6.setPositionAndOrientation(defaultVec);
    printf(" 6:  %.15g\n", Coil::computeMutualInductance(prim6, sec6));

    Coil prim7 = Coil(0.03, 0.0, 0.12, 120);
    Coil sec7 = Coil(0.03,  0.0, 0.12, 120);
    sec7.setPositionAndOrientation(defaultVec);
    printf(" 7:  %.15g\n", Coil::computeMutualInductance(prim7, sec7));

    Coil prim8 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec8 = Coil(0.03,  0.0, 0.0, 1);
    sec8.setPositionAndOrientation(defaultVec);
    printf(" 8:  %.15g\n", Coil::computeMutualInductance(prim8, sec8));

    Coil prim9 = Coil(0.03,  0.0, 0.0, 1);
    Coil sec9 = Coil(0.03, 0.03, 0.12, 3600);
    sec9.setPositionAndOrientation(defaultVec);
    printf(" 9:  %.15g\n", Coil::computeMutualInductance(prim9, sec9));

    Coil prim10 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec10 = Coil(0.03, 0.03, 0.12, 3600);
    sec10.setPositionAndOrientation(defaultVec);
    printf("10:  %.15g\n", Coil::computeMutualInductance(prim10, sec10));

    Coil prim11 = Coil(0.03, 0.03, 0.0, 30);
    Coil sec11 = Coil(0.03, 0.0, 0.12, 120);
    sec11.setPositionAndOrientation(defaultVec);
    printf("11:  %.15g\n", Coil::computeMutualInductance(prim11, sec11));

    Coil prim12 = Coil(0.03, 0.0, 0.12, 120);
    Coil sec12 = Coil(0.03, 0.03, 0.0, 30);
    sec12.setPositionAndOrientation(defaultVec);
    printf("12:  %.15g\n", Coil::computeMutualInductance(prim12, sec12));

    Coil prim13 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec13 = Coil(0.03, 0.0, 0.12, 120);
    sec13.setPositionAndOrientation(defaultVec);
    printf("13:  %.15g\n", Coil::computeMutualInductance(prim13, sec13));

    Coil prim14 = Coil(0.03, 0.0, 0.12, 120);
    Coil sec14 = Coil(0.03, 0.03, 0.12, 3600);
    sec14.setPositionAndOrientation(defaultVec);
    printf("14:  %.15g\n", Coil::computeMutualInductance(prim14, sec14));

    Coil prim15 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec15 = Coil(0.03, 0.03, 0.0, 30);
    sec15.setPositionAndOrientation(defaultVec);
    printf("15:  %.15g\n", Coil::computeMutualInductance(prim15, sec15));

    Coil prim16 = Coil(0.03, 0.03, 0.0, 30);
    Coil sec16 = Coil(0.03, 0.03, 0.12, 3600);
    sec16.setPositionAndOrientation(defaultVec);
    printf("16:  %.15g\n", Coil::computeMutualInductance(prim16, sec16));

    printf("\n");
}

void testMInductanceGeneralDifferentGeometries()
{
    vec3::Vector3 defaultVec = vec3::Vector3(1e-8, 0.0, 0.2);

    printf("Testing mutual inductance general in 16 different combinations\n\n");

    Coil prim1 = Coil(0.03, 0.0, 0.0, 1);
    Coil sec1 = Coil(0.03, 0.0, 0.0, 1);
    sec1.setPositionAndOrientation(defaultVec);
    printf(" 1:  %.15g\n", Coil::computeMutualInductance(prim1, sec1));

    Coil prim2 = Coil(0.03, 0.03, 0.0, 30);
    Coil sec2 = Coil(0.03, 0.0, 0.0, 1);
    sec2.setPositionAndOrientation(defaultVec);
    printf(" 2:  %.15g\n", Coil::computeMutualInductance(prim2, sec2));

    Coil prim3 = Coil(0.03, 0.0, 0.0, 1);
    Coil sec3 = Coil(0.03, 0.03, 0.0, 30);
    sec3.setPositionAndOrientation(defaultVec);
    printf(" 3:  %.15g\n", Coil::computeMutualInductance(prim3, sec3));

    Coil prim4 = Coil(0.03, 0.03, 0.0, 30);
    Coil sec4 = Coil(0.03, 0.03, 0.0, 30);
    sec4.setPositionAndOrientation(defaultVec);
    printf(" 4:  %.15g\n", Coil::computeMutualInductance(prim4, sec4));

    Coil prim5 = Coil(0.03, 0.0, 0.12, 120);
    Coil sec5 = Coil(0.03,  0.0, 0.0, 1);
    sec5.setPositionAndOrientation(defaultVec);
    printf(" 5:  %.15g\n", Coil::computeMutualInductance(prim5, sec5));

    Coil prim6 = Coil(0.03, 0.0, 0.0, 1);
    Coil sec6 = Coil(0.03,  0.0, 0.12, 120);
    sec6.setPositionAndOrientation(defaultVec);
    printf(" 6:  %.15g\n", Coil::computeMutualInductance(prim6, sec6));

    Coil prim7 = Coil(0.03, 0.0, 0.12, 120);
    Coil sec7 = Coil(0.03,  0.0, 0.12, 120);
    sec7.setPositionAndOrientation(defaultVec);
    printf(" 7:  %.15g\n", Coil::computeMutualInductance(prim7, sec7));

    Coil prim8 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec8 = Coil(0.03,  0.0, 0.0, 1);
    sec8.setPositionAndOrientation(defaultVec);
    printf(" 8:  %.15g\n", Coil::computeMutualInductance(prim8, sec8));

    Coil prim9 = Coil(0.03,  0.0, 0.0, 1);
    Coil sec9 = Coil(0.03, 0.03, 0.12, 3600);
    sec9.setPositionAndOrientation(defaultVec);
    printf(" 9:  %.15g\n", Coil::computeMutualInductance(prim9, sec9));

    Coil prim10 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec10 = Coil(0.03, 0.03, 0.12, 3600);
    sec10.setPositionAndOrientation(defaultVec);
    printf("10:  %.15g\n", Coil::computeMutualInductance(prim10, sec10));

    Coil prim11 = Coil(0.03, 0.03, 0.0, 30);
    Coil sec11 = Coil(0.03, 0.0, 0.12, 120);
    sec11.setPositionAndOrientation(defaultVec);
    printf("11:  %.15g\n", Coil::computeMutualInductance(prim11, sec11));

    Coil prim12 = Coil(0.03, 0.0, 0.12, 120);
    Coil sec12 = Coil(0.03, 0.03, 0.0, 30);
    sec12.setPositionAndOrientation(defaultVec);
    printf("12:  %.15g\n", Coil::computeMutualInductance(prim12, sec12));

    Coil prim13 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec13 = Coil(0.03, 0.0, 0.12, 120);
    sec13.setPositionAndOrientation(defaultVec);
    printf("13:  %.15g\n", Coil::computeMutualInductance(prim13, sec13));

    Coil prim14 = Coil(0.03, 0.0, 0.12, 120);
    Coil sec14 = Coil(0.03, 0.03, 0.12, 3600);
    sec14.setPositionAndOrientation(defaultVec);
    printf("14:  %.15g\n", Coil::computeMutualInductance(prim14, sec14));

    Coil prim15 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec15 = Coil(0.03, 0.03, 0.0, 30);
    sec15.setPositionAndOrientation(defaultVec);
    printf("15:  %.15g\n", Coil::computeMutualInductance(prim15, sec15));

    Coil prim16 = Coil(0.03, 0.03, 0.0, 30);
    Coil sec16 = Coil(0.03, 0.03, 0.12, 3600);
    sec16.setPositionAndOrientation(defaultVec);
    printf("16:  %.15g\n", Coil::computeMutualInductance(prim16, sec16));

    printf("\n");
}

