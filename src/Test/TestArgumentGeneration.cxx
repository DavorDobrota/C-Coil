#include "Test.h"
#include "Coil.h"

#include <cstdio>

void testCoilMutualInductanceZAxisDifferentGeometries()
{
    Coil prim1 = Coil(0.03, 0.0, 0.0, 1);
    Coil sec1 = Coil(0.03, 0.0, 0.0, 1);
    printf("%.15g\n", Coil::computeMutualInductance(prim1, sec1, 0.2));

    Coil prim2 = Coil(0.03, 0.03, 0.0, 30);
    Coil sec2 = Coil(0.03, 0.0, 0.0, 1);
    printf("%.15g\n", Coil::computeMutualInductance(prim2, sec2, 0.2));

    Coil prim3 = Coil(0.03, 0.0, 0.0, 1);
    Coil sec3 = Coil(0.03, 0.03, 0.0, 30);
    printf("%.15g\n", Coil::computeMutualInductance(prim3, sec3, 0.2));

    Coil prim4 = Coil(0.03, 0.03, 0.0, 30);
    Coil sec4 = Coil(0.03, 0.03, 0.0, 30);
    printf("%.15g\n", Coil::computeMutualInductance(prim4, sec4, 0.2));

    Coil prim5 = Coil(0.03, 0.0, 0.12, 120);
    Coil sec5 = Coil(0.03,  0.0, 0.0, 1);
    printf("%.15g\n", Coil::computeMutualInductance(prim5, sec5, 0.2));

    Coil prim6 = Coil(0.03, 0.0, 0.0, 1);
    Coil sec6 = Coil(0.03,  0.0, 0.12, 120);
    printf("%.15g\n", Coil::computeMutualInductance(prim6, sec6, 0.2));

    Coil prim7 = Coil(0.03, 0.0, 0.12, 120);
    Coil sec7 = Coil(0.03,  0.0, 0.12, 120);
    printf("%.15g\n", Coil::computeMutualInductance(prim7, sec7, 0.2));

    Coil prim8 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec8 = Coil(0.03,  0.0, 0.0, 1);
    printf("%.15g\n", Coil::computeMutualInductance(prim8, sec8, 0.2));

    Coil prim9 = Coil(0.03,  0.0, 0.0, 1);
    Coil sec9 = Coil(0.03, 0.03, 0.12, 3600);
    printf("%.15g\n", Coil::computeMutualInductance(prim9, sec9, 0.2));

    Coil prim10 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec10 = Coil(0.03, 0.03, 0.12, 3600);
    printf("%.15g\n", Coil::computeMutualInductance(prim10, sec10, 0.2));

    Coil prim11 = Coil(0.03, 0.03, 0.0, 30);
    Coil sec11 = Coil(0.03, 0.0, 0.12, 120);
    printf("%.15g\n", Coil::computeMutualInductance(prim11, sec11, 0.2));

    Coil prim12 = Coil(0.03, 0.0, 0.12, 120);
    Coil sec12 = Coil(0.03, 0.03, 0.0, 30);
    printf("%.15g\n", Coil::computeMutualInductance(prim12, sec12, 0.2));

    Coil prim13 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec13 = Coil(0.03, 0.0, 0.12, 120);
    printf("%.15g\n", Coil::computeMutualInductance(prim13, sec13, 0.2));

    Coil prim14 = Coil(0.03, 0.0, 0.12, 120);
    Coil sec14 = Coil(0.03, 0.03, 0.12, 3600);
    printf("%.15g\n", Coil::computeMutualInductance(prim14, sec14, 0.2));

    Coil prim15 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec15 = Coil(0.03, 0.03, 0.0, 30);
    printf("%.15g\n", Coil::computeMutualInductance(prim15, sec15, 0.2));

    Coil prim16 = Coil(0.03, 0.03, 0.0, 30);
    Coil sec16 = Coil(0.03, 0.03, 0.12, 3600);
    printf("%.15g\n", Coil::computeMutualInductance(prim16, sec16, 0.2));

    printf("\n");
}

void testCoilMutualInductanceGeneralDifferentGeometries()
{
    Coil prim1 = Coil(0.03, 0.0, 0.0, 1);
    Coil sec1 = Coil(0.03, 0.0, 0.0, 1);
    printf("%.12g\n", Coil::computeMutualInductance(prim1, sec1, 0.2, 0.0, 1e-15));

    Coil prim2 = Coil(0.03, 0.03, 0.0, 30);
    Coil sec2 = Coil(0.03, 0.0, 0.0, 1);
    printf("%.12g\n", Coil::computeMutualInductance(prim2, sec2, 0.2, 0.0, 1e-15));

    Coil prim3 = Coil(0.03, 0.0, 0.0, 1);
    Coil sec3 = Coil(0.03, 0.03, 0.0, 30);
    printf("%.12g\n", Coil::computeMutualInductance(prim3, sec3, 0.2, 0.0, 1e-15));

    Coil prim4 = Coil(0.03, 0.03, 0.0, 30);
    Coil sec4 = Coil(0.03, 0.03, 0.0, 30);
    printf("%.12g\n", Coil::computeMutualInductance(prim4, sec4, 0.2, 0.0, 1e-15));

    Coil prim5 = Coil(0.03, 0.0, 0.12, 120);
    Coil sec5 = Coil(0.03,  0.0, 0.0, 1);
    printf("%.12g\n", Coil::computeMutualInductance(prim5, sec5, 0.2, 0.0, 1e-15));

    Coil prim6 = Coil(0.03, 0.0, 0.0, 1);
    Coil sec6 = Coil(0.03,  0.0, 0.12, 120);
    printf("%.12g\n", Coil::computeMutualInductance(prim6, sec6, 0.2, 0.0, 1e-15));

    Coil prim7 = Coil(0.03, 0.0, 0.12, 120);
    Coil sec7 = Coil(0.03,  0.0, 0.12, 120);
    printf("%.12g\n", Coil::computeMutualInductance(prim7, sec7, 0.2, 0.0, 1e-15));

    Coil prim8 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec8 = Coil(0.03,  0.0, 0.0, 1);
    printf("%.12g\n", Coil::computeMutualInductance(prim8, sec8, 0.2, 0.0, 1e-15));

    Coil prim9 = Coil(0.03,  0.0, 0.0, 1);
    Coil sec9 = Coil(0.03, 0.03, 0.12, 3600);
    printf("%.12g\n", Coil::computeMutualInductance(prim9, sec9, 0.2, 0.0, 1e-15));

    Coil prim10 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec10 = Coil(0.03, 0.03, 0.12, 3600);
    printf("%.12g\n", Coil::computeMutualInductance(prim10, sec10, 0.2, 0.0, 1e-15));

    Coil prim11 = Coil(0.03, 0.03, 0.0, 30);
    Coil sec11 = Coil(0.03, 0.0, 0.12, 120);
    printf("%.12g\n", Coil::computeMutualInductance(prim11, sec11, 0.2, 0.0, 1e-15));

    Coil prim12 = Coil(0.03, 0.0, 0.12, 120);
    Coil sec12 = Coil(0.03, 0.03, 0.0, 30);
    printf("%.12g\n", Coil::computeMutualInductance(prim12, sec12, 0.2, 0.0, 1e-15));

    Coil prim13 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec13 = Coil(0.03, 0.0, 0.12, 120);
    printf("%.12g\n", Coil::computeMutualInductance(prim13, sec13, 0.2, 0.0, 1e-15));

    Coil prim14 = Coil(0.03, 0.0, 0.12, 120);
    Coil sec14 = Coil(0.03, 0.03, 0.12, 3600);
    printf("%.12g\n", Coil::computeMutualInductance(prim14, sec14, 0.2, 0.0, 1e-15));

    Coil prim15 = Coil(0.03, 0.03, 0.12, 3600);
    Coil sec15 = Coil(0.03, 0.03, 0.0, 30);
    printf("%.12g\n", Coil::computeMutualInductance(prim15, sec15, 0.2, 0.0, 1e-15));

    Coil prim16 = Coil(0.03, 0.03, 0.0, 30);
    Coil sec16 = Coil(0.03, 0.03, 0.12, 3600);
    printf("%.12g\n", Coil::computeMutualInductance(prim16, sec16, 0.2, 0.0, 1e-15));

    printf("\n");
}

