#include <cstdio>
#include <cmath>

#include "Coil.h"
#include "CoilGroup.h"
#include "Test.h"


int main()
{
    Coil coil = Coil(0.1, 0.1, 0.1, 10000);

    for (int i = 0; i <= 140; ++i)
    {
        printf("%.15g\n", coil.computeAndSetSelfInductance(PrecisionFactor(1.0 + i * 0.1 )));
    }

    return 0;
}
