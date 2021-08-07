#include <iostream>

#include "OldCoil.h"
#include "Coil.h"
#include "Test.h"

extern thread_pool tp;

int main()
{
    testCoilAmpereForceZAxisMTScaling(16);
    testCoilAmpereForceZGeneralMTScaling(16);

    return 0;
}
