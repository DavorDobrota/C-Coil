#include <iostream>
#include "OldCoil.h"
#include "Polynomial.h"
#include "Coil.h"
#include "Test.h"

extern thread_pool tp;

int main()
{

//    testCoilAmpereForceForFilamentsZAxis();
//    testCoilAmpereForceGeneralCase();

//testCoilAmpereForceThinCoils();

//testCoilMutualInductanceZAxis();

    testCoilMutualInductanceGeneralGraphs();

    return 0;
}
