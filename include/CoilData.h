
#ifndef GENERAL_COIL_PROGRAM_COILDATA_H
#define GENERAL_COIL_PROGRAM_COILDATA_H

#include "Coil.h"

#define TYPE double

// TODO - complete module, it is not supposed to be part of this project
const int arrSize = 12;

struct CoilData
{
    CoilData(const Coil &coil);

    int numOfTurns;
    TYPE currentDensity;

    TYPE innerRadius;
    TYPE thickness;
    TYPE length;

    int angularIterations;
    int lengthIterations;
    int thicknessIncrements;

    TYPE positionArray[arrSize];
    TYPE weightArray[arrSize];
};

#endif //GENERAL_COIL_PROGRAM_COILDATA_H
