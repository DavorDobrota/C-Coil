#include "CoilData.h"


CoilData::CoilData(const Coil &coil)
{
    numOfTurns = coil.getNumOfTurns();
    currentDensity = coil.getCurrentDensity();

    innerRadius = coil.getInnerRadius();
    thickness = coil.getThickness();
    length = coil.getLength();

    PrecisionArguments arguments = coil.getPrecisionSettings();

    angularIterations = arguments.angularBlockCount * arguments.angularIncrementCount;
    lengthIterations = arguments.lengthBlockCount * arguments.lengthIncrementCount;
    thicknessIncrements = arguments.thicknessBlockCount * arguments.thicknessIncrementCount;

}
