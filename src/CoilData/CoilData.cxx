#include "CoilData.h"

#include <sstream>


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

CoilData::operator std::string() const
{
    std::stringstream output;

    auto stringifyArray = [](auto ar, int size) -> std::string
    {
        std::stringstream output;

        output << "[";

        for(int i = 0; i < size; i++)
        {
            if(i != 0)
                output << ", ";
            output << ar[i];
        }

        output << "]";

        return output.str();
    };

    output << "CoilData("
        << "num_of_turns" << numOfTurns
        << ", current_density=" << currentDensity
        << ", inner_radius=" << innerRadius
        << ", thickness=" << thickness
        << ", length=" << length
        << ", angular_iterations=" << angularIterations
        << ", length_iterations=" << lengthIterations
        << ", thickness_increments=" << thicknessIncrements
        << ", position_array=" << stringifyArray(positionArray, arrSize)
        << ", weight_array=" << stringifyArray(weightArray, arrSize)
        << ")";

    return output.str();
}
