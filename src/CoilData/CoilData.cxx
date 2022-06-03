#include "CoilData.h"

#include <sstream>


//CoilData::CoilData()
//{
//    numOfTurns = 0;
//    currentDensity = 0.0;
//
//    innerRadius = 0.0;
//    thickness = 0.0;
//    length = 0.0;
//
//    angularIterations = 0;
//    lengthIterations = 0;
//    thicknessIncrements = 0;
//
//    for (int i = 0; i < GPU_INCREMENTS; ++i)
//    {
//        positionArray[i] = Legendre::positionMatrix[GPU_INCREMENTS - 1][i];
//        weightArray[i] = Legendre::weightsMatrix[GPU_INCREMENTS - 1][i];
//    }
//
//    positionVector{0.0, 0.0, 0.0};


//}

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
        << "const_factor=" << constFactor
        << ", inner_radius=" << innerRadius
        << ", thickness=" << thickness
        << ", length=" << length
        << ", angular_iterations=" << angularIncrements
        << ", length_iterations=" << lengthIncrements
        << ", thickness_increments=" << thicknessIncrements
        << ", position_array=" << stringifyArray(positionArray, arrSize)
        << ", weight_array=" << stringifyArray(weightArray, arrSize)
        << ")";

    return output.str();
}
