#include "Tensor.h"

#include <cmath>
#include <cstdio>

#include <sstream>


vec3::CoordVector3::CoordVector3() : CoordVector3(CARTESIAN, 0.0, 0.0, 0.0) {}

vec3::CoordVector3::CoordVector3(CoordinateSystem system, double elem1, double elem2, double elem3)
                                : coordinateSystem(system), elem1(elem1), elem2(elem2), elem3(elem3)
{
    switch (system)
    {
        case CYLINDRICAL:
            if (elem2 < 0)
            {
                fprintf(stderr, "CYLINDRICAL: Radius cannot be negative: %.15f\n", elem2);
                throw "Radius cannot be negative";
            }
            elem3 = std::fmod(elem3, 2*M_PI);
            if (elem3 < 0)
                elem3 += 2*M_PI;
            break;
        case SPHERICAL:
            if (elem1 < 0)
            {
                fprintf(stderr, "SPHERICAL: Radius cannot be negative: %.15f\n", elem1);
                throw "Radius cannot be negative";
            }

            if (elem2 < 0 || elem2 > M_PI)
            {
                fprintf(stderr, "SPHERICAL: Theta must be between 0 and PI: %.15f\n", elem2);
                throw "Theta must be between 0 and PI";
            }
            elem3 = std::fmod(elem3, 2*M_PI);
            if (elem3 < 0)
                elem3 += 2*M_PI;
            break;
        default:
            break;
    }
}

bool vec3::CoordVector3::isCartesian() const { return coordinateSystem == CARTESIAN; }

bool vec3::CoordVector3::isCylindrical() const { return coordinateSystem == CYLINDRICAL; }

bool vec3::CoordVector3::isSpherical() const { return coordinateSystem == SPHERICAL; }

void vec3::CoordVector3::convertToCartesian()
{
    switch (coordinateSystem)
    {
        case CYLINDRICAL:
            convertCylindricalToCartesian();
            break;
        case SPHERICAL:
            convertSphericalToCartesian();
            break;
        default:
            break;
    }
}

void vec3::CoordVector3::convertToCylindrical()
{
    switch (coordinateSystem)
    {
        case CARTESIAN:
            convertCartesianToCylindrical();
            break;
        case SPHERICAL:
            convertSphericalToCylindrical();
            break;
        default:
            break;
    }
}

void vec3::CoordVector3::convertToSpherical()
{
    switch (coordinateSystem)
    {
        case CARTESIAN:
            convertCartesianToSpherical();
            break;
        case CYLINDRICAL:
            convertCylindricalToSpherical();
            break;
        default:
            break;
    }
}

std::vector<vec3::CoordVector3> vec3::CoordVector3::convertAllToCartesian(std::vector<CoordVector3> &vector3Array)
{
    for (auto & i : vector3Array)
        i.convertToCartesian();
    return vector3Array;
}

std::vector<vec3::CoordVector3> vec3::CoordVector3::convertAllToCylindrical(std::vector<CoordVector3> &vector3Array)
{
    for (auto & i : vector3Array)
        i.convertToCylindrical();
    return vector3Array;
}

std::vector<vec3::CoordVector3> vec3::CoordVector3::convertAllToSpherical(std::vector<CoordVector3> &vector3Array)
{
    for (auto & i : vector3Array)
        i.convertToSpherical();
    return vector3Array;
}

void vec3::CoordVector3::convertCartesianToCylindrical()
{
    double newElem1 = elem3;
    double newElem2 = std::sqrt(elem1 * elem1 + elem2 * elem2);
    double newElem3 = std::atan2(elem2, elem1);

    elem1 = newElem1;
    elem2 = newElem2;
    elem3 = newElem3;

    if (elem3 < 0)
        elem3 += 2*M_PI;

    coordinateSystem = CYLINDRICAL;
}

void vec3::CoordVector3::convertCartesianToSpherical()
{
    double temp1 = elem1 * elem1;
    double temp2 = elem2 * elem2;

    double newElem1 = std::sqrt(temp1 + temp2 + elem3 * elem3);
    double newElem2 = std::atan2(std::sqrt(temp1 + temp2), elem3);
    double newElem3 = std::atan2(elem2, elem1);

    elem1 = newElem1;
    elem2 = newElem2;
    elem3 = newElem3;

    if (elem3 < 0)
        elem3 += 2*M_PI;

    coordinateSystem = SPHERICAL;
}

void vec3::CoordVector3::convertCylindricalToCartesian()
{
    double newElem1 = elem2 * std::cos(elem3);
    double newElem2 = elem2 * std::sin(elem3);
    double newElem3 = elem1;

    elem1 = newElem1;
    elem2 = newElem2;
    elem3 = newElem3;

    coordinateSystem = CARTESIAN;
}

void vec3::CoordVector3::convertCylindricalToSpherical()
{
    double newElem1 = std::sqrt(elem1 * elem1 + elem2 * elem2);
    double newElem2 = std::atan2(elem2, elem1);

    elem1 = newElem1;
    elem2 = newElem2;
    // the third component remains the same

    coordinateSystem = SPHERICAL;
}

void vec3::CoordVector3::convertSphericalToCartesian()
{
    double newElem1 = elem1 * std::sin(elem2) * std::cos(elem3);
    double newElem2 = elem1 * std::sin(elem2) * std::sin(elem3);
    double newElem3 = elem1 * std::cos(elem2);

    elem1 = newElem1;
    elem2 = newElem2;
    elem3 = newElem3;

    coordinateSystem = CARTESIAN;
}

void vec3::CoordVector3::convertSphericalToCylindrical()
{
    double newElem1 = elem1 * std::cos(elem2);
    double newElem2 = elem1 * std::sin(elem2);

    elem1 = newElem1;
    elem2 = newElem2;
    // the third component remains the same

    coordinateSystem = CYLINDRICAL;
}

vec3::FieldVector3 vec3::CoordVector3::convertToFieldVector(const CoordVector3 &vector)
{
    CoordVector3 vectorCopy = vector;
    vectorCopy.convertToCartesian();
    return vec3::FieldVector3(vectorCopy.elem1, vectorCopy.elem2, vectorCopy.elem3);
}

vec3::CoordVector3 vec3::CoordVector3::convertToCoordVector(const FieldVector3 &vector)
{
    return vec3::CoordVector3(CARTESIAN, vector.x, vector.y, vector.z);
}

vec3::CoordVector3::operator std::string() const
{
    std::stringstream output;

    output << "CoordVector3(";

    switch(getCoordinateSystem())
    {
        case vec3::CoordinateSystem::CARTESIAN:
            output << "x=" << elem1 << ", y=" << elem2
                   << ", z=" << elem3;
            break;
        case vec3::CoordinateSystem::CYLINDRICAL:
            output << "z=" << elem1 << ", r=" << elem2
                   << ", p=" << elem3;
            break;
        case vec3::CoordinateSystem::SPHERICAL:
            output << "r=" << elem1 << ", t=" << elem2
                   << ", p=" << elem3;
            break;
    };

    output << ")";

    return output.str();
}
