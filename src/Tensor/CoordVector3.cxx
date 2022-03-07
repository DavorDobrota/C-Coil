#include "Tensor.h"

#include <cmath>
#include <cstdio>

#include <sstream>


vec3::CoordVector3::CoordVector3() : CoordVector3(CARTESIAN, 0.0, 0.0, 0.0) {}

vec3::CoordVector3::CoordVector3(CoordinateSystem system, double comp1, double comp2, double comp3)
                                : coordinateSystem(system), comp1(comp1), comp2(comp2), comp3(comp3)
{
    switch (system)
    {
        case CYLINDRICAL:
            if (comp2 < 0)
            {
                fprintf(stderr, "CYLINDRICAL: Radius cannot be negative: %.15f\n", comp2);
                throw "Radius cannot be negative";
            }
            comp3 = std::fmod(comp3, 2*M_PI);
            if (comp3 < 0)
                comp3 += 2*M_PI;
            break;
        case SPHERICAL:
            if (comp1 < 0)
            {
                fprintf(stderr, "SPHERICAL: Radius cannot be negative: %.15f\n", comp1);
                throw "Radius cannot be negative";
            }

            if (comp2 < 0 || comp2 > M_PI)
            {
                fprintf(stderr, "SPHERICAL: Theta must be between 0 and PI: %.15f\n", comp2);
                throw "Theta must be between 0 and PI";
            }
            comp3 = std::fmod(comp3, 2*M_PI);
            if (comp3 < 0)
                comp3 += 2*M_PI;
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
    double newElem1 = comp3;
    double newElem2 = std::sqrt(comp1 * comp1 + comp2 * comp2);
    double newElem3 = std::atan2(comp2, comp1);

    comp1 = newElem1;
    comp2 = newElem2;
    comp3 = newElem3;

    if (comp3 < 0)
        comp3 += 2*M_PI;

    coordinateSystem = CYLINDRICAL;
}

void vec3::CoordVector3::convertCartesianToSpherical()
{
    double temp1 = comp1 * comp1;
    double temp2 = comp2 * comp2;

    double newElem1 = std::sqrt(temp1 + temp2 + comp3 * comp3);
    double newElem2 = std::atan2(std::sqrt(temp1 + temp2), comp3);
    double newElem3 = std::atan2(comp2, comp1);

    comp1 = newElem1;
    comp2 = newElem2;
    comp3 = newElem3;

    if (comp3 < 0)
        comp3 += 2*M_PI;

    coordinateSystem = SPHERICAL;
}

void vec3::CoordVector3::convertCylindricalToCartesian()
{
    double newElem1 = comp2 * std::cos(comp3);
    double newElem2 = comp2 * std::sin(comp3);
    double newElem3 = comp1;

    comp1 = newElem1;
    comp2 = newElem2;
    comp3 = newElem3;

    coordinateSystem = CARTESIAN;
}

void vec3::CoordVector3::convertCylindricalToSpherical()
{
    double newElem1 = std::sqrt(comp1 * comp1 + comp2 * comp2);
    double newElem2 = std::atan2(comp2, comp1);

    comp1 = newElem1;
    comp2 = newElem2;
    // the third component remains the same

    coordinateSystem = SPHERICAL;
}

void vec3::CoordVector3::convertSphericalToCartesian()
{
    double newElem1 = comp1 * std::sin(comp2) * std::cos(comp3);
    double newElem2 = comp1 * std::sin(comp2) * std::sin(comp3);
    double newElem3 = comp1 * std::cos(comp2);

    comp1 = newElem1;
    comp2 = newElem2;
    comp3 = newElem3;

    coordinateSystem = CARTESIAN;
}

void vec3::CoordVector3::convertSphericalToCylindrical()
{
    double newElem1 = comp1 * std::cos(comp2);
    double newElem2 = comp1 * std::sin(comp2);

    comp1 = newElem1;
    comp2 = newElem2;
    // the third component remains the same

    coordinateSystem = CYLINDRICAL;
}

vec3::FieldVector3 vec3::CoordVector3::convertToFieldVector(const CoordVector3 &vector)
{
    CoordVector3 vectorCopy = vector;
    vectorCopy.convertToCartesian();
    return vec3::FieldVector3(vectorCopy.comp1, vectorCopy.comp2, vectorCopy.comp3);
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
            output << "x=" << comp1 << ", y=" << comp2
                   << ", z=" << comp3;
            break;
        case vec3::CoordinateSystem::CYLINDRICAL:
            output << "z=" << comp1 << ", r=" << comp2
                   << ", p=" << comp3;
            break;
        case vec3::CoordinateSystem::SPHERICAL:
            output << "r=" << comp1 << ", t=" << comp2
                   << ", p=" << comp3;
            break;
    };

    output << ")";

    return output.str();
}
