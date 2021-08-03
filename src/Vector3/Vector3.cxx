#include "Vector3.h"


vec3::Vector3::Vector3(CoordinateSystem system, double comp1, double comp2, double comp3)
{
    coordinateSystem = system;
    component1 = comp1;
    component2 = comp2;
    component3 = comp3;
}

bool vec3::Vector3::isCartesian() { return coordinateSystem == CARTESIAN; }

bool vec3::Vector3::isCylindrical() { return coordinateSystem == CYLINDRICAL; }

bool vec3::Vector3::isSpherical() { return coordinateSystem == SPHERICAL; }

void vec3::Vector3::convertToCartesian()
{
    switch (coordinateSystem)
    {
        case CARTESIAN:
            break;
        case CYLINDRICAL:
            convertCylindricalToCartesian();
            break;
        case SPHERICAL:
            convertSphericalToCartesian();
            break;
    }
}

void vec3::Vector3::convertToCylindrical()
{
    switch (coordinateSystem)
    {
        case CARTESIAN:
            convertCartesianToCylindrical();
            break;
        case CYLINDRICAL:
            break;
        case SPHERICAL:
            convertSphericalToCylindrical();
            break;
    }
}

void vec3::Vector3::convertToSpherical()
{
    switch (coordinateSystem)
    {
        case CARTESIAN:
            convertCartesianToSpherical();
            break;
        case CYLINDRICAL:
            convertCylindricalToSpherical();
            break;
        case SPHERICAL:
            break;
    }
}

void vec3::Vector3::convertAllToCartesian(std::vector<Vector3> &vector3Array)
{
    for (auto & i : vector3Array)
        i.convertToCartesian();
}

void vec3::Vector3::convertAllToCylindrical(std::vector<Vector3> &vector3Array)
{
    for (auto & i : vector3Array)
        i.convertToCylindrical();
}

void vec3::Vector3::convertAllToSpherical(std::vector<Vector3> &vector3Array)
{
    for (auto & i : vector3Array)
        i.convertToSpherical();
}
