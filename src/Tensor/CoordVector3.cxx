#include "Tensor.h"

#include <cmath>
#include <cstdio>


vec3::CoordVector3::CoordVector3() : CoordVector3(CARTESIAN, 0.0, 0.0, 0.0) {}

vec3::CoordVector3::CoordVector3(CoordinateSystem system, double comp1, double comp2, double comp3)
                                : coordinateSystem(system), component1(comp1), component2(comp2), component3(comp3)
{
    switch (system)
    {
        case CYLINDRICAL:
            if (component2 < 0)
            {
                fprintf(stderr, "CYLINDRICAL: Radius cannot be negative: %.15f\n", comp2);
                throw "Radius cannot be negative";
            }
            component3 = std::fmod(comp3, 2*M_PI);
            if (component3 < 0)
                component3 += 2*M_PI;
            break;
        case SPHERICAL:
            if (component1 < 0)
            {
                fprintf(stderr, "SPHERICAL: Radius cannot be negative: %.15f\n", comp1);
                throw "Radius cannot be negative";
            }

            if (component2 < 0 || component2 > M_PI)
            {
                fprintf(stderr, "SPHERICAL: Theta must be between 0 and PI: %.15f\n", comp2);
                throw "Theta must be between 0 and PI";
            }
            component3 = std::fmod(comp3, 2*M_PI);
            if (component3 < 0)
                component3 += 2*M_PI;
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
    double newComp1 = component3;
    double newComp2 = std::sqrt(component1 * component1 + component2 * component2);
    double newComp3 = std::atan2(component2, component1);

    component1 = newComp1;
    component2 = newComp2;
    component3 = newComp3;

    if (component3 < 0)
        component3 += 2*M_PI;

    coordinateSystem = CYLINDRICAL;
}

void vec3::CoordVector3::convertCartesianToSpherical()
{
    double temp1 = component1 * component1;
    double temp2 = component2 * component2;

    double newComp1 = std::sqrt(temp1 + temp2 + component3 * component3);
    double newComp2 = std::atan2(std::sqrt(temp1 + temp2), component3);
    double newComp3 = std::atan2(component2, component1);

    component1 = newComp1;
    component2 = newComp2;
    component3 = newComp3;

    if (component3 < 0)
        component3 += 2*M_PI;

    coordinateSystem = SPHERICAL;
}

void vec3::CoordVector3::convertCylindricalToCartesian()
{
    double newComp1 = component2 * std::cos(component3);
    double newComp2 = component2 * std::sin(component3);
    double newComp3 = component1;

    component1 = newComp1;
    component2 = newComp2;
    component3 = newComp3;

    coordinateSystem = CARTESIAN;
}

void vec3::CoordVector3::convertCylindricalToSpherical()
{
    double newComp1 = std::sqrt(component1 * component1 + component2 * component2);
    double newComp2 = std::atan2(component2, component1);

    component1 = newComp1;
    component2 = newComp2;
    // the third component remains the same

    coordinateSystem = SPHERICAL;
}

void vec3::CoordVector3::convertSphericalToCartesian()
{
    double newComp1 = component1 * std::sin(component2) * std::cos(component3);
    double newComp2 = component1 * std::sin(component2) * std::sin(component3);
    double newComp3 = component1 * std::cos(component2);

    component1 = newComp1;
    component2 = newComp2;
    component3 = newComp3;

    coordinateSystem = CARTESIAN;
}

void vec3::CoordVector3::convertSphericalToCylindrical()
{
    double newComp1 = component1 * std::cos(component2);
    double newComp2 = component1 * std::sin(component2);

    component1 = newComp1;
    component2 = newComp2;
    // the third component remains the same

    coordinateSystem = CYLINDRICAL;
}

vec3::FieldVector3 vec3::CoordVector3::convertToFieldVector(CoordVector3 inputVector)
{
    inputVector.convertToCartesian();
    return vec3::FieldVector3(inputVector.component1, inputVector.component2, inputVector.component3);
}

vec3::CoordVector3 vec3::CoordVector3::convertToCoordVector(const FieldVector3 &inputVector)
{
    return vec3::CoordVector3(CARTESIAN, inputVector.xComponent, inputVector.yComponent, inputVector.zComponent);
}
