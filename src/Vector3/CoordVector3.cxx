#include "Vector3.h"

#include <cmath>

vec3::CoordVector3::CoordVector3() : Vector3(CARTESIAN, 0.0, 0.0, 0.0){}

vec3::CoordVector3::CoordVector3(CoordinateSystem system, double comp1, double comp2, double comp3): Vector3(system, comp1, comp2, comp3)
{
    switch (system)
    {
        case CARTESIAN:
            break;
        case CYLINDRICAL:
            if (component2 < 0)
                throw "Radius cannot be negative";
            component3 = std::fmod(comp3, 2*M_PI);
            if (component3 < 0)
                component3 += 2*M_PI;
        case SPHERICAL:
            if (component2 < 0 || component2 > M_PI)
                throw "Theta must be between 0 and PI";
            component3 = std::fmod(comp3, 2*M_PI);
            if (component3 < 0)
                component3 += 2*M_PI;
    }
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
