#include "Vector3.h"

#include <cmath>

namespace
{
    const short cartesianIndex = 1;
    const short cylindricalIndex = 2;
    const short sphericalIndex = 3;
}

Vector3::Vector3() : Vector3(0, 0.0, 0.0, 0.0) {}

Vector3::Vector3(short index, double comp1, double comp2, double comp3)
{
    coordinateSystemIndex = index;
    component1 = comp1;
    component2 = comp2;
    component3 = comp3;
}

Vector3 Vector3::getVectorCartesian(double x, double y, double z)
{
    return Vector3(cartesianIndex, x, y, z);
}

Vector3 Vector3::getVectorCylindrical(double z, double r, double phi)
{
    return Vector3(cylindricalIndex, z, r, std::fmod(phi, 2*M_PI));
}

Vector3 Vector3::getVectorSpherical(double r, double theta, double phi)
{
    if (theta > M_PI)
        throw "Theta must be between 0 and PI";
    return Vector3(sphericalIndex, r, theta, std::fmod(phi, 2*M_PI));
}

bool Vector3::isCartesian() { return coordinateSystemIndex == cartesianIndex; }

bool Vector3::isCylindrical() { return coordinateSystemIndex == cylindricalIndex; }

bool Vector3::isSpherical() { return coordinateSystemIndex == sphericalIndex; }

void Vector3::convertCartesianToCylindrical()
{
    double newComp1 = component3;
    double newComp2 = std::sqrt(component1 * component1 + component2 * component2);
    double newComp3 = std::atan2(component2, component1);

    component1 = newComp1;
    component2 = newComp2;
    component3 = newComp3;

    if (component3 < 0)
        component3 += 2*M_PI;

    coordinateSystemIndex = cylindricalIndex;
}

void Vector3::convertCartesianToSpherical()
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

    coordinateSystemIndex = sphericalIndex;
}

void Vector3::convertCylindricalToCartesian()
{
    double newComp1 = component2 * std::cos(component3);
    double newComp2 = component2 * std::sin(component3);
    double newComp3 = component1;

    component1 = newComp1;
    component2 = newComp2;
    component3 = newComp3;

    coordinateSystemIndex = cartesianIndex;
}

void Vector3::convertCylindricalToSpherical()
{
    double newComp1 = std::sqrt(component1 * component1 + component2 * component2);
    double newComp2 = std::atan2(component2, component1);

    component1 = newComp1;
    component2 = newComp2;
    // the third component remains the same

    coordinateSystemIndex = sphericalIndex;
}

void Vector3::convertSphericalToCartesian()
{
    double newComp1 = component1 * std::sin(component2) * std::cos(component3);
    double newComp2 = component1 * std::sin(component2) * std::sin(component3);
    double newComp3 = component1 * std::cos(component2);

    component1 = newComp1;
    component2 = newComp2;
    component3 = newComp3;

    coordinateSystemIndex = cartesianIndex;
}

void Vector3::convertSphericalToCylindrical()
{
    double newComp1 = component1 * std::cos(component2);
    double newComp2 = component1 * std::sin(component2);

    component1 = newComp1;
    component2 = newComp2;
    // the third component remains the same

    coordinateSystemIndex = cylindricalIndex;
}

void Vector3::convertToCartesian()
{
    switch (coordinateSystemIndex)
    {
        case cylindricalIndex:
            convertCylindricalToCartesian();
            break;
        case sphericalIndex:
            convertSphericalToCartesian();
            break;
    }
}

void Vector3::convertToCylindrical()
{
    switch (coordinateSystemIndex)
    {
        case cartesianIndex:
            convertCartesianToCylindrical();
            break;
        case sphericalIndex:
            convertSphericalToCylindrical();
            break;
    }
}

void Vector3::convertToSpherical()
{
    switch (coordinateSystemIndex)
    {
        case cartesianIndex:
            convertCartesianToSpherical();
            break;
        case cylindricalIndex:
            convertCylindricalToSpherical();
            break;
    }
}

double Vector3::scalarProduct(Vector3 vec1, Vector3 vec2)
{
    vec1.convertToCartesian();
    vec2.convertToCartesian();

    return vec1.component1 * vec2.component1 + vec1.component2 * vec2.component2 + vec1.component3 * vec2.component3;
}

Vector3 Vector3::crossProduct(Vector3 vec1, Vector3 vec2)
{
    vec1.convertToCartesian();
    vec2.convertToCartesian();

    double xComponent = vec1.component2 * vec2.component3 - vec1.component3 * vec2.component2;
    double yComponent = vec1.component3 * vec2.component1 - vec1.component1 * vec2.component3;
    double zComponent = vec1.component1 * vec2.component2 - vec1.component2 * vec2.component1;

    return Vector3(cartesianIndex, xComponent, yComponent, zComponent);
}


