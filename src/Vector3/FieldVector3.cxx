#include "Vector3.h"

#include <cmath>
#include <utility>

vec3::FieldVector3::FieldVector3() : FieldVector3(CARTESIAN, 0.0, 0.0, 0.0, CoordVector3()) {}

vec3::FieldVector3::FieldVector3(CoordinateSystem system, double comp1, double comp2, double comp3, CoordVector3 position)
        : Vector3(system, comp1, comp2, comp3)
{
    positionVector = std::move(position);
}

double vec3::FieldVector3::scalarProduct(FieldVector3 vec1, FieldVector3 vec2)
{
    vec1.convertToCartesian();
    vec2.convertToCartesian();

    return vec1.component1 * vec2.component1 + vec1.component2 * vec2.component2 + vec1.component3 * vec2.component3;
}

vec3::FieldVector3 vec3::FieldVector3::crossProduct(FieldVector3 vec1, FieldVector3 vec2)
{
    vec1.convertToCartesian();
    vec2.convertToCartesian();

    double xComponent = vec1.component2 * vec2.component3 - vec1.component3 * vec2.component2;
    double yComponent = vec1.component3 * vec2.component1 - vec1.component1 * vec2.component3;
    double zComponent = vec1.component1 * vec2.component2 - vec1.component2 * vec2.component1;

    return FieldVector3(CARTESIAN, xComponent, yComponent, zComponent, CoordVector3());
}

void vec3::FieldVector3::convertCartesianToCylindrical()
{
    positionVector.convertToCylindrical();

    double newComp1 = component3;
    double newComp2
}

void vec3::FieldVector3::convertCartesianToSpherical()
{

}

void vec3::FieldVector3::convertCylindricalToCartesian()
{

}

void vec3::FieldVector3::convertCylindricalToSpherical()
{

}

void vec3::FieldVector3::convertSphericalToCartesian()
{

}

void vec3::FieldVector3::convertSphericalToCylindrical()
{

}
