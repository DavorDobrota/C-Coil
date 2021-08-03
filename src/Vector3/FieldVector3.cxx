#include "Vector3.h"


vec3::FieldVector3::FieldVector3() : FieldVector3(0.0, 0.0, 0.0) {}

vec3::FieldVector3::FieldVector3(double x, double y, double z) : xComponent(x), yComponent(y), zComponent(z) {}

double vec3::FieldVector3::scalarProduct(FieldVector3 vec1, FieldVector3 vec2)
{
    return vec1.xComponent * vec2.xComponent + vec1.yComponent * vec2.yComponent + vec1.zComponent * vec2.zComponent;
}

vec3::FieldVector3 vec3::FieldVector3::crossProduct(FieldVector3 vec1, FieldVector3 vec2)
{

    double xOutput = vec1.yComponent * vec2.zComponent - vec1.zComponent * vec2.yComponent;
    double yOutput = vec1.zComponent * vec2.xComponent - vec1.xComponent * vec2.zComponent;
    double zOutput = vec1.xComponent * vec2.yComponent - vec1.yComponent * vec2.xComponent;

    return FieldVector3(xOutput, yOutput, zOutput);
}

