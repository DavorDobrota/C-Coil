#include "Tensor.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <sstream>


vec3::FieldVector3::FieldVector3() : FieldVector3(0.0, 0.0, 0.0) {}

vec3::FieldVector3::FieldVector3(double x, double y, double z) : x(x), y(y), z(z) {}


vec3::FieldVector3 vec3::FieldVector3::operator+(const FieldVector3 &otherVec) const
{
    return FieldVector3(this->x + otherVec.x,
                        this->y + otherVec.y,
                        this->z + otherVec.z);
}

vec3::FieldVector3 vec3::FieldVector3::operator+=(const FieldVector3 &otherVec)
{
    this->x += otherVec.x;
    this->y += otherVec.y;
    this->z += otherVec.z;
    return *this;
}

vec3::FieldVector3 vec3::FieldVector3::operator-(const FieldVector3 &otherVec) const
{
    return FieldVector3(this->x - otherVec.x,
                        this->y - otherVec.y,
                        this->z - otherVec.z);
}

vec3::FieldVector3 vec3::FieldVector3::operator-=(const FieldVector3 &otherVec)
{
    this->x -= otherVec.x;
    this->y -= otherVec.y;
    this->z -= otherVec.z;
    return *this;
}

vec3::FieldVector3 vec3::FieldVector3::operator*(double multiplier) const
{
    return FieldVector3(x * multiplier,
                        y * multiplier,
                        z * multiplier);
}

vec3::FieldVector3 vec3::FieldVector3::operator*=(double multiplier)
{
    x *= multiplier;
    y *= multiplier;
    z *= multiplier;
    return *this;
}

double vec3::FieldVector3::magnitude() const
{
    return std::sqrt(x * x + y * y + z * z);
}


double vec3::FieldVector3::scalarProduct(FieldVector3 vector1, FieldVector3 vector2)
{
    return vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z;
}

vec3::FieldVector3 vec3::FieldVector3::crossProduct(FieldVector3 vector1, FieldVector3 vector2)
{

    double xOutput = vector1.y * vector2.z - vector1.z * vector2.y;
    double yOutput = vector1.z * vector2.x - vector1.x * vector2.z;
    double zOutput = vector1.x * vector2.y - vector1.y * vector2.x;

    return FieldVector3(xOutput, yOutput, zOutput);
}

vec3::FieldVector3::operator std::string() const
{
    std::stringstream output;

    output << "FieldVector3(" << "x=" << x << ", y=" << y << ", z=" << z << ")";

    return output.str();
}
