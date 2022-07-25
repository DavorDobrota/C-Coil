#include "Tensor.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <sstream>


vec3::Vector3::Vector3() : Vector3(0.0, 0.0, 0.0) {}

vec3::Vector3::Vector3(double x, double y, double z) : x(x), y(y), z(z) {}


vec3::Triplet::Triplet() : Triplet(0.0, 0.0, 0.0) {}

vec3::Triplet::Triplet(double first, double second, double third) : first(first), second(second), third(third) {}


vec3::Vector3 vec3::Vector3::operator+(const Vector3 &otherVec) const
{
    return Vector3(this->x + otherVec.x,
                        this->y + otherVec.y,
                        this->z + otherVec.z);
}

vec3::Vector3 vec3::Vector3::operator+=(const Vector3 &otherVec)
{
    this->x += otherVec.x;
    this->y += otherVec.y;
    this->z += otherVec.z;
    return *this;
}

vec3::Vector3 vec3::Vector3::operator-(const Vector3 &otherVec) const
{
    return Vector3(this->x - otherVec.x,
                        this->y - otherVec.y,
                        this->z - otherVec.z);
}

vec3::Vector3 vec3::Vector3::operator-=(const Vector3 &otherVec)
{
    this->x -= otherVec.x;
    this->y -= otherVec.y;
    this->z -= otherVec.z;
    return *this;
}

vec3::Vector3 vec3::Vector3::operator*(double multiplier) const
{
    return Vector3(x * multiplier,
                        y * multiplier,
                        z * multiplier);
}

vec3::Vector3 vec3::Vector3::operator*=(double multiplier)
{
    x *= multiplier;
    y *= multiplier;
    z *= multiplier;
    return *this;
}

double vec3::Vector3::abs() const
{
    return std::sqrt(x * x + y * y + z * z);
}


double vec3::Vector3::scalarProduct(Vector3 vector1, Vector3 vector2)
{
    return vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z;
}

vec3::Vector3 vec3::Vector3::crossProduct(Vector3 vector1, Vector3 vector2)
{

    double xOutput = vector1.y * vector2.z - vector1.z * vector2.y;
    double yOutput = vector1.z * vector2.x - vector1.x * vector2.z;
    double zOutput = vector1.x * vector2.y - vector1.y * vector2.x;

    return Vector3(xOutput, yOutput, zOutput);
}


vec3::Vector3 vec3::Vector3::getFromCylindricalCoords(double z, double r, double phi)
{
    return vec3::Vector3(r * std::cos(phi), r * std::sin(phi), z);
}

vec3::Vector3 vec3::Vector3::getFromSphericalCoords(double r, double theta, double phi)
{
    double rCoord = r * std::sin(theta);
    return vec3::Vector3(rCoord * std::cos(phi), rCoord * std::sin(phi), r * std::cos(phi));
}

vec3::Triplet vec3::Vector3::getAsCylindricalCoords() const
{
    return vec3::Triplet(z, std::sqrt(x * x + y * y), std::atan2(y, x));
}

vec3::Triplet vec3::Vector3::getAsSphericalCoords() const
{
    double temp = x * x + y * y;

    return vec3::Triplet(std::sqrt(temp + z * z), std::atan2(temp, z), std::atan2(y, x));
}


vec3::Vector3::operator std::string() const
{
    std::stringstream output;

    output << "Vector3(" << "x=" << x << ", y=" << y << ", z=" << z << ")";

    return output.str();
}

vec3::Triplet::operator std::string() const
{
    std::stringstream output;

    output << "Triplet(" << "first=" << first << ", second=" << second << ", third=" << third << ")";

    return output.str();
}
