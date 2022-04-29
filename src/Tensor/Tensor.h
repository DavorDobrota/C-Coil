#ifndef GENERAL_COIL_PROGRAM_VECTOR3_H
#define GENERAL_COIL_PROGRAM_VECTOR3_H

#include <vector>
#include <iostream>


namespace vec3
{
    enum CoordinateSystem {CARTESIAN, CYLINDRICAL, SPHERICAL};

    class CoordVector3;
    class FieldVector3;
    class Matrix3;

    class CoordVector3
    {
        public:

            double comp1;
            double comp2;
            double comp3;

            CoordVector3();
            explicit CoordVector3(CoordinateSystem system, double comp1, double comp2, double comp3);

            [[nodiscard]] bool isCartesian() const;
            [[nodiscard]] bool isCylindrical() const;
            [[nodiscard]] bool isSpherical() const;

            void convertToCartesian();
            void convertToCylindrical();
            void convertToSpherical();

            static std::vector<CoordVector3> convertAllToCartesian(std::vector<CoordVector3> &Vector3Array);
            static std::vector<CoordVector3> convertAllToCylindrical(std::vector<CoordVector3> &Vector3Array);
            static std::vector<CoordVector3> convertAllToSpherical(std::vector<CoordVector3> &Vector3Array);

            static FieldVector3 convertToFieldVector(const CoordVector3 &vector);
            static CoordVector3 convertToCoordVector(const FieldVector3 &vector);

            [[nodiscard]] CoordinateSystem getCoordinateSystem() const { return coordinateSystem; }

            explicit operator std::string() const;

        private:

            CoordinateSystem coordinateSystem;

            void convertCartesianToCylindrical();
            void convertCartesianToSpherical();
            void convertCylindricalToCartesian();
            void convertCylindricalToSpherical();
            void convertSphericalToCartesian();
            void convertSphericalToCylindrical();

    };

    class FieldVector3
    {
        public:

            double x;
            double y;
            double z;

            FieldVector3();
            explicit FieldVector3(double x, double y, double z);

            FieldVector3 operator+(const FieldVector3 &otherVec) const;
            FieldVector3 operator+=(const FieldVector3 &otherVec);
            FieldVector3 operator-(const FieldVector3 &otherVec) const;
            FieldVector3 operator-=(const FieldVector3 &otherVec);
            FieldVector3 operator*(double multiplier) const;
            FieldVector3 operator*=(double multiplier);

            [[nodiscard]] double magnitude() const;

            static double scalarProduct(FieldVector3 vector1, FieldVector3 vector2);
            static FieldVector3 crossProduct(FieldVector3 vector1, FieldVector3 vector2);

            explicit operator std::string() const;
    };

    class Matrix3
    {
        public:

            double xx;
            double xy;
            double xz;

            double yx;
            double yy;
            double yz;

            double zx;
            double zy;
            double zz;

            Matrix3();
            explicit Matrix3(double xx, double xy, double xz, double yx, double yy, double yz, double zx, double zy, double zz);

            Matrix3 operator+(const Matrix3 &mat) const;
            Matrix3 operator+=(const Matrix3 &mat);
            Matrix3 operator*(const Matrix3 &mat) const;
            FieldVector3 operator*(const FieldVector3 &vec) const;

            operator std::string() const;
    };
}


#endif //GENERAL_COIL_PROGRAM_Vector3_H
