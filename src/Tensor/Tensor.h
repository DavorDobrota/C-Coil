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

    class Triplet;

    class Vector3Array;
    class Matrix3Array;

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

            static FieldVector3 getFromCylindricalCoords(double z, double r, double phi);
            static FieldVector3 getFromSphericalCoords(double r, double theta, double phi);

            [[nodiscard]] Triplet getAsCylindricalCoords() const;
            [[nodiscard]] Triplet getAsSphericalCoords() const;

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

            [[nodiscard]] double det() const;

            Matrix3 operator+(const Matrix3 &mat) const;
            Matrix3 operator+=(const Matrix3 &mat);

            Matrix3 operator*(double multiplier) const;
            void operator*=(double multiplier);
            Matrix3 operator*(const Matrix3 &mat) const;
            FieldVector3 operator*(const FieldVector3 &vec) const;

            explicit operator std::string() const;
    };

    class Triplet
    {
        public:

            double first;
            double second;
            double third;

            Triplet();
            explicit Triplet(double first, double second, double third);

            explicit operator std::string() const;
    };

    class Vector3Array
    {
        private:

            std::vector<FieldVector3> vectorArray;

        public:

            Vector3Array();
            explicit Vector3Array(const std::vector<FieldVector3> &vectorArray);

            void append(const FieldVector3 &appendedVector3);
            void append(double x, double y, double z);
            void reserve(size_t reserveSize);
            void resize(size_t newSize);
            [[nodiscard]] size_t size() const;

            std::vector<FieldVector3> & getStdVectorRef();

            [[nodiscard]] std::vector<double> x() const;
            [[nodiscard]] std::vector<double> y() const;
            [[nodiscard]] std::vector<double> z() const;
            [[nodiscard]] std::vector<double> abs() const;

            FieldVector3 operator[](int index);
            void operator+=(const FieldVector3 &appendedVector3);

            explicit operator std::string() const;
    };

    class Matrix3Array
    {
        private:

        std::vector<Matrix3> matrixArray;

        public:

            Matrix3Array();
            explicit Matrix3Array(const std::vector<Matrix3> &matrixArray);

            void append(const Matrix3 &appendedMatrix3);
            void append(double xx, double xy, double xz, double yx, double yy, double yz, double zx, double zy, double zz);
            void reserve(size_t reserveSize);
            void resize(size_t newSize);
            [[nodiscard]] size_t size() const;

            std::vector<Matrix3> & getStdVectorRef();

            [[nodiscard]] std::vector<double> xx() const;
            [[nodiscard]] std::vector<double> xy() const;
            [[nodiscard]] std::vector<double> xz() const;
            [[nodiscard]] std::vector<double> yx() const;
            [[nodiscard]] std::vector<double> yy() const;
            [[nodiscard]] std::vector<double> yz() const;
            [[nodiscard]] std::vector<double> zx() const;
            [[nodiscard]] std::vector<double> zy() const;
            [[nodiscard]] std::vector<double> zz() const;
            [[nodiscard]] std::vector<double> det() const;

            Matrix3 operator[](int index);
            void operator+=(const Matrix3 &appendedMatrix3);

            explicit operator std::string() const;
    };
}


#endif //GENERAL_COIL_PROGRAM_Vector3_H
