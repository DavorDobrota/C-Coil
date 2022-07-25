#ifndef GENERAL_COIL_PROGRAM_VECTOR3_H
#define GENERAL_COIL_PROGRAM_VECTOR3_H

#include <vector>
#include <iostream>


namespace vec3
{

    class Vector3;
    class Matrix3;

    class Triplet;

    class Vector3Array;
    class Matrix3Array;


    class Vector3
    {
        public:

            double x;
            double y;
            double z;

            Vector3();
            explicit Vector3(double x, double y, double z);

            Vector3 operator+(const Vector3 &otherVec) const;
            Vector3 operator+=(const Vector3 &otherVec);
            Vector3 operator-(const Vector3 &otherVec) const;
            Vector3 operator-=(const Vector3 &otherVec);
            Vector3 operator*(double multiplier) const;
            Vector3 operator*=(double multiplier);

            [[nodiscard]] double abs() const;

            static double scalarProduct(Vector3 vector1, Vector3 vector2);
            static Vector3 crossProduct(Vector3 vector1, Vector3 vector2);

            static Vector3 getFromCylindricalCoords(double z, double r, double phi);
            static Vector3 getFromSphericalCoords(double r, double theta, double phi);

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
            Vector3 operator*(const Vector3 &vec) const;

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

            std::vector<Vector3> vectorArray;

        public:

            Vector3Array();
            explicit Vector3Array(size_t initSize);
            explicit Vector3Array(const std::vector<Vector3> &vectorArray);

            void append(const Vector3 &appendedVector3);
            void append(double x, double y, double z);
            void reserve(size_t reserveSize);
            void resize(size_t newSize);
            void clear();
            [[nodiscard]] size_t size() const;

            std::vector<Vector3> & getStdVectorRef();

            [[nodiscard]] std::vector<double> x() const;
            [[nodiscard]] std::vector<double> y() const;
            [[nodiscard]] std::vector<double> z() const;
            [[nodiscard]] std::vector<double> abs() const;

            Vector3 & operator[](size_t index);
            const Vector3 & operator[](size_t index) const;
            void operator+=(const Vector3 &appendedVector3);

            explicit operator std::string() const;
    };

    class Matrix3Array
    {
        private:

        std::vector<Matrix3> matrixArray;

        public:

            Matrix3Array();
            explicit Matrix3Array(size_t initSize);
            explicit Matrix3Array(const std::vector<Matrix3> &matrixArray);

            void append(const Matrix3 &appendedMatrix3);
            void append(double xx, double xy, double xz, double yx, double yy, double yz, double zx, double zy, double zz);
            void reserve(size_t reserveSize);
            void resize(size_t newSize);
            void clear();
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

            Matrix3 & operator[](size_t index);
            const Matrix3 & operator[](size_t index) const;
            void operator+=(const Matrix3 &appendedMatrix3);

            explicit operator std::string() const;
    };
}


#endif //GENERAL_COIL_PROGRAM_Vector3_H
