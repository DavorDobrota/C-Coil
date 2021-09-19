#ifndef GENERAL_COIL_PROGRAM_VECTOR3_H
#define GENERAL_COIL_PROGRAM_VECTOR3_H

#include <vector>

namespace vec3
{
    enum CoordinateSystem {CARTESIAN, CYLINDRICAL, SPHERICAL};

    class CoordVector3;
    class FieldVector3;
    class Matrix3;

    class CoordVector3
    {
        public:

            double component1;
            double component2;
            double component3;

            CoordVector3();
            explicit CoordVector3(CoordinateSystem system, double comp1, double comp2, double comp3);

            [[nodiscard]] bool isCartesian() const;
            [[nodiscard]] bool isCylindrical() const;
            [[nodiscard]] bool isSpherical() const;

            void convertToCartesian();
            void convertToCylindrical();
            void convertToSpherical();

            static void convertAllToCartesian(std::vector<CoordVector3> &Vector3Array);
            static void convertAllToCylindrical(std::vector<CoordVector3> &Vector3Array);
            static void convertAllToSpherical(std::vector<CoordVector3> &Vector3Array);

            static FieldVector3 convertToFieldVector(CoordVector3 inputVector);
            static CoordVector3 convertToCoordVector(const FieldVector3 &inputVector);

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

            double xComponent;
            double yComponent;
            double zComponent;

            FieldVector3();
            explicit FieldVector3(double x, double y, double z);

            FieldVector3 operator+(const FieldVector3 &otherVec) const;
            void operator+=(const FieldVector3 &otherVec);
            FieldVector3 operator-(const FieldVector3 &otherVec) const;
            void operator-=(const FieldVector3 &otherVec);
            FieldVector3 operator*(double multiplier) const;
            void operator*=(double multiplier);

            double magnitude() const;

            static double scalarProduct(FieldVector3 vec1, FieldVector3 vec2);
            static FieldVector3 crossProduct(FieldVector3 vec1, FieldVector3 vec2);
    };

    class Matrix3
    {
        public:

            double xxElement;
            double xyElement;
            double xzElement;

            double yxElement;
            double yyElement;
            double yzElement;

            double zxElement;
            double zyElement;
            double zzElement;

            Matrix3();
            explicit Matrix3(double xx, double xy, double xz, double yx, double yy, double yz, double zx, double zy, double zz);

            Matrix3 operator+(const Matrix3 &mat) const;
            void operator+=(const Matrix3 &mat);
            Matrix3 operator*(const Matrix3 &mat) const;
            FieldVector3 operator*(const FieldVector3 &vec) const;
    };
}


#endif //GENERAL_COIL_PROGRAM_Vector3_H
