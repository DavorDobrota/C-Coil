#ifndef GENERAL_COIL_PROGRAM_VECTOR3_H
#define GENERAL_COIL_PROGRAM_VECTOR3_H

#include <vector>

namespace vec3
{
    enum CoordinateSystem {CARTESIAN, CYLINDRICAL, SPHERICAL};

    class Vector3
    {
        public:

            double component1;
            double component2;
            double component3;

            bool isCartesian();
            bool isCylindrical();
            bool isSpherical();

            void convertToCartesian();
            void convertToCylindrical();
            void convertToSpherical();

            static void convertAllToCartesian(std::vector<Vector3> &Vector3Array);
            static void convertAllToCylindrical(std::vector<Vector3> &Vector3Array);
            static void convertAllToSpherical(std::vector<Vector3> &Vector3Array);

        protected:

            CoordinateSystem coordinateSystem;

            explicit Vector3(CoordinateSystem system, double comp1, double comp2, double comp3);

            virtual void convertCartesianToCylindrical() = 0;
            virtual void convertCartesianToSpherical() = 0;
            virtual void convertCylindricalToCartesian() = 0;
            virtual void convertCylindricalToSpherical() = 0;
            virtual void convertSphericalToCartesian() = 0;
            virtual void convertSphericalToCylindrical() = 0;

    };

    class CoordVector3: public Vector3
    {
        public:

            CoordVector3();
            explicit CoordVector3(CoordinateSystem system, double comp1, double comp2, double comp3);

        private:

            void convertCartesianToCylindrical() override;
            void convertCartesianToSpherical() override;
            void convertCylindricalToCartesian() override;
            void convertCylindricalToSpherical() override;
            void convertSphericalToCartesian() override;
            void convertSphericalToCylindrical() override;

    };

    class FieldVector3: public Vector3
    {
        public:

            CoordVector3 positionVector;

            FieldVector3();
            explicit FieldVector3(CoordinateSystem system, double comp1, double comp2, double comp3, CoordVector3 position);

            static double scalarProduct(FieldVector3 vec1, FieldVector3 vec2);
            static FieldVector3 crossProduct(FieldVector3 vec1, FieldVector3 vec2);

        private:

            void convertCartesianToCylindrical() override;
            void convertCartesianToSpherical() override;
            void convertCylindricalToCartesian() override;
            void convertCylindricalToSpherical() override;
            void convertSphericalToCartesian() override;
            void convertSphericalToCylindrical() override;
    };
}


#endif //GENERAL_COIL_PROGRAM_Vector3_H
