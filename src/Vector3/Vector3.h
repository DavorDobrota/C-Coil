#ifndef GENERAL_COIL_PROGRAM_VECTOR3_H
#define GENERAL_COIL_PROGRAM_VECTOR3_H

struct Vector3
{
    double component1;
    double component2;
    double component3;

    Vector3();

    static Vector3 getVectorCartesian(double x, double y, double z);
    static Vector3 getVectorCylindrical(double z, double r, double phi);
    static Vector3 getVectorSpherical(double r, double theta, double phi);

    bool isCartesian();
    bool isCylindrical();
    bool isSpherical();

    void convertToCartesian();
    void convertToCylindrical();
    void convertToSpherical();

    static double scalarProduct(Vector3 vec1, Vector3 vec2);
    static Vector3 crossProduct(Vector3 vec1, Vector3 vec2);

    private:

        short coordinateSystemIndex;

        explicit Vector3(short index, double comp1, double comp2, double comp3);

        void convertCartesianToCylindrical();
        void convertCartesianToSpherical();
        void convertCylindricalToCartesian();
        void convertCylindricalToSpherical();
        void convertSphericalToCartesian();
        void convertSphericalToCylindrical();

};

#endif //GENERAL_COIL_PROGRAM_VECTOR3_H
