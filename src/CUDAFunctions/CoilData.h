#ifndef GENERAL_COIL_PROGRAM_COILDATA_H
#define GENERAL_COIL_PROGRAM_COILDATA_H

#define TYPE float


struct CoilData
{
    CoilData() :

        constFactor{0.0},
        useFastMethod{false},

        innerRadius{0.0},
        thickness{0.0},
        length{0.0},

        angularIncrements{0},
        lengthIncrements{0},
        thicknessIncrements{0},

        thicknessPositionArray{},
        thicknessWeightArray{},

        cosPrecomputeArray{},
        angularWeightArray{},

        positionVector{},
        transformArray{},
        invTransformArray{}
    {}

    TYPE constFactor;
    bool useFastMethod;

    TYPE innerRadius;
    TYPE thickness;
    TYPE length;

    int angularIncrements;
    int lengthIncrements;
    int thicknessIncrements;

    TYPE thicknessPositionArray[GPU_INCREMENTS];
    TYPE thicknessWeightArray[GPU_INCREMENTS];

    TYPE cosPrecomputeArray[GPU_INCREMENTS];
    TYPE angularWeightArray[GPU_INCREMENTS];

    TYPE positionVector[3];
    TYPE transformArray[9];
    TYPE invTransformArray[9];
};

struct DataVector
{
    DataVector() :
        x{0.0},
        y{0.0},
        z{0.0}
    {}

    TYPE x;
    TYPE y;
    TYPE z;
};

struct DataMatrix
{
    DataMatrix() :
        xx{0.0},
        xy{0.0},
        xz{0.0},

        yx{0.0},
        yy{0.0},
        yz{0.0},

        zx{0.0},
        zy{0.0},
        zz{0.0}
    {}

    TYPE xx;
    TYPE xy;
    TYPE xz;

    TYPE yx;
    TYPE yy;
    TYPE yz;

    TYPE zx;
    TYPE zy;
    TYPE zz;
};

#endif //GENERAL_COIL_PROGRAM_COILDATA_H
