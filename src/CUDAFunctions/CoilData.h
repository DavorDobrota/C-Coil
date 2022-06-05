#ifndef GENERAL_COIL_PROGRAM_COILDATA_H
#define GENERAL_COIL_PROGRAM_COILDATA_H

#define TYPE float
#define GPU_INCREMENTS 20


struct CoilData
{
    CoilData() :

        constFactor{0.0},

        innerRadius{0.0},
        thickness{0.0},
        length{0.0},

        angularIncrements{0},
        lengthIncrements{0},
        thicknessIncrements{0},

        angularBlocks{0},
        lengthBlocks{0},
        thicknessBlocks{0},

        positionArray{},
        weightArray{},
        cosPrecomputeArray{},

        positionVector{},
        transformArray{},
        invTransformArray{}
    {}

    TYPE constFactor;

    TYPE innerRadius;
    TYPE thickness;
    TYPE length;

    int angularIncrements;
    int lengthIncrements;
    int thicknessIncrements;

    int angularBlocks;
    int lengthBlocks;
    int thicknessBlocks;

    TYPE positionArray[GPU_INCREMENTS];
    TYPE weightArray[GPU_INCREMENTS];
    TYPE cosPrecomputeArray[GPU_INCREMENTS];

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

#endif //GENERAL_COIL_PROGRAM_COILDATA_H
