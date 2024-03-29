#ifndef GENERAL_COIL_PROGRAM_COILDATA_H
#define GENERAL_COIL_PROGRAM_COILDATA_H


struct CoilData
{
    CoilData() :

        constFactor{0.0},
        current{0.0},
        frequency{0.0},

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
    TYPE current;
    TYPE frequency;

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

struct VectorData
{
    VectorData() :

        x{0.0},
        y{0.0},
        z{0.0}
    {}

    TYPE x;
    TYPE y;
    TYPE z;
};

struct MatrixData
{
    MatrixData() :

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

struct ForceTorqueData

{
    ForceTorqueData() :

        forceX{0.0},
        forceY{0.0},
        forceZ{0.0},

        torqueX{0.0},
        torqueY{0.0},
        torqueZ{0.0}
    {}

    TYPE forceX;
    TYPE forceY;
    TYPE forceZ;

    TYPE torqueX;
    TYPE torqueY;
    TYPE torqueZ;
};

struct SecondaryCoilData
{
    SecondaryCoilData() :

        innerRadius{0.0},
        thickness{0.0},
        length{0.0},

        correctionFactor{0.0},

        angularIncrements{0},
        thicknessIncrements{0},
        lengthIncrements{0},

        angularPositionArray{},
        angularWeightArray{},
        thicknessPositionArray{},
        thicknessWeightArray{},
        lengthPositionArray{},
        lengthWeightArray{}
    {}

    TYPE innerRadius;
    TYPE thickness;
    TYPE length;

    TYPE correctionFactor;

    int angularIncrements;
    int lengthIncrements;
    int thicknessIncrements;

    TYPE angularPositionArray[GPU_INCREMENTS];
    TYPE angularWeightArray[GPU_INCREMENTS];
    TYPE thicknessPositionArray[GPU_INCREMENTS];
    TYPE thicknessWeightArray[GPU_INCREMENTS];
    TYPE lengthPositionArray[GPU_INCREMENTS];
    TYPE lengthWeightArray[GPU_INCREMENTS];
};

struct SecondaryCoilPositionData
{
    SecondaryCoilPositionData() :

        positionVector{},
        alphaAngle{0.0},
        betaAngle{0.0}
    {}

    TYPE positionVector[3];
    TYPE alphaAngle;
    TYPE betaAngle;
};

struct CoilPairPositionData
{
    CoilPairPositionData() :

        primPositionVector{},
        primAlphaAngle{0.0},
        primBetaAngle{0.0},

        secPositionVector{},
        secAlphaAngle{0.0},
        secBetaAngle{0.0}
    {}

    TYPE primPositionVector[3];
    TYPE primAlphaAngle;
    TYPE primBetaAngle;

    TYPE secPositionVector[3];
    TYPE secAlphaAngle;
    TYPE secBetaAngle;
};

struct CoilPairArgumentsData
{
    CoilPairArgumentsData() :

        constFactor{0.0},
        correctionFactor{0.0},
        useFastMethod{false},

        primInnerRadius{0.0},
        primThickness{0.0},
        primLength{0.0},

        primAngularIncrements{0},
        primThicknessIncrements{0},
        primLengthIncrements{0},

        secInnerRadius{0.0},
        secThickness{0.0},
        secLength{0.0},

        secAngularIncrements{0},
        secThicknessIncrements{0},
        secLengthIncrements{0},

        primThicknessPositionArray{},
        primThicknessWeightArray{},
        primCosPrecomputeArray{},
        primAngularWeightArray{},

        secAngularPositionArray{},
        secAngularWeightArray{},
        secThicknessPositionArray{},
        secThicknessWeightArray{},
        secLengthPositionArray{},
        secLengthWeightArray{}
    {}

    TYPE constFactor;
    TYPE correctionFactor;
    bool useFastMethod;

    TYPE primInnerRadius;
    TYPE primThickness;
    TYPE primLength;

    int primAngularIncrements;
    int primLengthIncrements;
    int primThicknessIncrements;

    TYPE secInnerRadius;
    TYPE secThickness;
    TYPE secLength;

    int secAngularIncrements;
    int secLengthIncrements;
    int secThicknessIncrements;

    TYPE primThicknessPositionArray[GPU_INCREMENTS];
    TYPE primThicknessWeightArray[GPU_INCREMENTS];
    TYPE primCosPrecomputeArray[GPU_INCREMENTS];
    TYPE primAngularWeightArray[GPU_INCREMENTS];

    TYPE secAngularPositionArray[GPU_INCREMENTS];
    TYPE secAngularWeightArray[GPU_INCREMENTS];
    TYPE secThicknessPositionArray[GPU_INCREMENTS];
    TYPE secThicknessWeightArray[GPU_INCREMENTS];
    TYPE secLengthPositionArray[GPU_INCREMENTS];
    TYPE secLengthWeightArray[GPU_INCREMENTS];
};

#endif //GENERAL_COIL_PROGRAM_COILDATA_H
