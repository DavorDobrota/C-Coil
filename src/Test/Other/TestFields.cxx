#include "Test.h"
#include "CoilGroup.h"


void Test::testCoilGroupFieldsMTD()
{
    CoilGroup coilGroup = CoilGroup();

    const int pointCount = 20;
    const int coilCount = 10;
    auto precision = PrecisionFactor(4.0);

    vec3::Vector3Array positions(pointCount);

    for (int i = 0; i < pointCount; ++i)
        positions[i] = vec3::Vector3(0.5, 0.0, 0.1 * double(i));

    vec3::Vector3Array refValues;
    vec3::Vector3Array MTDValues;

    for (int i = 0; i < coilCount; ++i)
    {
        coilGroup.addCoil(
            0.1, 0.1, 0.1, 10000, 10,
            PrecisionFactor(), 8,
            vec3::Vector3(0.0, 0.0, 0.2 * double(i))
        );
    }

    coilGroup.setDefaultPrecisionFactor(precision);

    printf("Testing results of CoilGroup compute fields methods for CPU_ST and CPU_MT\n");
    printf("The values should be very close and the error 1e-5 or less\n\n");

    refValues = coilGroup.computeAllAPotentialVectors(positions, CPU_ST);
    MTDValues = coilGroup.computeAllAPotentialVectors(positions, CPU_MT);

    for (int i = 0; i < pointCount; ++i)
        printf("%18.15g | %18.15g : %.3g\n",
               refValues[i].abs(), MTDValues[i].abs(),
               std::abs((MTDValues[i].abs() - refValues[i].abs()) / refValues[i].abs())
        );
    printf("\n");

    refValues = coilGroup.computeAllBFieldVectors(positions, CPU_ST);
    MTDValues = coilGroup.computeAllBFieldVectors(positions, CPU_MT);

    for (int i = 0; i < pointCount; ++i)
        printf("%18.15g | %18.15g : %.3g\n",
               refValues[i].abs(), MTDValues[i].abs(),
               std::abs((MTDValues[i].abs() - refValues[i].abs()) / refValues[i].abs())
        );
    printf("\n");

    refValues = coilGroup.computeAllEFieldVectors(positions, CPU_ST);
    MTDValues = coilGroup.computeAllEFieldVectors(positions, CPU_MT);

    for (int i = 0; i < pointCount; ++i)
        printf("%18.15g | %18.15g : %.3g\n",
               refValues[i].abs(), MTDValues[i].abs(),
               std::abs((MTDValues[i].abs() - refValues[i].abs()) / refValues[i].abs())
        );
    printf("\n");

    vec3::Matrix3Array refMatrices;
    vec3::Matrix3Array MTDMatrices;

    refMatrices = coilGroup.computeAllBGradientMatrices(positions, CPU_ST);
    MTDMatrices = coilGroup.computeAllBGradientMatrices(positions, CPU_MT);

    for (int i = 0; i < pointCount; ++i)
        printf("%18.15g | %18.15g : %.3g\n",
               refMatrices[i].det(), MTDMatrices[i].det(),
               std::abs((MTDMatrices[i].det() - refMatrices[i].det()) / refMatrices[i].det())
        );
    printf("\n");
}

void Test::testCoilGroupFieldsGPU()
{
    CoilGroup coilGroup = CoilGroup();

    const int pointCount = 20;
    const int coilCount = 10;
    auto precision = PrecisionFactor(4.0);

    vec3::Vector3Array positions(pointCount);

    for (int i = 0; i < pointCount; ++i)
        positions[i] = vec3::Vector3(0.5, 0.0, 0.1 * double(i));

    vec3::Vector3Array refValues;
    vec3::Vector3Array GPUValues;

    for (int i = 0; i < coilCount; ++i)
    {
        coilGroup.addCoil(
                0.1, 0.1, 0.1, 10000, 10,
                PrecisionFactor(), 8,
                vec3::Vector3(0.0, 0.0, 0.2 * double(i))
        );
    }

    coilGroup.setDefaultPrecisionFactor(precision);

    printf("Testing results of CoilGroup compute fields methods for CPU_ST and CPU_MT\n");
    printf("The values should be very close and the error 1e-5 or less\n\n");

    refValues = coilGroup.computeAllAPotentialVectors(positions, CPU_ST);
    GPUValues = coilGroup.computeAllAPotentialVectors(positions, GPU);

    for (int i = 0; i < pointCount; ++i)
        printf("%18.15g | %18.15g : %.3g\n",
               refValues[i].abs(), GPUValues[i].abs(),
               std::abs((GPUValues[i].abs() - refValues[i].abs()) / refValues[i].abs())
        );
    printf("\n");

    refValues = coilGroup.computeAllBFieldVectors(positions, CPU_ST);
    GPUValues = coilGroup.computeAllBFieldVectors(positions, GPU);

    for (int i = 0; i < pointCount; ++i)
        printf("%18.15g | %18.15g : %.3g\n",
               refValues[i].abs(), GPUValues[i].abs(),
               std::abs((GPUValues[i].abs() - refValues[i].abs()) / refValues[i].abs())
        );
    printf("\n");

    refValues = coilGroup.computeAllEFieldVectors(positions, CPU_ST);
    GPUValues = coilGroup.computeAllEFieldVectors(positions, GPU);

    for (int i = 0; i < pointCount; ++i)
        printf("%18.15g | %18.15g : %.3g\n",
               refValues[i].abs(), GPUValues[i].abs(),
               std::abs((GPUValues[i].abs() - refValues[i].abs()) / refValues[i].abs())
        );
    printf("\n");

    vec3::Matrix3Array refMatrices;
    vec3::Matrix3Array GPUMatrices;

    refMatrices = coilGroup.computeAllBGradientMatrices(positions, CPU_ST);
    GPUMatrices = coilGroup.computeAllBGradientMatrices(positions, GPU);

    for (int i = 0; i < pointCount; ++i)
        printf("%18.15g | %18.15g : %.3g\n",
               refMatrices[i].det(), GPUMatrices[i].det(),
               std::abs((GPUMatrices[i].det() - refMatrices[i].det()) / refMatrices[i].det())
        );
    printf("\n");
}
