#include "Test.h"
#include "Coil.h"
#include "CoilGroup.h"


void testCoilMInductanceArrangements()
{
    Coil prim = Coil(0.1, 0.1, 0.1, 10000);
    Coil sec = Coil(0.1, 0.1, 0.1, 10000);

    const int configCount = 20;
    auto precision = PrecisionFactor(5.0);

    vec3::Vector3Array primPositions(configCount);
    vec3::Vector3Array secPositions(configCount);
    std::vector<double> primYAxisAngle(configCount);
    std::vector<double> primZAxisAngle(configCount);
    std::vector<double> secYAxisAngle(configCount);
    std::vector<double> secZAxisAngle(configCount);

    for (int i = 0; i < configCount; ++i)
    {
        primPositions[i] = vec3::Vector3(0.0, 0.0, 0.0);
        secPositions[i] = vec3::Vector3(0.0, 0.1, 0.2 + double(i) * 0.01);
        primYAxisAngle[i] = 0.5;
        primZAxisAngle[i] = 0.5;
        secYAxisAngle[i] = 0.6;
        secZAxisAngle[i] = 0.6;
    }

    printf("Testing results of Coil::computeAllMutualInductanceArrangements for CPU and GPU\n");
    printf("The values should be very close and the error 1e-5 or less\n\n");

    std::vector<double> mutualInductanceCPU(configCount);
    std::vector<double> mutualInductanceGPU(configCount);

    mutualInductanceCPU = Coil::computeAllMutualInductanceArrangements(
        prim, sec, primPositions,secPositions,
        primYAxisAngle, primZAxisAngle,secYAxisAngle, secZAxisAngle,
        precision, CPU_ST
    );

    mutualInductanceGPU = Coil::computeAllMutualInductanceArrangements(
            prim, sec, primPositions,secPositions,
            primYAxisAngle, primZAxisAngle,secYAxisAngle, secZAxisAngle,
            precision, GPU
    );

    for (int i = 0; i < configCount; ++i)
        printf("%18.15g | %18.15g : %.3g\n",
               mutualInductanceCPU[i], mutualInductanceGPU[i],
               std::abs((mutualInductanceGPU[i] - mutualInductanceCPU[i]) / mutualInductanceCPU[i])
        );
}

void testCoilForceArrangements()
{
    Coil prim = Coil(0.1, 0.1, 0.1, 10000);
    Coil sec = Coil(0.1, 0.1, 0.1, 10000);

    int configCount = 20;
    auto precision = PrecisionFactor(5.0);

    vec3::Vector3Array primPositions(configCount);
    vec3::Vector3Array secPositions(configCount);
    std::vector<double> primYAxisAngle(configCount);
    std::vector<double> primZAxisAngle(configCount);
    std::vector<double> secYAxisAngle(configCount);
    std::vector<double> secZAxisAngle(configCount);

    for (int i = 0; i < configCount; ++i)
    {
        primPositions[i] = vec3::Vector3(0.0, 0.0, 0.0);
        secPositions[i] = vec3::Vector3(0.0, 0.1, 0.2 + double(i) * 0.01);
        primYAxisAngle[i] = 0.5;
        primZAxisAngle[i] = 0.5;
        secYAxisAngle[i] = 0.6;
        secZAxisAngle[i] = 0.6;
    }

    printf("Testing results of Coil::computeAllAmpereForceArrangements for CPU and GPU\n");
    printf("The values should be very close and the error 1e-5 or less\n\n");

    std::vector<std::pair<vec3::Vector3, vec3::Vector3>> forceTorqueCPU(configCount);
    std::vector<std::pair<vec3::Vector3, vec3::Vector3>> forceTorqueGPU(configCount);

    forceTorqueCPU = Coil::computeAllAmpereForceArrangements(
            prim, sec, primPositions,secPositions,
            primYAxisAngle, primZAxisAngle,secYAxisAngle, secZAxisAngle,
            precision, CPU_ST
    );

    forceTorqueGPU = Coil::computeAllAmpereForceArrangements(
            prim, sec, primPositions,secPositions,
            primYAxisAngle, primZAxisAngle,secYAxisAngle, secZAxisAngle,
            precision, GPU
    );

    for (int i = 0; i < configCount; ++i)
        printf("%18.15g | %18.15g : %.3g\n",
               forceTorqueCPU[i].first.abs(), forceTorqueGPU[i].first.abs(),
               std::abs((forceTorqueGPU[i].first.abs() - forceTorqueCPU[i].first.abs()) / forceTorqueCPU[i].first.abs())
        );
}

void testGroupMInductanceArrangements()
{
    CoilGroup coilGroup = CoilGroup();

    const int configCount = 20;
    const int coilCount = 10;
    auto precision = PrecisionFactor(4.0);

    vec3::Vector3Array secPositions(configCount);
    std::vector<double> secYAxisAngle(configCount);
    std::vector<double> secZAxisAngle(configCount);

    for (int i = 0; i < coilCount; ++i)
    {
        Coil tempCoil = Coil(0.1, 0.1, 0.1, 10000, 10);
        tempCoil.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.2 * double(i)));

        coilGroup.addCoil(tempCoil);
    }

    coilGroup.setDefaultPrecisionFactor(PrecisionFactor(5.0));

    Coil secondary = Coil(0.1, 0.1, 0.1, 10000, 5);

    for (int i = 0; i < configCount; ++i)
    {
        secPositions[i] = vec3::Vector3(0.5, 0.0, 0.1 * double(i));
        secYAxisAngle[i] = 0.5;
        secZAxisAngle[i] = 0.6;
    }

    printf("Testing results of CoilGroup::computeAllMutualInductanceArrangements for CPU and GPU\n");
    printf("The values should be very close and the error 1e-5 or less\n\n");

    std::vector<double> mutualInductanceCPU(configCount);
    std::vector<double> mutualInductanceGPU(configCount);

    mutualInductanceCPU = coilGroup.computeAllMutualInductanceArrangements(
            secondary, secPositions,secYAxisAngle, secZAxisAngle,
            precision, CPU_ST
    );

    mutualInductanceGPU = coilGroup.computeAllMutualInductanceArrangements(
            secondary, secPositions,secYAxisAngle, secZAxisAngle,
            precision, GPU
    );

    for (int i = 0; i < configCount; ++i)
        printf("%18.15g | %18.15g : %.3g\n",
               mutualInductanceCPU[i], mutualInductanceGPU[i],
               std::abs((mutualInductanceGPU[i] - mutualInductanceCPU[i]) / mutualInductanceCPU[i])
        );
}

void testGroupForceArrangements()
{
    CoilGroup coilGroup = CoilGroup();

    const int configCount = 20;
    const int coilCount = 10;
    auto precision = PrecisionFactor(4.0);

    vec3::Vector3Array secPositions(configCount);
    std::vector<double> secYAxisAngle(configCount);
    std::vector<double> secZAxisAngle(configCount);

    for (int i = 0; i < coilCount; ++i)
    {
        Coil tempCoil = Coil(0.1, 0.1, 0.1, 10000, 10);
        tempCoil.setPositionAndOrientation(vec3::Vector3(0.0, 0.0, 0.2 * double(i)));

        coilGroup.addCoil(tempCoil);
    }

    coilGroup.setDefaultPrecisionFactor(PrecisionFactor(5.0));

    Coil secondary = Coil(0.1, 0.1, 0.1, 10000, 5);

    for (int i = 0; i < configCount; ++i)
    {
        secPositions[i] = vec3::Vector3(0.5, 0.0, 0.1 * double(i));
        secYAxisAngle[i] = 0.5;
        secZAxisAngle[i] = 0.6;
    }

    printf("Testing results of CoilGroup::computeAllMutualInductanceArrangements for CPU and GPU\n");
    printf("The values should be very close and the error 1e-5 or less\n\n");

    std::vector<std::pair<vec3::Vector3, vec3::Vector3>> forceTorqueCPU(configCount);
    std::vector<std::pair<vec3::Vector3, vec3::Vector3>> forceTorqueGPU(configCount);

    forceTorqueCPU = coilGroup.computeAllAmpereForceArrangements(
            secondary, secPositions,secYAxisAngle, secZAxisAngle,
            precision, CPU_ST
    );

    forceTorqueGPU = coilGroup.computeAllAmpereForceArrangements(
            secondary, secPositions,secYAxisAngle, secZAxisAngle,
            precision, GPU
    );

    for (int i = 0; i < configCount; ++i)
        printf("%18.15g | %18.15g : %.3g\n",
               forceTorqueCPU[i].first.abs(), forceTorqueGPU[i].first.abs(),
               std::abs((forceTorqueGPU[i].first.abs() - forceTorqueCPU[i].first.abs()) / forceTorqueCPU[i].first.abs())
        );
}