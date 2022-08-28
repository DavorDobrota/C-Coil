#include "CoilGroup.h"
#include "ThreadPool.h"

#include <functional>
#include <atomic>


namespace
{
    threadPool::ThreadPoolControl g_threadPool;
}


vec3::Vector3Array CoilGroup::calculateAllAPotentialMTD(const vec3::Vector3Array &pointVectors) const
{
    g_threadPool.setTaskCount(memberCoils.size());
    g_threadPool.getCompletedTasks().store(0ull);

    std::vector<std::atomic<double>> resultArrX(pointVectors.size());
    std::vector<std::atomic<double>> resultArrY(pointVectors.size());
    std::vector<std::atomic<double>> resultArrZ(pointVectors.size());

    auto calcThread = []
    (
            int idx,
            const Coil &coil,
            const vec3::Vector3Array &pointVectors,
            std::vector<std::atomic<double>> &valuesX,
            std::vector<std::atomic<double>> &valuesY,
            std::vector<std::atomic<double>> &valuesZ
    ){
        for (int i = 0; i < pointVectors.size(); ++i)
        {
            vec3::Vector3 temp = coil.computeAPotentialVector(pointVectors[i]);
            valuesX[i].store(valuesX[i].load() + temp.x);
            valuesY[i].store(valuesY[i].load() + temp.y);
            valuesZ[i].store(valuesZ[i].load() + temp.z);
        }
        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (const auto & memberCoil : memberCoils)
    {
        g_threadPool.push
        (
            calcThread,
            std::ref(*memberCoil),
            std::ref(pointVectors),
            std::ref(resultArrX),
            std::ref(resultArrY),
            std::ref(resultArrZ)
        );
    }

    g_threadPool.synchronizeThreads();

    vec3::Vector3Array output;
    output.resize(pointVectors.size());

    for (int i = 0; i < pointVectors.size(); ++i)
        output.append(resultArrX[i], resultArrY[i], resultArrZ[i]);

    return output;
}

vec3::Vector3Array CoilGroup::calculateAllBFieldMTD(const vec3::Vector3Array &pointVectors) const
{
    g_threadPool.setTaskCount(memberCoils.size());
    g_threadPool.getCompletedTasks().store(0ull);

    std::vector<std::atomic<double>> resultArrX(pointVectors.size());
    std::vector<std::atomic<double>> resultArrY(pointVectors.size());
    std::vector<std::atomic<double>> resultArrZ(pointVectors.size());

    auto calcThread = []
    (
            int idx,
            const Coil &coil,
            const vec3::Vector3Array &pointVectors,
            std::vector<std::atomic<double>> &valuesX,
            std::vector<std::atomic<double>> &valuesY,
            std::vector<std::atomic<double>> &valuesZ
    ){
        for (int i = 0; i < pointVectors.size(); ++i)
        {
            vec3::Vector3 temp = coil.computeBFieldVector(pointVectors[i]);
            valuesX[i].store(valuesX[i].load() + temp.x);
            valuesY[i].store(valuesY[i].load() + temp.y);
            valuesZ[i].store(valuesZ[i].load() + temp.z);
        }
        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (const auto & memberCoil : memberCoils)
    {
        g_threadPool.push
        (
            calcThread,
            std::ref(*memberCoil),
            std::ref(pointVectors),
            std::ref(resultArrX),
            std::ref(resultArrY),
            std::ref(resultArrZ)
        );
    }
    g_threadPool.synchronizeThreads();

    vec3::Vector3Array output;
    output.resize(pointVectors.size());

    for (int i = 0; i < pointVectors.size(); ++i)
        output.append(resultArrX[i], resultArrY[i], resultArrZ[i]);

    return output;
}

vec3::Vector3Array CoilGroup::calculateAllEFieldMTD(const vec3::Vector3Array &pointVectors) const
{
    g_threadPool.setTaskCount(memberCoils.size());
    g_threadPool.getCompletedTasks().store(0ull);

    std::vector<std::atomic<double>> resultArrX(pointVectors.size());
    std::vector<std::atomic<double>> resultArrY(pointVectors.size());
    std::vector<std::atomic<double>> resultArrZ(pointVectors.size());

    auto calcThread = []
    (
            int idx,
            const Coil &coil,
            const vec3::Vector3Array &pointVectors,
            std::vector<std::atomic<double>> &valuesX,
            std::vector<std::atomic<double>> &valuesY,
            std::vector<std::atomic<double>> &valuesZ
    ){
        for (int i = 0; i < pointVectors.size(); ++i)
        {
            vec3::Vector3 temp = coil.computeEFieldVector(pointVectors[i]);
            valuesX[i].store(valuesX[i].load() + temp.x);
            valuesY[i].store(valuesY[i].load() + temp.y);
            valuesZ[i].store(valuesZ[i].load() + temp.z);
        }
        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (const auto & memberCoil : memberCoils)
    {
        g_threadPool.push
        (
            calcThread,
            std::ref(*memberCoil),
            std::ref(pointVectors),
            std::ref(resultArrX),
            std::ref(resultArrY),
            std::ref(resultArrZ)
        );
    }
    g_threadPool.synchronizeThreads();

    vec3::Vector3Array output;
    output.resize(pointVectors.size());

    for (int i = 0; i < pointVectors.size(); ++i)
        output.append(resultArrX[i], resultArrY[i], resultArrZ[i]);

    return output;
}

vec3::Matrix3Array CoilGroup::calculateAllBGradientMTD(const vec3::Vector3Array &pointVectors) const
{
    g_threadPool.setTaskCount(memberCoils.size());
    g_threadPool.getCompletedTasks().store(0ull);

    std::vector<std::atomic<double>> resultArrXX(pointVectors.size());
    std::vector<std::atomic<double>> resultArrXY(pointVectors.size());
    std::vector<std::atomic<double>> resultArrXZ(pointVectors.size());

    std::vector<std::atomic<double>> resultArrYX(pointVectors.size());
    std::vector<std::atomic<double>> resultArrYY(pointVectors.size());
    std::vector<std::atomic<double>> resultArrYZ(pointVectors.size());

    std::vector<std::atomic<double>> resultArrZX(pointVectors.size());
    std::vector<std::atomic<double>> resultArrZY(pointVectors.size());
    std::vector<std::atomic<double>> resultArrZZ(pointVectors.size());

    auto calcThread = []
    (
            int idx,
            const Coil &coil,
            const vec3::Vector3Array &pointVectors,
            std::vector<std::atomic<double>> &valuesXX,
            std::vector<std::atomic<double>> &valuesXY,
            std::vector<std::atomic<double>> &valuesXZ,
            std::vector<std::atomic<double>> &valuesYX,
            std::vector<std::atomic<double>> &valuesYY,
            std::vector<std::atomic<double>> &valuesYZ,
            std::vector<std::atomic<double>> &valuesZX,
            std::vector<std::atomic<double>> &valuesZY,
            std::vector<std::atomic<double>> &valuesZZ
    ){
        for (int i = 0; i < pointVectors.size(); ++i)
        {
            vec3::Matrix3 temp = coil.computeBGradientMatrix(pointVectors[i]);

            valuesXX[i].store(valuesXX[i].load() + temp.xx);
            valuesXY[i].store(valuesXY[i].load() + temp.xy);
            valuesXZ[i].store(valuesXZ[i].load() + temp.xz);

            valuesYX[i].store(valuesYX[i].load() + temp.yx);
            valuesYY[i].store(valuesYY[i].load() + temp.yy);
            valuesYZ[i].store(valuesYZ[i].load() + temp.yz);

            valuesZX[i].store(valuesZX[i].load() + temp.zx);
            valuesZY[i].store(valuesZY[i].load() + temp.zy);
            valuesZZ[i].store(valuesZZ[i].load() + temp.zz);
        }
        g_threadPool.getCompletedTasks().fetch_add(1ull);
    };

    for (const auto & memberCoil : memberCoils)
    {
        g_threadPool.push
            (
                calcThread,
                std::ref(*memberCoil),
                std::ref(pointVectors),
                std::ref(resultArrXX),
                std::ref(resultArrXY),
                std::ref(resultArrXZ),
                std::ref(resultArrYX),
                std::ref(resultArrYY),
                std::ref(resultArrYZ),
                std::ref(resultArrZX),
                std::ref(resultArrZY),
                std::ref(resultArrZZ)
            );
    }
    g_threadPool.synchronizeThreads();

    vec3::Matrix3Array output;
    output.resize(pointVectors.size());

    for (int i = 0; i < pointVectors.size(); ++i)
        output.append(resultArrXX[i], resultArrXY[i], resultArrXZ[i],
                      resultArrYX[i], resultArrYY[i], resultArrYZ[i],
                      resultArrZX[i], resultArrZY[i], resultArrZZ[i]);

    return output;
}
