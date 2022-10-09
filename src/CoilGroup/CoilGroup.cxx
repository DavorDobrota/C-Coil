#include "CoilGroup.h"
#include "ThreadPool.h"

#include <sstream>

#define _USE_MATH_DEFINES
#include <math.h>
#include <utility>
#include <functional>

using namespace std::placeholders;


namespace
{
    threadPool::ThreadPoolControl g_threadPool;
}


CoilGroup::CoilGroup(std::vector<std::shared_ptr<Coil>> memberCoils, PrecisionFactor precisionFactor, int threadCount) :
            memberCoils(std::move(memberCoils)), defaultPrecisionFactor(precisionFactor), threadCount(threadCount) {}


PrecisionFactor CoilGroup::getDefaultPrecisionFactor() const { return defaultPrecisionFactor; }

int CoilGroup::getThreadCount() const { return threadCount; }

const std::vector<std::shared_ptr<Coil>> &CoilGroup::getMemberCoils() const { return memberCoils; }


void CoilGroup::setDefaultPrecisionFactor(PrecisionFactor precisionFactor)
{
    defaultPrecisionFactor = precisionFactor;

    for (auto& memberCoil : memberCoils)
        memberCoil->setDefaultPrecision(precisionFactor);
}

void CoilGroup::setThreadCount(int threadCount)
{
    this->threadCount = threadCount;

    for (auto& memberCoil : memberCoils)
        memberCoil->setThreadCount(threadCount);
}

void CoilGroup::addCoil(double innerRadius, double thickness, double length, int numOfTurns, double current,
                        PrecisionFactor precisionFactor, int coilThreads, vec3::Vector3 coordinatePosition,
                        double yAxisAngle, double zAxisAngle)
{
    this->memberCoils.push_back(
        std::make_shared<Coil>(
            innerRadius, thickness, length, numOfTurns, current,
            precisionFactor, coilThreads, coordinatePosition, yAxisAngle, zAxisAngle
        )
    );
}

void CoilGroup::removeCoil(size_t index)
{
    if(index >= memberCoils.size())
        throw std::out_of_range("Coil index out of range!");

    memberCoils.erase(memberCoils.begin() + index);
}

Coil& CoilGroup::operator[](size_t index) const
{
    if(index >= memberCoils.size())
        throw std::out_of_range("Coil index out of range!");

    return static_cast<Coil &>(*memberCoils[index]);
}


bool CoilGroup::isPointInside(vec3::Vector3 pointVector)
{
    for (auto& memberCoil : memberCoils)
        if (memberCoil->isPointInside(pointVector))
            return true;
    return false;
}


CoilGroup::operator std::string() const
{
    std::stringstream output;

    auto stringifyVector = [](auto &ar) -> std::string
    {
        std::stringstream output;

        output << "[";

        for(int i = 0; i < ar.size(); i++)
        {
            if(i != 0)
                output << ", ";
            output << std::string(*ar[i]);
        }

        output << "]";

        return output.str();
    };

    output << "CoilGroup("
           << "member_coils=" << stringifyVector(memberCoils)
           << ", default_precision_factor=" << std::string(defaultPrecisionFactor)
           << ", thread_count=" << threadCount
           << ")";

    return output.str();
}
