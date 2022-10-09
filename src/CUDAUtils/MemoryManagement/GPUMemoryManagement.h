#ifndef COIL_EVOLUTION_GPUMEMORYMANAGEMENT_H
#define COIL_EVOLUTION_GPUMEMORYMANAGEMENT_H

#include <vector>


namespace GPUMem
{
    void freeMemory();

    std::vector<void*> getBuffers(std::vector<long long> bufferSizes, long long elementSize = 1ll);

    template<typename T>
    std::vector<T*> getBuffers(std::vector<long long> bufferSizes)
    {
        std::vector<T*> ret;

        for(void *ptr : getBuffers(std::move(bufferSizes), sizeof(T))) {
            ret.push_back(static_cast<T*>(ptr));
        }

        return ret;
    }
}

#endif //COIL_EVOLUTION_GPUMEMORYMANAGEMENT_H
