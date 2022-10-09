#include "GPUMemoryManagement.h"
#include "CUDAUtils/ErrorCheck/CUDAErrorCheck.h"


struct GPUBuffer {
    long long size;
    void *buffer;
};

namespace
{
    std::vector<GPUBuffer> buffers;
}

void generateNewBuffers(std::vector<long long> bufferSizes, long long elementSize)
{
    GPUMem::freeMemory();
    buffers.resize(bufferSizes.size());

    for(long long i = 0; i < buffers.size(); i++) {
        gpuErrchk(cudaMalloc(&buffers[i].buffer, bufferSizes[i] * elementSize));
        buffers[i].size = bufferSizes[i] * elementSize;
    }
}


namespace GPUMem
{
    void freeMemory()
    {
        for(auto &buff : buffers) {
            gpuErrchk(cudaFree(buff.buffer));
        }
        buffers.resize(0);
    }

    std::vector<void*> getBuffers(std::vector<long long> bufferSizes, long long elementSize)
    {
        std::vector<void*> ret;

        if(buffers.size() < bufferSizes.size()) {
            generateNewBuffers(bufferSizes, elementSize);
        } else {
            for(long long i = 0; i < buffers.size(); i++) {
                if(buffers[i].size < bufferSizes[i] * elementSize) {
                    generateNewBuffers(bufferSizes, elementSize);
                    break;
                }
            }
        }

        for(auto &buff : buffers) {
            ret.push_back(buff.buffer);
        }

        return ret;
    }
}
