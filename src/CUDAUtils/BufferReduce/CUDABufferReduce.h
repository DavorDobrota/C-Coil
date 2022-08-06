#ifndef COIL_EVOLUTION_CUDABUFFERREDUCE_H
#define COIL_EVOLUTION_CUDABUFFERREDUCE_H


namespace
{
    __device__
    void warpReduce(volatile TYPE *buffer, unsigned int idx)
    {
        if (NTHREADS >= 64) { if (idx < 32) { buffer[idx] += buffer[idx + 32]; } }
        if (NTHREADS >= 32) buffer[idx] += buffer[idx + 16];
        if (NTHREADS >= 16) buffer[idx] += buffer[idx + 8];
        if (NTHREADS >= 8) buffer[idx] += buffer[idx + 4];
        if (NTHREADS >= 4) buffer[idx] += buffer[idx + 2];
        if (NTHREADS >= 2) buffer[idx] += buffer[idx + 1];
    }

    __device__
    void bufferReduce(volatile TYPE *buffer, unsigned int idx)
    {
        if (NTHREADS >= 2048) { if (idx < 1024) { buffer[idx] += buffer[idx + 1024]; } __syncthreads(); }
        if (NTHREADS >= 1024) { if (idx < 512) { buffer[idx] += buffer[idx + 512]; } __syncthreads(); }
        if (NTHREADS >= 512) { if (idx < 256) { buffer[idx] += buffer[idx + 256]; } __syncthreads(); }
        if (NTHREADS >= 256) { if (idx < 128) { buffer[idx] += buffer[idx + 128]; } __syncthreads(); }
        if (NTHREADS >= 128) { if (idx < 64) { buffer[idx] += buffer[idx + 64]; } __syncthreads(); }
        warpReduce(buffer, idx);
    }
}

#endif //COIL_EVOLUTION_CUDABUFFERREDUCE_H
