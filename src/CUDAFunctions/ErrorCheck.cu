#include "CUDAErrorCheck.h"

#include <cstdio>

void gpuAssert(cudaError_t code, const char *file, int line, bool abort)
{
    if(code != cudaSuccess)
    {
        fprintf
        (
            stderr,
            "GPUassert: %s %s %d\n",
            cudaGetErrorString(code), file, line
        );
        if (abort) exit(code);
    }
}

// TODO: fix this
//void printMemStats(bool pause)
//{
//    size_t free, total;
//    
//    cudaMemGetInfo(&free, &total);
//    printf
//    (
//        "%lf GB free of total %lf GB\n",
//        (double)free / 1000000000, (double)total / 1000000000
//    );
//    
//    if(pause)
//        consolePause();
//}
