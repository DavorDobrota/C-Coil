#ifndef COPPER_DEBUG_ERROR_CHECK
#define COPPER_DEBUG_ERROR_CHECK

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true);

#endif // COPPER_DEBUG_ERROR_CHECK