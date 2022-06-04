#include <PyMain.h>


PYBIND11_MODULE(coil_evolution, m)
{
    initTensor(m);
    initCoil(m);
    initCoilGroup(m);
    initBenchmark(m);
    initCompare(m);
}
