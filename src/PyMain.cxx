#include <PyMain.h>


PYBIND11_MODULE(c_coil, m)
{
    initTensor(m);
    initCoil(m);
    initCoilGroup(m);
    initBenchmark(m);
    initCompare(m);
}
