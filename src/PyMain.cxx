#include <PyMain.h>


PYBIND11_MODULE(coil_evolution, m)
{
    initTensor(m);
    initCoilData(m);
    initCoilGroup(m);
}
