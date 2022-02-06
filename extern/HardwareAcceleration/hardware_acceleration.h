#ifndef COPPER_HARDWARE_ACCELERATED_FUNCTIONS
#define COPPER_HARDWARE_ACCELERATED_FUNCTIONS

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    #ifdef DLL_BUILD
        #define DLL __declspec(dllexport)
    #else
        #define DLL __declspec(dllimport)
    #endif
#else
    #define DLL
#endif // Win32

extern "C"
{
    DLL void Calculate_hardware_accelerated_b
    (
        long long num_ops,
        float *theta,
        float *point_distance,
        float j,
        float inner_radius,
        float length,
        float thickness,
        float inc_inner_radius,
        float inc_length,
        float inc_fi,
        float *Bh = nullptr,
        float *Bz = nullptr,
        float *B = nullptr
    );

    DLL void Calculate_hardware_accelerated_a
    (
        long long num_ops,
        float *theta,
        float *point_distance,
        float j,
        float inner_radius,
        float length,
        float thickness,
        float inc_inner_radius,
        float inc_length,
        float inc_fi,
        float *Ax = nullptr,
        float *Ay = nullptr,
        float *A = nullptr
    );
}

#endif // COPPER_HARDWARE_ACCELERATED_FUNCTIONS