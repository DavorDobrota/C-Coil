#include "Benchmark.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>
#include <chrono>


void benchMathFunctions()
{
    using namespace std::chrono;

    int nOps = 200'000'000;
    double temp, interval;
    high_resolution_clock::time_point begin_time;

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= nOps; ++i)
    {
        temp += std::sqrt(i);
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("sqrt  : %.1f MOps/s\n", 1e-6 * nOps / interval);

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= nOps; ++i)
    {
        temp += std::log10(400000.0 / i);
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("log10 : %.1f MOps/s\n", 1e-6 * nOps / interval);

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= nOps; ++i)
    {
        temp += std::log10((float) (i));
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("log10f: %.1f MOps/s\n", 1e-6 * nOps / interval);

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= nOps; ++i)
    {
        temp += std::log(i);
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("log   : %.1f MOps/s\n", 1e-6 * nOps / interval);

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= nOps; ++i)
    {
        temp += std::cos( 1.0 / i);
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("cos   : %.1f MOps/s\n", 1e-6 * nOps / interval);
}
