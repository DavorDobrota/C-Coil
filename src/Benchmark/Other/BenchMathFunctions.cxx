#include "Benchmark.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>
#include <chrono>


void benchMathFunctions()
{
    using namespace std::chrono;

    int opCount = 200'000'000;
    double temp, interval;
    high_resolution_clock::time_point begin_time;

    printf("Performance estimate of std math functions\n\n");

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= opCount; ++i)
    {
        temp += i;
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("sum   : %.1f MOps/s\n", 1e-6 * opCount / interval);


    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= opCount; ++i)
    {
        temp += std::sqrt(i);
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("sqrt  : %.1f MOps/s\n", 1e-6 * opCount / interval);

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= opCount; ++i)
    {
        temp += std::log(i);
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("log   : %.1f MOps/s\n", 1e-6 * opCount / interval);

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= opCount; ++i)
    {
        temp += std::log10(i);
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("log10 : %.1f MOps/s\n", 1e-6 * opCount / interval);

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= opCount; ++i)
    {
        temp += std::log2((i));
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("log2: %.1f MOps/s\n", 1e-6 * opCount / interval);

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= opCount; ++i)
    {
        temp += std::sin( 1.0 / i);
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("sin   : %.1f MOps/s\n", 1e-6 * opCount / interval);

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= opCount; ++i)
    {
        temp += std::cos(1.0 / i);
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("cos   : %.1f MOps/s\n", 1e-6 * opCount / interval);

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= opCount; ++i)
    {
        temp += std::atan2(1.0, i);
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("atan2 : %.1f MOps/s\n", 1e-6 * opCount / interval);

    temp = 0.0;
    begin_time = high_resolution_clock::now();
    for (int i = 1; i <= opCount; ++i)
    {
        temp += std::exp(1.0 / i);
    }
    printf("%.15f\n", temp);
    interval = duration_cast<duration<double>>(high_resolution_clock::now() - begin_time).count();
    printf("exp   : %.1f MOps/s\n", 1e-6 * opCount / interval);
}
