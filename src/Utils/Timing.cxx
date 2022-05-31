#include "Timing.h"

#include <chrono>
#include <vector>

using namespace std::chrono;


namespace
{
    std::vector<high_resolution_clock::time_point> g_start_vec;
    high_resolution_clock::time_point g_stop;
    duration<double> g_interval_duration;
}

void recordStartPoint()
{
    g_start_vec.push_back(high_resolution_clock::now());
}

double getIntervalDuration()
{
    g_stop = high_resolution_clock::now();
    g_interval_duration = duration_cast<duration<double>>
    (
        g_stop - g_start_vec.back()
    );
    
    g_start_vec.pop_back();
    
    return g_interval_duration.count();
}
