#ifndef GENERAL_COIL_PROGRAM_THREADPOOL_H
#define GENERAL_COIL_PROGRAM_THREADPOOL_H

#include "ctpl_stl.h"

#include <ctype.h>


/// @brief Contains a wrapper for the CTPL library (https://github.com/vit-vit/CTPL) which is used for multithreading.
namespace threadPool
{
    extern ctpl::thread_pool threadPool;

    class ThreadPoolControl
    {
    public:
        size_t getTaskCount() const;
        void setTaskCount(size_t taskCount);

        std::atomic_int_fast64_t &getCompletedTasks() const;
        void setCompletedTasks(int completedTasks);

        void synchronizeThreads() const;

        int getSize();
        void setSize(int threadCount);

        template<class ... Args>
        void push(Args&& ...args) { threadPool.push(args...); }

    private:
        std::atomic_int_fast64_t completedTasks;
        size_t taskCount;
    };
}

#endif //GENERAL_COIL_PROGRAM_THREADPOOL_H
