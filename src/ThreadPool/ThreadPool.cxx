#include "ThreadPool.h"

#include <ctype.h>

#include "ctpl_stl.h"


namespace threadPool 
{
    ctpl::thread_pool threadPool;

    size_t ThreadPoolControl::getTaskCount() const { return taskCount; }
    void ThreadPoolControl::setTaskCount(size_t taskCount){ this->taskCount = taskCount; }

    std::atomic_int_fast64_t &ThreadPoolControl::getCompletedTasks() const { return const_cast<std::atomic_int_fast64_t&>(completedTasks); }
    void ThreadPoolControl::setCompletedTasks(int completedTasks){ this->completedTasks = completedTasks; }

    void ThreadPoolControl::synchronizeThreads() const { while(completedTasks.load() < taskCount){}; }

    int ThreadPoolControl::getSize(){ return threadPool.size(); }
    void ThreadPoolControl::setSize(int threadCount){ threadPool.resize(threadCount); }
}
