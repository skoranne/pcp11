////////////////////////////////////////////////////////////////////////////////
// File   : threadpool.h
// Author : Sandeep Koranne (C) 2020. All rights reserved.
// Purpose: Simple thread pool for futures
//
////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <functional>
#include <atomic>
#include <queue>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <future>

#pragma once

namespace THREAD_POOL {

  class ThreadPool {
  private:
    const unsigned int MAX_THREADS = 8;
    std::queue< std::function<void()> > task_queue;
    std::mutex mutex;
    std::condition_variable condition_variable;
    std::atomic<bool> work;
  public:
    ThreadPool(unsigned int MaxThreads): MAX_THREADS{MaxThreads},work{true} {}
    ThreadPool( const ThreadPool& ) = delete;
    ThreadPool& operator=( const ThreadPool& ) = delete;
    ThreadPool( ThreadPool&& ) = delete;

    void add( std::function<void()> f) {
      assert( work && "Queue is shut down." );
      {
	std::unique_lock<std::mutex> lock{ mutex };
	task_queue.emplace( f );
      }
      condition_variable.notify_one();
    }
      
    void complete() {
      {
	std::unique_lock<std::mutex> lock{ mutex };
	work = false;
      }
      condition_variable.notify_all();
    }

    void run() {
      while( true ) {
	std::function<void()> func;
	{
	  std::unique_lock<std::mutex> lock{mutex};
	  condition_variable.wait(lock, [this]() {return !task_queue.empty() || !work; });
	  if (!work && task_queue.empty())
	    {
	      return;
	    }
	  func = task_queue.front();
	  task_queue.pop();
	}
	func();
      }
    }
      
  };

} // end of namespace THREAD_POOL
