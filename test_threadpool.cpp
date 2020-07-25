////////////////////////////////////////////////////////////////////////////////
// File   : threadpool.cpp
// Author : Sandeep Koranne (C) 2020. All rights reserved.
// Purpose: Simple thread pool for futures
//
////////////////////////////////////////////////////////////////////////////////

#include "threadpool.h"
#include <iostream>
#include <vector>
#include <memory>
#include <gmpxx.h>

using namespace THREAD_POOL;

int collatz_count( int N ) {
  mpz_class GN = N;
  int ORIG_N = N;
  int retval = 1;
  while( GN != 1 ) {
    if( GN < 0 ) std::cout << "ORIG N overflows: " << ORIG_N << "\t" << GN << std::endl;
    assert( GN > 0 ); // must be overflowing
    if( ( GN % 2 ) == 0 ) GN = GN /2;
    else GN = 3*GN+1;
    retval++;
    if( GN == 0 ) { std::cout << ORIG_N << std::endl; break;}
  }
  return retval;
}

struct CollatzOperator
{
  unsigned int NUM;
  unsigned int count;
  explicit CollatzOperator(unsigned int i=1): NUM{i}, count{0} {}
  void operator()() { count = collatz_count( NUM ); }
  std::function<void()> getLambda() { return [this]()->void { this->operator()(); }; }
  //~CollatzOperator() { std::cout << __FUNCTION__ << std::endl; }
};


static void TestThreadPool(unsigned int N, unsigned int M)
{
  const unsigned int NUM_THREADS = std::thread::hardware_concurrency();
  std::cout << "System supports " << NUM_THREADS << " threads.\n";
  unsigned int ACTUAL_THREAD_COUNT = std::min( N, NUM_THREADS );
  ThreadPool mypool{ACTUAL_THREAD_COUNT};
  std::vector< std::thread > my_threads;
  for( unsigned int i=0; i < ACTUAL_THREAD_COUNT; ++i ) {
    my_threads.push_back( std::thread( &ThreadPool::run, &mypool ) );
  }
  std::vector< CollatzOperator > jobs(M);
  for( unsigned int i=0; i < M; ++i ) {
    jobs[i].NUM = i+1;
    mypool.add( jobs[i].getLambda() );
    // Bind is not preferable as then we have to use std::shared_ptr
    //mypool.add( std::bind( &CollatzOperator::operator(), jobs[i] ) );
  }
  mypool.complete();
  for( unsigned int i=0; i < ACTUAL_THREAD_COUNT; ++i ) {
    my_threads[i].join();
  }
  std::pair<unsigned int, unsigned int> MAX_COLLATZ;
  MAX_COLLATZ.first = jobs[0].NUM;
  MAX_COLLATZ.second= jobs[0].count;
  for( const auto& x : jobs ) {
    if( ( x.NUM % 10000 ) == 0 ) std::cout << x.NUM << "\t" << x.count << "\n";
    if( MAX_COLLATZ.second < x.count ) {
      MAX_COLLATZ.first = x.NUM;
      MAX_COLLATZ.second= x.count;
    }
  }
  std::cout << "MAX COLLATZ = " << MAX_COLLATZ.first << "\t" << MAX_COLLATZ.second << std::endl;
}


static void Usage( const char* progName )
{
  std::cout << progName << " N (threads) M (number-trials) " << std::endl;
  exit(-1);
}
int main(int argc, char* argv[])
{
  int N = 10, M=10000;
  if( argc == 3 ) 
    N = atoi( argv[1] ), M = atoi( argv[2] );
  else
    Usage( argv[0] );
  TestThreadPool(N,M);
  return (0);
}
