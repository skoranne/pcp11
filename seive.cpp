////////////////////////////////////////////////////////////////////////////////
// File   : seive.cpp
// Author : Sandeep Koranne (C) 2017 All rights reserved.
// Purpose: Prime counter using Seive-of-Erastothenes
//
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <vector>
#include <cmath>

void seive( std::vector<bool>& ans )
{
  const size_t E = ans.size();
  //for( size_t i=0, e=ans.size(); i < e; ++i ) ans[i] = true; // all prime
  for( size_t i=2;  i < E; ++i ) {
    if( ans[i] ) {
      for( size_t j=2; i*j < E; ++j ) ans[i*j] = false;
    }
  }
}

unsigned long countPrimes( unsigned long limit )
{
  int N = limit;
  std::vector<bool> ans( N, true );
  seive( ans );
  unsigned long sum = 0;
  for( size_t i=2, e=ans.size(); i < e; ++i )
    if( ans[i] ) {
      //std::cout << i << " is a prime.\n";
      sum++;
    }
  std::cout << "There are " << sum << " primes upto " << limit << std::endl;
}

int main( int argc, char* argv[] )
{
  unsigned long limit = 2000;
  if( argc > 1 ) limit = atoi( argv[1] );
  unsigned long P = countPrimes( limit );
  return 0;
}

  
