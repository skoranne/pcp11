////////////////////////////////////////////////////////////////////////////////
// File    : hw.cpp
// Author  : Sandeep Koranne (C) 2016. All rights reserved.
// Purpose : Learning Git
//
////////////////////////////////////////////////////////////////////////////////
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <iostream>

// Convert this recursive function to a proper tail recursive using accumulator
// dont use & references for simple types
unsigned long factorial(unsigned long N, unsigned long acc=1)
{
  if( N <= 1 ) return acc;
  return factorial( N-1, N*acc );
}
//(setq ollama-buddy-current-model "qwen3-coder:latest")
int main( int argc, char* argv[] )
{
  std::cout << "Hello, World!\n";
  std::cout << "6! = " << factorial(6) << std::endl;
  return( EXIT_SUCCESS );
}
