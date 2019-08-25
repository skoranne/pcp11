////////////////////////////////////////////////////////////////////////////////
// File   : computing_with_array.cpp
// Author : Sandeep Koranne (C) 2019. All rights reserved.
//
// Purpose: Test of valarray
//
////////////////////////////////////////////////////////////////////////////////
#include <valarray>
#include <cassert>
#include <iostream>

static void TestValArray()
{
  using VALUES = std::valarray<double>;
  constexpr int N = 10;
  VALUES x(N);
  for( int i=0; i < N; ++i ) {
    x[i] = i-5;
  }
  VALUES z = 2.0*x;
  VALUES y=sin(x);
  x = 2.0*x+3.0*y;
  for( auto ix: z ) { std::cout << ix << " "; }
}

int main()
{
  TestValArray();

  return 0;
}

