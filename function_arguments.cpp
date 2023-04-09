////////////////////////////////////////////////////////////////////////////////
// File   : function_arguments.cpp
// Author : Sandeep Koranne (C) 2023. All rights reserved.
// Purpose: In CommonLisp we have map over hash to print it. Do this in C++
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cassert>
#include <functional>
#include <map>

namespace {
  using M = std::map<int,int>;
  template <typename Map>
  // Note the use of the function argument
  void PrintMap( const Map& M, std::function<void(typename Map::key_type,typename Map::mapped_type)> F )
  {
    for( auto const& [key,val] : M ) { // decomposition declaration
      F(key,val);
      std::cout << std::endl;
    }
  }
  static void Test1()
  {
    M localMap{ {1,1},{2,1},{3,2},{8,21} }; // initializer list ctor
    // local lambda for printing
    PrintMap( localMap, [](int a, int b) { std::cout << a << " " << b << " "; } );
  }
}


int main()
{
  Test1();
  return 0;
}
