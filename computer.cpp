////////////////////////////////////////////////////////////////////////////////
// File   : computer.cpp
// Author : Sandeep Koranne (C) 2021. All rights reserved.
// Purpose: Use [Optional], Any and Variant for interpreter.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <complex>
#include <ratio>
#include <any>
#include <optional>
#include <functional>
#include <variant>
#include <tuple>

namespace Computer
{
  using namespace std::complex_literals; // this allows us to do 1.+2i
  struct Overflow {};
  using Complex = std::complex<double>;
  using ComputeTypes = std::variant<long,double,Complex>;
  using Types = std::optional< ComputeTypes >;
  constexpr ComputeTypes ZEROL{long{0}};
  constexpr ComputeTypes ZEROC{Complex{0}};
  Complex getValue( ComputeTypes a )
  {
    Complex x{0};
    std::visit( [&x,a]( auto&& arg ) {
		  using T = std::decay_t<decltype(arg)>;
		  if constexpr (std::is_same_v<T,Complex>) x = std::get<Complex>(a);
                  if constexpr (std::is_same_v<T, long>) x = std::get<long>(a);
                  if constexpr (std::is_same_v<T,double>) x = std::get<double>(a);
		}, a );
    return x;
  }
  ComputeTypes operator+( ComputeTypes a, ComputeTypes b )
  {
    Complex x{getValue(a)}, y{getValue(b)};
    return ComputeTypes{Complex{x+y}};
  }
  ComputeTypes& operator+=( ComputeTypes& lhs, ComputeTypes rhs )
  {
    lhs = lhs + rhs;
    return lhs;
  }
  std::ostream& operator<<( std::ostream& os, ComputeTypes ct )
  {
    std::visit( [&os](auto&& arg) { os << arg; }, ct );
    return os;
  }
}

using namespace Computer;

// In this function we have the barebones of an interpreter lifting numeric types
// to Complex and performing computation

static Types SimpleTest()
{
  constexpr long VALUE{5};
  constexpr ComputeTypes x = long{VALUE};
  static_assert( std::get<long>(x) == VALUE, "Elementary check" );
  static_assert( std::get<0>(x) == VALUE, "Elementary check" );  
  auto retval = Types(x);
  retval = retval;
  std::cout << "sizeof(ComputeTypes) = " << sizeof(ComputeTypes) << std::endl;
  std::cout << "sizeof(Types) = " << sizeof(Types) << std::endl;
  ComputeTypes y{Complex{1.+2i}};
  ComputeTypes z{Complex{2.+1i}};
  y = z+y;
  y += ComputeTypes{long{2}};
  z = y; // at this time z has type Complex
  std::cout << "y = " << y << std::endl;
  if( std::get<Complex>(z) == Complex{0} )
    return std::nullopt;
  else
    return z;
}

auto TupleTest(void)  //-> std::tuple<int,std::string,double>
{
  return std::make_tuple(1,"one",1.0d);
}

int main()
{
  auto x = SimpleTest();
  if( x.has_value() )
    std::cout << "We got " << std::get<Complex>(x.value_or(ZEROC)) << std::endl;
  else
    std::cout << "Computation resulted in Complex{0}" << std::endl;
  
  // This is an example of structured binding
  auto [ix,sx,dx] = TupleTest(); 
  std::cout << ix << "\t" << sx << "\t" << dx << std::endl;
  return 0;
}

#if 0
sizeof(ComputeTypes) = 24
sizeof(Types) = 32
y = (5,3)
We got (5,3)
1	one	1
#endif
