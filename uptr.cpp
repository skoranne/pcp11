////////////////////////////////////////////////////////////////////////////////
// File   : uptr.cpp
// Author : Sandeep Koranne, (C) 2020. All rights reserved.
// Purpose: Example program for checking C++-11 stuff.
// Revision control on github
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cassert>
#include <memory> // for unique_ptr
#include <cstdlib>
#include <future>
#include <mutex>
#include <thread>
#include <typeinfo> // for typeinfo(A).name()
#include <type_traits> // for POD stuff

namespace Example {
  enum class Operations : int { NONE=0, ADD, SUB, MUL, DIV, MOD, EXP, LOG, LSH, RSH };
int OperationsToInt( Operations );
Operations IntToOperations( int );
class Foo {
private:
  int x;
public:
  explicit Foo(int _x): x{_x} {}
  int getX() const { return x; }
  void IncrementNumber() { ++x; }
  friend std::ostream& operator<<( std::ostream&, const Foo& );
  ~Foo() { std::cout << __PRETTY_FUNCTION__ << " " << x << std::endl; }
  Foo( const Foo& ) = delete;
  Foo& operator=( const Foo& ) = delete;
};

  std::ostream& operator<<( std::ostream&, const Foo& );
} // end of namespace Example

std::ostream& Example::operator<<( std::ostream& os, const Foo& foo ) {
  return os << foo.x << " ";
}

////////////////////////////////////////////////////////////////////////////////
// Conversion operations to int and vice-versa
//
////////////////////////////////////////////////////////////////////////////////
static void PrettyPrint( const char* name, int value )
{
  std::cout << name << "\t" << value << std::endl;
}

[[nodiscard]] int Example::OperationsToInt( Example::Operations op )
{
  switch( op ) {
  case Operations::NONE: return 0;
  case Operations::ADD : return 1;
  case Operations::SUB : return 2;
  case Operations::MUL : return 3;
  case Operations::DIV : return 4;
  case Operations::MOD : return 5;
  case Operations::EXP : return 6;
  case Operations::LOG : return 7;
  case Operations::LSH : return 8;
  case Operations::RSH : return 9;
  default  : return (-1);
  }
  return (-1);
}

[[nodiscard]] Example::Operations Example::IntToOperations(int op)
{
  switch(op) {
  case 0: return Operations::NONE;
  case 1: return Operations::ADD;
  case 2: return Operations::SUB;
  case 3: return Operations::MUL;
  case 4: return Operations::DIV;
  case 5: return Operations::MOD;
  case 6: return Operations::EXP;
  case 7: return Operations::LOG;
  case 8: return Operations::LSH;
  case 9: return Operations::RSH;
  default: return Operations::NONE;
  }
}

#define PRINTER(name) PrettyPrint(#name, OperationsToInt(name))
static void DescribeOperations()
{
  PRINTER( Example::Operations::NONE );
  PRINTER( Example::Operations::ADD );
  PRINTER( Example::Operations::SUB );
  PRINTER( Example::Operations::MUL );
  PRINTER( Example::Operations::DIV );
  PRINTER( Example::Operations::MOD );
  PRINTER( Example::Operations::EXP );
  PRINTER( Example::Operations::LOG );
  PRINTER( Example::Operations::LSH );
  PRINTER( Example::Operations::RSH );
}
#undef PRINTER

void Test3()
{
  DescribeOperations();
}


std::mutex COUT_MUTEX;
using LCK_GRD = std::lock_guard< std::mutex >;
using namespace Example;
using FUP = std::unique_ptr< Foo >;
[[nodiscard]] FUP OperateOnFUP( FUP x )
{
  {
    LCK_GRD lck{ COUT_MUTEX };
    std::cout << __PRETTY_FUNCTION__ << " " << x->getX() 
	      << " running on " << std::this_thread::get_id() << std::endl;
  }
  x->IncrementNumber();
  return std::move( x );
}

void OperateOnFoo( const Foo& x )
{
  int y = x.getX();
  std::cout << " y = " << y << std::endl;
}

int GLOBAL_CLOSURE_VAR = 0;
auto ClosureExample( int& x )
{
  return [=](int y) { return (x+y+GLOBAL_CLOSURE_VAR);};
}


void Test2()
{
  GLOBAL_CLOSURE_VAR = 10;
  auto C10 = ClosureExample( GLOBAL_CLOSURE_VAR ); // a function which adds 10 to its arg
  GLOBAL_CLOSURE_VAR = 20;
  auto C20 = ClosureExample( GLOBAL_CLOSURE_VAR ); // a function which adds 20 to its arg
  GLOBAL_CLOSURE_VAR = 3;
  std::cout << "Ans10 = " << C10(0) << std::endl; // now it gives 13
  std::cout << "Ans20 = " << C20(0) << std::endl;
  // C++ also does lexical closure, so answer is still based on old value
}


void Test1()
{
  Foo f1(10);
  OperateOnFoo( f1 );
  int y{};
  auto fnl = [=,&f1](const Foo& x)->int{ return ( y+f1.getX() + x.getX() ); };
  y = fnl( f1 );
  std::cout << "As a result of lambda = " << y << std::endl;
  static_assert( std::is_copy_constructible<Foo>::value == false, "No COPY CTOR" );
  static_assert( std::is_copy_assignable<Foo>::value == false, "No operator=" );
  static_assert( std::is_trivial<Foo>::value == false, "Trivial" );
  assert( y == 20 );
  {
    LCK_GRD lck{ COUT_MUTEX };
    std::cout << __PRETTY_FUNCTION__ << " " << f1 
	      << " running on " << std::this_thread::get_id() << std::endl;
  }
  FUP a{ new Foo{5} };
  FUP b{ new Foo{2} };
  y = a->getX() + b->getX();
  std::cout << " y = " << y << " = " << *a << " + " << *b << std::endl;
  auto f0 = async( OperateOnFUP, std::move(a) );
  {
    LCK_GRD lck{ COUT_MUTEX };
    std::cout << "typeof(f0) = " << typeid( f0 ).name() << std::endl;
  }
  using PT = std::packaged_task<FUP(FUP)>; // declare task type
  PT task{ OperateOnFUP };                 // create task
  std::future<FUP> result = task.get_future(); // create place holder for result
  //task( std::move( b ) ); // this is launched on main thread
  std::thread tid0(std::move( task ), std::move(b) ); // now thread is launched with b
  a = f0.get();
  tid0.join(); // without this, the program terminates
  b = result.get();
  b.reset();
  std::cout << "Done with call\n";
}

int main()
{
  Test3();
  Test2();
  Test1();
  return( EXIT_SUCCESS );
}


