#include <cassert>
#include <vector>
#include <iostream>
#include <valarray>
#include <algorithm>
#include <functional>

namespace BoundaryIntegration {
  struct Point {
    double x,y;
  };
  struct Line {
    Point u,v; // u->v line, as a vector
  };
  using PolyLine = std::vector<Line>;
  template <typename FUNC, template <typename FUNC2> class ALGORITHM>  
  double CalculateLineIntegral(ALGORITHM<FUNC> f, Line L);
  template <typename FUNC, template <typename FUNC2> class ALGORITHM>
  double CalculateLineIntegral(ALGORITHM<FUNC> f, const PolyLine& PL);

  template <typename FUNC>
  struct MIDPOINT
  {
    MIDPOINT( FUNC f ): m_f( f ) {}
    FUNC m_f;
    double CalculateLineIntegral( Line L );
  };
};

template<typename FUNC>
double BoundaryIntegration::MIDPOINT<FUNC>::CalculateLineIntegral( Line L )
{
  Point mid_point{0.5*(L.u.x+L.v.x), 0.5*(L.u.y+L.v.y) };
  std::cout << "Evaluating at (" << mid_point.x << "," << mid_point.y << ")\n";
  return m_f( mid_point );
}

template <typename FUNC, template <typename FUNC2> class ALGORITHM>
double BoundaryIntegration::CalculateLineIntegral(ALGORITHM<FUNC> f, const PolyLine& PL)
{
  double retval = 0.0;
  for( auto L : PL ) {
    retval += CalculateLineIntegral( f, L );
  }
  return retval;
}

template <typename FUNC, template <typename FUNC2> class ALGORITHM>
double BoundaryIntegration::CalculateLineIntegral(ALGORITHM<FUNC> f, const Line L)
{
  return f.CalculateLineIntegral( L );
}

using namespace BoundaryIntegration;

static void TestSimpleLine()
{
  using FUNC = std::function< double(Point) >;
  using MIDPOINT = BoundaryIntegration::MIDPOINT<FUNC>;  
 
  FUNC f = [](Point p){ return p.x; };
  MIDPOINT mp(f);
  PolyLine PL{ {0,0,1,0},{1,0,1,1}, {1,1,0,1}, {0,1,0,0} };
  double ans = CalculateLineIntegral( mp, PL );
  std::cout << "Ans = " << ans << std::endl;
}

// Generate C++ function for factorial using unsigned long

unsigned long Factorial(unsigned int n) {
    if (n == 0) return 1;
    else return n * Factorial(n-1);
}

int main()
{
  TestSimpleLine();
  // TODO: make a call to factorial to test it and std::cout the result
  std::cout << "Factorial(5) = " << Factorial(5) << std::endl;

  return 0;
}

