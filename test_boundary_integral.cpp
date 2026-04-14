// test_boundary_integral.cpp
// Unit test for BoundaryIntegration::CalculateLineIntegral using the MIDPOINT algorithm.
// The test constructs a unit square polyline and integrates the function f(x,y)=x.
// Expected result (given the current midpoint implementation which returns only f(midpoint))
// is the sum of x‑coordinates of the four edge midpoints: 0.5 + 1 + 0.5 + 0 = 2.0.

#include "boundary_integral.cpp"
#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    using FUNC = std::function<double(Point)>;
    using MIDPOINT = BoundaryIntegration::MIDPOINT<FUNC>;

    // Scalar field: f(x,y) = x
    FUNC f = [](Point p){ return p.x; };
    MIDPOINT mp(f);

    // Define a closed unit square as a polyline
    BoundaryIntegration::PolyLine PL{
        {{0,0},{1,0}},   // bottom edge
        {{1,0},{1,1}},   // right edge
        {{1,1},{0,1}},   // top edge
        {{0,1},{0,0}}    // left edge
    };

    double result = BoundaryIntegration::CalculateLineIntegral(mp, PL);
    double expected = 2.0; // sum of x‑coordinates of midpoints
    const double eps = 1e-9;
    assert(std::abs(result - expected) < eps && "CalculateLineIntegral returned unexpected value");

    std::cout << "CalculateLineIntegral unit test passed. Result = " << result << std::endl;
    return 0;
}
