#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>

class Polynomial {
public:
    std::vector<double> coefficients;

    double evaluate(double x) const;
};

#endif  // POLYNOMIAL_H
