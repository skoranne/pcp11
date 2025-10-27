class Polynomial {
public:
    std::vector<double> coefficients;

    double operator()(double x) const {
        return std::accumulate(coefficients.begin(), coefficients.end(), 0.0,
            [x](double accumulator, double coefficient) {
                return accumulator * x + coefficient;
            });
    }
};
