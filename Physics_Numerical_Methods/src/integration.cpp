#include "integration.h"
#include <vector>

namespace numerical {

double trapezoidal_rule(
    std::function<double(double)> f,
    double a,
    double b,
    int n
) {
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));
    
    for (int i = 1; i < n; i++) {
        sum += f(a + i * h);
    }
    
    return h * sum;
}

double simpson_rule(
    std::function<double(double)> f,
    double a,
    double b,
    int n
) {
    if (n % 2 != 0) n++;
    
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    
    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        sum += (i % 2 == 0) ? 2 * f(x) : 4 * f(x);
    }
    
    return (h / 3.0) * sum;
}

double gauss_quadrature(
    std::function<double(double)> f,
    double a,
    double b
) {
    // Puntos y pesos de Gauss-Legendre (n=5)
    std::vector<double> points = {
        -0.9061798459, -0.5384693101, 0.0, 
         0.5384693101,  0.9061798459
    };
    std::vector<double> weights = {
        0.2369268851, 0.4786286705, 0.5688888889,
        0.4786286705, 0.2369268851
    };
    
    double sum = 0.0;
    double c = (b - a) / 2.0;
    double d = (b + a) / 2.0;
    
    for (size_t i = 0; i < points.size(); i++) {
        sum += weights[i] * f(c * points[i] + d);
    }
    
    return c * sum;
}

} // namespace numerical
