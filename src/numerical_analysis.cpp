#include "numerical_analysis.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace numerical {

double newton_raphson(
    std::function<double(double)> f,
    std::function<double(double)> df,
    double x0,
    double tolerance,
    int max_iter
) {
    double x = x0;
    
    for (int i = 0; i < max_iter; i++) {
        double fx = f(x);
        double dfx = df(x);
        
        if (std::abs(dfx) < 1e-10) {
            throw std::runtime_error("Derivada muy pequeña en Newton-Raphson");
        }
        
        double x_new = x - fx / dfx;
        
        if (std::abs(x_new - x) < tolerance) {
            std::cout << "Newton-Raphson convergió en " << i+1 << " iteraciones\n";
            return x_new;
        }
        
        x = x_new;
    }
    
    std::cerr << "Newton-Raphson: No convergió en " << max_iter << " iteraciones\n";
    return x;
}

double bisection_method(
    std::function<double(double)> f,
    double a,
    double b,
    double tolerance
) {
    if (f(a) * f(b) > 0) {
        throw std::runtime_error("f(a) y f(b) deben tener signos opuestos");
    }
    
    int iter = 0;
    while (std::abs(b - a) > tolerance) {
        double c = (a + b) / 2.0;
        
        if (f(c) == 0.0) {
            return c;
        }
        
        if (f(a) * f(c) < 0) {
            b = c;
        } else {
            a = c;
        }
        
        iter++;
    }
    
    std::cout << "Bisección convergió en " << iter << " iteraciones\n";
    return (a + b) / 2.0;
}

double numerical_derivative(
    std::function<double(double)> f,
    double x,
    double h
) {
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

} // namespace numerical
