#include "ode_solvers.h"
#include "integration.h"
#include "numerical_analysis.h"
#include <iostream>
#include <cmath>
#include <cassert>

void test_euler_method() {
    std::cout << "Testing Euler Method... ";
    
    // dy/dt = -2y, y(0) = 1
    // Solución exacta: y(t) = e^(-2t)
    auto f = [](double t, double y) { return -2.0 * y; };
    
    auto result = numerical::euler_method(f, 1.0, 0.0, 1.0, 1000);
    double exact = std::exp(-2.0);
    double error = std::abs(result.back() - exact);
    
    assert(error < 0.01);  // Error tolerable para Euler
    std::cout << "✓ PASSED (error = " << error << ")\n";
}

void test_runge_kutta_4() {
    std::cout << "Testing Runge-Kutta 4... ";
    
    auto f = [](double t, double y) { return -2.0 * y; };
    
    auto result = numerical::runge_kutta_4(f, 1.0, 0.0, 1.0, 100);
    double exact = std::exp(-2.0);
    double error = std::abs(result.back() - exact);
    
    assert(error < 1e-6);  // RK4 es mucho más preciso
    std::cout << "✓ PASSED (error = " << error << ")\n";
}

void test_trapezoidal_rule() {
    std::cout << "Testing Trapezoidal Rule... ";
    
    // Integrar x² de 0 a 1 (resultado exacto = 1/3)
    auto f = [](double x) { return x * x; };
    
    double result = numerical::trapezoidal_rule(f, 0.0, 1.0, 1000);
    double exact = 1.0 / 3.0;
    double error = std::abs(result - exact);
    
    assert(error < 1e-5);
    std::cout << "✓ PASSED (error = " << error << ")\n";
}

void test_simpson_rule() {
    std::cout << "Testing Simpson Rule... ";
    
    auto f = [](double x) { return x * x; };
    
    double result = numerical::simpson_rule(f, 0.0, 1.0, 1000);
    double exact = 1.0 / 3.0;
    double error = std::abs(result - exact);
    
    assert(error < 1e-10);
    std::cout << "✓ PASSED (error = " << error << ")\n";
}

void test_newton_raphson() {
    std::cout << "Testing Newton-Raphson... ";
    
    // Encontrar √2 (raíz de x² - 2 = 0)
    auto f = [](double x) { return x*x - 2.0; };
    auto df = [](double x) { return 2.0 * x; };
    
    double result = numerical::newton_raphson(f, df, 1.0, 1e-10);
    double exact = std::sqrt(2.0);
    double error = std::abs(result - exact);
    
    assert(error < 1e-10);
    std::cout << "✓ PASSED (error = " << error << ")\n";
}

void test_bisection_method() {
    std::cout << "Testing Bisection Method... ";
    
    // Encontrar raíz de x³ - x - 2 = 0 en [1, 2]
    auto f = [](double x) { return x*x*x - x - 2.0; };
    
    double result = numerical::bisection_method(f, 1.0, 2.0, 1e-6);
    double exact = 1.5213797068;  // Raíz conocida
    double error = std::abs(result - exact);
    
    assert(error < 1e-5);
    std::cout << "✓ PASSED (error = " << error << ")\n";
}

void test_gauss_quadrature() {
    std::cout << "Testing Gauss Quadrature... ";
    
    // Integrar sin(x) de 0 a π (resultado = 2)
    auto f = [](double x) { return std::sin(x); };
    
    double result = numerical::gauss_quadrature(f, 0.0, std::acos(-1.0));
    double exact = 2.0;
    double error = std::abs(result - exact);
    
    assert(error < 1e-4);
    std::cout << "✓ PASSED (error = " << error << ")\n";
}

void test_numerical_derivative() {
    std::cout << "Testing Numerical Derivative... ";
    
    // Derivada de x² en x=3 (resultado = 6)
    auto f = [](double x) { return x * x; };
    
    double result = numerical::numerical_derivative(f, 3.0);
    double exact = 6.0;
    double error = std::abs(result - exact);
    
    assert(error < 1e-5);
    std::cout << "✓ PASSED (error = " << error << ")\n";
}

int main() {
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════╗\n";
    std::cout << "║  Physics Numerical Methods - Unit Tests ║\n";
    std::cout << "╚══════════════════════════════════════════╝\n";
    std::cout << "\n";
    
    try {
        test_euler_method();
        test_runge_kutta_4();
        test_trapezoidal_rule();
        test_simpson_rule();
        test_newton_raphson();
        test_bisection_method();
        test_gauss_quadrature();
        test_numerical_derivative();
        
        std::cout << "\n";
        std::cout << "════════════════════════════════════════════\n";
        std::cout << "  ✓ ALL TESTS PASSED (" << 8 << "/8)\n";
        std::cout << "════════════════════════════════════════════\n";
        std::cout << "\n";
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\n✗ TEST FAILED: " << e.what() << "\n";
        return 1;
    }
}
