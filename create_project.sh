#!/bin/bash

# Script para crear la estructura completa del proyecto
# Physics_Numerical_Methods

echo "🚀 Creando estructura del proyecto Physics_Numerical_Methods..."

# Crear directorios
mkdir -p Physics_Numerical_Methods/{include,src,examples,tests,docs,build}

cd Physics_Numerical_Methods

# ==================== README.md ====================
cat > README.md << 'EOF'
# 🔬 Physics Numerical Methods

[![C++17](https://img.shields.io/badge/C++-17-blue.svg)](https://isocpp.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Implementaciones en C++ de métodos numéricos para resolver problemas físicos.

## 📦 Instalación Rápida

```bash
git clone https://github.com/juareymal-ui/Physics_Numerical_Methods.git
cd Physics_Numerical_Methods
mkdir build && cd build
cmake ..
make
```

## 🚀 Uso Rápido

```bash
# Ejecutar simulación de péndulo
./bin/pendulum_simulation

# Ejecutar oscilador armónico
./bin/harmonic_oscillator
```

## 📚 Documentación

Ver [docs/theory.md](docs/theory.md) para teoría detallada.

## 👤 Autor

**Juan Maldonado** - [GitHub](https://github.com/juareymal-ui)
EOF

# ==================== CMakeLists.txt ====================
cat > CMakeLists.txt << 'EOF'
cmake_minimum_required(VERSION 3.15)
project(PhysicsNumericalMethods VERSION 1.0 LANGUAGES CXX)

# C++17 estándar
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Opciones de compilación
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra")

# Directorios de include
include_directories(${PROJECT_SOURCE_DIR}/include)

# Librería principal
add_library(numerical_methods
    src/ode_solvers.cpp
    src/integration.cpp
    src/physics_simulations.cpp
    src/numerical_analysis.cpp
)

# Ejemplos
add_executable(pendulum_simulation examples/pendulum_simulation.cpp)
target_link_libraries(pendulum_simulation numerical_methods)

add_executable(harmonic_oscillator examples/harmonic_oscillator.cpp)
target_link_libraries(harmonic_oscillator numerical_methods)

add_executable(planetary_motion examples/planetary_motion.cpp)
target_link_libraries(planetary_motion numerical_methods)

add_executable(heat_equation examples/heat_equation.cpp)
target_link_libraries(heat_equation numerical_methods)

# Tests
add_executable(unit_tests tests/unit_tests.cpp)
target_link_libraries(unit_tests numerical_methods)

# Directorios de salida
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
EOF

# ==================== include/ode_solvers.h ====================
cat > include/ode_solvers.h << 'EOF'
#ifndef ODE_SOLVERS_H
#define ODE_SOLVERS_H

#include <vector>
#include <functional>

namespace numerical {

/**
 * @brief Método de Euler para resolver EDOs
 * @param f Función dy/dt = f(t, y)
 * @param y0 Condición inicial
 * @param t0 Tiempo inicial
 * @param tf Tiempo final
 * @param n_steps Número de pasos
 * @return Vector con la solución
 */
std::vector<double> euler_method(
    std::function<double(double, double)> f,
    double y0,
    double t0,
    double tf,
    int n_steps
);

/**
 * @brief Método de Runge-Kutta de 4to orden
 * @param f Función dy/dt = f(t, y)
 * @param y0 Condición inicial
 * @param t0 Tiempo inicial
 * @param tf Tiempo final
 * @param n_steps Número de pasos
 * @return Vector con la solución
 */
std::vector<double> runge_kutta_4(
    std::function<double(double, double)> f,
    double y0,
    double t0,
    double tf,
    int n_steps
);

/**
 * @brief Estructura para sistemas de EDOs
 */
struct StateVector {
    std::vector<double> values;
    
    StateVector(size_t size = 0) : values(size, 0.0) {}
    StateVector(std::initializer_list<double> init) : values(init) {}
    
    double& operator[](size_t i) { return values[i]; }
    const double& operator[](size_t i) const { return values[i]; }
    size_t size() const { return values.size(); }
};

/**
 * @brief RK4 para sistemas de EDOs
 */
std::vector<StateVector> runge_kutta_4_system(
    std::function<StateVector(double, const StateVector&)> f,
    const StateVector& y0,
    double t0,
    double tf,
    int n_steps
);

} // namespace numerical

#endif // ODE_SOLVERS_H
EOF

# ==================== include/integration.h ====================
cat > include/integration.h << 'EOF'
#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <functional>

namespace numerical {

/**
 * @brief Regla del trapecio para integración numérica
 * @param f Función a integrar
 * @param a Límite inferior
 * @param b Límite superior
 * @param n Número de subdivisiones
 * @return Aproximación de la integral
 */
double trapezoidal_rule(
    std::function<double(double)> f,
    double a,
    double b,
    int n
);

/**
 * @brief Regla de Simpson para integración numérica
 * @param f Función a integrar
 * @param a Límite inferior
 * @param b Límite superior
 * @param n Número de subdivisiones (debe ser par)
 * @return Aproximación de la integral
 */
double simpson_rule(
    std::function<double(double)> f,
    double a,
    double b,
    int n
);

/**
 * @brief Cuadratura de Gauss (5 puntos)
 * @param f Función a integrar
 * @param a Límite inferior
 * @param b Límite superior
 * @return Aproximación de la integral
 */
double gauss_quadrature(
    std::function<double(double)> f,
    double a,
    double b
);

} // namespace numerical

#endif // INTEGRATION_H
EOF

# ==================== include/physics_simulations.h ====================
cat > include/physics_simulations.h << 'EOF'
#ifndef PHYSICS_SIMULATIONS_H
#define PHYSICS_SIMULATIONS_H

#include <string>
#include <utility>

namespace physics {

/**
 * @brief Simulador de péndulo simple
 */
class SimplePendulum {
private:
    double g_;  // Gravedad
    double L_;  // Longitud
    
public:
    SimplePendulum(double gravity = 9.81, double length = 1.0);
    
    /**
     * @brief Ecuación del movimiento del péndulo
     * @return (d_theta, d_omega)
     */
    std::pair<double, double> equation(double theta, double omega) const;
    
    /**
     * @brief Ejecutar simulación
     * @param theta0 Ángulo inicial (rad)
     * @param omega0 Velocidad angular inicial (rad/s)
     * @param t_max Tiempo máximo (s)
     * @param steps Número de pasos
     * @param filename Archivo de salida
     */
    void simulate(double theta0, double omega0, double t_max, 
                  int steps, const std::string& filename) const;
};

/**
 * @brief Simulador de oscilador armónico amortiguado
 */
class DampedOscillator {
private:
    double omega_;  // Frecuencia natural
    double gamma_;  // Coeficiente de amortiguamiento
    
public:
    DampedOscillator(double omega, double gamma);
    
    std::pair<double, double> equation(double x, double v) const;
    
    void simulate(double x0, double v0, double t_max,
                  int steps, const std::string& filename) const;
};

/**
 * @brief Ecuación del calor 1D
 */
class HeatEquation1D {
private:
    double alpha_;  // Difusividad térmica
    double L_;      // Longitud del dominio
    int nx_;        // Puntos espaciales
    
public:
    HeatEquation1D(double alpha, double length, int n_points);
    
    void solve(double t_max, int nt, const std::string& filename);
};

} // namespace physics

#endif // PHYSICS_SIMULATIONS_H
EOF

# ==================== include/numerical_analysis.h ====================
cat > include/numerical_analysis.h << 'EOF'
#ifndef NUMERICAL_ANALYSIS_H
#define NUMERICAL_ANALYSIS_H

#include <functional>

namespace numerical {

/**
 * @brief Método de Newton-Raphson para encontrar raíces
 * @param f Función f(x)
 * @param df Derivada f'(x)
 * @param x0 Estimación inicial
 * @param tolerance Tolerancia de convergencia
 * @param max_iter Máximo número de iteraciones
 * @return Raíz aproximada
 */
double newton_raphson(
    std::function<double(double)> f,
    std::function<double(double)> df,
    double x0,
    double tolerance = 1e-6,
    int max_iter = 100
);

/**
 * @brief Método de bisección para encontrar raíces
 * @param f Función f(x)
 * @param a Límite inferior del intervalo
 * @param b Límite superior del intervalo
 * @param tolerance Tolerancia
 * @return Raíz aproximada
 */
double bisection_method(
    std::function<double(double)> f,
    double a,
    double b,
    double tolerance = 1e-6
);

/**
 * @brief Calcular derivada numérica (diferencias centrales)
 * @param f Función
 * @param x Punto de evaluación
 * @param h Tamaño del paso
 * @return Aproximación de f'(x)
 */
double numerical_derivative(
    std::function<double(double)> f,
    double x,
    double h = 1e-5
);

} // namespace numerical

#endif // NUMERICAL_ANALYSIS_H
EOF

echo "✅ Headers creados"

# ==================== src/ode_solvers.cpp ====================
cat > src/ode_solvers.cpp << 'EOF'
#include "ode_solvers.h"

namespace numerical {

std::vector<double> euler_method(
    std::function<double(double, double)> f,
    double y0,
    double t0,
    double tf,
    int n_steps
) {
    double h = (tf - t0) / n_steps;
    std::vector<double> y(n_steps + 1);
    
    y[0] = y0;
    double t = t0;
    
    for (int i = 0; i < n_steps; i++) {
        y[i + 1] = y[i] + h * f(t, y[i]);
        t += h;
    }
    
    return y;
}

std::vector<double> runge_kutta_4(
    std::function<double(double, double)> f,
    double y0,
    double t0,
    double tf,
    int n_steps
) {
    double h = (tf - t0) / n_steps;
    std::vector<double> y(n_steps + 1);
    
    y[0] = y0;
    double t = t0;
    
    for (int i = 0; i < n_steps; i++) {
        double k1 = h * f(t, y[i]);
        double k2 = h * f(t + h/2, y[i] + k1/2);
        double k3 = h * f(t + h/2, y[i] + k2/2);
        double k4 = h * f(t + h, y[i] + k3);
        
        y[i + 1] = y[i] + (k1 + 2*k2 + 2*k3 + k4) / 6.0;
        t += h;
    }
    
    return y;
}

std::vector<StateVector> runge_kutta_4_system(
    std::function<StateVector(double, const StateVector&)> f,
    const StateVector& y0,
    double t0,
    double tf,
    int n_steps
) {
    double h = (tf - t0) / n_steps;
    std::vector<StateVector> result(n_steps + 1, StateVector(y0.size()));
    
    result[0] = y0;
    double t = t0;
    
    for (int i = 0; i < n_steps; i++) {
        StateVector k1 = f(t, result[i]);
        
        StateVector temp1(y0.size());
        for (size_t j = 0; j < y0.size(); j++) {
            temp1[j] = result[i][j] + h * k1[j] / 2.0;
        }
        StateVector k2 = f(t + h/2, temp1);
        
        StateVector temp2(y0.size());
        for (size_t j = 0; j < y0.size(); j++) {
            temp2[j] = result[i][j] + h * k2[j] / 2.0;
        }
        StateVector k3 = f(t + h/2, temp2);
        
        StateVector temp3(y0.size());
        for (size_t j = 0; j < y0.size(); j++) {
            temp3[j] = result[i][j] + h * k3[j];
        }
        StateVector k4 = f(t + h, temp3);
        
        for (size_t j = 0; j < y0.size(); j++) {
            result[i + 1][j] = result[i][j] + 
                (h / 6.0) * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
        }
        
        t += h;
    }
    
    return result;
}

} // namespace numerical
EOF

# ==================== src/integration.cpp ====================
cat > src/integration.cpp << 'EOF'
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
EOF

# ==================== src/physics_simulations.cpp ====================
cat > src/physics_simulations.cpp << 'EOF'
#include "physics_simulations.h"
#include "ode_solvers.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>

namespace physics {

// ========== SimplePendulum ==========
SimplePendulum::SimplePendulum(double gravity, double length)
    : g_(gravity), L_(length) {}

std::pair<double, double> SimplePendulum::equation(double theta, double omega) const {
    double d_theta = omega;
    double d_omega = -(g_ / L_) * std::sin(theta);
    return {d_theta, d_omega};
}

void SimplePendulum::simulate(double theta0, double omega0, double t_max,
                               int steps, const std::string& filename) const {
    double dt = t_max / steps;
    double theta = theta0;
    double omega = omega0;
    
    std::ofstream file(filename);
    file << "time,theta,omega,energy\n";
    
    for (int i = 0; i <= steps; i++) {
        double t = i * dt;
        double energy = 0.5 * L_ * L_ * omega * omega + 
                       g_ * L_ * (1 - std::cos(theta));
        
        file << t << "," << theta << "," << omega << "," << energy << "\n";
        
        // RK4
        auto [k1_theta, k1_omega] = equation(theta, omega);
        auto [k2_theta, k2_omega] = equation(theta + 0.5*dt*k1_theta, 
                                             omega + 0.5*dt*k1_omega);
        auto [k3_theta, k3_omega] = equation(theta + 0.5*dt*k2_theta, 
                                             omega + 0.5*dt*k2_omega);
        auto [k4_theta, k4_omega] = equation(theta + dt*k3_theta, 
                                             omega + dt*k3_omega);
        
        theta += (dt/6.0) * (k1_theta + 2*k2_theta + 2*k3_theta + k4_theta);
        omega += (dt/6.0) * (k1_omega + 2*k2_omega + 2*k3_omega + k4_omega);
    }
    
    file.close();
    std::cout << "Simulación completada: " << filename << "\n";
}

// ========== DampedOscillator ==========
DampedOscillator::DampedOscillator(double omega, double gamma)
    : omega_(omega), gamma_(gamma) {}

std::pair<double, double> DampedOscillator::equation(double x, double v) const {
    return {v, -2*gamma_*v - omega_*omega_*x};
}

void DampedOscillator::simulate(double x0, double v0, double t_max,
                                 int steps, const std::string& filename) const {
    double dt = t_max / steps;
    double x = x0, v = v0;
    
    std::ofstream file(filename);
    file << "time,position,velocity,energy\n";
    
    for (int i = 0; i <= steps; i++) {
        double t = i * dt;
        double energy = 0.5 * v*v + 0.5 * omega_*omega_ * x*x;
        
        file << t << "," << x << "," << v << "," << energy << "\n";
        
        // RK4
        auto [k1x, k1v] = equation(x, v);
        auto [k2x, k2v] = equation(x + 0.5*dt*k1x, v + 0.5*dt*k1v);
        auto [k3x, k3v] = equation(x + 0.5*dt*k2x, v + 0.5*dt*k2v);
        auto [k4x, k4v] = equation(x + dt*k3x, v + dt*k3v);
        
        x += (dt/6.0) * (k1x + 2*k2x + 2*k3x + k4x);
        v += (dt/6.0) * (k1v + 2*k2v + 2*k3v + k4v);
    }
    
    file.close();
    std::cout << "Simulación completada: " << filename << "\n";
}

// ========== HeatEquation1D ==========
HeatEquation1D::HeatEquation1D(double alpha, double length, int n_points)
    : alpha_(alpha), L_(length), nx_(n_points) {}

void HeatEquation1D::solve(double t_max, int nt, const std::string& filename) {
    double dx = L_ / (nx_ - 1);
    double dt = t_max / nt;
    
    std::vector<double> u(nx_), u_new(nx_);
    
    // Condición inicial: pulso gaussiano
    for (int i = 0; i < nx_; i++) {
        double x = i * dx;
        u[i] = std::exp(-50 * std::pow(x - L_/2, 2));
    }
    
    std::ofstream file(filename);
    file << "x";
    for (int n = 0; n <= nt; n += (nt/10)) {
        file << ",t" << n;
    }
    file << "\n";
    
    for (int i = 0; i < nx_; i++) {
        file << i * dx << "," << u[i];
    }
    
    for (int n = 0; n < nt; n++) {
        for (int i = 1; i < nx_ - 1; i++) {
            u_new[i] = u[i] + alpha_ * dt / (dx * dx) * 
                      (u[i+1] - 2*u[i] + u[i-1]);
        }
        
        u_new[0] = 0;
        u_new[nx_-1] = 0;
        u = u_new;
        
        if ((n+1) % (nt/10) == 0) {
            for (int i = 0; i < nx_; i++) {
                file << "," << u[i];
            }
        }
    }
    
    file << "\n";
    file.close();
    std::cout << "Ecuación del calor resuelta: " << filename << "\n";
}

} // namespace physics
EOF

# ==================== src/numerical_analysis.cpp ====================
cat > src/numerical_analysis.cpp << 'EOF'
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
EOF

echo "✅ Source files creados"
