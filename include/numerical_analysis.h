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
