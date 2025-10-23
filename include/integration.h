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
