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
