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
