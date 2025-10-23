#define _USE_MATH_DEFINES
#include "physics_simulations.h"
#include <iostream>
#include <cmath>

int main() {
    std::cout << "=== Simulación de Péndulo Simple ===\n\n";
    
    // Crear péndulo: g=9.81 m/s², L=1.0 m
    physics::SimplePendulum pendulum(9.81, 1.0);
    
    // Condiciones iniciales
    double theta0 = M_PI / 4.0;  // 45 grados
    double omega0 = 0.0;          // Velocidad inicial = 0
    
    std::cout << "Parámetros:\n";
    std::cout << "  g = 9.81 m/s²\n";
    std::cout << "  L = 1.0 m\n";
    std::cout << "  θ₀ = " << theta0 * 180 / M_PI << "°\n";
    std::cout << "  ω₀ = " << omega0 << " rad/s\n\n";
    
    // Simular 10 segundos
    pendulum.simulate(theta0, omega0, 10.0, 1000, "pendulum_data.csv");
    
    std::cout << "\n✓ Datos guardados en pendulum_data.csv\n";
    std::cout << "\nPara visualizar, ejecuta:\n";
    std::cout << "  python3 visualize_pendulum.py\n";
    
    return 0;
}
