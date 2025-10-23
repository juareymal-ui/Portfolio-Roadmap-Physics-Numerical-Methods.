#define _USE_MATH_DEFINES
#include "physics_simulations.h"
#include <iostream>
#include <cmath>

int main() {
    std::cout << "=== Oscilador Armónico Amortiguado ===\n\n";
    
    // Parámetros
    double omega = 2.0 * M_PI;  // ω = 2π rad/s (periodo de 1 segundo)
    double gamma = 0.1;          // γ = 0.1 s⁻¹ (amortiguamiento débil)
    
    physics::DampedOscillator oscillator(omega, gamma);
    
    // Condiciones iniciales
    double x0 = 1.0;   // Posición inicial: 1 metro
    double v0 = 0.0;   // Velocidad inicial: 0
    
    std::cout << "Parámetros:\n";
    std::cout << "  ω = " << omega << " rad/s\n";
    std::cout << "  γ = " << gamma << " s⁻¹\n";
    std::cout << "  x₀ = " << x0 << " m\n";
    std::cout << "  v₀ = " << v0 << " m/s\n\n";
    
    // Tipo de amortiguamiento
    double discriminant = gamma * gamma - omega * omega;
    if (discriminant < 0) {
        std::cout << "Tipo: Subamortiguado (oscila)\n\n";
    } else if (discriminant > 0) {
        std::cout << "Tipo: Sobreamortiguado (no oscila)\n\n";
    } else {
        std::cout << "Tipo: Críticamente amortiguado\n\n";
    }
    
    // Simular 10 segundos
    oscillator.simulate(x0, v0, 10.0, 1000, "oscillator_data.csv");
    
    std::cout << "✓ Datos guardados en oscillator_data.csv\n";
    std::cout << "\nPara visualizar, ejecuta:\n";
    std::cout << "  python3 visualize_oscillator.py\n";
    
    return 0;
}
