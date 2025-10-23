#include "physics_simulations.h"
#include <iostream>

int main() {
    std::cout << "=== Ecuación del Calor 1D ===\n\n";
    
    // Parámetros
    double alpha = 0.01;  // Difusividad térmica (m²/s)
    double L = 1.0;       // Longitud del dominio (m)
    int nx = 100;         // Número de puntos espaciales
    
    std::cout << "Parámetros:\n";
    std::cout << "  α (difusividad) = " << alpha << " m²/s\n";
    std::cout << "  L (longitud) = " << L << " m\n";
    std::cout << "  Puntos espaciales = " << nx << "\n";
    std::cout << "  Condición inicial: Pulso gaussiano\n";
    std::cout << "  Condiciones de frontera: T = 0 en los bordes\n\n";
    
    physics::HeatEquation1D heat(alpha, L, nx);
    
    // Resolver para 1 segundo
    double t_max = 1.0;
    int nt = 1000;
    
    std::cout << "Resolviendo ecuación del calor...\n";
    std::cout << "Tiempo de simulación: " << t_max << " s\n";
    std::cout << "Pasos temporales: " << nt << "\n\n";
    
    heat.solve(t_max, nt, "heat_data.csv");
    
    std::cout << "\n✓ Datos guardados en heat_data.csv\n";
    std::cout << "\nPara visualizar la evolución térmica:\n";
    std::cout << "  python3 visualize_heat.py\n";
    
    return 0;
}
