#include "ode_solvers.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <array>

class PlanetaryMotion {
private:
    double G;  // Constante gravitacional
    double M;  // Masa del sol
    
public:
    PlanetaryMotion(double mass_sun) : G(6.67430e-11), M(mass_sun) {}
    
    // Estado: [x, y, vx, vy]
    numerical::StateVector derivatives(double t, const numerical::StateVector& state) {
        double x = state[0], y = state[1];
        double vx = state[2], vy = state[3];
        
        double r = std::sqrt(x*x + y*y);
        double r3 = r * r * r;
        
        double ax = -G * M * x / r3;
        double ay = -G * M * y / r3;
        
        return {vx, vy, ax, ay};
    }
    
    void simulate(const numerical::StateVector& initial, 
                  double t_max, int steps, const std::string& filename) {
        
        auto f = [this](double t, const numerical::StateVector& s) {
            return this->derivatives(t, s);
        };
        
        auto result = numerical::runge_kutta_4_system(f, initial, 0.0, t_max, steps);
        
        std::ofstream file(filename);
        file << "time,x,y,vx,vy,r,energy\n";
        
        double dt = t_max / steps;
        for (int i = 0; i <= steps; i++) {
            double t = i * dt;
            double x = result[i][0];
            double y = result[i][1];
            double vx = result[i][2];
            double vy = result[i][3];
            
            double r = std::sqrt(x*x + y*y);
            double v_sq = vx*vx + vy*vy;
            double energy = 0.5 * v_sq - G * M / r;
            
            file << t << "," << x << "," << y << "," 
                 << vx << "," << vy << "," << r << "," << energy << "\n";
        }
        
        file.close();
    }
};

int main() {
    std::cout << "=== Movimiento Planetario (Tierra alrededor del Sol) ===\n\n";
    
    // Masa del Sol
    double M_sun = 1.989e30;  // kg
    
    PlanetaryMotion simulation(M_sun);
    
    // Condiciones iniciales de la Tierra
    double AU = 1.496e11;           // 1 Unidad Astronómica en metros
    double v_earth = 29783;         // Velocidad orbital en m/s
    
    numerical::StateVector initial = {AU, 0.0, 0.0, v_earth};
    
    std::cout << "Parámetros:\n";
    std::cout << "  Masa del Sol: " << M_sun << " kg\n";
    std::cout << "  Distancia inicial: " << AU/1e9 << " Gm\n";
    std::cout << "  Velocidad orbital: " << v_earth << " m/s\n\n";
    
    // Simular 1 año (31,536,000 segundos)
    double year = 365.25 * 24 * 3600;
    
    std::cout << "Simulando 1 año terrestre...\n";
    simulation.simulate(initial, year, 10000, "orbit_data.csv");
    
    std::cout << "\n✓ Datos guardados en orbit_data.csv\n";
    std::cout << "\nPara visualizar la órbita:\n";
    std::cout << "  python3 visualize_orbit.py\n";
    
    return 0;
}
