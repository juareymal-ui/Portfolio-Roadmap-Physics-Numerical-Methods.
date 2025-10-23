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
