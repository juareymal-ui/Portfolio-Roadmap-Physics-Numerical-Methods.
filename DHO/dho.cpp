// dho.cpp
// Damped Harmonic Oscillator solver using RK4
// Equation: x'' + 2ζω x' + ω² x = 0
//
// Output: CSV file with time, x(t), v(t)

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>

// State = (x, v)
struct State { double x, v; };

// RHS of the system
// dx/dt = v
// dv/dt = -2ζω v - ω² x
void rhs(const State& s, State& dsdt, double omega, double zeta) {
    dsdt.x = s.v;
    dsdt.v = -2.0 * zeta * omega * s.v - (omega * omega) * s.x;
}

// One RK4 step
State rk4_step(const State& s, double dt, double omega, double zeta) {
    State k1, k2, k3, k4, tmp;

    rhs(s, k1, omega, zeta);

    tmp.x = s.x + 0.5 * dt * k1.x;
    tmp.v = s.v + 0.5 * dt * k1.v;
    rhs(tmp, k2, omega, zeta);

    tmp.x = s.x + 0.5 * dt * k2.x;
    tmp.v = s.v + 0.5 * dt * k2.v;
    rhs(tmp, k3, omega, zeta);

    tmp.x = s.x + dt * k3.x;
    tmp.v = s.v + dt * k3.v;
    rhs(tmp, k4, omega, zeta);

    State snew;
    snew.x = s.x + (dt / 6.0) * (k1.x + 2*k2.x + 2*k3.x + k4.x);
    snew.v = s.v + (dt / 6.0) * (k1.v + 2*k2.v + 2*k3.v + k4.v);
    return snew;
}

int main() {
    // Parameters
    double omega = 2.0 * M_PI;   // natural frequency (rad/s)
    double zeta  = 2.0;          // damping ratio (0 = no damping)
    double x0    = 1.0;          // initial displacement
    double v0    = 0.0;          // initial velocity
    double dt    = 0.001;        // timestep (s)
    double tmax  = 10.0;         // simulation duration (s)

    // Initial state
    State s{ x0, v0 };
    double t = 0.0;

    // Open CSV output
    std::ofstream ofs("dho_output.csv");
    ofs << std::setprecision(12);
    ofs << "t,x,v\n";

    // Time loop
    while (t <= tmax + 1e-12) {
        ofs << t << "," << s.x << "," << s.v << "\n";
        s = rk4_step(s, dt, omega, zeta);
        t += dt;
    }

    ofs.close();
    std::cout << "Done. Output written to dho_output.csv\n";
    return 0;
}
