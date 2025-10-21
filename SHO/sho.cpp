// sho.cpp
// Simple Harmonic Oscillator solver (RK4 and Velocity-Verlet)
// Outputs CSV: t,x_num,v_num,x_exact,energy_num
//
// Compile (Linux / MinGW):
//   g++ -std=c++17 -O2 sho.cpp -o sho
// Or with MSVC (Developer Command Prompt):
//   cl /EHsc /O2 sho.cpp

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>

struct State { double x, v; };

double exact_x(double x0, double v0, double omega, double t) {
    return x0 * std::cos(omega * t) + (v0 / omega) * std::sin(omega * t);
}

double energy(double x, double v, double m, double k) {
    // For SHO: potential = 1/2 k x^2, kinetic = 1/2 m v^2
    return 0.5 * m * v * v + 0.5 * k * x * x;
}

// RHS of system: dx/dt = v, dv/dt = -omega^2 * x
void rhs(const State& s, State& dsdt, double omega) {
    dsdt.x = s.v;
    dsdt.v = - (omega * omega) * s.x;
}

// One RK4 step for the system
State rk4_step(const State& s, double dt, double omega) {
    State k1, k2, k3, k4, tmp;

    rhs(s, k1, omega);

    tmp.x = s.x + 0.5 * dt * k1.x;
    tmp.v = s.v + 0.5 * dt * k1.v;
    rhs(tmp, k2, omega);

    tmp.x = s.x + 0.5 * dt * k2.x;
    tmp.v = s.v + 0.5 * dt * k2.v;
    rhs(tmp, k3, omega);

    tmp.x = s.x + dt * k3.x;
    tmp.v = s.v + dt * k3.v;
    rhs(tmp, k4, omega);

    State snew;
    snew.x = s.x + (dt / 6.0) * (k1.x + 2.0 * k2.x + 2.0 * k3.x + k4.x);
    snew.v = s.v + (dt / 6.0) * (k1.v + 2.0 * k2.v + 2.0 * k3.v + k4.v);

    return snew;
}

// One Velocity-Verlet step (symplectic)
State verlet_step(const State& s, double dt, double omega) {
    // acceleration a = -omega^2 x
    double a0 = - (omega * omega) * s.x;

    State s_half;
    s_half.x = s.x + s.v * dt + 0.5 * a0 * dt * dt; // position at t+dt
    double a1 = - (omega * omega) * s_half.x;       // acceleration at t+dt
    State snew;
    snew.x = s_half.x;
    snew.v = s.v + 0.5 * (a0 + a1) * dt;
    return snew;
}

int main(int argc, char** argv) {
    // Parameters (you can edit these or parse from argv)
    const double omega = 2.0 * M_PI; // angular frequency (rad/s) -> period T = 1 s
    const double m = 1.0;
    const double k = m * omega * omega;

    const double x0 = 1.0;    // initial position
    const double v0 = 0.0;    // initial velocity

    const double t0 = 0.0;
    const double tmax = 10.0; // seconds
    const double dt = 0.001;  // time step (try dt=0.001 or dt = T/200 => T=1 => dt=0.005)

    std::string method = "rk4";
    if (argc >= 2) method = argv[1]; // "rk4" or "verlet"

    std::cout << "SHO solver: method=" << method << " dt=" << dt << " omega=" << omega << "\n";

    std::ofstream ofs("sho_output.csv");
    ofs << std::setprecision(12);
    ofs << "t,x_num,v_num,x_exact,energy_num\n";

    State s{ x0, v0 };
    double t = t0;
    while (t <= tmax + 1e-12) {
        double x_ex = exact_x(x0, v0, omega, t);
        double E = energy(s.x, s.v, m, k);
        ofs << t << "," << s.x << "," << s.v << "," << x_ex << "," << E << "\n";

        // Step
        if (method == "rk4") {
            s = rk4_step(s, dt, omega);
        } else {
            s = verlet_step(s, dt, omega);
        }

        t += dt;
    }

    ofs.close();
    std::cout << "Done. Output written to sho_output.csv\n";
    return 0;
}
