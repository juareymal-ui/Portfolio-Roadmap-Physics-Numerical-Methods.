#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <fstream>
#include <iomanip>

// ===================== BASE CLASSES =====================
class ODESystem {
public:
    virtual std::vector<double> computeDerivatives(const std::vector<double>& state, double t) const = 0;
    virtual double energy(const std::vector<double>& state) const { return 0.0; }
    virtual std::vector<std::string> getStateNames() const = 0;
    virtual ~ODESystem() = default;
};

class Integrator {
public:
    virtual void step(ODESystem& system, std::vector<double>& state, double& t, double dt) = 0;
    virtual ~Integrator() = default;
};

// ===================== PHYSICAL SYSTEMS =====================
class HarmonicOscillator : public ODESystem {
    double m, k;
public:
    HarmonicOscillator(double mass, double springConstant) : m(mass), k(springConstant) {}
    
    std::vector<double> computeDerivatives(const std::vector<double>& state, double t) const override {
        return {state[1], -k/m * state[0]}; // [dx/dt, dv/dt]
    }
    
    double energy(const std::vector<double>& state) const override {
        return 0.5 * m * state[1]*state[1] + 0.5 * k * state[0]*state[0];
    }
    
    std::vector<std::string> getStateNames() const override {
        return {"position", "velocity"};
    }
};

class Pendulum : public ODESystem {
    double g, L;
public:
    Pendulum(double gravity, double length) : g(gravity), L(length) {}
    
    std::vector<double> computeDerivatives(const std::vector<double>& state, double t) const override {
        return {state[1], -g/L * std::sin(state[0])}; // [dθ/dt, dω/dt]
    }
    
    double energy(const std::vector<double>& state) const override {
        return 0.5 * L*L * state[1]*state[1] + g * L * (1 - std::cos(state[0]));
    }
    
    std::vector<std::string> getStateNames() const override {
        return {"angle", "angular_velocity"};
    }
};

class NBodyGravitation : public ODESystem {
    std::vector<double> masses;
    double G;
    int n_bodies;
public:
    NBodyGravitation(const std::vector<double>& masses, double gravitationalConstant) 
        : masses(masses), G(gravitationalConstant), n_bodies(masses.size()) {}
    
    std::vector<double> computeDerivatives(const std::vector<double>& state, double t) const override {
        std::vector<double> derivatives(6*n_bodies, 0.0);
        
        // Copy velocities to position derivatives
        for (int i = 0; i < n_bodies; ++i) {
            derivatives[6*i] = state[6*i+3];   // dx/dt = vx
            derivatives[6*i+1] = state[6*i+4]; // dy/dt = vy
            derivatives[6*i+2] = state[6*i+5]; // dz/dt = vz
        }
        
        // Calculate accelerations
        for (int i = 0; i < n_bodies; ++i) {
            for (int j = i+1; j < n_bodies; ++j) {
                double dx = state[6*j] - state[6*i];
                double dy = state[6*j+1] - state[6*i+1];
                double dz = state[6*j+2] - state[6*i+2];
                double r2 = dx*dx + dy*dy + dz*dz;
                double r = std::sqrt(r2);
                double r3 = r2 * r;
                
                double force_mag = G / r3;
                double fx = force_mag * dx;
                double fy = force_mag * dy;
                double fz = force_mag * dz;
                
                derivatives[6*i+3] += fx * masses[j];
                derivatives[6*i+4] += fy * masses[j];
                derivatives[6*i+5] += fz * masses[j];
                
                derivatives[6*j+3] -= fx * masses[i];
                derivatives[6*j+4] -= fy * masses[i];
                derivatives[6*j+5] -= fz * masses[i];
            }
        }
        
        // Convert forces to accelerations
        for (int i = 0; i < n_bodies; ++i) {
            derivatives[6*i+3] /= masses[i];
            derivatives[6*i+4] /= masses[i];
            derivatives[6*i+5] /= masses[i];
        }
        
        return derivatives;
    }
    
    double energy(const std::vector<double>& state) const override {
        double kinetic = 0.0, potential = 0.0;
        for (int i = 0; i < n_bodies; ++i) {
            kinetic += 0.5 * masses[i] * (state[6*i+3]*state[6*i+3] + 
                                         state[6*i+4]*state[6*i+4] + 
                                         state[6*i+5]*state[6*i+5]);
            for (int j = i+1; j < n_bodies; ++j) {
                double dx = state[6*j] - state[6*i];
                double dy = state[6*j+1] - state[6*i+1];
                double dz = state[6*j+2] - state[6*i+2];
                double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                potential -= G * masses[i] * masses[j] / r;
            }
        }
        return kinetic + potential;
    }
    
    std::vector<std::string> getStateNames() const override {
        std::vector<std::string> names;
        for (int i = 0; i < n_bodies; ++i) {
            names.push_back("body" + std::to_string(i+1) + "_x");
            names.push_back("body" + std::to_string(i+1) + "_y");
            names.push_back("body" + std::to_string(i+1) + "_z");
            names.push_back("body" + std::to_string(i+1) + "_vx");
            names.push_back("body" + std::to_string(i+1) + "_vy");
            names.push_back("body" + std::to_string(i+1) + "_vz");
        }
        return names;
    }
};

// ===================== INTEGRATORS =====================
class RK4Integrator : public Integrator {
public:
    void step(ODESystem& system, std::vector<double>& state, double& t, double dt) override {
        auto k1 = system.computeDerivatives(state, t);
        std::vector<double> state2(state.size());
        for (size_t i = 0; i < state.size(); ++i) state2[i] = state[i] + 0.5*dt*k1[i];
        
        auto k2 = system.computeDerivatives(state2, t + 0.5*dt);
        std::vector<double> state3(state.size());
        for (size_t i = 0; i < state.size(); ++i) state3[i] = state[i] + 0.5*dt*k2[i];
        
        auto k3 = system.computeDerivatives(state3, t + 0.5*dt);
        std::vector<double> state4(state.size());
        for (size_t i = 0; i < state.size(); ++i) state4[i] = state[i] + dt*k3[i];
        
        auto k4 = system.computeDerivatives(state4, t + dt);
        
        for (size_t i = 0; i < state.size(); ++i) {
            state[i] += dt/6.0 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
        }
        t += dt;
    }
};

class AdaptiveRK45Integrator : public Integrator {
    double tol;
    double min_dt;
    int max_iterations;
public:
    AdaptiveRK45Integrator(double tolerance = 1e-6, double min_step = 1e-10, int max_iters = 1000) 
        : tol(tolerance), min_dt(min_step), max_iterations(max_iters) {}
    
    void step(ODESystem& system, std::vector<double>& state, double& t, double dt) override {
        std::vector<double> original_state = state;
        double original_t = t;
        
        for (int iter = 0; iter < max_iterations; ++iter) {
            state = original_state;
            t = original_t;
            
            // Dormand-Prince RK45 method
            auto k1 = system.computeDerivatives(state, t);
            std::vector<double> temp(state.size());
            for (size_t i = 0; i < state.size(); ++i) 
                temp[i] = state[i] + dt * (1.0/5.0) * k1[i];
            
            auto k2 = system.computeDerivatives(temp, t + dt/5.0);
            for (size_t i = 0; i < state.size(); ++i) 
                temp[i] = state[i] + dt * (3.0/40.0) * k1[i] + dt * (9.0/40.0) * k2[i];
            
            auto k3 = system.computeDerivatives(temp, t + 3.0*dt/10.0);
            for (size_t i = 0; i < state.size(); ++i) 
                temp[i] = state[i] + dt * (44.0/45.0) * k1[i] - dt * (56.0/15.0) * k2[i] + 
                         dt * (32.0/9.0) * k3[i];
            
            auto k4 = system.computeDerivatives(temp, t + 4.0*dt/5.0);
            for (size_t i = 0; i < state.size(); ++i) 
                temp[i] = state[i] + dt * (19372.0/6561.0) * k1[i] - dt * (25360.0/2187.0) * k2[i] + 
                         dt * (64448.0/6561.0) * k3[i] - dt * (212.0/729.0) * k4[i];
            
            auto k5 = system.computeDerivatives(temp, t + 8.0*dt/9.0);
            for (size_t i = 0; i < state.size(); ++i) 
                temp[i] = state[i] + dt * (9017.0/3168.0) * k1[i] - dt * (355.0/33.0) * k2[i] + 
                         dt * (46732.0/5247.0) * k3[i] + dt * (49.0/176.0) * k4[i] - 
                         dt * (5103.0/18656.0) * k5[i];
            
            auto k6 = system.computeDerivatives(temp, t + dt);
            
            // 5th order solution
            std::vector<double> state5(state.size());
            for (size_t i = 0; i < state.size(); ++i) {
                state5[i] = state[i] + dt * (35.0/384.0) * k1[i] + dt * (500.0/1113.0) * k3[i] + 
                           dt * (125.0/192.0) * k4[i] - dt * (2187.0/6784.0) * k5[i] + 
                           dt * (11.0/84.0) * k6[i];
            }
            
            // 4th order solution for error estimation
            std::vector<double> state4(state.size());
            for (size_t i = 0; i < state.size(); ++i) {
                state4[i] = state[i] + dt * (5179.0/57600.0) * k1[i] + dt * (7571.0/16695.0) * k3[i] + 
                           dt * (393.0/640.0) * k4[i] - dt * (92097.0/339200.0) * k5[i] + 
                           dt * (187.0/2100.0) * k6[i] + dt * (1.0/40.0) * k2[i];
            }
            
            // Error estimation
            double error = 0.0;
            for (size_t i = 0; i < state.size(); ++i) {
                double diff = std::abs(state5[i] - state4[i]);
                error += diff * diff;
            }
            error = std::sqrt(error / state.size());
            
            if (error <= tol || dt <= min_dt) {
                state = state5;
                t += dt;
                break;
            }
            
            // Adjust step size
            double scale = 0.9 * std::pow(tol/error, 0.2);
            scale = std::max(0.1, std::min(5.0, scale));
            dt *= scale;
            
            if (dt < min_dt) dt = min_dt;
        }
    }
};

class VelocityVerletIntegrator : public Integrator {
public:
    void step(ODESystem& system, std::vector<double>& state, double& t, double dt) override {
        size_t n = state.size() / 2; // Assume [positions, velocities]
        
        // Get current acceleration
        auto accel = system.computeDerivatives(state, t);
        std::vector<double> current_accel(n, 0.0);
        for (size_t i = 0; i < n; ++i) {
            current_accel[i] = accel[i + n]; // acceleration is second half of derivatives
        }
        
        // Update positions using current velocities and acceleration
        for (size_t i = 0; i < n; ++i) {
            state[i] += state[i + n] * dt + 0.5 * current_accel[i] * dt * dt;
        }
        
        // Get new acceleration
        auto new_accel_derivs = system.computeDerivatives(state, t + dt);
        std::vector<double> new_accel(n, 0.0);
        for (size_t i = 0; i < n; ++i) {
            new_accel[i] = new_accel_derivs[i + n];
        }
        
        // Update velocities using average of old and new acceleration
        for (size_t i = 0; i < n; ++i) {
            state[i + n] += 0.5 * (current_accel[i] + new_accel[i]) * dt;
        }
        
        t += dt;
    }
};

// ===================== SIMULATION & OUTPUT =====================
class DataLogger {
private:
    std::ofstream file;
    int precision;
    
public:
    DataLogger(const std::string& filename, int prec = 6) : precision(prec) {
        file.open(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
        }
    }
    
    ~DataLogger() {
        if (file.is_open()) {
            file.close();
        }
    }
    
    void writeHeader(const std::vector<std::string>& column_names) {
        if (!file.is_open()) return;
        
        file << std::setprecision(precision);
        for (size_t i = 0; i < column_names.size(); ++i) {
            file << column_names[i];
            if (i < column_names.size() - 1) file << ",";
        }
        file << std::endl;
    }
    
    void writeRow(const std::vector<double>& values) {
        if (!file.is_open()) return;
        
        file << std::setprecision(precision);
        for (size_t i = 0; i < values.size(); ++i) {
            file << values[i];
            if (i < values.size() - 1) file << ",";
        }
        file << std::endl;
    }
    
    bool isOpen() const { return file.is_open(); }
};

void simulate(const std::string& filename, const std::string& systemName, 
              ODESystem& system, Integrator& integrator, 
              std::vector<double> initialState, double tStart, double tEnd, 
              double dt, int maxSteps = 10000) {
    
    // Create data logger
    DataLogger logger(filename);
    if (!logger.isOpen()) {
        std::cerr << "Failed to create data file: " << filename << std::endl;
        return;
    }
    
    // Prepare column names
    std::vector<std::string> column_names = {"time", "energy", "energy_error"};
    auto state_names = system.getStateNames();
    column_names.insert(column_names.end(), state_names.begin(), state_names.end());
    
    // Write header
    logger.writeHeader(column_names);
    
    std::vector<double> state = initialState;
    double t = tStart;
    double initialEnergy = system.energy(state);
    
    std::cout << "Simulating " << systemName << " from t=" << tStart << " to t=" << tEnd << std::endl;
    std::cout << "Saving data to: " << filename << std::endl;
    std::cout << "Initial energy: " << initialEnergy << std::endl;
    
    int stepCount = 0;
    int outputFrequency = 1;
    
    // Determine output frequency based on total steps
    int totalExpectedSteps = (tEnd - tStart) / dt;
    if (totalExpectedSteps > 1000) {
        outputFrequency = totalExpectedSteps / 1000;
    }
    
    while (t < tEnd && stepCount < maxSteps) {
        double stepSize = std::min(dt, tEnd - t);
        integrator.step(system, state, t, stepSize);
        stepCount++;
        
        // Write data at specified frequency
        if (stepCount % outputFrequency == 0) {
            double currentEnergy = system.energy(state);
            std::vector<double> row_data = {t, currentEnergy, currentEnergy - initialEnergy};
            row_data.insert(row_data.end(), state.begin(), state.end());
            logger.writeRow(row_data);
        }
        
        // Progress indicator
        if (stepCount % 1000 == 0) {
            std::cout << "Progress: " << (t - tStart) / (tEnd - tStart) * 100 << "%" << std::endl;
        }
    }
    
    // Write final state
    double finalEnergy = system.energy(state);
    std::vector<double> final_row = {t, finalEnergy, finalEnergy - initialEnergy};
    final_row.insert(final_row.end(), state.begin(), state.end());
    logger.writeRow(final_row);
    
    std::cout << "Final energy: " << finalEnergy << ", Energy error: " << finalEnergy - initialEnergy << std::endl;
    std::cout << "Completed in " << stepCount << " steps" << std::endl;
    std::cout << "Data saved to: " << filename << "\n" << std::endl;
}

int main() {
    // Harmonic Oscillator: mass=1.0, spring constant=1.0
    HarmonicOscillator ho(1.0, 1.0);
    std::vector<double> hoState = {1.0, 0.0}; // x=1.0, v=0.0

    // Pendulum: gravity=9.8, length=1.0
    Pendulum pendulum(9.8, 1.0);
    std::vector<double> pendulumState = {M_PI/4, 0.0}; // θ=π/4, ω=0.0

    // N-Body: 2 bodies (Sun and Earth) - simplified for demonstration
    std::vector<double> masses = {1.0, 0.000003}; // Normalized masses
    NBodyGravitation nbody(masses, 1.0); // Normalized G
    std::vector<double> nbodyState(12);
    // Body 1 at origin, body 2 in circular orbit
    nbodyState[0] = 0; nbodyState[1] = 0; nbodyState[2] = 0; // Body 1 position
    nbodyState[6] = 1.0; nbodyState[7] = 0; nbodyState[8] = 0; // Body 2 position
    nbodyState[3] = 0; nbodyState[4] = 0; nbodyState[5] = 0; // Body 1 velocity
    nbodyState[9] = 0; nbodyState[10] = 1.0; nbodyState[11] = 0; // Body 2 velocity

    // Integrators
    RK4Integrator rk4;
    AdaptiveRK45Integrator rk45(1e-6, 1e-10, 1000);
    VelocityVerletIntegrator vv;

    // Run simulations with CSV output
    std::cout << "=== HARMONIC OSCILLATOR ===" << std::endl;
    simulate("harmonic_oscillator_rk4.csv", "Harmonic Oscillator (RK4)", ho, rk4, hoState, 0.0, 10.0, 0.01, 1000);
    
    std::cout << "=== PENDULUM ===" << std::endl;
    simulate("pendulum_adaptive_rk45.csv", "Pendulum (Adaptive RK45)", pendulum, rk45, pendulumState, 0.0, 10.0, 0.1, 1000);
    
    std::cout << "=== N-BODY GRAVITATION ===" << std::endl;
    simulate("nbody_velocity_verlet.csv", "N-Body (Velocity Verlet)", nbody, vv, nbodyState, 0.0, 10.0, 0.01, 1000);

    return 0;
}