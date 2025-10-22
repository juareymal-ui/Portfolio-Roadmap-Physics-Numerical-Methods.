#include "ode_solvers.h"

namespace numerical {

std::vector<double> euler_method(
    std::function<double(double, double)> f,
    double y0,
    double t0,
    double tf,
    int n_steps
) {
    double h = (tf - t0) / n_steps;
    std::vector<double> y(n_steps + 1);
    
    y[0] = y0;
    double t = t0;
    
    for (int i = 0; i < n_steps; i++) {
        y[i + 1] = y[i] + h * f(t, y[i]);
        t += h;
    }
    
    return y;
}

std::vector<double> runge_kutta_4(
    std::function<double(double, double)> f,
    double y0,
    double t0,
    double tf,
    int n_steps
) {
    double h = (tf - t0) / n_steps;
    std::vector<double> y(n_steps + 1);
    
    y[0] = y0;
    double t = t0;
    
    for (int i = 0; i < n_steps; i++) {
        double k1 = h * f(t, y[i]);
        double k2 = h * f(t + h/2, y[i] + k1/2);
        double k3 = h * f(t + h/2, y[i] + k2/2);
        double k4 = h * f(t + h, y[i] + k3);
        
        y[i + 1] = y[i] + (k1 + 2*k2 + 2*k3 + k4) / 6.0;
        t += h;
    }
    
    return y;
}

std::vector<StateVector> runge_kutta_4_system(
    std::function<StateVector(double, const StateVector&)> f,
    const StateVector& y0,
    double t0,
    double tf,
    int n_steps
) {
    double h = (tf - t0) / n_steps;
    std::vector<StateVector> result(n_steps + 1, StateVector(y0.size()));
    
    result[0] = y0;
    double t = t0;
    
    for (int i = 0; i < n_steps; i++) {
        StateVector k1 = f(t, result[i]);
        
        StateVector temp1(y0.size());
        for (size_t j = 0; j < y0.size(); j++) {
            temp1[j] = result[i][j] + h * k1[j] / 2.0;
        }
        StateVector k2 = f(t + h/2, temp1);
        
        StateVector temp2(y0.size());
        for (size_t j = 0; j < y0.size(); j++) {
            temp2[j] = result[i][j] + h * k2[j] / 2.0;
        }
        StateVector k3 = f(t + h/2, temp2);
        
        StateVector temp3(y0.size());
        for (size_t j = 0; j < y0.size(); j++) {
            temp3[j] = result[i][j] + h * k3[j];
        }
        StateVector k4 = f(t + h, temp3);
        
        for (size_t j = 0; j < y0.size(); j++) {
            result[i + 1][j] = result[i][j] + 
                (h / 6.0) * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
        }
        
        t += h;
    }
    
    return result;
}

} // namespace numerical
