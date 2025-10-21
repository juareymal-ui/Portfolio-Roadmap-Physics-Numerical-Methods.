#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

int main() {
    // Parameters
    const double c = 1.0;          // wave speed
    const double R = 1.0;          // max radius
    const int Nr = 200;            // number of spatial points
    const double dr = R / Nr;
    const double dt = 0.5 * dr / c; // CFL condition: dt <= dr/c
    const int Nt = 500;            // number of time steps

    // Spatial grid
    std::vector<double> r(Nr+1);
    for(int j=0;j<=Nr;++j) r[j] = j*dr;

    // Initialize v = r*u
    std::vector<double> v_prev(Nr+1,0.0);
    std::vector<double> v_curr(Nr+1,0.0);
    std::vector<double> v_next(Nr+1,0.0);

    // Example initial condition: Gaussian pulse
    for(int j=0;j<=Nr;++j) {
        double x = r[j] - 0.5*R;
        double u0 = std::exp(-100*x*x);   // u(r,0)
        v_curr[j] = r[j]*u0;
        v_prev[j] = v_curr[j];            // initial velocity = 0
    }

    // Open CSV
    std::ofstream ofs("wave_spherical.csv");
    ofs << std::setprecision(12);
    ofs << "t";
    for(int j=0;j<=Nr;++j) ofs << ",u" << j;
    ofs << "\n";

    // Time stepping
    for(int n=0;n<Nt;++n) {
        double t = n*dt;

        // Output current solution u = v/r
        ofs << t;
        for(int j=0;j<=Nr;++j) {
            double u = (r[j] != 0.0) ? v_curr[j]/r[j] : 0.0;
            ofs << "," << u;
        }
        ofs << "\n";

        // Leapfrog scheme
        for(int j=1;j<Nr;++j) {
            double lambda = c*dt/dr;
            v_next[j] = 2*v_curr[j] - v_prev[j] + lambda*lambda*(v_curr[j+1]-2*v_curr[j]+v_curr[j-1]);
        }

        // Boundary conditions (fixed ends)
        v_next[0] = 0.0;
        v_next[Nr] = 0.0;

        // Advance time steps
        v_prev = v_curr;
        v_curr = v_next;
    }

    ofs.close();
    std::cout << "Simulation complete. Output saved to wave_spherical.csv\n";
    return 0;
}
