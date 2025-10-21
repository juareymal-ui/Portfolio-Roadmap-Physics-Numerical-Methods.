#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>

using namespace std;
using cd = complex<double>;

int main() {
    // Physical constants
    const double hbar = 1.0;
    const double m = 1.0;
    const double rmax = 10.0;

    // Numerical parameters
    const int Nr = 500;
    const double dr = rmax / Nr;
    const double dt = 0.005;
    const int Nt = 1000;

    // Radial grid
    vector<double> r(Nr+1);
    for(int j=0;j<=Nr;++j) r[j] = j*dr;

    // Potential: Hydrogen-like V = -1/r (avoid r=0)
    vector<double> V(Nr+1,0.0);
    for(int j=1;j<=Nr;++j) V[j] = -1.0 / r[j];
    V[0] = V[1]; // avoid division by zero

    // Initial wavefunction: Gaussian packet
    vector<cd> u(Nr+1,0.0);
    vector<cd> u_new(Nr+1,0.0);
    for(int j=1;j<=Nr;++j){
        double x = r[j]-5.0;
        u[j] = r[j]*exp(-10.0*x*x);
    }

    // Tridiagonal matrix elements for Crank-Nicolson
    cd im(0.0,1.0);
    vector<cd> a(Nr-1), b(Nr-1), c(Nr-1);
    cd alpha = im*dt/(2.0*hbar) * (hbar*hbar/(2.0*m*dr*dr));

    for(int j=0;j<Nr-1;++j){
        a[j] = -alpha;
        c[j] = -alpha;
        b[j] = 1.0 + 2.0*alpha + im*dt/(2.0*hbar)*V[j+1];
    }

    // Thomas algorithm helper (solve tridiagonal)
    auto thomas = [&](vector<cd>& a, vector<cd>& b, vector<cd>& c, vector<cd>& d){
        vector<cd> c_star(Nr-1), d_star(Nr-1);
        c_star[0] = c[0]/b[0];
        d_star[0] = d[0]/b[0];
        for(int i=1;i<Nr-1;++i){
            cd m = b[i]-a[i]*c_star[i-1];
            c_star[i] = c[i]/m;
            d_star[i] = (d[i]-a[i]*d_star[i-1])/m;
        }
        vector<cd> x(Nr-1);
        x[Nr-2] = d_star[Nr-2];
        for(int i=Nr-3;i>=0;--i) x[i] = d_star[i]-c_star[i]*x[i+1];
        return x;
    };

    // Time evolution
    ofstream ofs("tdse_radial.csv");
    ofs << setprecision(12);
    ofs << "t";
    for(int j=0;j<=Nr;++j) ofs << ",u_real_"<<j<<",u_imag_"<<j;
    ofs << "\n";

    for(int n=0;n<Nt;++n){
        double t = n*dt;

        // Output
        ofs << t;
        for(int j=0;j<=Nr;++j) ofs << "," << real(u[j]) << "," << imag(u[j]);
        ofs << "\n";

        // RHS vector d = (1 - i dt/2 H) u^n
        vector<cd> d(Nr-1);
        for(int j=1;j<Nr;++j)
            d[j-1] = (1.0 - 2.0*alpha - im*dt/(2.0*hbar)*V[j])*u[j] + alpha*(u[j+1]+u[j-1]);

        // Solve tridiagonal system
        vector<cd> u_inner = thomas(a,b,c,d);

        // Update u_new
        for(int j=1;j<Nr;++j) u[j] = u_inner[j-1];
        u[0] = 0.0; u[Nr] = 0.0;
    }

    ofs.close();
    cout << "TDSE radial simulation complete.\n";
    return 0;
}
