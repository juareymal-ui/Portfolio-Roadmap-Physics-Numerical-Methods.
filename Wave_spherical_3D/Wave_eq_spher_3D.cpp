#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

int main() {
    // Physical parameters
    const double c = 1.0;        // wave speed
    const double r_max = 1.0;
    const double theta_max = M_PI;
    const double phi_max = 2.0*M_PI;

    // Grid size
    const int Nr = 20;
    const int Ntheta = 20;
    const int Nphi = 20;
    const double dr = r_max / Nr;
    const double dtheta = theta_max / Ntheta;
    const double dphi = phi_max / Nphi;

    const double dt = 0.3 * std::min({dr, r_max*dtheta, r_max*dphi})/c; // CFL condition
    const int Nt = 50;

    // 3D arrays: use vector of vectors of vectors
    std::vector<std::vector<std::vector<double>>> u_prev(Nr+1,std::vector<std::vector<double>>(Ntheta+1,std::vector<double>(Nphi+1,0.0)));
    std::vector<std::vector<std::vector<double>>> u_curr = u_prev;
    std::vector<std::vector<std::vector<double>>> u_next = u_prev;

    // Initial condition: small Gaussian pulse at center
    for(int j=1;j<Nr;++j)
        for(int k=1;k<Ntheta;++k)
            for(int l=0;l<Nphi;++l){
                double r = j*dr;
                double theta = k*dtheta;
                double phi = l*dphi;
                double x = r*sin(theta)*cos(phi);
                double y = r*sin(theta)*sin(phi);
                double z = r*cos(theta);
                double dist2 = x*x + y*y + z*z;
                u_curr[j][k][l] = std::exp(-100*dist2);
                u_prev[j][k][l] = u_curr[j][k][l]; // zero initial velocity
            }

    // Open CSV for a single snapshot output (flattened)
    std::ofstream ofs("wave3d_snapshot.csv");
    ofs << std::setprecision(12);

    // Time stepping
    for(int n=0;n<Nt;++n){
        for(int j=1;j<Nr;++j){
            double r = j*dr;
            for(int k=1;k<Ntheta;++k){
                double theta = k*dtheta;
                for(int l=0;l<Nphi;++l){
                    int lp = (l+1)%Nphi;
                    int lm = (l-1+Nphi)%Nphi;

                    double dr_term = ( (r+0.5*dr)*(r+0.5*dr)*(u_curr[j+1][k][l]-u_curr[j][k][l])
                                     - (r-0.5*dr)*(r-0.5*dr)*(u_curr[j][k][l]-u_curr[j-1][k][l]) ) / (r*r*dr*dr);

                    double dtheta_term = ( sin(theta+0.5*dtheta)*(u_curr[j][k+1][l]-u_curr[j][k][l])
                                        - sin(theta-0.5*dtheta)*(u_curr[j][k][l]-u_curr[j][k-1][l]) ) / (r*r*sin(theta)*dtheta*dtheta);

                    double dphi_term = (u_curr[j][k][lp]-2*u_curr[j][k][l]+u_curr[j][k][lm]) / (r*r*sin(theta)*sin(theta)*dphi*dphi);

                    u_next[j][k][l] = 2*u_curr[j][k][l] - u_prev[j][k][l] + c*c*dt*dt*(dr_term + dtheta_term + dphi_term);
                }
            }
        }

        // Update time steps
        u_prev = u_curr;
        u_curr = u_next;
    }

    // Save a single snapshot (r-theta-phi flattened)
    for(int j=0;j<=Nr;++j)
        for(int k=0;k<=Ntheta;++k)
            for(int l=0;l<=Nphi;++l)
                ofs << j*dr << "," << k*dtheta << "," << l*dphi << "," << u_curr[j][k][l] << "\n";

    ofs.close();
    std::cout << "3D wave simulation done. Snapshot saved.\n";
    return 0;
}
