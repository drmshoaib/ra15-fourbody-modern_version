#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "../../Integrator/include/Force.h"
#include "../../Integrator/include/RA15.h"

static double calcEnergy(const double* x, const double* v, const double* m)
{
    auto dist = [&](int i, int j) {
        const double dx = x[2*i] - x[2*j];
        const double dy = x[2*i+1] - x[2*j+1];
        return std::sqrt(dx*dx + dy*dy);
    };

    const double U = m[0]*m[1]/dist(0,1) + m[0]*m[2]/dist(0,2) + m[0]*m[3]/dist(0,3)
                   + m[1]*m[2]/dist(1,2) + m[1]*m[3]/dist(1,3) + m[2]*m[3]/dist(2,3);

    double K = 0.0;
    for (int i = 0; i < 4; ++i)
        K += 0.5 * m[i] * (v[2*i]*v[2*i] + v[2*i+1]*v[2*i+1]);

    return K - U;
}

int main(int argc, char* argv[])
{
    const std::string inPath  = (argc > 1) ? argv[1] : "data/initialcond.txt";
    const std::string outPath = (argc > 2) ? argv[2] : "result1.txt";

    double X[8], V[8], M[4], T;

    std::ifstream in(inPath);
    if (!in) {
        std::cerr << "Cannot open input file: " << inPath << "\n";
        return 1;
    }
    for (int i = 0; i < 8; ++i) in >> X[i];
    for (int i = 0; i < 8; ++i) in >> V[i];
    for (int i = 0; i < 4; ++i) in >> M[i];
    in >> T;
    in.close();

    std::ofstream out(outPath);
    if (!out) {
        std::cerr << "Cannot create output file: " << outPath << "\n";
        return 1;
    }
    out << std::setprecision(20) << std::scientific;

    // ll=8  → adaptive step size, accuracy ~10^-8 per Radau sequence
    // outputInterval → time between recorded snapshots
    const int    ll             = 8;
    const double outputInterval = 0.1;
    const int    outputSteps    = static_cast<int>(T / outputInterval);

    const double E0      = calcEnergy(X, V, M);
    double       time    = 0.0;
    double       maxDrift = 0.0;

    out << "# time x0 y0 x1 y1 x2 y2 x3 y3 vx0 vy0 vx1 vy1 vx2 vy2 vx3 vy3 energy\n";

    for (int i = 0; i < outputSteps; ++i)
    {
        const double E     = calcEnergy(X, V, M);
        const double drift = std::abs((E - E0) / E0);
        if (drift > maxDrift) maxDrift = drift;

        out << time
            << " " << X[0] << " " << X[1] << " " << X[2] << " " << X[3]
            << " " << X[4] << " " << X[5] << " " << X[6] << " " << X[7]
            << " " << V[0] << " " << V[1] << " " << V[2] << " " << V[3]
            << " " << V[4] << " " << V[5] << " " << V[6] << " " << V[7]
            << " " << E << "\n";

        // Integrate one output interval with adaptive step size
        ra15(X, V, outputInterval, 0.0, ll, 8, -2, computeForce, M);
        time += outputInterval;
    }

    std::cout << "Done. " << outputSteps << " snapshots written to " << outPath << "\n";
    std::cout << "Max |dE/E0| over full run: " << maxDrift << "\n";
    return 0;
}
