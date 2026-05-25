#include <cstdio>
#include <cmath>
#include <cstdlib>
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
    const char* inPath  = (argc > 1) ? argv[1] : "data/initialcond.txt";
    const char* outPath = (argc > 2) ? argv[2] : "result1.txt";

    double X[8], V[8], M[4], T;

    FILE* in = fopen(inPath, "r");
    if (!in) {
        fprintf(stderr, "Cannot open input file: %s\n", inPath);
        return 1;
    }
    for (int i = 0; i < 8; ++i) fscanf(in, "%lf", &X[i]);
    for (int i = 0; i < 8; ++i) fscanf(in, "%lf", &V[i]);
    for (int i = 0; i < 4; ++i) fscanf(in, "%lf", &M[i]);
    fscanf(in, "%lf", &T);
    fclose(in);

    FILE* out = fopen(outPath, "w");
    if (!out) {
        fprintf(stderr, "Cannot create output file: %s\n", outPath);
        return 1;
    }

    // Integration parameters matching original: fixed step 0.001, 1,000,000 steps
    const double timestep  = 0.001;
    const int    maxSteps  = 1000000;
    const double energyTol = 1.0e-9;

    const double E0 = calcEnergy(X, V, M);
    double time = 0.0;

    // Header
    fprintf(out, "# time x0 y0 x1 y1 x2 y2 x3 y3 vx0 vy0 vx1 vy1 vx2 vy2 vx3 vy3 energy\n");

    for (int i = 0; i < maxSteps; ++i)
    {
        fprintf(out,
            "%.20f %.20f %.20f %.20f %.20f "
            "%.20f %.20f %.20f %.20f "
            "%.20f %.20f %.20f %.20f "
            "%.20f %.20f %.20f %.20f %.20f\n",
            time,
            X[0], X[1], X[2], X[3], X[4], X[5], X[6], X[7],
            V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7],
            calcEnergy(X, V, M));

        // ll = -1 → constant step size XL = timestep
        ra15(X, V, timestep, timestep, -1, 8, -2, computeForce, M);

        const double E = calcEnergy(X, V, M);
        if (std::abs((E - E0) / E0) > energyTol) {
            fprintf(out, "# Energy conservation violated at t=%.6f (dE/E0=%.3e)\n",
                    time, (E - E0) / E0);
            fclose(out);
            return 1;
        }

        time += timestep;
    }

    fclose(out);
    printf("Done. Output written to %s\n", outPath);
    return 0;
}
