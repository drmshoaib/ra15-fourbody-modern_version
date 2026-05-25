#include "Force.h"
#include <cmath>

void computeForce(const double* X, const double* V, double t, double* F, const double* M)
{
    (void)V; (void)t;

    double rh[4][4] = {};

    // Inverse-cube distances for every pair (i < j)
    for (int n = 0; n < 4; ++n) {
        const int j = n * 2;
        for (int l = n + 1; l < 4; ++l) {
            const int k  = l * 2;
            const double dx = X[j] - X[k];
            const double dy = X[j + 1] - X[k + 1];
            const double r2 = dx * dx + dy * dy;
            rh[n][l] = 1.0 / (r2 * std::sqrt(r2));
            rh[l][n] = rh[n][l];
        }
    }

    // Accelerations: a_i = sum_{j!=i} M_j * (r_j - r_i) / |r_ij|^3
    for (int n = 0; n < 4; ++n) {
        const int j = n * 2;
        F[j]     = 0.0;
        F[j + 1] = 0.0;
        for (int l = 0; l < 4; ++l) {
            if (l == n) continue;
            const int k = l * 2;
            F[j]     += M[l] * (X[k]     - X[j])     * rh[n][l];
            F[j + 1] += M[l] * (X[k + 1] - X[j + 1]) * rh[n][l];
        }
    }
}
