#pragma once

// Function signature for the force/acceleration evaluator:
//   X[nv]   current positions
//   V[nv]   current velocities (may be ignored for conservative systems)
//   t       current time
//   F[nv]   output: accelerations
//   masses  extra parameter block passed through unchanged
using ForceFunc = void(*)(const double* X, const double* V, double t,
                          double* F, const double* masses);

// Everhart's 15th-order Gauss-Radau integrator (RA15).
//
// x[nv], v[nv]  in/out: position and velocity
// tf            integration span (negative for backward integration)
// xl            initial sequence size hint (0 = auto)
// ll            accuracy control: ss = 10^(-ll), typical range 6-12;
//               if ll < 0, xl is used as a fixed (constant) step size
// nv            number of simultaneous second-order equations
// nclass        -2 = y'' = F(y,t)     [conservative, most common]
//                2 = y'' = F(y',y,t)  [general second-order]
//                1 = y'  = F(y,t)     [first-order]
// force         acceleration evaluator
// masses        passed through to every force() call
void ra15(double* x, double* v, double tf, double xl, int ll,
          int nv, int nclass, ForceFunc force, const double* masses);
