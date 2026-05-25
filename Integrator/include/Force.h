#pragma once

// Gravitational accelerations for 4 coplanar bodies in barycentric coordinates.
// X[8]: positions  {x0,y0, x1,y1, x2,y2, x3,y3}
// V[8]: velocities (unused; reserved for velocity-dependent extensions)
// t   : current time (unused)
// F[8]: output accelerations
// M[4]: masses
void computeForce(const double* X, const double* V, double t, double* F, const double* M);
