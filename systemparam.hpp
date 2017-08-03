#pragma once
//------------------------------------------------------------------------
#include <math.h>
//------------------------------------------------------------------------
extern double L;
const double dt = 0.005;
const double CUTOFF = 3.0;
const double MARGIN = 0.5;
const double ML2 = (CUTOFF + MARGIN) * (CUTOFF + MARGIN);
const double CL2 = (CUTOFF*CUTOFF);
const double RC2 = 1.0 / CL2;
const double RC6 = RC2 * RC2 * RC2;
const double RC12 = RC6 * RC6;
const double C0 = - 4.0 * (RC12 - RC6);
const double WCA_CL2 = pow(2.0, 1.0 / 3.0);


inline void adjust_periodic(double &dx, double &dy, double &dz) {
  const double LH = L * 0.5;
  if (dx < -LH)dx += L;
  if (dx > LH) dx -= L;
  if (dy < -LH)dy += L;
  if (dy > LH) dy -= L;
  if (dz < -LH)dz += L;
  if (dz > LH) dz -= L;
}
//------------------------------------------------------------------------
