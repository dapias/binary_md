//------------------------------------------------------------------------
#include <iostream>
#include <assert.h>
#include "systemparam.hpp"
#include "observer.hpp"
//------------------------------------------------------------------------
void
Observer::two_temperatures(double &t1, double &t2, std::shared_ptr<Variables> &vars) {
  int na = 0;
  int nb = 0;
  double ka = 0.0;
  double kb = 0.0;
  for (auto &a : vars->atoms) {
    double v2 = a.px * a.px + a.py * a.py + a.pz * a.pz;
    if (a.type == 0) {
      ++na;
      ka += v2;
    } else {
      ++nb;
      kb += v2;
    }
  }
  t1 = ka / na / 3.0;
  t2 = kb / nb / 3.0;
}
//------------------------------------------------------------------------
void
Observer::two_temperatures_droplet(double r, double &t1, double &t2, std::shared_ptr<Variables> &vars) {
  int na = 0;
  int nb = 0;
  double ka = 0.0;
  double kb = 0.0;
  const double LH = L * 0.5;
  const double r2 = r * r;
  for (auto &a : vars->atoms) {
    double x = a.qx - LH;
    double y = a.qy - LH;
    double z = a.qz - LH;
    double v2 = a.px * a.px + a.py * a.py + a.pz * a.pz;
    if ((x * x + y * y + z * z) < r2) {
      ++na;
      ka += v2;
    } else {
      ++nb;
      kb += v2;
    }
  }
  t1 = ka / na / 3.0;
  t2 = kb / nb / 3.0;
}
//------------------------------------------------------------------------
double
Observer::kinetic_energy(std::shared_ptr<Variables> &vars) {
  double k = 0;
  for (auto &a : vars->atoms) {
    k += a.px * a.px;
    k += a.py * a.py;
    k += a.pz * a.pz;
  }
  k /= static_cast<double>(vars->number_of_atoms());
  return k * 0.5;
};
//------------------------------------------------------------------------
double
Observer::potential_energy(std::shared_ptr<Variables> &vars, std::vector<Pair> &pairs) {
  double v = 0.0;
  const int pp = pairs.size();
  const int pn = vars->number_of_atoms();
  Atom *atoms = vars->atoms.data();
  for (int k = 0; k < pp; k++) {
    const int i = pairs[k].i;
    const int j = pairs[k].j;
    double dx = atoms[j].qx - atoms[i].qx;
    double dy = atoms[j].qy - atoms[i].qy;
    double dz = atoms[j].qz - atoms[i].qz;
    adjust_periodic(dx, dy, dz);
    double r2 = (dx * dx + dy * dy + dz * dz);
    if (atoms[i].type != atoms[j].type) {
      if (r2 > WCA_CL2)continue;
      double r6 = r2 * r2 * r2;
      double r12 = r6 * r6;
      v += 4.0 * (1.0 / r12 - 1.0 / r6) + 1.0;
    } else {
      if (r2 > CL2)continue;
      double r6 = r2 * r2 * r2;
      double r12 = r6 * r6;
      v += 4.0 * (1.0 / r12 - 1.0 / r6) + C0;
    }
  }
  v /= static_cast<double>(pn);
  return v;
}
//------------------------------------------------------------------------
void
Observer::liquidgas_temperatures(double &tl, double &tg, std::shared_ptr<Variables> &vars) {
  int na = 0;
  int nb = 0;
  double ka = 0.0;
  double kb = 0.0;
  const double LH = L * 0.5;
  for (auto &a : vars->atoms) {
    double v2 = a.px * a.px + a.py * a.py + a.pz * a.pz;
    if (a.qy < LH) {
      ++na;
      ka += v2;
    } else {
      ++nb;
      kb += v2;
    }
  }
  tl = ka / na / 3.0;
  tg = kb / nb / 3.0;
}
//------------------------------------------------------------------------
void
Observer::local_energy(std::vector<double> &local_energy, std::shared_ptr<Variables> &vars, std::vector<Pair> &pairs) {
  const int pp = pairs.size();
  const int pn = vars->number_of_atoms();
  local_energy.resize(pn);
  std::fill(local_energy.begin(), local_energy.end(), 0.0);
  Atom *atoms = vars->atoms.data();
  for (int k = 0; k < pp; k++) {
    const int i = pairs[k].i;
    const int j = pairs[k].j;
    double dx = atoms[j].qx - atoms[i].qx;
    double dy = atoms[j].qy - atoms[i].qy;
    double dz = atoms[j].qz - atoms[i].qz;
    adjust_periodic(dx, dy, dz);
    double r2 = (dx * dx + dy * dy + dz * dz);
    if (atoms[i].type != atoms[j].type) {
      if (r2 > WCA_CL2)continue;
      double r6 = r2 * r2 * r2;
      double r12 = r6 * r6;
      double v =  4.0 * (1.0 / r12 - 1.0 / r6) + 1.0;
      local_energy[i] += v;
      local_energy[j] += v;
    } else {
      if (r2 > CL2)continue;
      double r6 = r2 * r2 * r2;
      double r12 = r6 * r6;
      double v = 4.0 * (1.0 / r12 - 1.0 / r6) + C0;
      local_energy[i] += v;
      local_energy[j] += v;
    }
  }
  for (int i = 0; i < pn; i++) {
    double px = atoms[i].px;
    double py = atoms[i].py;
    double pz = atoms[i].pz;
    local_energy[i] += 0.5 * (px * px + py * py + pz * pz);
  }
}
//------------------------------------------------------------------------
