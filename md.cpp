//------------------------------------------------------------------------
#include <iostream>
#include <assert.h>
#include <math.h>
#include <random>
#include <fstream>
#include "md.hpp"
#include "confmaker.hpp"
#include "systemparam.hpp"
#include "observer.hpp"
#include "variables.hpp"
#include "parameter.hpp"
#include "thermostats.hpp"
//------------------------------------------------------------------------
MD::MD(void) : vars(new Variables()), mesh(new MeshList()), obs(new Observer()) {
  margin_length = 0.0;
  aimed_temperature = 1.0;
  runmap["MonoThermalize"] = &MD::run_mono_thermalize;
  runmap["MonoRestart"] = &MD::run_mono_restart;
  runmap["Binary"] = &MD::run_binary;
  thermostat_map["None"] = std::make_shared<None>();
  thermostat_map["Langevin"] = std::make_shared<Langevin>();
  thermostat_map["VelocityScaling"] = std::make_shared<VelocityScaling>();
  thermostat_map["NoseHoover"] = std::make_shared<NoseHoover>();
  thermostat_map["KineticMoments"] = std::make_shared<KineticMoments>();
  thermostat_map["NoseHooverLangevin"] = std::make_shared<KineticMoments>();
}
//------------------------------------------------------------------------
void
MD::make_pair(void) {
  pairs.clear();
  Atom *atoms = vars->atoms.data();
  const int pn = vars->number_of_atoms();
  for (int i = 0; i < pn - 1; i++) {
    for (int j = i + 1; j < pn; j++) {
      double dx = atoms[j].qx - atoms[i].qx;
      double dy = atoms[j].qy - atoms[i].qy;
      double dz = atoms[j].qz - atoms[i].qz;
      adjust_periodic(dx, dy, dz);
      double r2 = (dx * dx + dy * dy + dz * dz);
      if (r2 > ML2)continue;
      Pair p;
      p.i = i;
      p.j = j;
      pairs.push_back(p);
    }
  }
}
//------------------------------------------------------------------------
void
MD::check_pairlist(void) {
  double vmax2 = 0.0;
  for (auto &a : vars->atoms) {
    double v2 = a.px * a.px + a.py * a.py + a.pz * a.pz;
    if (vmax2 < v2) vmax2 = v2;
  }
  double vmax = sqrt(vmax2);
  margin_length -= vmax * 2.0 * dt;
  if (margin_length < 0.0) {
    margin_length = MARGIN;
    mesh->make_pair(vars, pairs);
  }
}
//------------------------------------------------------------------------
void
MD::update_position(void) {
  const double dt2 = dt * 0.5;
  for (auto &a : vars->atoms) {
    a.qx += a.px * dt2;
    a.qy += a.py * dt2;
    a.qz += a.pz * dt2;
  }
}
//------------------------------------------------------------------------
void
MD::calculate_force(void) {
  const int pn = vars->number_of_atoms();
  Atom *atoms = vars->atoms.data();
  for (int i = 0; i < pn - 1; i++) {
    for (int j = i + 1; j < pn; j++) {
      double dx = atoms[j].qx - atoms[i].qx;
      double dy = atoms[j].qy - atoms[i].qy;
      double dz = atoms[j].qz - atoms[i].qz;
      adjust_periodic(dx, dy, dz);
      double r2 = (dx * dx + dy * dy + dz * dz);
      if (r2 > CL2)continue;
      double r6 = r2 * r2 * r2;
      double df = (24.0 * r6 - 48.0) / (r6 * r6 * r2) * dt;
      atoms[i].px += df * dx;
      atoms[i].py += df * dy;
      atoms[i].pz += df * dz;
      atoms[j].px -= df * dx;
      atoms[j].py -= df * dy;
      atoms[j].pz -= df * dz;
    }
  }
}
//------------------------------------------------------------------------
void
MD::calculate_force_pair(void) {
  const int pp = pairs.size();
  Atom *atoms = vars->atoms.data();
  for (int k = 0; k < pp; k++) {
    const int i = pairs[k].i;
    const int j = pairs[k].j;
    double dx = atoms[j].qx - atoms[i].qx;
    double dy = atoms[j].qy - atoms[i].qy;
    double dz = atoms[j].qz - atoms[i].qz;
    adjust_periodic(dx, dy, dz);
    double r2 = (dx * dx + dy * dy + dz * dz);
    if (r2 > CL2)continue;
    double r6 = r2 * r2 * r2;
    double df = (24.0 * r6 - 48.0) / (r6 * r6 * r2) * dt;
    atoms[i].px += df * dx;
    atoms[i].py += df * dy;
    atoms[i].pz += df * dz;
    atoms[j].px -= df * dx;
    atoms[j].py -= df * dy;
    atoms[j].pz -= df * dz;
  }
}
//------------------------------------------------------------------------
void
MD::calculate_force_list(void) {
  Atom *atoms = vars->atoms.data();
  const int pn = vars->number_of_atoms();
  int *neighbor_list = vars->neighbor_list.data();
  int *i_position = vars->i_position.data();
  int *j_count = vars->j_count.data();
  for (int i = 0; i < pn; i++) {
    const double qix = atoms[i].qx;
    const double qiy = atoms[i].qy;
    const double qiz = atoms[i].qz;
    double pix = atoms[i].px;
    double piy = atoms[i].py;
    double piz = atoms[i].pz;
    const int ip = i_position[i];
    for (int k = 0; k < j_count[i]; k++) {
      const int j = neighbor_list[ip + k];
      double dx = atoms[j].qx - qix;
      double dy = atoms[j].qy - qiy;
      double dz = atoms[j].qz - qiz;
      adjust_periodic(dx, dy, dz);
      double r2 = (dx * dx + dy * dy + dz * dz);
      if (atoms[i].type != atoms[j].type) {
        if (r2 > WCA_CL2)continue;
      } else {
        if (r2 > CL2)continue;
      }
      double r6 = r2 * r2 * r2;
      double df = (24.0 * r6 - 48.0) / (r6 * r6 * r2) * dt;
      pix += df * dx;
      piy += df * dy;
      piz += df * dz;
      atoms[j].px -= df * dx;
      atoms[j].py -= df * dy;
      atoms[j].pz -= df * dz;
    }
    atoms[i].px = pix;
    atoms[i].py = piy;
    atoms[i].pz = piz;
  }
}
//------------------------------------------------------------------------
void
MD::periodic(void) {
  for (auto &a : vars->atoms) {
    if (a.qx < 0.0) a.qx += L;
    if (a.qy < 0.0) a.qy += L;
    if (a.qz < 0.0) a.qz += L;
    if (a.qx > L) a.qx -= L;
    if (a.qy > L) a.qy -= L;
    if (a.qz > L) a.qz -= L;
    assert(a.qx < L);
    assert(a.qy < L);
    assert(a.qz < L);
  }
}
//------------------------------------------------------------------------
void
MD::calculate(void) {
  update_position();
  check_pairlist();
  calculate_force_list();
  update_position();
  thermostat->update(this);
  periodic();
  vars->time += dt;
}
//------------------------------------------------------------------------
bool
MD::init(Parameter &param) {
  L = param.GetIntegerDef("SystemSize", 30);
  mesh->init();
  TotalSteps = param.GetIntegerDef("TotalSteps", 20000);
  ObserveInterval = param.GetIntegerDef("ObserveInterval", 100);
  aimed_temperature = param.GetDoubleDef("AimedTemperature", 0.8);
  hbtype = param.GetStringDef("Thermostat", "NoseHoover");
  if (thermostat_map.count(hbtype) == 0) {
    std::cerr << "Error: Unknown Thermstat " << hbtype << std::endl;
    return false;
  }
  thermostat = thermostat_map[hbtype];
  return true;
}
//------------------------------------------------------------------------
void
MD::showinfo(void) {
  const int N = vars->number_of_atoms();
  const double density = static_cast<double>(N) / L / L / L;
  std::cout << "# N = " << vars->number_of_atoms() << std::endl;
  std::cout << "# L = " << L << std::endl;
  std::cout << "# Density = " << density << std::endl;
  std::cout << "# CUTOFF = " << CUTOFF << std::endl;
  std::cout << "# TotalSteps = " << TotalSteps << std::endl;
  std::cout << "# ObserveInterval = " << ObserveInterval << std::endl;
  std::cout << "# dt = " << dt << std::endl;
  std::cout << "# Thermostat = " << hbtype << std::endl;
  std::cout << "# aimed_temperature = " << aimed_temperature << std::endl;

}
//------------------------------------------------------------------------
void
MD::run_mono_thermalize(Parameter &param) {
  const double radius = param.GetDoubleDef("DropletRadius", 20);
  double density = param.GetDoubleDef("Density", 0.5);
  ConfMaker::make_droplet(radius, density, vars);
  vars->set_initial_temperature(aimed_temperature);
  showinfo();
  for (int i = 0; i < TotalSteps; i++) {
    if ( (i % ObserveInterval) == 0) {
      std::cout << vars->time << " ";
      std::cout << obs->temperature(vars) << " ";
      std::cout << obs->kinetic_energy(vars) << " ";
      std::cout << obs->potential_energy(vars, pairs) << std::endl;
      vars->export_cdview();
    }
    calculate();
  }
  vars->save_as_cdview("droplet.cdv");
}
//------------------------------------------------------------------------
void
MD::run_mono_restart(Parameter &param) {
  std::string import_file = param.GetStringDef("ImportFile", "droplet.cdv");
  const double radius = param.GetDoubleDef("DropletRadius", 20);
  vars->import_cdview(import_file.c_str());
  showinfo();
  const double LH = L * 0.5;
  const double scale = sqrt(0.9);
  const double iscale = 1.0 / scale;
  for (auto &a : vars->atoms) {
    double x = a.qx - LH;
    double y = a.qy - LH;
    double z = a.qz - LH;
    if (x * x + y * y + z * z < radius * radius) {
      a.px *= scale;
      a.py *= scale;
      a.pz *= scale;
    } else {
      a.px *= iscale;
      a.py *= iscale;
      a.pz *= iscale;
    }
  }
  double t1s = 0.0;
  double t2s = 0.0;
  double ts = 0.0;
  for (int i = 0; i < TotalSteps; i++) {
    calculate();
    double t1, t2;
    obs->two_temperatures_droplet(radius, t1, t2, vars);
    t1s += t1;
    t2s += t2;
    ts += obs->temperature(vars);
    if (i != 0 &&  (i % ObserveInterval) == 0) {
      double k = obs->kinetic_energy(vars);
      double v = obs->potential_energy(vars, pairs);
      const double di = static_cast<double>(ObserveInterval);
      ts /= di;
      t1s /= di;
      t2s /= di;
      std::cout << vars->time << " ";
      std::cout << ts << " ";
      std::cout << t1s << " ";
      std::cout << t2s << " ";
      std::cout << v << " ";
      std::cout << k + v << " ";
      std::cout << std::endl;
      ts = 0.0;
      t1s = 0.0;
      t2s = 0.0;
    }
  }
}
//------------------------------------------------------------------------
void
MD::run_binary(Parameter & param) {
  const double density = param.GetDoubleDef("Density", 0.5);
  const std::string boundary = param.GetStringDef("Boundary", "Flat");
  if (boundary == "Droplet") {
    const double radius = param.GetDoubleDef("DropletRadius", 10);
    ConfMaker::make_binarydroplet(radius, density, vars);
    std::cout << "# Boundary = Droplet" << std::endl;
    std::cout << "# DropletRadius = " << radius << std::endl;
  } else {
    std::cout << "# Boundary = Flat" << std::endl;
    ConfMaker::make_binaryflat(density, vars);
  }
  showinfo();
  vars->set_initial_temperature(aimed_temperature);
  int n1 = 0;
  int n2 = 0;
  for (auto &a : vars->atoms) {
    if (a.type == 0)n1++;
    else n2++;
  }
  std::cout << "# A-atoms = " << n1 << std::endl;
  std::cout << "# B-atoms = " << n2 << std::endl;

  double t1s = 0.0;
  double t2s = 0.0;
  double ts = 0.0;
  for (int i = 0; i < TotalSteps; i++) {
    double t1, t2;
    obs->two_temperatures(t1, t2, vars);
    t1s += t1;
    t2s += t2;
    ts += obs->temperature(vars);
    if ( (i % ObserveInterval) == 0) {
      double k = obs->kinetic_energy(vars);
      double v = obs->potential_energy(vars, pairs);
      const double di = static_cast<double>(ObserveInterval);
      ts /= di;
      t1s /= di;
      t2s /= di;
      std::cout << vars->time << " ";
      std::cout << ts << " ";
      std::cout << t1s << " ";
      std::cout << t2s << " ";
      std::cout << v << " ";
      std::cout << k + v << " ";
      std::cout << std::endl;
      ts = 0.0;
      t1s = 0.0;
      t2s = 0.0;
      vars->export_cdview();
    }
    calculate();
  }
}
//------------------------------------------------------------------------
void
MD::run(Parameter & param) {
  if (!init(param))return;
  std::string mode = param.GetStringDef("Mode", "");
  if (mode == "") {
    std::cout << "Mode is not specified." << std::endl;
    return;
  } else if (runmap.count(mode) == 0) {
    std::cout << "Mode " << mode << " is not defined. " << std::endl;
    return;
  }
  std::cout << "# Mode = " << mode << std::endl;
  (this->*runmap[mode])(param);
}
//------------------------------------------------------------------------
