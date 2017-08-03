#pragma once
#include <memory>
#include "variables.hpp"
//------------------------------------------------------------------------
class Observer {
public:
  double kinetic_energy(std::shared_ptr<Variables> &vars);
  double potential_energy(std::shared_ptr<Variables> &vars, std::vector<Pair> & pairs);
  void local_energy(std::vector<double> &local_energy, std::shared_ptr<Variables> &vars, std::vector<Pair> & pairs);
  double temperature(std::shared_ptr<Variables> &vars) {return kinetic_energy(vars) / 1.5;}
  double total_energy(std::shared_ptr<Variables> &vars, std::vector<Pair> &pairs) {return kinetic_energy(vars) + potential_energy(vars, pairs);}
  void two_temperatures(double &t1, double &t2, std::shared_ptr<Variables> &vars);
  void two_temperatures_droplet(double r, double &t1, double &t2, std::shared_ptr<Variables> &vars);
  void liquidgas_temperatures(double &t1, double &t2, std::shared_ptr<Variables> &vars);
};
//------------------------------------------------------------------------
