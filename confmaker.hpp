//------------------------------------------------------------------------
#pragma once
//------------------------------------------------------------------------
#include <memory>
#include "variables.hpp"
//------------------------------------------------------------------------
class AtomAdder {
public:
  virtual void add(double x, double y, double z, std::shared_ptr<Variables> &vars) {
    Atom a(x, y, z);
    vars->add_atom(a);
  };
};
//------------------------------------------------------------------------
class ConfMaker {
  static void make_FCC(double density, AtomAdder &adder, std::shared_ptr<Variables> &vars);
public:
  static void make_monoflat(double leftdensity, double rightdensity, std::shared_ptr<Variables> &vars);
  static void make_droplet(double radius, double density, std::shared_ptr<Variables> &vars);
  static void make_binarydroplet(double radius, const double density, std::shared_ptr<Variables> &vars);
  static void make_binaryflat(const double density, std::shared_ptr<Variables> &vars);
};
//------------------------------------------------------------------------
