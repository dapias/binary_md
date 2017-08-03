//------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "systemparam.hpp"
#include "confmaker.hpp"
//------------------------------------------------------------------------
void
ConfMaker::make_FCC(double density, AtomAdder &adder, std::shared_ptr<Variables> &vars) {
  const double s0 = 1.0 / pow(density * 0.25, 1.0 / 3.0);
  const int is = static_cast<int>(L / s0);
  const double s = static_cast<double>(L) / is;
  const double hs = s * 0.5;
  const double ds = s * 0.375;
  for (int iz = 0; iz < is; iz++) {
    for (int iy = 0; iy < is; iy++) {
      for (int ix = 0; ix < is; ix++) {
        const double bx = ix * s + ds;
        const double by = iy * s + ds;
        const double bz = iz * s + ds;
        adder.add(bx, by, bz, vars);
        adder.add(bx + hs, by, bz, vars);
        adder.add(bx, by + hs, bz, vars);
        adder.add(bx, by, bz + hs, vars);
      }
    }
  }
}
//------------------------------------------------------------------------
void
ConfMaker::make_droplet(double radius, double density, std::shared_ptr<Variables> &vars) {
  class DropletAdder : public AtomAdder {
  private:
    const double radius;
  public:
    DropletAdder(const double r) : radius(r) {}
    void add(double x, double y, double z, std::shared_ptr<Variables> &vars) {
      const double r2 = radius * radius;
      const double LH = L * 0.5;
      if ((x - LH) * (x - LH) + (y - LH) * (y - LH) + (z - LH) * (z - LH) < r2) {
        Atom a(x, y, z);
        vars->add_atom(a);
      }
    }
  };
  DropletAdder adder(radius);
  make_FCC(density, adder, vars);

}
//------------------------------------------------------------------------
void
ConfMaker::make_binarydroplet(double radius, double density, std::shared_ptr<Variables> &vars) {
  class BinaryDropletAdder : public AtomAdder {
  private:
    const double radius;
  public:
    BinaryDropletAdder(const double r) : radius(r) {}
    void add(double x, double y, double z, std::shared_ptr<Variables> &vars) {
      const double r2 = radius * radius;
      Atom a(x, y, z);
      const double LH = L * 0.5;
      if ((x - LH) * (x - LH) + (y - LH) * (y - LH) + (z - LH) * (z - LH) < r2) {
        a.type = 0;
      } else {
        a.type = 1;
      }
      vars->add_atom(a);
    }
  };
  BinaryDropletAdder adder(radius);
  make_FCC(density, adder, vars);
}
//------------------------------------------------------------------------
void
ConfMaker::make_binaryflat(double density, std::shared_ptr<Variables> &vars) {
  class BinaryFlatAdder : public AtomAdder {
  public:
    void add(double x, double y, double z, std::shared_ptr<Variables> &vars) {
      Atom a(x, y, z);
      const double LH = L * 0.5;
      if (y < LH) {
        a.type = 0;
      } else {
        a.type = 1;
      }
      vars->add_atom(a);
    }
  };
  BinaryFlatAdder adder;
  make_FCC(density, adder, vars);
}
//------------------------------------------------------------------------
void
ConfMaker::make_monoflat(double leftdensity, double rightdensity, std::shared_ptr<Variables> &vars) {
  class MonoFlatAdder : public AtomAdder {
  public:
    bool left;
    void add(double x, double y, double z, std::shared_ptr<Variables> &vars) {
      Atom a(x, y, z);
      const double LH = L * 0.5;
      if ((left && y < LH) || (!left && y >= LH)) {
        vars->add_atom(a);
      }
    }
  };
  MonoFlatAdder adder;
  adder.left = true;
  make_FCC(leftdensity, adder, vars);
  adder.left = false;
  make_FCC(rightdensity, adder, vars);
}
//------------------------------------------------------------------------
