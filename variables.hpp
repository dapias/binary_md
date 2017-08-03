#pragma once
#include <vector>
//------------------------------------------------------------------------
struct Atom {
  double qx, qy, qz;
  double px, py, pz;
  int type;
  Atom(void) {
    px = 0.0;
    py = 0.0;
    pz = 0.0;
    type = 0;
  }
  Atom(double x, double y, double z) {
    qx = x;
    qy = y;
    qz = z;
    px = 0.0;
    py = 0.0;
    pz = 0.0;
    type = 0;
  }
  Atom(int t, double x, double y, double z, double vx, double vy, double vz) {
    type = t;
    qx = x;
    qy = y;
    qz = z;
    px = vx;
    py = vy;
    pz = vz;
  }
};
//------------------------------------------------------------------------
struct Pair {
  int i, j;
};
//------------------------------------------------------------------------
class Variables {
public:
  std::vector<Atom> atoms;
  std::vector<int> neighbor_list;
  std::vector<int> i_position;
  std::vector<int> j_count;
  double time;
  double zeta; //For Nose-Hoover
  double eta; //For Kinetic-Moments
  Variables(void) {time = 0.0; zeta = 0.0;}
  void add_atom(Atom &a) {atoms.push_back(a);}
  void save_as_cdview(const char *filename);
  void export_cdview(void);
  void import_cdview(const char *filename);
  int number_of_atoms(void) {return static_cast<int>(atoms.size());}
  void set_initial_velocity(const double);
  void set_initial_temperature(const double);
  void make_neighbor_list(std::vector<Pair> &pair);
};
//------------------------------------------------------------------------
