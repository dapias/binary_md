//------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <assert.h>
#include <sstream>
#include "systemparam.hpp"
#include "variables.hpp"
//------------------------------------------------------------------------
void
Variables::set_initial_velocity(const double V0) {
  std::mt19937 mt(2);
  std::uniform_real_distribution<double> ud(0.0, 1.0);
  double avx = 0.0;
  double avy = 0.0;
  double avz = 0.0;
  for (auto &a : atoms) {
    double z = ud(mt) * 2.0 - 1.0;
    double phi = 2.0 * ud(mt) * M_PI;
    double vx = V0 * sqrt(1 - z * z) * cos(phi);
    double vy = V0 * sqrt(1 - z * z) * sin(phi);
    double vz = V0 * z;
    a.px = vx;
    a.py = vy;
    a.pz = vz;
    avx += vx;
    avy += vy;
    avz += vz;
  }
  const int pn = atoms.size();
  avx /= static_cast<double>(pn);
  avy /= static_cast<double>(pn);
  avz /= static_cast<double>(pn);
  for (auto &a : atoms) {
    a.px -= avx;
    a.py -= avy;
    a.pz -= avz;
  }
}
//------------------------------------------------------------------------
void
Variables::set_initial_temperature(const double T0) {
  set_initial_velocity(sqrt(3.0 * T0));
}
//------------------------------------------------------------------------
void
Variables::save_as_cdview(const char *filename) {
  std::ofstream ofs(filename);
  int i = 0;
  for (auto &a : atoms) {
    ofs << i << " " << a.type << " ";
    ofs << a.qx << " ";
    ofs << a.qy << " ";
    ofs << a.qz << " ";
    ofs << a.px << " ";
    ofs << a.py << " ";
    ofs << a.pz << " ";
    ofs << std::endl;
    ++i;
  }
}
//------------------------------------------------------------------------
void
Variables::export_cdview(void) {
  static int count = 0;
  char filename[256];
  sprintf(filename, "conf%04d.cdv", count);
  ++count;
  save_as_cdview(filename);
}
//------------------------------------------------------------------------
void
Variables::import_cdview(const char *filename) {
  std::ifstream ifs(filename);
  if (!ifs) {
    std::cerr << "Cannot open file " << filename << "." << std::endl;
    return;
  }
  std::string str;
  while (getline(ifs, str)) {
    std::istringstream is(str);
    std::string str2;
    std::vector<double> v;
    getline(is, str2, ' ');
    getline(is, str2, ' ');
    Atom a;
    a.type = std::stoi(str2);
    getline(is, str2, ' ');
    a.qx = std::stof(str2);
    getline(is, str2, ' ');
    a.qy = std::stof(str2);
    getline(is, str2, ' ');
    a.qz = std::stof(str2);

    getline(is, str2, ' ');
    a.px = std::stof(str2);
    getline(is, str2, ' ');
    a.py = std::stof(str2);
    getline(is, str2, ' ');
    a.pz = std::stof(str2);
    add_atom(a);
  };
}
//------------------------------------------------------------------------
void
Variables::make_neighbor_list(std::vector<Pair> &pairs) {
  const int pn = atoms.size();
  const int pp = pairs.size();
  neighbor_list.clear();
  neighbor_list.resize(pp);
  i_position.resize(pn);
  j_count.resize(pn);
  std::fill(j_count.begin(), j_count.end(), 0);
  for (auto &p : pairs) {
    j_count[p.i]++;
  }
  i_position[0] = 0;
  int sum = 0;
  for (int i = 0; i < pn - 1; i++) {
    sum += j_count[i];
    i_position[i + 1] = sum;
  }
  std::vector<int> pointer(pn);
  std::fill(pointer.begin(), pointer.end(), 0);
  for (auto &p : pairs) {
    int pos = i_position[p.i] + pointer[p.i];
    neighbor_list[pos] = p.j;
    pointer[p.i]++;
  }
}
//------------------------------------------------------------------------
