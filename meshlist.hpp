#pragma once
//------------------------------------------------------------------------
#include <vector>
#include <memory>
#include "variables.hpp"
//------------------------------------------------------------------------
class MeshList {
private:
  double mesh_size;
  int m;
  int number_of_mesh;
  bool initialized;
  std::vector<int> count;
  std::vector<int> indexes;
  std::vector<int> sorted_buffer;
  void search(int index, std::shared_ptr<Variables> &vars, std::vector<Pair> &pairs);
  void search_other(int id, int ix, int iy, int iz, std::shared_ptr<Variables> &vars, std::vector<Pair> &pairs);
public:
  MeshList(void);
  void init(void);
  void make_pair(std::shared_ptr<Variables> &vars, std::vector<Pair> &pairs);
};
//------------------------------------------------------------------------
