//------------------------------------------------------------------------
#include <iostream>
#include <string>
#include "parameter.hpp"
#include "md.hpp"
//------------------------------------------------------------------------
int
main(int argc, char **argv) {
  std::string filename;
  if (argc > 1) {
    filename = argv[1];
    std::cout << "# Input file is " << filename << std::endl;
  } else {
    std::cout << "# Input file is not specified. input.cfg is used." << std::endl;
    filename = "input.cfg";
  }
  Parameter param(filename.c_str());
  MD *md = new MD();
  if (param.is_valid()) {
    md->run(param);
  }
  delete md;
}
//------------------------------------------------------------------------
