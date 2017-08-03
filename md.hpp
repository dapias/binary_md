#pragma once
#include <string>
#include <map>
#include <memory>
#include "parameter.hpp"
#include "variables.hpp"
#include "observer.hpp"
#include "meshlist.hpp"
#include "thermostats.hpp"
//------------------------------------------------------------------------
class MD {
private:
  typedef void(MD::*RUNFUNC)(Parameter &param);
  std::shared_ptr<Variables> vars;
  std::unique_ptr<MeshList> mesh;
  std::shared_ptr<Observer> obs;
  std::vector<Pair> pairs;
  std::map<std::string, RUNFUNC> runmap;
  std::shared_ptr<Thermostat> thermostat;
  std::map<std::string, std::shared_ptr<Thermostat>> thermostat_map;
  double margin_length;
  void makeconf(void);
  void update_position(void);
  void calculate_force(void);
  void calculate_force_pair(void);
  void calculate_force_list(void);
  void periodic(void);
  void calculate(void);
  void make_pair(void);
  void check_pairlist(void);
  bool init(Parameter &param);
  void showinfo(void);
  //For Simulation
  std::string hbtype;
  double aimed_temperature;
  int TotalSteps;
  int ObserveInterval;
  // Run Funtions
  void run_rt(Parameter &param);
  void run_mono_thermalize(Parameter &param);
  void run_mono_restart(Parameter &param);
  void run_binary(Parameter &param);

public:
  MD(void);
  std::shared_ptr<Observer> get_observer(void){return obs;}
  std::shared_ptr<Variables> get_variables(void){return vars;}
  double get_aimed_temperature(void){return aimed_temperature;};
  void run(Parameter &param);
};
//------------------------------------------------------------------------
