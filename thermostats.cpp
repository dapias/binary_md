//------------------------------------------------------------------------
#include <random>
#include "systemparam.hpp"
#include "md.hpp"
//------------------------------------------------------------------------
void
VelocityScaling::update(MD *md){
  std::shared_ptr<Variables> vars = md->get_variables();
  double t = md->get_observer()->temperature(vars);
  double ratio = sqrt(md->get_aimed_temperature() / t);
  for (auto &a : vars->atoms) {
    a.px *= ratio;
    a.py *= ratio;
    a.pz *= ratio;
  }
}
//------------------------------------------------------------------------
void
Langevin::update(MD *md){
	static std::mt19937 mt(1);
  const double gamma = 1.0;
  const double T = md->get_aimed_temperature();
  const double D = sqrt(2.0 * gamma * T / dt);
  std::normal_distribution<double> nd(0.0, D);
  std::shared_ptr<Variables> vars = md->get_variables();
  for (auto &a : vars->atoms) {
    a.px += (-gamma * a.px + nd(mt)) * dt;
    a.py += (-gamma * a.py + nd(mt)) * dt;
    a.pz += (-gamma * a.pz + nd(mt)) * dt;
  }
}
//------------------------------------------------------------------------
void
NoseHoover::update(MD *md){
  std::shared_ptr<Variables> vars = md->get_variables();
	std::shared_ptr<Observer> obs = md->get_observer();
  double t = obs->temperature(vars);
  double at = md->get_aimed_temperature();
  double Q = 0.01;
  vars->zeta += (t - at) * dt;
	const double zq = vars->zeta / Q;
  for (auto &a : vars->atoms) {
    a.px -= a.px * zq * dt;
    a.py -= a.py * zq * dt;
    a.pz -= a.pz * zq * dt;
  }
}
//------------------------------------------------------------------------
void
NoseHooverLangevin::update(MD *md){
	static std::mt19937 mt(1);
  double Q = 0.01;
  const double gamma = 1.0;
  const double T = md->get_aimed_temperature();
  const double D = sqrt(2.0 * gamma * T * Q/ dt);
	std::normal_distribution<double> nd(0.0,D);
  std::shared_ptr<Variables> vars = md->get_variables();
	std::shared_ptr<Observer> obs = md->get_observer();
  double t = obs->temperature(vars);
  double at = md->get_aimed_temperature();
  vars->zeta += (t - at - gamma*vars->zeta + nd(mt)) * dt;
	const double zq = vars->zeta / Q;
  for (auto &a : vars->atoms) {
    a.px -= a.px * zq * dt;
    a.py -= a.py * zq * dt;
    a.pz -= a.pz * zq * dt;
  }
}
//------------------------------------------------------------------------
void
KineticMoments::update(MD *md){
  std::shared_ptr<Variables> vars = md->get_variables();
  double at = md->get_aimed_temperature();
	const double Qzeta = 0.01;
	const double Qeta = 1.0;
  double p4 = 0.0;
  double p2 = 0.0;
  const int N = vars->number_of_atoms();
  for (auto &a : vars->atoms) {
    p4 += a.px * a.px * a.px * a.px;
    p4 += a.py * a.py * a.py * a.py;
    p4 += a.pz * a.pz * a.pz * a.pz;
    p2 += a.px * a.px;
    p2 += a.py * a.py;
    p2 += a.pz * a.pz;
  }
  p4 /= (3.0 * static_cast<double>(N));
  p2 /= (3.0 * static_cast<double>(N));
  vars->zeta += (p2 - at) * dt;
  vars->eta += (p4 - 3.0 * at * p2) * dt;
	const double zq = vars->zeta / Qzeta;
	const double eq = vars->eta / Qeta;
  for (auto &a : vars->atoms) {
    a.px -= (a.px * zq + a.px * a.px * a.px * eq) * dt;
    a.py -= (a.py * zq + a.py * a.py * a.py * eq) * dt;
    a.pz -= (a.pz * zq + a.pz * a.pz * a.pz * eq) * dt;
  }
}
//------------------------------------------------------------------------

