#pragma once
class MD;
//------------------------------------------------------------------------
class Thermostat {
  public:
    virtual void update(MD *) = 0; 
};
//------------------------------------------------------------------------
class None : public Thermostat {
  public:
    void update(MD *){}; 
};
//------------------------------------------------------------------------
class VelocityScaling : public Thermostat {
  public:
    void update(MD *); 
};
//------------------------------------------------------------------------
class Langevin : public Thermostat {
  public:
    void update(MD *md); 
};
//------------------------------------------------------------------------
class NoseHoover : public Thermostat {
  public:
    void update(MD *md); 
};
//------------------------------------------------------------------------
class NoseHooverLangevin : public Thermostat {
  public:
    void update(MD *md); 
};
//------------------------------------------------------------------------
class KineticMoments : public Thermostat {
  public:
    void update(MD *md); 
};
//------------------------------------------------------------------------
