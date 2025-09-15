#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include "core/types.hpp"

namespace fmx::atm {

struct Indices {
  double F10_7{120.0};
  double F10_7A{120.0};
  int Kp{3};
  bool has_ap_array{false};
  double Ap[7] = {0,0,0,0,0,0,0}; // daily + six 3hr values as per MSIS when legacy switch used
  double Ap_daily{0.0};
  double Ap_now{0.0};
};

struct SpeciesState {
  double rho;   // mass density [kg/m^3]
  double mass;  // molecular mass [kg]
};

struct AtmosphereState {
  double T_K;                // temperature [K]
  fmx::Vec3 wind_ms;         // neutral wind [m/s]
  std::vector<SpeciesState> species; // species list
};

class Atmosphere {
public:
  virtual ~Atmosphere() = default;
  virtual AtmosphereState evaluate(double alt_km, double lat_deg, double lon_deg,
                                   const std::string& utc_iso,
                                   const Indices& idx) const = 0;
};

class StubAtmosphere : public Atmosphere {
public:
  AtmosphereState evaluate(double alt_km, double, double,
                           const std::string&, const Indices& idx) const override;
};

} // namespace fmx::atm
