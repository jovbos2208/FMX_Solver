#pragma once

#include "atm/Atmosphere.hpp"

namespace fmx::atm {

// Thin wrapper around NRLMSIS 2.0. If FMX_WITH_NRLMSIS2 is not defined,
// falls back to a plausible stub and warns once.
class NRLMSIS2Atmosphere : public Atmosphere {
public:
  AtmosphereState evaluate(double alt_km, double lat_deg, double lon_deg,
                           const std::string& utc_iso,
                           const Indices& idx) const override;
};

} // namespace fmx::atm

