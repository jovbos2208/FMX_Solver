#pragma once

#include "core/types.hpp"
#include <string>

namespace fmx::atm {

// Thin wind-only wrapper for HWM14.
struct HWM14 {
  // Returns neutral wind vector [m/s] at given state.
  static fmx::Vec3 evaluate_wind(double alt_km, double lat_deg, double lon_deg,
                                 const std::string& utc_iso, double ap3hr);
};

} // namespace fmx::atm
