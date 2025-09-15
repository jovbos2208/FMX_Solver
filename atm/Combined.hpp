#pragma once

#include "atm/Atmosphere.hpp"
#include "atm/NRLMSIS2.hpp"
#include "atm/HWM14.hpp"

namespace fmx::atm {

// Combines NRLMSIS2 (T, species) with HWM14 (wind)
class CombinedAtmosphere : public Atmosphere {
public:
  AtmosphereState evaluate(double alt_km, double lat_deg, double lon_deg,
                           const std::string& utc_iso,
                           const Indices& idx) const override {
    NRLMSIS2Atmosphere msis;
    AtmosphereState st = msis.evaluate(alt_km, lat_deg, lon_deg, utc_iso, idx);
    // Choose ap3hr for wind
    double ap3hr = 0.0;
    if (idx.has_ap_array) ap3hr = idx.Ap[1];
    else if (idx.Ap_now > 0.0) ap3hr = idx.Ap_now;
    else ap3hr = 10.0; // fallback
    st.wind_ms = HWM14::evaluate_wind(alt_km, lat_deg, lon_deg, utc_iso, ap3hr);
    return st;
  }
};

} // namespace fmx::atm
