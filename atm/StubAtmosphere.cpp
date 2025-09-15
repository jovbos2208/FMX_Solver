#include "atm/Atmosphere.hpp"
#include "core/units.hpp"
#include <algorithm>

namespace fmx::atm {

static inline double clamp(double x, double lo, double hi) { return std::max(lo, std::min(hi, x)); }

AtmosphereState StubAtmosphere::evaluate(double alt_km, double, double,
                                         const std::string&, const Indices& idx) const {
  AtmosphereState st{};
  // Simple parametric temperature model: 700–1000 K in thermosphere
  double baseT = 700.0 + 3.0 * clamp(idx.F10_7 - 70.0, 0.0, 200.0);
  st.T_K = clamp(baseT, 600.0, 1200.0);
  st.wind_ms = {0.0, 0.0, 0.0};

  // Crude species mass densities vs altitude: O dominates ~200–500 km, with
  // exponential decay; small fractions of He and H at high altitudes.
  double h = clamp(alt_km, 100.0, 600.0);
  double rho_O  = 5e-11 * std::exp(-(h - 200.0)/60.0); // kg/m^3
  double rho_He = 2e-12 * std::exp(-(h - 300.0)/80.0);
  double rho_H  = 1e-12 * std::exp(-(h - 400.0)/120.0);
  if (h < 180.0) { rho_O *= 2.0; }
  if (idx.Kp >= 5) { rho_O *= 1.5; }

  st.species.push_back({rho_O,  fmx::units::m_O});
  st.species.push_back({rho_He, fmx::units::m_He});
  st.species.push_back({rho_H,  fmx::units::m_H});
  return st;
}

} // namespace fmx::atm

