#include "solver/RegimeAdapter.hpp"
#include "core/units.hpp"
#include "solver/PanelSolver.hpp"
#include <algorithm>
#include <cmath>

namespace fmx::solver {

static double eff_diameter_m(const std::vector<Species>& species) {
  // Approximate effective collision diameter by mole-fraction weighted average of hard-sphere diameters.
  // Typical kinetic diameters [m]: O: 3.0e-10, N2: 3.7e-10, O2: 3.5e-10, He: 2.6e-10, H: 2.2e-10
  struct SD { double m; double d; } refs[] = {
    { units::m_O,  3.0e-10 },
    { units::m_N2, 3.7e-10 },
    { units::m_O2, 3.5e-10 },
    { units::m_He, 2.6e-10 },
    { units::m_H,  2.2e-10 },
  };
  double n_tot = 0.0;
  double sum = 0.0;
  for (const auto& sp : species) {
    // Convert mass density to number density assuming uniform T not needed here; we only need mole fractions by mass
    // Use mass weighting fallback: weight by rho/m
    for (const auto& r : refs) {
      if (std::abs(sp.mass - r.m) < 1e-30) {
        double n = (sp.rho > 0.0 && sp.mass > 0.0) ? (sp.rho / sp.mass) : 0.0;
        sum += n * r.d;
        n_tot += n;
        break;
      }
    }
  }
  if (n_tot <= 0.0) return 3.5e-10; // default
  return sum / n_tot;
}

RegimeDiagnostics apply_regime_blend(const fmx::atm::AtmosphereState& st,
                                     const std::vector<Species>& species,
                                     double T_K,
                                     const RegimeConfig& cfg,
                                     Output& io) {
  RegimeDiagnostics diag{};
  if (!cfg.enabled) return diag;

  // Estimate total number density n and effective diameter d
  double n = 0.0; // [1/m^3]
  for (const auto& sp : species) {
    if (sp.mass > 0.0) n += sp.rho / sp.mass;
  }
  double d = eff_diameter_m(species);
  // Mean free path lambda = 1 / (sqrt(2) * pi * d^2 * n)
  double lambda = (n > 0.0 && d > 0.0) ? (1.0 / (std::sqrt(2.0) * M_PI * d * d * n)) : 1e9;
  double L = std::max(1e-6, cfg.L_char_m);
  double Kn = lambda / L;
  diag.Kn = Kn;
  double beta = std::exp(-cfg.gamma * Kn);
  diag.beta = beta;

  // Surrogate correction: scalar scaling of FM force: s(Kn) = 1/(1 + a*Kn^b)
  // Blend: F = (1-beta)*F_FM + beta*s*F_FM = ((1-beta) + beta*s) * F_FM
  double s = 1.0;
  if (cfg.corr_mode == RegimeConfig::CorrMode::Scalar) {
    double a = std::max(0.0, cfg.corr_a);
    double b = std::max(0.0, cfg.corr_b);
    s = 1.0 / (1.0 + a * std::pow(std::max(1e-12, Kn), b));
  }
  double eff = (1.0 - beta) + beta * s;
  io.F *= eff;
  io.M *= eff;
  (void)st; (void)T_K;
  return diag;
}

} // namespace fmx::solver
