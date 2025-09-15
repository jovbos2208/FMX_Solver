#pragma once

#include <vector>
#include "core/types.hpp"
#include "atm/Atmosphere.hpp"

namespace fmx::solver {

struct RegimeConfig {
  bool enabled{false};
  double L_char_m{1.0};
  double gamma{1.0};
  // Correction surrogate options
  enum class CorrMode { None, Scalar, PerFacet };
  CorrMode corr_mode{CorrMode::None};
  double corr_a{0.0}; // scale = 1/(1 + a*Kn^b)
  double corr_b{1.0};
  // Per-facet angle-sensitive scaling of CN/CT
  double aN{0.2}, bN{1.0};
  double aT{0.8}, bT{1.0};
};

struct RegimeDiagnostics {
  double Kn{0.0};
  double beta{0.0};
};

// Estimate Kn and apply blend (placeholder: no change to forces yet; returns diagnostics)
struct Species; // fwd
struct Output;  // fwd

RegimeDiagnostics apply_regime_blend(const fmx::atm::AtmosphereState& st,
                                     const std::vector<Species>& species,
                                     double T_K,
                                     const RegimeConfig& cfg,
                                     Output& io);

} // namespace fmx::solver
