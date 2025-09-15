#include "gsi/CLL.hpp"
#include "gsi/Sentman.hpp"
#include <algorithm>
#include <cmath>

namespace fmx::gsi {

// Placeholder implementation: fall back to Sentman coefficients and apply a simple
// scaling to approximate reduced tangential slip for alpha_t<1.0.
// This is a structural stub; a proper CLL evaluation will be added next (tabulated or closed form).
std::tuple<double,double> coefficients(double theta, double Ma, double tau,
                                       const CLLParams& p) {
  const auto [CN_s, CT_s] = fmx::gsi::coefficients(theta, Ma, tau, SentmanParams{1.0});
  // Heuristic: normal roughly follows Sentman; tangential reduced by (alpha_t)
  double CN = CN_s * (0.5 + 0.5*std::clamp(p.alpha_n, 0.0, 1.0));
  double CT = CT_s * std::clamp(p.alpha_t, 0.0, 1.0);
  return {CN, CT};
}

} // namespace fmx::gsi

