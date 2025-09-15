// Sentman free-molecular coefficients (baseline stub)
#pragma once

#include <tuple>

namespace fmx::gsi {

struct SentmanParams {
  // Energy accommodation coefficient (alpha_E)
  double alpha_E{1.0};
};

// Compute dimensionless normal/tangential coefficients (C_N, C_T)
// theta: angle between -c_hat and panel normal [rad]
// Ma: species Mach number |c|/sqrt(kT/m)
// tau: Tw/T
// Returns (C_N, C_T)
std::tuple<double,double> coefficients(double theta, double Ma, double tau,
                                       const SentmanParams& p = {});

} // namespace fmx::gsi
