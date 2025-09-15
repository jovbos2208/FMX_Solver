// Cercignani–Lampis–Lord (CLL) gas–surface interaction kernel (interface)
#pragma once

#include <tuple>

namespace fmx::gsi {

struct CLLParams {
  double alpha_n{1.0}; // normal energy accommodation (0..1)
  double alpha_t{1.0}; // tangential momentum accommodation (0..1)
};

// Compute dimensionless normal/tangential coefficients (C_N, C_T) for CLL
// theta: angle between -c_hat and panel normal [rad]
// Ma: species Mach number |c|/sqrt(kT/m)
// tau: Tw/T
// Returns (C_N, C_T)
std::tuple<double,double> coefficients(double theta, double Ma, double tau,
                                       const CLLParams& p);

} // namespace fmx::gsi

