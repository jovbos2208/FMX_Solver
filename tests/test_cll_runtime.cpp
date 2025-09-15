#include <iostream>
#include <cmath>
#include "gsi/CLLRuntime.hpp"

int main() {
  using fmx::gsi::CLLRuntime; using fmx::gsi::CLLQuantization;
  CLLQuantization q; q.gh_order = 16; q.theta_deg_step = 0.5; q.tau_step = 0.05; q.alpha_step = 0.1; q.workers = 1;
  CLLRuntime rt(q);
  // Reproducibility: same key twice yields same result
  auto r1 = rt.query(0.0, 8.0, 1.0, 1.0, 1.0);
  auto r2 = rt.query(0.0, 8.0, 1.0, 1.0, 1.0);
  if (std::abs(r1.first - r2.first) > 1e-12 || std::abs(r1.second - r2.second) > 1e-12) {
    std::cerr << "CLLRuntime not reproducible for identical keys\n"; return 1;
  }
  // Edge cases: normal incidence -> CT~0, CN>0
  auto rn = rt.query(0.0, 8.0, 1.0, 1.0, 1.0);
  if (!(std::abs(rn.second) < 1e-6 && rn.first > 0.0)) {
    std::cerr << "CLL CN/CT unexpected at normal incidence: CN="<<rn.first<<" CT="<<rn.second<<"\n"; return 1;
  }
  // Grazing trend: CT decreases with Ma at grazing angle
  double thg = 89.0 * M_PI/180.0;
  auto rg8  = rt.query(thg, 8.0, 1.0, 1.0, 1.0);
  auto rg16 = rt.query(thg,16.0, 1.0, 1.0, 1.0);
  if (!(std::abs(rg16.second) < std::abs(rg8.second))) {
    std::cerr << "CLL grazing CT did not decrease with Ma: CT8="<<rg8.second<<" CT16="<<rg16.second<<"\n"; return 1;
  }
  return 0;
}

