#include <cmath>
#include <iostream>
#include "gsi/Sentman.hpp"

int main() {
  using fmx::gsi::coefficients;
  using fmx::gsi::SentmanParams;

  // Normal incidence: theta=0 -> no tangential shear, positive normal load
  {
    auto [CN, CT] = coefficients(0.0, 8.0, 1.0, SentmanParams{});
    if (std::abs(CT) > 1e-6) {
      std::cerr << "CT not ~0 at normal incidence: " << CT << "\n";
      return 1;
    }
    if (!(CN > 0.0)) {
      std::cerr << "CN should be >0 at normal incidence: " << CN << "\n";
      return 1;
    }
  }

  // Grazing incidence: theta=pi/2 -> CN ~ 0; CT decreases with Ma
  {
    auto [CN8, CT8] = coefficients(0.5*M_PI, 8.0, 1.0, SentmanParams{});
    auto [CN16, CT16] = coefficients(0.5*M_PI, 16.0, 1.0, SentmanParams{});
    if (std::abs(CN8) > 5e-3) {
      std::cerr << "CN not ~0 at grazing: CN8=" << CN8 << "\n";
      return 1;
    }
    if (!(std::abs(CT16) < std::abs(CT8))) {
      std::cerr << "CT did not decrease with Ma at grazing: CT8=" << CT8 << " CT16=" << CT16 << "\n";
      return 1;
    }
  }

  // Low speed ratio: coefficients remain finite
  {
    auto [CN, CT] = coefficients(0.3, 0.05, 1.0, SentmanParams{});
    if (!std::isfinite(CN) || !std::isfinite(CT)) {
      std::cerr << "Coeffs not finite at low Ma\n";
      return 1;
    }
  }

  return 0;
}
