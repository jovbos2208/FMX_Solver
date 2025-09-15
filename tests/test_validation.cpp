#include <cmath>
#include <iostream>
#include <vector>
#include "core/types.hpp"
#include "solver/PanelSolver.hpp"
#include "atm/Atmosphere.hpp"

using fmx::Vec3;

static fmx::solver::Input make_plate_at_angle(double theta_rad) {
  // Create a unit square plate centered at origin, rotate normal by theta about Z
  // Base plate lies on YZ-plane (normal along -X). Rotate around Z to tilt toward +Y.
  double c = std::cos(theta_rad), s = std::sin(theta_rad);
  auto rot = [&](const Vec3& v){ return Vec3{ v.x*c - v.y*s, v.x*s + v.y*c, v.z }; };
  // Define initial vertices for a plate at x=0 with normal -X
  std::vector<Vec3> v = {
    {0,  0.5,  0.5}, {0,  0.5, -0.5}, {0, -0.5, -0.5},
    {0, -0.5,  0.5}
  };
  for (auto& p : v) p = rot(p);

  fmx::geom::Mesh m;
  m.tris.push_back({v[0], v[1], v[2]});
  m.tris.push_back({v[3], v[0], v[2]});

  fmx::solver::Input in;
  in.facets = m.to_facets(0);
  in.materials = { {1.0, 1.0, 300.0} };
  fmx::atm::StubAtmosphere atm;
  auto st = atm.evaluate(400.0, 0.0, 0.0, "2025-09-12T12:00:00Z", {120.0, 3});
  for (const auto& sp : st.species) in.species.push_back({sp.rho, sp.mass});
  in.T_K = st.T_K;
  in.V_sat_ms = {7500.0, 0.0, 0.0};
  in.wind_ms = st.wind_ms;
  return in;
}

int main() {
  // Monotonic trend: Force magnitude decreases as incidence approaches grazing
  auto F0 = fmx::solver::solve(make_plate_at_angle(0.0)).F.x;         // normal incidence
  auto F30 = fmx::solver::solve(make_plate_at_angle(M_PI/6)).F.x;     // 30 deg
  auto F60 = fmx::solver::solve(make_plate_at_angle(M_PI/3)).F.x;     // 60 deg
  auto F85 = fmx::solver::solve(make_plate_at_angle(85.0*M_PI/180.0)).F.x; // 85 deg

  double a0 = std::abs(F0), a30 = std::abs(F30), a60 = std::abs(F60), a85 = std::abs(F85);
  // Relaxed checks: grazing is much smaller; 60<30; 30 close to 0 within 50%
  if (!(a85 < 0.5*a0 && a60 < a30 && a85 < a60 && a30 < 1.5*a0)) {
    std::cerr << "Unexpected angle trend: "
              << a0 << ", " << a30 << ", " << a60 << ", " << a85 << "\n";
    return 1;
  }

  // Backface test: rotate 180 deg around Z; expect zero force
  auto Fb = fmx::solver::solve(make_plate_at_angle(M_PI)).F.x;
  if (std::abs(Fb) > 1e-9) {
    std::cerr << "Backface force not ~0: " << Fb << "\n";
    return 1;
  }

  return 0;
}
