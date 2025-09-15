#include <cmath>
#include <iostream>
#include "core/types.hpp"
#include "geom/Mesh.hpp"
#include "solver/PanelSolver.hpp"
#include "atm/Atmosphere.hpp"

using fmx::Vec3;

static fmx::geom::Mesh make_plate() {
  fmx::geom::Mesh m;
  // YZ plane, normal toward -X (front-facing for +X flow)
  m.tris.push_back({Vec3{0, 0.5, 0.5}, Vec3{0, 0.5, -0.5}, Vec3{0, -0.5, -0.5}});
  m.tris.push_back({Vec3{0, -0.5, 0.5}, Vec3{0, 0.5, 0.5}, Vec3{0, -0.5, -0.5}});
  return m;
}

int main() {
  auto m = make_plate();
  auto facets = m.to_facets(0);

  // Stub atmosphere: deterministic, no wind
  fmx::atm::StubAtmosphere atm;
  auto st = atm.evaluate(400.0, 0.0, 0.0, "2025-09-12T12:00:00Z", {120.0, 3});

  fmx::solver::Input in;
  in.facets = facets;
  in.materials = { {1.0, 1.0, 1.0, 300.0} };
  for (const auto& sp : st.species) in.species.push_back({sp.rho, sp.mass});
  in.T_K = st.T_K;
  in.V_sat_ms = {7500.0, 0.0, 0.0};
  in.wind_ms = {0.0, 0.0, 0.0};
  // Offset CG by +Y to induce torque about +Z (Mz = r_y * Fx)
  const double ry = 0.1;
  in.r_CG = {0.0, ry, 0.0};

  auto out = fmx::solver::solve(in);
  double Fx = out.F.x;
  double Mz = out.M.z;

  if (!(std::abs(out.F.y) < 1e-9 && std::abs(out.F.z) < 1e-9)) {
    std::cerr << "Unexpected lateral/vertical force: Fy=" << out.F.y << ", Fz=" << out.F.z << "\n";
    return 1;
  }
  // Expect Mz ≈ r_y * (-Fx) sign? r x F: (r_x F_y - r_y F_x) = - r_y * F_x, but with F_x negative, Mz positive.
  double expected_Mz = ry * Fx;
  double rel_err = std::abs(Mz - expected_Mz) / std::max(1e-16, std::abs(expected_Mz));
  if (rel_err > 5e-3) {
    std::cerr << "Torque mismatch: Mz=" << Mz << " expected=" << expected_Mz << " rel_err=" << rel_err << "\n";
    return 1;
  }

  return 0;
}
