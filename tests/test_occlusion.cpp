#include <cmath>
#include <iostream>
#include <vector>
#include "core/types.hpp"
#include "geom/Mesh.hpp"
#include "geom/BVH.hpp"
#include "atm/Atmosphere.hpp"
#include "solver/PanelSolver.hpp"

using fmx::Vec3;

static fmx::solver::Input make_input_with_mesh(const fmx::geom::Mesh& m) {
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
  // Single front-facing plate (YZ plane at x=0, normal toward -X)
  fmx::geom::Mesh m1;
  m1.tris.push_back({Vec3{0, 0.5, 0.5}, Vec3{0, 0.5, -0.5}, Vec3{0, -0.5, -0.5}});
  m1.tris.push_back({Vec3{0, -0.5, 0.5}, Vec3{0, 0.5, 0.5}, Vec3{0, -0.5, -0.5}});

  // Two inline plates: same as m1 plus another at x=+0.2
  fmx::geom::Mesh m2 = m1;
  m2.tris.push_back({Vec3{0.2, 0.5, 0.5}, Vec3{0.2, 0.5, -0.5}, Vec3{0.2, -0.5, -0.5}});
  m2.tris.push_back({Vec3{0.2, -0.5, 0.5}, Vec3{0.2, 0.5, 0.5}, Vec3{0.2, -0.5, -0.5}});

  // Compute forces
  auto in1 = make_input_with_mesh(m1);
  fmx::geom::BVHOccluder occ1(m1.tris);
  in1.occluder = &occ1;
  auto out1 = fmx::solver::solve(in1);

  auto in2_noocc = make_input_with_mesh(m2);
  in2_noocc.occluder = nullptr; // no occlusion
  auto out2_noocc = fmx::solver::solve(in2_noocc);

  auto in2_occ = make_input_with_mesh(m2);
  fmx::geom::BVHOccluder occ2(m2.tris);
  in2_occ.occluder = &occ2;
  auto out2_occ = fmx::solver::solve(in2_occ);

  // Expectations:
  // - No occlusion: two plates contribute ~2x single plate force (same orientation/area)
  // - With occlusion: rear plate shadowed -> total ≈ single plate
  const double F1 = std::abs(out1.F.x);
  const double F2_no = std::abs(out2_noocc.F.x);
  const double F2_occ = std::abs(out2_occ.F.x);

  if (!(F2_no > 1.8 * F1)) {
    std::cerr << "Two-plate force without occlusion not ~2x single (F2_no=" << F2_no << ", F1=" << F1 << ")\n";
    return 1;
  }
  if (!(std::abs(F2_occ - F1) / F1 < 0.05)) {
    std::cerr << "Occluded two-plate force not ~single (F2_occ=" << F2_occ << ", F1=" << F1 << ")\n";
    return 1;
  }

  return 0;
}

