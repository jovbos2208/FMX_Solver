#include <cmath>
#include <iostream>
#include "core/types.hpp"
#include "geom/Mesh.hpp"
#include "solver/PanelSolver.hpp"
#include "atm/Atmosphere.hpp"

using fmx::Vec3;

static fmx::geom::Mesh make_plate() {
  fmx::geom::Mesh m;
  m.tris.push_back({Vec3{0, 0.5, 0.5}, Vec3{0, 0.5, -0.5}, Vec3{0, -0.5, -0.5}});
  m.tris.push_back({Vec3{0, -0.5, 0.5}, Vec3{0, 0.5, 0.5}, Vec3{0, -0.5, -0.5}});
  return m;
}

static fmx::geom::Mesh make_cube(double a=1.0) {
  fmx::geom::Mesh m;
  double h=a*0.5;
  auto add_quad = [&](Vec3 a, Vec3 b, Vec3 c, Vec3 d){ m.tris.push_back({a,b,c}); m.tris.push_back({d,a,c}); };
  // +X
  add_quad({ h,-h,-h},{ h, h,-h},{ h, h, h},{ h,-h, h});
  // -X (front-facing for +X flow)
  add_quad({-h, h, h},{-h, h,-h},{-h,-h,-h},{-h,-h, h});
  // +Y
  add_quad({-h, h, h},{ h, h, h},{ h, h,-h},{-h, h,-h});
  // -Y
  add_quad({-h,-h,-h},{ h,-h,-h},{ h,-h, h},{-h,-h, h});
  // +Z
  add_quad({-h,-h, h},{ h,-h, h},{ h, h, h},{-h, h, h});
  // -Z
  add_quad({-h, h,-h},{ h, h,-h},{ h,-h,-h},{-h,-h,-h});
  return m;
}

static fmx::solver::Output solve_mesh(const fmx::geom::Mesh& m) {
  auto facets = m.to_facets(0);
  fmx::atm::StubAtmosphere atm;
  auto st = atm.evaluate(400.0, 0.0, 0.0, "2025-09-12T12:00:00Z", {120.0, 3});
  fmx::solver::Input in;
  in.facets = facets;
  in.materials = { {1.0, 1.0, 1.0, 300.0} };
  for (const auto& sp : st.species) in.species.push_back({sp.rho, sp.mass});
  in.T_K = st.T_K;
  in.V_sat_ms = {7500.0, 0.0, 0.0};
  in.wind_ms = {0.0, 0.0, 0.0};
  in.r_CG = {0.0, 0.0, 0.0};
  return fmx::solver::solve(in);
}

int main() {
  auto plate_out = solve_mesh(make_plate());
  auto cube_out  = solve_mesh(make_cube());

  // Symmetry checks: lateral/vertical forces and all moments should be ~0
  auto near0 = [](double v){ return std::abs(v) < 1e-9; };
  if (!near0(plate_out.F.y) || !near0(plate_out.F.z)) {
    std::cerr << "Plate Fy/Fz not ~0: " << plate_out.F.y << ", " << plate_out.F.z << "\n";
    return 1;
  }
  if (!near0(cube_out.F.y) || !near0(cube_out.F.z)) {
    std::cerr << "Cube Fy/Fz not ~0: " << cube_out.F.y << ", " << cube_out.F.z << "\n";
    return 1;
  }
  if (!(near0(cube_out.M.x) && near0(cube_out.M.y) && near0(cube_out.M.z))) {
    std::cerr << "Cube torque not ~0: M=[" << cube_out.M.x << ", " << cube_out.M.y << ", " << cube_out.M.z << "]\n";
    return 1;
  }

  // Direction: drag dominates in -X direction
  if (!(plate_out.F.x < 0.0 && cube_out.F.x < 0.0)) {
    std::cerr << "Expected negative Fx drag for plate/cube\n";
    return 1;
  }
  // Expect cube drag magnitude >= plate drag (grazing side faces add or equal drag)
  if (!(std::abs(cube_out.F.x) + 1e-12 >= std::abs(plate_out.F.x))) {
    std::cerr << "Cube drag not larger than plate: |Fx_cube|=" << std::abs(cube_out.F.x)
              << " |Fx_plate|=" << std::abs(plate_out.F.x) << "\n";
    return 1;
  }

  return 0;
}
