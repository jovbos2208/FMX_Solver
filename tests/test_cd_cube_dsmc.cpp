#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include "core/types.hpp"
#include "geom/Mesh.hpp"
#include "solver/PanelSolver.hpp"
#include "core/units.hpp"

using fmx::Vec3;

struct DBRow { double rho_tot; double n_N2, n_O2, n_O, n_He, n_H, n_Ar, n_N, n_AO, n_NO; double T_K; };

static bool read_db_row_by_line(const std::string& path, int target_line, DBRow& row) {
  std::ifstream in(path);
  if (!in) return false;
  std::string line;
  // header
  if (!std::getline(in, line)) return false;
  int line_no = 0;
  while (std::getline(in, line)) {
    if (line_no == target_line) {
      std::vector<std::string> cols; cols.reserve(16);
      size_t pos=0; while (pos<line.size()) {
        size_t comma = line.find(',', pos);
        if (comma==std::string::npos) comma = line.size();
        cols.push_back(line.substr(pos, comma-pos));
        pos = comma+1;
      }
      if (cols.size() < 13) return false;
      try {
        row.rho_tot = std::stod(cols[0]);
        row.n_N2 = std::stod(cols[1]);
        row.n_O2 = std::stod(cols[2]);
        row.n_O  = std::stod(cols[3]);
        row.n_He = std::stod(cols[4]);
        row.n_H  = std::stod(cols[5]);
        row.n_Ar = std::stod(cols[6]);
        row.n_N  = std::stod(cols[7]);
        row.n_AO = std::stod(cols[8]);
        row.n_NO = std::stod(cols[9]);
        row.T_K  = std::stod(cols[10]);
        return true;
      } catch (...) { return false; }
    }
    ++line_no;
  }
  return false;
}

static fmx::geom::Mesh make_cube(double a=1.0) {
  fmx::geom::Mesh m; double h=a*0.5;
  auto q=[&](Vec3 a,Vec3 b,Vec3 c,Vec3 d){ m.tris.push_back({a,b,c}); m.tris.push_back({d,a,c}); };
  q({ h,-h,-h},{ h, h,-h},{ h, h, h},{ h,-h, h});
  q({-h, h, h},{-h, h,-h},{-h,-h,-h},{-h,-h, h});
  q({-h, h, h},{ h, h, h},{ h, h,-h},{-h, h,-h});
  q({-h,-h,-h},{ h,-h,-h},{ h,-h, h},{-h,-h, h});
  q({-h,-h, h},{ h,-h, h},{ h, h, h},{-h, h, h});
  q({-h, h,-h},{ h, h,-h},{ h,-h,-h},{-h,-h,-h});
  return m;
}

int main() {
  const std::string csv = "database_200km.csv";
  // According to user: use row index 4654 (0-based after header) for atmospheric input
  const int target_line = 4654;
  DBRow row; if (!read_db_row_by_line(csv, target_line, row)) {
    std::cout << "[cd_cube_dsmc] database_200km.csv not found or row missing; skipping test.\n";
    return 0; // skip
  }
  // AoA (deg) and Cd reference per user note
  double AoA_deg = -3.604146e+00;
  double cd_ref = 9.843571e-01;
  // Build cube and rotate by AoA about Z
  auto mesh = make_cube(1.0);
  auto facets = mesh.to_facets(0);
  double th = AoA_deg * M_PI/180.0; double c=std::cos(th), s=std::sin(th);
  for (auto& f : facets) {
    auto rotv = [&](Vec3 v){ return Vec3{ v.x*c - v.y*s, v.x*s + v.y*c, v.z }; };
    f.n = rotv(f.n).normalized();
    f.r_center = rotv(f.r_center);
  }
  // Assemble atmosphere from DB row
  fmx::solver::Input in;
  in.facets = facets;
  in.materials = { {1.0, 1.0, 1.0, 300.0} };
  // Convert number densities [1/m^3] to mass densities [kg/m^3]
  using namespace fmx::units;
  auto add_sp = [&](double n, double m){ in.species.push_back({ n * m, m }); };
  add_sp(row.n_O,  m_O);
  add_sp(row.n_N2, m_N2);
  add_sp(row.n_O2, m_O2);
  add_sp(row.n_He, m_He);
  add_sp(row.n_H,  m_H);
  add_sp(row.n_Ar, m_Ar);
  add_sp(row.n_N,  m_N);
  add_sp(row.n_AO, m_O);
  add_sp(row.n_NO, m_NO);
  in.T_K = row.T_K;
  // Orbital velocity at 200 km altitude (circular): V = sqrt(mu / (R+h))
  const double mu_E = 3.986004418e14; // m^3/s^2
  const double R_E = 6371000.0;       // m
  const double h = 200000.0;          // m
  double V_orbit = std::sqrt(mu_E / (R_E + h));
  in.V_sat_ms = {V_orbit, 0.0, 0.0};
  in.wind_ms = {0.0, 0.0, 0.0};
  in.r_CG = {0.0, 0.0, 0.0};
  auto out = fmx::solver::solve(in);
  // Compute Cd = -Fx / (0.5 rho V^2 Sref), Sref = 1 m^2
  double rho_tot = 0.0; for (const auto& sp : in.species) rho_tot += sp.rho;
  double V = (in.V_sat_ms - in.wind_ms).norm();
  double q = 0.5 * rho_tot * V * V;
  double Sref = 1.0;
  double Cd = (-out.F.x) / (q * Sref);
  double rel_err = std::abs(Cd - cd_ref) / std::max(1e-12, std::abs(cd_ref));
  // Allow generous tolerance (15%) because of model differences (Stub atmosphere, FM model vs DSMC specifics)
  if (rel_err > 0.15) {
    std::cerr << "Cd mismatch: got="<<Cd<<" ref="<<cd_ref<<" rel_err="<<rel_err<<"\n";
    return 1;
  }
  std::cout << "[cd_cube_dsmc] OK: Cd="<<Cd<<" ref="<<cd_ref<<" rel_err="<<rel_err<<"\n";
  return 0;
}
