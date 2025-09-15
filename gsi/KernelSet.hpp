// GSI KernelSet: lookup for tabulated (theta, Ma, tau, alpha_n, alpha_t) -> (C_N, C_T)
#pragma once

#include <string>
#include <vector>
#include <tuple>

namespace fmx::gsi {

struct KernelRow {
  double theta;   // [rad]
  double Ma;      // speed ratio (|c|/sqrt(kT/m))
  double tau;     // Tw/T
  double alpha_n; // normal accommodation
  double alpha_t; // tangential accommodation
  double CN;      // normal coefficient
  double CT;      // tangential coefficient
};

class KernelSet {
public:
  bool load_csv(const std::string& path);
  // Nearest-neighbor lookup (placeholder); later replaced by interpolation
  std::tuple<double,double> query(double theta, double Ma, double tau, double an, double at) const;
  bool valid() const { return !rows.empty(); }
private:
  std::vector<KernelRow> rows;
  // Optional grid representation for multilinear interpolation
  std::vector<double> ax_theta, ax_Ma, ax_tau, ax_an, ax_at;
  std::vector<double> grid_CN; // size = NT * NM * Nk * Na * Nat
  std::vector<double> grid_CT;
  bool is_grid{false};
  inline size_t idx5(size_t it, size_t im, size_t ik, size_t ia, size_t ib) const {
    size_t NT=ax_theta.size(), NM=ax_Ma.size(), NK=ax_tau.size(), NA=ax_an.size(), NB=ax_at.size();
    return (((ib*NA + ia)*NK + ik)*NM + im)*NT + it;
  }
  std::tuple<double,double> query_grid(double theta, double Ma, double tau, double an, double at) const;
};

} // namespace fmx::gsi
