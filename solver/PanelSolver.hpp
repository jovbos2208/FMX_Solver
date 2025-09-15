// Serial panel solver (no occlusion) computing forces and moments
#pragma once

#include <vector>
#include "core/types.hpp"
#include "gsi/Sentman.hpp"
#include "gsi/CLL.hpp"
#include "gsi/KernelSet.hpp"
#include "gsi/CLLRuntime.hpp"
#include "solver/RegimeAdapter.hpp"
#include "geom/Occluder.hpp"

namespace fmx::solver {

enum class GsiModel { Sentman, CLL };

struct Material {
  double alpha_n{1.0};
  double alpha_t{1.0};
  double alpha_E{1.0};
  double Tw_K{300.0};
};

struct Species {
  double rho;   // species mass density [kg/m^3]
  double mass;  // species molecular mass [kg]
};

struct Input {
  std::vector<fmx::Facet> facets;
  std::vector<Material> materials; // indexed by facet.material_id
  std::vector<Species> species;    // species list (e.g., O, N2, O2, He, H)
  double T_K{800.0};               // ambient temperature [K]
  fmx::Vec3 V_sat_ms{0,0,0};       // satellite velocity in ECEF/whatever frame
  fmx::Vec3 wind_ms{0,0,0};        // atmospheric wind in same frame
  fmx::Vec3 r_CG{0,0,0};           // center of gravity for moments
  const fmx::geom::Occluder* occluder{nullptr}; // optional occlusion
  GsiModel gsi_model{GsiModel::Sentman};
  const fmx::gsi::KernelSet* cll_kernel{nullptr}; // optional CLL table
  fmx::gsi::CLLRuntime* cll_runtime{nullptr};     // optional CLL runtime service
  // Regime adapter (optional, per-facet blending)
  const RegimeConfig* regime{nullptr};
  double regime_Kn{0.0};
  double regime_beta{0.0};
};

struct Output {
  fmx::Vec3 F;  // total force [N]
  fmx::Vec3 M;  // total moment about CG [N*m]
};

Output solve_serial(const Input& in);
Output solve(const Input& in);

} // namespace fmx::solver
