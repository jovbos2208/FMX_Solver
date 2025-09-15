#include "solver/PanelSolver.hpp"
#include "core/units.hpp"
#include <cmath>

namespace fmx::solver {

using fmx::Vec3;

static inline double clamp(double x, double lo, double hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

Output solve_serial(const Input& in) {
  Output out{};
  const Vec3 c = in.V_sat_ms - in.wind_ms; // relative velocity
  const double c_norm = c.norm();
  if (c_norm == 0.0) return out;
  const Vec3 chat = c / c_norm;

  for (const auto& f : in.facets) {
    // Occlusion test: cast along -chat from facet center
    if (in.occluder) {
      fmx::geom::Ray ray{f.r_center, (-chat)};
      if (in.occluder->any_hit(ray, 1e9)) continue; // occluded -> no contribution
    }
    // Incidence cosine: mu = -c_hat · n; if <= 0, no flux on this facet
    double mu = -Vec3::dot(chat, f.n);
    if (mu <= 0.0 || f.area <= 0.0) continue;

    // Tangential direction: projection of -c_hat onto facet plane
    Vec3 tvec = (-chat) - (mu) * f.n;
    double tnorm = tvec.norm();
    Vec3 that = (tnorm > 0.0) ? (tvec / tnorm) : Vec3{0,0,0};

    const auto& mat = (f.material_id < in.materials.size()) ? in.materials[f.material_id] : Material{};
    const double tau = (in.T_K > 0.0) ? (mat.Tw_K / in.T_K) : 1.0;

    Vec3 Fi{0,0,0};
    for (const auto& sp : in.species) {
      const double Ma = c_norm / std::sqrt(fmx::units::k_B * in.T_K / sp.mass);
      double theta = std::acos(clamp(mu, 0.0, 1.0));
      double CN=0.0, CT=0.0;
      if (in.gsi_model == GsiModel::Sentman) {
        std::tie(CN, CT) = fmx::gsi::coefficients(theta, Ma, tau, fmx::gsi::SentmanParams{mat.alpha_E});
      } else {
        if (in.cll_runtime) {
          auto res = in.cll_runtime->query(theta, Ma, tau, mat.alpha_n, mat.alpha_t);
          CN = res.first; CT = res.second;
        } else if (in.cll_kernel && in.cll_kernel->valid()) {
          std::tie(CN, CT) = in.cll_kernel->query(theta, Ma, tau, mat.alpha_n, mat.alpha_t);
        } else {
          std::tie(CN, CT) = fmx::gsi::coefficients(theta, Ma, tau, fmx::gsi::CLLParams{mat.alpha_n, mat.alpha_t});
        }
      }
      // Regime per-facet blending (optional)
      if (in.regime && in.regime->enabled && in.regime->corr_mode == RegimeConfig::CorrMode::PerFacet) {
        double theta = std::acos(clamp(mu, 0.0, 1.0));
        double sN = 1.0 / (1.0 + in.regime->aN * std::pow(std::max(1e-12, in.regime_Kn), in.regime->bN) * std::sin(theta));
        double sT = 1.0 / (1.0 + in.regime->aT * std::pow(std::max(1e-12, in.regime_Kn), in.regime->bT) * std::sin(theta));
        double effN = (1.0 - in.regime_beta) + in.regime_beta * sN;
        double effT = (1.0 - in.regime_beta) + in.regime_beta * sT;
        CN *= effN; CT *= effT;
      }
      const double p_inf = sp.rho * c_norm * c_norm;
      Vec3 dF = (f.n * (CN) + that * (CT)) * (p_inf * f.area);
      Fi += dF;
    }

    out.F += Fi;
    Vec3 r = f.r_center - in.r_CG;
    out.M += Vec3::cross(r, Fi);
  }

  return out;
}

Output solve(const Input& in) {
  const Vec3 c = in.V_sat_ms - in.wind_ms; // relative velocity
  const double c_norm = c.norm();
  if (c_norm == 0.0) return {};
  const Vec3 chat = c / c_norm;

#if defined(FMX_USE_OPENMP)
  double Fx=0, Fy=0, Fz=0;
  double Mx=0, My=0, Mz=0;
  const std::size_t N = in.facets.size();
  #pragma omp parallel for reduction(+:Fx,Fy,Fz,Mx,My,Mz)
  for (long long i = 0; i < static_cast<long long>(N); ++i) {
    const auto& f = in.facets[static_cast<std::size_t>(i)];
    if (in.occluder) {
      fmx::geom::Ray ray{f.r_center, (-chat)};
      if (in.occluder->any_hit(ray, 1e9)) continue;
    }

    double mu = -Vec3::dot(chat, f.n);
    if (mu <= 0.0 || f.area <= 0.0) continue;

    Vec3 tvec = (-chat) - (mu) * f.n;
    double tnorm = tvec.norm();
    Vec3 that = (tnorm > 0.0) ? (tvec / tnorm) : Vec3{0,0,0};

    const auto& mat = (f.material_id < in.materials.size()) ? in.materials[f.material_id] : Material{};
    const double tau = (in.T_K > 0.0) ? (mat.Tw_K / in.T_K) : 1.0;

    double dFx=0, dFy=0, dFz=0;
    for (const auto& sp : in.species) {
      const double Ma = c_norm / std::sqrt(fmx::units::k_B * in.T_K / sp.mass);
      double theta = std::acos(clamp(mu, 0.0, 1.0));
      double CN=0.0, CT=0.0;
      if (in.gsi_model == GsiModel::Sentman) {
        std::tie(CN, CT) = fmx::gsi::coefficients(theta, Ma, tau, fmx::gsi::SentmanParams{mat.alpha_E});
      } else {
        if (in.cll_runtime) {
          auto res = in.cll_runtime->query(theta, Ma, tau, mat.alpha_n, mat.alpha_t);
          CN = res.first; CT = res.second;
        } else if (in.cll_kernel && in.cll_kernel->valid()) {
          std::tie(CN, CT) = in.cll_kernel->query(theta, Ma, tau, mat.alpha_n, mat.alpha_t);
        } else {
          std::tie(CN, CT) = fmx::gsi::coefficients(theta, Ma, tau, fmx::gsi::CLLParams{mat.alpha_n, mat.alpha_t});
        }
      }
      const double p_inf = sp.rho * c_norm * c_norm;
      Vec3 dF = (f.n * (CN) + that * (CT)) * (p_inf * f.area);
      dFx += dF.x; dFy += dF.y; dFz += dF.z;
    }

    Fx += dFx; Fy += dFy; Fz += dFz;
    Vec3 r = f.r_center - in.r_CG;
    Vec3 Mi = Vec3::cross(r, {dFx, dFy, dFz});
    Mx += Mi.x; My += Mi.y; Mz += Mi.z;
  }
  return { {Fx,Fy,Fz}, {Mx,My,Mz} };
#else
  return solve_serial(in);
#endif
}

} // namespace fmx::solver
