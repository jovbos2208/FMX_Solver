#include "gsi/Sentman.hpp"
#include <cmath>
#include <algorithm>

namespace fmx::gsi {

static inline double clamp(double x, double lo, double hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

// 8-point Gauss-Hermite nodes and weights for exp(-x^2) on (-inf, inf)
// Source: standard tables; nodes symmetric, weights positive.
static inline void gh8(const double*& x, const double*& w, int& n) {
  static const double nodes[8] = {
    -2.930637420257244019223, -1.981656756695842925855, -1.157193712446780194721,
    -0.381186990207322116854,  0.381186990207322116854,  1.157193712446780194721,
     1.981656756695842925855,  2.930637420257244019223
  };
  static const double weights[8] = {
    0.000199604072211367619206, 0.017077983007413475456, 0.20780232581489187954,
    0.66114701255824129103,     0.66114701255824129103,  0.20780232581489187954,
    0.017077983007413475456,    0.000199604072211367619206
  };
  x = nodes; w = weights; n = 8;
}

// Compute dimensionless incoming moments and flux for drift a = (ax, 0, az)
// over half-space w_z < 0 under weight exp(-|w-a|^2).
struct MomentsIn {
  double Mxz{0}, Mzz{0};
  double flux{0}; // number flux: ∫_{wz<0} (-wz) e^{-|w-a|^2} dw
  double eflux{0}; // energy flux moment: ∫_{wz<0} (-wz) (wx^2+wy^2+wz^2) e^{-|w-a|^2} dw
};

static inline MomentsIn incoming_moments(double ax, double az) {
  const double *x, *w; int n; gh8(x, w, n);
  MomentsIn mi{};
  const double ex_ax = std::exp(-ax*ax);
  const double ex_az = std::exp(-az*az);
  for (int ix = 0; ix < n; ++ix) {
    const double wx = x[ix];
    const double wx_w = w[ix];
    const double fx = std::exp(2.0*ax*wx) * ex_ax;
    for (int iy = 0; iy < n; ++iy) {
      const double wy = x[iy]; (void)wy; // symmetry, a_y=0
      const double wy_w = w[iy];
      // a_y = 0 => factor for y is 1
      for (int iz = 0; iz < n; ++iz) {
        const double wz = x[iz];
        if (wz >= 0.0) continue; // incoming half-space w_z < 0
        const double wz_w = w[iz];
        const double fz = std::exp(2.0*az*wz) * ex_az;
        const double weight = wx_w * wy_w * wz_w;
        const double common = fx * fz * weight;
        mi.Mxz += wx * wz * common;
        mi.Mzz += wz * wz * common;
        mi.flux += (-wz) * common;
        mi.eflux += (-wz) * (wx*wx + wy*wy + wz*wz) * common;
      }
    }
  }
  // Include the π^{-3/2} factor
  const double inv_pi32 = 1.0 / (std::pow(M_PI, 1.5));
  mi.Mxz *= inv_pi32;
  mi.Mzz *= inv_pi32;
  mi.flux *= inv_pi32;
  return mi;
}

// Outgoing (diffuse, a=0) moments over half-space w_z > 0 under exp(-|w|^2)
struct MomentsOut0 {
  double Mxz{0}, Mzz{0};
  double flux{0};
};

static inline MomentsOut0 outgoing_moments_zero_drift() {
  static bool inited = false;
  static MomentsOut0 mo;
  if (inited) return mo;
  const double *x, *w; int n; gh8(x, w, n);
  for (int ix = 0; ix < n; ++ix) {
    const double wx = x[ix];
    const double wx_w = w[ix];
    for (int iy = 0; iy < n; ++iy) {
      const double wy_w = w[iy]; (void)wy_w; // symmetry in y
      for (int iz = 0; iz < n; ++iz) {
        const double wz = x[iz];
        if (wz <= 0.0) continue; // outgoing half-space
        const double wz_w = w[iz];
        const double weight = wx_w * wy_w * wz_w;
        mo.Mxz += wx * wz * weight;
        mo.Mzz += wz * wz * weight;
        mo.flux += (wz) * weight;
      }
    }
  }
  const double inv_pi32 = 1.0 / (std::pow(M_PI, 1.5));
  mo.Mxz *= inv_pi32;
  mo.Mzz *= inv_pi32;
  mo.flux *= inv_pi32;
  inited = true;
  return mo;
}

// Numerical Sentman via Gauss-Hermite quadrature with full accommodation and Tw=T.
static std::tuple<double,double> numeric_coefficients(double theta, double Ma, double tau_in, double alpha_E) {
  theta = clamp(theta, 0.0, M_PI_2);
  Ma = std::max(1e-8, Ma);
  // Convert to speed ratio S = c / sqrt(2 kT/m) from Ma = c / sqrt(kT/m)
  const double S = Ma / std::sqrt(2.0);
  // Geometry: mu = cos(theta); local basis e_z = n, e_x = t_hat (projection of -chat)
  const double mu = std::cos(theta);
  const double st = std::sqrt(std::max(0.0, 1.0 - mu*mu));
  // Drift vector a = S * chat, where chat·n = -mu, chat tangential component along -e_x with magnitude st.
  const double ax = -S * st;  // along +e_x is opposite to chat's tangential
  const double az = -S * mu;  // negative if mu>0 (front-facing)

  // Incoming moments (dimensionless)
  const MomentsIn mi = incoming_moments(ax, az);
  // Compute E_i/(k T_i) = mi.eflux / mi.flux
  const double Ei_over_kTi = (mi.flux > 0.0) ? (mi.eflux / mi.flux) : 3.0; // fallback 3 for equilibrium
  // Temperature ratio tau = T_r/T_i from energy accommodation (Tuttas et al.)
  const double tau = (1.0 - alpha_E) * 0.5 * Ei_over_kTi + alpha_E * tau_in;

  // Outgoing zero-drift moments (dimensionless at T_i); scale with tau
  MomentsOut0 mo = outgoing_moments_zero_drift();
  const double sqrt_tau = std::sqrt(std::max(0.0, tau));
  mo.flux *= sqrt_tau;
  mo.Mzz  *= tau;
  mo.Mxz   = 0.0; // symmetry
  // Number flux matching -> scale for outgoing density
  const double ratio = (mo.flux > 0.0) ? (mi.flux / mo.flux) : 0.0;

  // Traction components along local axes (per-unit ρ c^2 factor left out)
  const double tx_hat = mi.Mxz - ratio * mo.Mxz;
  const double tz_hat = mi.Mzz - ratio * mo.Mzz;

  // Normalize by ρ c^2: factor = v_th^2 / c^2 = 1 / S^2
  const double inv_S2 = 1.0 / (S * S);
  const double C_T = tx_hat * inv_S2; // along e_x = t_hat
  const double C_N = tz_hat * inv_S2; // along e_z = n

  return {C_N, C_T};
}

// Placeholder closed-form wrapper: delegates to numeric implementation for now.
static std::tuple<double,double> closed_form_coefficients(double theta, double Ma, double tau_in, double alpha_E) {
  // Closed-form Sentman (diffuse, full accommodation kernel), using
  // half-space Gaussian integrals with error functions. Includes
  // generalized reflected temperature via energy accommodation (tau).
  theta = clamp(theta, 0.0, M_PI_2);
  Ma = std::max(1e-8, Ma);
  const double S = Ma / std::sqrt(2.0);
  const double mu = std::cos(theta);
  const double st = std::sqrt(std::max(0.0, 1.0 - mu*mu));
  const double ax = -S * st;
  const double az = -S * mu; // equals S*cos(delta) with cos(delta)=chat·n=-mu

  // Helper lambdas
  auto H0 = [&](double a){
    // 1D flux integral: ∫_{-∞}^{-a} (-(u+a)) e^{-u^2} du
    // Result: 0.5 e^{-a^2} - 0.5 sqrt(pi) a erfc(a)
    return 0.5*std::exp(-a*a) - 0.5*std::sqrt(M_PI)*a*std::erfc(a);
  };
  auto K = [&](double a){
    // 1D moment: ∫_{-∞}^{-a} (-(u+a))(u+a)^2 e^{-u^2} du
    // Result: 0.5 e^{-a^2} - 0.5 sqrt(pi) (2a + a^3) erfc(a)
    return 0.5*std::exp(-a*a) - 0.5*std::sqrt(M_PI)*(2.0*a + a*a*a)*std::erfc(a);
  };

  // Dimensionless incoming flux and momentum moments (normalized by pi^{3/2})
  const double flux_in = (1.0/std::sqrt(M_PI)) * H0(az);
  // Mzz_in
  const double Mzz_in = (1.0/std::sqrt(M_PI)) * K(az);
  // Mxz_in
  const double Izz = 0.5*std::sqrt(M_PI)*std::erfc(az) * (az*az - 1.0) - 0.5*az*std::exp(-az*az);
  const double Mxz_in = (1.0/std::sqrt(M_PI)) * (ax * Izz);

  // Incoming energy flux -> Ei/(k Ti)
  // eflux_in = (1/sqrt(pi)) [ (ax^2 + 1) H0(az) + K(az) ]
  const double eflux_in = (1.0/std::sqrt(M_PI)) * ( (ax*ax + 1.0) * H0(az) + K(az) );
  const double Ei_over_kTi = (flux_in > 0.0) ? (eflux_in / flux_in) : 3.0;

  // Generalized temperature ratio tau via energy accommodation
  const double tau = (1.0 - alpha_E) * 0.5 * Ei_over_kTi + alpha_E * tau_in;
  const double sqrt_tau = std::sqrt(std::max(0.0, tau));

  // Outgoing diffuse zero-drift moments at Ti: flux0 = Mzz0 = 0.5/sqrt(pi)
  const double flux0 = 0.5 / std::sqrt(M_PI);
  const double Mzz0  = 0.5 / std::sqrt(M_PI);
  // Scale by tau
  const double flux_out = flux0 * sqrt_tau;
  const double Mzz_out  = Mzz0 * tau;

  // Density scale by matching number flux
  const double ratio = (flux_out > 0.0) ? (flux_in / flux_out) : 0.0;

  // Traction components (dimensionless, before dividing by S^2)
  const double tx_hat = Mxz_in;            // no outgoing shear
  const double tz_hat = Mzz_in - ratio * Mzz_out;

  const double inv_S2 = 1.0 / (S*S);
  const double C_T = tx_hat * inv_S2;
  const double C_N = tz_hat * inv_S2;
  return {C_N, C_T};
}

std::tuple<double,double> coefficients(double theta, double Ma, double tau_in,
                                       const SentmanParams& p) {
#if defined(FMX_USE_SENTMAN_CLOSED_FORM)
  return closed_form_coefficients(theta, Ma, tau_in, p.alpha_E);
#else
  return numeric_coefficients(theta, Ma, tau_in, p.alpha_E);
#endif
}

} // namespace fmx::gsi
