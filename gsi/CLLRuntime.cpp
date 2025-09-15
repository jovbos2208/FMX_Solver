#include "gsi/CLLRuntime.hpp"
#include <algorithm>
#include <cmath>

namespace fmx::gsi {

static double round_to_step(double x, double step) {
  if (step <= 0.0) return x;
  return std::round(x / step) * step;
}

static double nearest_in(const std::vector<double>& ax, double x) {
  if (ax.empty()) return x;
  double best = ax.front();
  double bestd = std::abs(x - best);
  for (double v : ax) {
    double d = std::abs(x - v);
    if (d < bestd) { best = v; bestd = d; }
  }
  return best;
}

CLLRuntime::CLLRuntime(const CLLQuantization& q) : q_(q) {
  int nw = std::max(1, q_.workers);
  for (int i=0;i<nw;++i) {
    workers_.emplace_back([this]{ worker_loop(); });
  }
}

CLLRuntime::~CLLRuntime() {
  stop_.store(true);
  q_cv_.notify_all();
  for (auto& t : workers_) { if (t.joinable()) t.join(); }
}

CLLRuntime::Key CLLRuntime::quantize(double theta, double Ma, double tau, double an, double at) const {
  Key k{};
  k.theta = round_to_step(theta, q_.theta_deg_step * M_PI / 180.0);
  k.Ma    = nearest_in(q_.Ma_anchors, Ma);
  k.tau   = round_to_step(tau, q_.tau_step);
  k.an    = round_to_step(an, q_.alpha_step);
  k.at    = round_to_step(at, q_.alpha_step);
  return k;
}

std::pair<double,double> CLLRuntime::query(double theta, double Ma, double tau, double an, double at) const {
  Key k = quantize(theta, Ma, tau, an, at);
  // Check cache
  {
    std::lock_guard<std::mutex> lk(mtx_);
    auto it = cache_.find(k);
    if (it != cache_.end()) {
      return it->second.fut.get();
    }
  }
  // Insert in-flight entry
  auto prom = std::make_shared<Promise>();
  Future fut = prom->get_future().share();
  {
    std::lock_guard<std::mutex> lk(mtx_);
    cache_.insert({k, Entry{fut}});
  }
  // Enqueue task
  {
    std::lock_guard<std::mutex> lk(q_mtx_);
    queue_.emplace_back(k, prom);
  }
  q_cv_.notify_one();
  // Wait for result
  return fut.get();
}

void CLLRuntime::worker_loop() {
  while (!stop_.load()) {
    std::pair<Key, std::shared_ptr<Promise>> task;
    {
      std::unique_lock<std::mutex> lk(q_mtx_);
      q_cv_.wait(lk, [&]{ return stop_.load() || !queue_.empty(); });
      if (stop_.load()) break;
      task = std::move(queue_.front());
      queue_.pop_front();
    }
    // Compute
    auto [CN, CT] = compute_cll(task.first.theta, task.first.Ma, task.first.tau, task.first.an, task.first.at, q_.gh_order);
    task.second->set_value({CN, CT});
  }
}

static void gh_nodes(int order, const double*& x, const double*& w, int& n) {
  static const double X8[8] = {
    -2.930637420257244019223, -1.981656756695842925855, -1.157193712446780194721,
    -0.381186990207322116854,  0.381186990207322116854,  1.157193712446780194721,
     1.981656756695842925855,  2.930637420257244019223
  };
  static const double W8[8] = {
    0.000199604072211367619206, 0.017077983007413475456, 0.20780232581489187954,
    0.66114701255824129103,     0.66114701255824129103,  0.20780232581489187954,
    0.017077983007413475456,    0.000199604072211367619206
  };
  static const double X16[16] = {
    -4.688738939305818364688, -3.869447904860122698719, -3.176999161979956026813, -2.546202157847481362159,
    -1.951787990916253977435, -1.382825855449042263714, -0.832978146885495870749, -0.297320658949765406372,
     0.297320658949765406372,  0.832978146885495870749,  1.382825855449042263714,  1.951787990916253977435,
     2.546202157847481362159,  3.176999161979956026813,  3.869447904860122698719,  4.688738939305818364688
  };
  static const double W16[16] = {
    2.26484121477513627e-10, 1.56739878838904129e-07, 2.48585858557841925e-05, 0.00104982845311528136,
    0.015083977521474187,     0.081274388361574411,    0.181081356895535475,    0.234217974906362647,
    0.234217974906362647,     0.181081356895535475,    0.081274388361574411,    0.015083977521474187,
    0.00104982845311528136,   2.48585858557841925e-05, 1.56739878838904129e-07, 2.26484121477513627e-10
  };
  if (order >= 16) { x = X16; w = W16; n = 16; }
  else { x = X8; w = W8; n = 8; }
}

std::pair<double,double> CLLRuntime::compute_cll(double theta, double Ma, double tau, double an, double at, int gh_order) {
  // Deterministic Gauss–Hermite integration (8-point) of incoming moments
  // and analytical mixture for outgoing CLL (diffuse + specular proportions).
  const double mu = std::cos(theta);
  const double st = std::sqrt(std::max(0.0, 1.0 - mu*mu));
  const double S = Ma / std::sqrt(2.0);
  const double ax = -S * st;
  const double az = -S * mu;
  const double sqrt_tau = std::sqrt(std::max(0.0, tau));

  const double *X, *W; int N;
  gh_nodes(gh_order, X, W, N);
  const double inv_pi32 = 1.0 / (std::pow(M_PI, 1.5));
  const double ex_ax = std::exp(-ax*ax);
  const double ex_az = std::exp(-az*az);

  double flux_in=0.0, Mxz_in=0.0, Mzz_in=0.0;
  double Tz3_in=0.0, Xz2_in=0.0;
  for (int ix=0; ix<N; ++ix) {
    double wx = X[ix]; double wx_w = W[ix]; double fx = std::exp(2.0*ax*wx) * ex_ax;
    for (int iy=0; iy<N; ++iy) {
      double wy = X[iy]; (void)wy; double wy_w = W[iy];
      for (int iz=0; iz<N; ++iz) {
        double wz = X[iz]; if (wz >= 0.0) continue; double wz_w = W[iz];
        double fz = std::exp(2.0*az*wz) * ex_az;
        double weight = wx_w * wy_w * wz_w;
        double common = fx * fz * weight * inv_pi32;
        flux_in += (-wz) * common;
        Mxz_in  += (wx * wz) * common;
        Mzz_in  += (wz * wz) * common;
        Tz3_in  += ((-wz)*(-wz)*(-wz)) * common;
        Xz2_in  += (wx * wz * wz) * common * (-1.0);
      }
    }
  }

  const double flux0_tau = 0.5 / std::sqrt(M_PI) * sqrt_tau;
  const double Mzz0_tau  = 0.5 / std::sqrt(M_PI) * tau;
  double flux_out = (1.0 - an) * flux_in + an * flux0_tau;
  double Mzz_out  = (1.0 - an) * Tz3_in  + an * Mzz0_tau;
  // Tangential shear: CLL tangential kernel yields mean tangential momentum reduced by (1-αt);
  // diffuse contribution has zero mean in tangential component, so only the specular-correlated part contributes.
  // Use (1-αt) times the incident cross-moment Xz2_in; the flux scaling via 'ratio' accounts for different flux_out.
  double Mxz_out  = (1.0 - at) * Xz2_in;

  double ratio = (flux_out > 0.0) ? (flux_in / flux_out) : 0.0;
  double tx_hat = Mxz_in - ratio * Mxz_out;
  double tz_hat = Mzz_in - ratio * Mzz_out;
  const double inv_S2 = 1.0 / (S*S);
  double CT = tx_hat * inv_S2;
  double CN = tz_hat * inv_S2;
  return {CN, CT};
}

} // namespace fmx::gsi
