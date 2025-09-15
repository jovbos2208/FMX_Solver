#include "gsi/KernelSet.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <algorithm>
#include <cmath>

namespace fmx::gsi {

static void split_csv_line(const std::string& s, std::vector<double>& out) {
  out.clear();
  std::istringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    try { out.push_back(std::stod(tok)); } catch (...) {}
  }
}

bool KernelSet::load_csv(const std::string& path) {
  std::ifstream in(path);
  if (!in) return false;
  rows.clear();
  ax_theta.clear(); ax_Ma.clear(); ax_tau.clear(); ax_an.clear(); ax_at.clear();
  grid_CN.clear(); grid_CT.clear(); is_grid = false;
  std::string line;
  // Detect grid header: a line starting with "dims: Ntheta,NMa,Ntau,Nan,Nat"
  std::streampos start_pos = in.tellg();
  while (std::getline(in, line)) {
    if (line.rfind("#", 0) == 0 || line.empty()) continue;
    if (line.rfind("dims:", 0) == 0) {
      // Grid format
      is_grid = true;
      int NT=0,NM=0,NK=0,NA=0,NB=0;
      {
        auto s = line.substr(5);
        std::vector<double> tmp; split_csv_line(s, tmp);
        if (tmp.size() >= 5) { NT=(int)tmp[0]; NM=(int)tmp[1]; NK=(int)tmp[2]; NA=(int)tmp[3]; NB=(int)tmp[4]; }
      }
      auto expect_axis = [&](const char* prefix, std::vector<double>& ax){
        std::string l; if (!std::getline(in, l)) return false; auto pos = l.find(':'); if (pos==std::string::npos) return false; std::string vals = l.substr(pos+1); split_csv_line(vals, ax); return !ax.empty(); };
      if (!expect_axis("theta:", ax_theta)) return false;
      if (!expect_axis("Ma:", ax_Ma)) return false;
      if (!expect_axis("tau:", ax_tau)) return false;
      if (!expect_axis("alpha_n:", ax_an)) return false;
      if (!expect_axis("alpha_t:", ax_at)) return false;
      if ((int)ax_theta.size()!=NT || (int)ax_Ma.size()!=NM || (int)ax_tau.size()!=NK || (int)ax_an.size()!=NA || (int)ax_at.size()!=NB) return false;
      size_t total = (size_t)NT*NM*NK*NA*NB;
      grid_CN.resize(total);
      grid_CT.resize(total);
      // Read CN and CT grids: one line per element: CN,CT
      size_t idx=0;
      while (idx<total && std::getline(in, line)) {
        if (line.empty() || line[0]=='#') continue;
        std::vector<double> vals; split_csv_line(line, vals);
        if (vals.size()>=2) { grid_CN[idx]=vals[0]; grid_CT[idx]=vals[1]; ++idx; }
      }
      return idx==total;
    } else {
      // Not grid: rewind and parse as point CSV
      in.clear(); in.seekg(start_pos);
      break;
    }
  }
  in.clear(); in.seekg(start_pos);
  // Point list format
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;
    std::istringstream ss(line);
    KernelRow r{};
    char comma;
    // Expect: theta_rad,Ma,tau,alpha_n,alpha_t,CN,CT
    if (!(ss >> r.theta)) continue; ss >> comma;
    if (!(ss >> r.Ma)) continue; ss >> comma;
    if (!(ss >> r.tau)) continue; ss >> comma;
    if (!(ss >> r.alpha_n)) continue; ss >> comma;
    if (!(ss >> r.alpha_t)) continue; ss >> comma;
    if (!(ss >> r.CN)) continue; ss >> comma;
    if (!(ss >> r.CT)) continue;
    rows.push_back(r);
  }
  return !rows.empty();
}

std::tuple<double,double> KernelSet::query(double theta, double Ma, double tau, double an, double at) const {
  if (is_grid && !grid_CN.empty()) return query_grid(theta, Ma, tau, an, at);
  if (rows.empty()) return {0.0, 0.0};
  // Simple nearest neighbor in normalized space
  auto norm = [](double x, double lo, double hi){ return (hi>lo) ? ( (x-lo)/(hi-lo) ) : 0.0; };
  double tmin=0.0, tmax=3.141592653589793/2.0; // 0..pi/2
  // Scan to find rough bounds for normalization
  double Ma_min=std::numeric_limits<double>::infinity(), Ma_max=0.0;
  double tau_min=std::numeric_limits<double>::infinity(), tau_max=0.0;
  for (const auto& r : rows) { Ma_min=std::min(Ma_min,r.Ma); Ma_max=std::max(Ma_max,r.Ma); tau_min=std::min(tau_min,r.tau); tau_max=std::max(tau_max,r.tau);} 
  double bestd = std::numeric_limits<double>::infinity();
  const KernelRow* best = &rows.front();
  for (const auto& r : rows) {
    double dtheta = norm(theta, tmin, tmax) - norm(r.theta, tmin, tmax);
    double dMa    = norm(Ma, Ma_min, Ma_max) - norm(r.Ma, Ma_min, Ma_max);
    double dtau   = norm(tau, tau_min, tau_max) - norm(r.tau, tau_min, tau_max);
    double dan    = an - r.alpha_n;
    double dat    = at - r.alpha_t;
    double d = dtheta*dtheta + dMa*dMa + dtau*dtau + dan*dan + dat*dat;
    if (d < bestd) { bestd = d; best = &r; }
  }
  return {best->CN, best->CT};
}

static size_t find_lower(const std::vector<double>& ax, double x) {
  if (x <= ax.front()) return 0;
  if (x >= ax.back()) return ax.size()-2;
  auto it = std::upper_bound(ax.begin(), ax.end(), x);
  size_t i = (size_t)std::distance(ax.begin(), it);
  return i>0? i-1 : 0;
}

std::tuple<double,double> KernelSet::query_grid(double theta, double Ma, double tau, double an, double at) const {
  // Clamp and find lower indices for each axis
  auto clampv = [](double x,double lo,double hi){ return x<lo?lo:(x>hi?hi:x); };
  theta = clampv(theta, ax_theta.front(), ax_theta.back());
  Ma    = clampv(Ma,    ax_Ma.front(),    ax_Ma.back());
  tau   = clampv(tau,   ax_tau.front(),   ax_tau.back());
  an    = clampv(an,    ax_an.front(),    ax_an.back());
  at    = clampv(at,    ax_at.front(),    ax_at.back());
  size_t iT = find_lower(ax_theta, theta);
  size_t iM = find_lower(ax_Ma,    Ma);
  size_t iK = find_lower(ax_tau,   tau);
  size_t iA = find_lower(ax_an,    an);
  size_t iB = find_lower(ax_at,    at);
  double t = (theta - ax_theta[iT]) / std::max(1e-12, ax_theta[iT+1]-ax_theta[iT]);
  double m = (Ma    - ax_Ma[iM])   / std::max(1e-12, ax_Ma[iM+1]-ax_Ma[iM]);
  double k = (tau   - ax_tau[iK])  / std::max(1e-12, ax_tau[iK+1]-ax_tau[iK]);
  double a = (an    - ax_an[iA])   / std::max(1e-12, ax_an[iA+1]-ax_an[iA]);
  double b = (at    - ax_at[iB])   / std::max(1e-12, ax_at[iB+1]-ax_at[iB]);
  // 5D multilinear interpolation (32 corners)
  double CN=0.0, CT=0.0;
  for (int db=0; db<=1; ++db) {
    for (int da=0; da<=1; ++da) {
      for (int dk=0; dk<=1; ++dk) {
        for (int dm=0; dm<=1; ++dm) {
          for (int dt=0; dt<=1; ++dt) {
            double w = (dt? t:1-t) * (dm? m:1-m) * (dk? k:1-k) * (da? a:1-a) * (db? b:1-b);
            size_t it = iT + dt;
            size_t im = iM + dm;
            size_t ik = iK + dk;
            size_t ia = iA + da;
            size_t ib = iB + db;
            size_t idx = idx5(it, im, ik, ia, ib);
            CN += w * grid_CN[idx];
            CT += w * grid_CT[idx];
          }
        }
      }
    }
  }
  return {CN, CT};
}

} // namespace fmx::gsi
