#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <optional>
#include <algorithm>
#include <random>
#include <chrono>
#include "core/types.hpp"
#include "gsi/Sentman.hpp"
#include "geom/Mesh.hpp"
#include "geom/BVH.hpp"
#include "solver/PanelSolver.hpp"
#include "core/units.hpp"
#include "atm/Atmosphere.hpp"
#include "atm/NRLMSIS2.hpp"
#include "atm/HWM14.hpp"
#include "atm/Combined.hpp"
#include "gsi/KernelSet.hpp"
#include "gsi/CLLRuntime.hpp"
#include "solver/RegimeAdapter.hpp"

namespace {
using fmx::Vec3;

struct CliConfig {
  std::string geometry;
  Vec3 cg{0,0,0};
  double alpha_E{1.0};
  double alpha_n{1.0};
  double alpha_t{1.0};
  double Tw_K{300.0};
  double alt_km{400.0};
  double lat_deg{0.0};
  double lon_deg{0.0};
  std::string utc{"2025-09-12T12:00:00Z"};
  Vec3 V_sat{7500.0, 0.0, 0.0};
  double F10_7{120.0};
  double F10_7A{120.0};
  int Kp{3};
  std::string atm_model{"Stub"};
  bool has_ap7{false};
  double Ap7[7] = {0,0,0,0,0,0,0};
  double Ap_daily{0.0};
  double Ap_now{0.0};
  std::string gsi_model{"Sentman"};
  std::string gsi_table_path;
  // CLL runtime
  double rt_theta_deg_step{2.0};
  double rt_tau_step{0.1};
  double rt_alpha_step{0.25};
  std::vector<double> rt_Ma_anchors{0.5,1.0,2.0,4.0,8.0,12.0,16.0};
  int rt_sample_count{4000};
  int rt_workers{1};
  // Regime adapter
  bool regime_enabled{false};
  double regime_L_char_m{0.0}; // 0 -> auto bbox
  double regime_gamma{1.0};
  std::string regime_corr_mode{"none"};
  double regime_corr_a{0.0};
  double regime_corr_b{1.0};
  // UQ
  int uq_samples{0};
  double uq_sigma_rho{0.3}; // log-normal sigma for densities
  double uq_alpha_spread{0.05}; // uniform +/- on alpha_n/alpha_t
  double uq_Tw_spread{5.0}; // +/- K on Tw
};

static std::string slurp(const std::string& path) {
  std::ifstream in(path);
  std::ostringstream ss; ss << in.rdbuf();
  return ss.str();
}

static std::vector<double> parse_list(const std::string& s) {
  std::vector<double> v; std::string tok; for (size_t i=0,j=0; i<=s.size(); ++i) {
    if (i==s.size() || s[i]==',') { tok = s.substr(j, i-j); try{ v.push_back(std::stod(tok)); }catch(...){} j=i+1; }
  } return v;
}

static bool find_string(const std::string& s, const std::string& key, std::string& out) {
  auto pos = s.find("\"" + key + "\""); if (pos == std::string::npos) return false;
  pos = s.find(':', pos); if (pos == std::string::npos) return false;
  pos = s.find('"', pos); if (pos == std::string::npos) return false;
  auto end = s.find('"', pos+1); if (end == std::string::npos) return false;
  out = s.substr(pos+1, end-pos-1); return true;
}

static bool find_number(const std::string& s, const std::string& key, double& out) {
  auto pos = s.find("\"" + key + "\""); if (pos == std::string::npos) return false;
  pos = s.find(':', pos); if (pos == std::string::npos) return false;
  auto end = s.find_first_of(",}\n\r", pos+1);
  std::string t = s.substr(pos+1, end-pos-1);
  try { out = std::stod(t); return true; } catch (...) { return false; }
}

static bool find_int(const std::string& s, const std::string& key, int& out) {
  double d; if (!find_number(s, key, d)) return false; out = static_cast<int>(d); return true;
}

static bool find_arrayD(const std::string& s, const std::string& key, std::vector<double>& out) {
  auto pos = s.find("\"" + key + "\""); if (pos == std::string::npos) return false;
  pos = s.find('[', pos); if (pos == std::string::npos) return false;
  auto end = s.find(']', pos); if (end == std::string::npos) return false;
  std::string arr = s.substr(pos+1, end-pos-1);
  out.clear();
  std::istringstream ss(arr);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    try { out.push_back(std::stod(tok)); } catch (...) {}
  }
  return !out.empty();
}

static bool find_array3(const std::string& s, const std::string& key, Vec3& out) {
  auto pos = s.find("\"" + key + "\""); if (pos == std::string::npos) return false;
  pos = s.find('[', pos); if (pos == std::string::npos) return false;
  auto end = s.find(']', pos); if (end == std::string::npos) return false;
  std::string arr = s.substr(pos+1, end-pos-1);
  std::replace(arr.begin(), arr.end(), ',', ' ');
  std::istringstream ss(arr);
  if (!(ss >> out.x >> out.y >> out.z)) return false;
  return true;
}

static CliConfig parse_config(const std::string& json) {
  CliConfig c;
  // Simple top-level keys
  find_string(json, "geometry", c.geometry);
  find_array3(json, "cg", c.cg);
  // Nested materials.default
  auto mpos = json.find("\"materials\"");
  if (mpos != std::string::npos) {
    auto dpos = json.find("\"default\"", mpos);
    if (dpos != std::string::npos) {
      std::string sub = json.substr(dpos, std::min<size_t>(json.size()-dpos, 2000));
      double aE; if (find_number(sub, "alpha_E", aE)) c.alpha_E = aE;
      double an; if (find_number(sub, "alpha_n", an)) c.alpha_n = an;
      double at; if (find_number(sub, "alpha_t", at)) c.alpha_t = at;
      double Tw; if (find_number(sub, "Tw_K", Tw)) c.Tw_K = Tw;
    }
  }
  // Atmosphere.indices
  auto apos = json.find("\"atmosphere\"");
  if (apos != std::string::npos) {
    std::string sub = json.substr(apos, std::min<size_t>(json.size()-apos, 2000));
    std::string model; if (find_string(sub, "model", model)) c.atm_model = model;
    double F10; if (find_number(sub, "F10_7", F10)) c.F10_7 = F10;
    double F10A; if (find_number(sub, "F10_7A", F10A)) c.F10_7A = F10A;
    double apd; if (find_number(sub, "Ap_daily", apd)) c.Ap_daily = apd;
    double apnow; if (find_number(sub, "Ap_now", apnow)) c.Ap_now = apnow;
    int kp; if (find_int(sub, "Kp", kp)) c.Kp = kp;
    std::vector<double> apv; if (find_arrayD(sub, "Ap", apv) && apv.size() >= 2) {
      c.has_ap7 = false;
      if (apv.size() >= 7) { for (int i=0;i<7;i++) c.Ap7[i] = apv[i]; c.has_ap7 = true; }
      else { c.Ap_daily = apv[0]; c.Ap_now = apv[1]; }
    }
  }
  // State
  auto spos = json.find("\"state\"");
  if (spos != std::string::npos) {
    std::string sub = json.substr(spos, std::min<size_t>(json.size()-spos, 4000));
    double alt; if (find_number(sub, "alt_km", alt)) c.alt_km = alt;
    double lat; if (find_number(sub, "lat_deg", lat)) c.lat_deg = lat;
    double lon; if (find_number(sub, "lon_deg", lon)) c.lon_deg = lon;
    std::string utc; if (find_string(sub, "utc", utc)) c.utc = utc;
    Vec3 V; if (find_array3(sub, "V_sat_mps", V)) c.V_sat = V;
  }
  // GSI model
  auto gpos = json.find("\"gsi\"");
  if (gpos != std::string::npos) {
    std::string sub = json.substr(gpos, std::min<size_t>(json.size()-gpos, 1000));
    std::string model; if (find_string(sub, "model", model)) c.gsi_model = model;
    std::string table; if (find_string(sub, "table_path", table)) c.gsi_table_path = table;
    // runtime block
    auto rpos = sub.find("\"runtime\"");
    if (rpos != std::string::npos) {
      std::string rsub = sub.substr(rpos, std::min<size_t>(sub.size()-rpos, 1000));
      double v;
      if (find_number(rsub, "theta_deg_step", v)) c.rt_theta_deg_step = v;
      if (find_number(rsub, "tau_step", v)) c.rt_tau_step = v;
      if (find_number(rsub, "alpha_step", v)) c.rt_alpha_step = v;
      int w;
      if (find_int(rsub, "workers", w)) c.rt_workers = w;
      if (find_int(rsub, "sample_count", w)) c.rt_sample_count = w;
      std::string mas; if (find_string(rsub, "Ma_anchors", mas)) c.rt_Ma_anchors = parse_list(mas);
    }
  }
  // Regime
  auto rpos = json.find("\"regime\"");
  if (rpos != std::string::npos) {
    std::string sub = json.substr(rpos, std::min<size_t>(json.size()-rpos, 1000));
    double v; if (find_number(sub, "L_char_m", v)) c.regime_L_char_m = v;
    if (find_number(sub, "gamma", v)) c.regime_gamma = v;
    std::string on; if (find_string(sub, "enabled", on)) c.regime_enabled = (on=="true" || on=="1");
    std::string mode; if (find_string(sub, "corr_mode", mode)) c.regime_corr_mode = mode;
    if (find_number(sub, "corr_a", v)) c.regime_corr_a = v;
    if (find_number(sub, "corr_b", v)) c.regime_corr_b = v;
  }
  // UQ
  auto upos = json.find("\"uq\"");
  if (upos != std::string::npos) {
    std::string sub = json.substr(upos, std::min<size_t>(json.size()-upos, 1000));
    double v; int n;
    if (find_int(sub, "samples", n)) c.uq_samples = n;
    if (find_number(sub, "sigma_rho", v)) c.uq_sigma_rho = v;
    if (find_number(sub, "alpha_spread", v)) c.uq_alpha_spread = v;
    if (find_number(sub, "Tw_spread", v)) c.uq_Tw_spread = v;
  }
  return c;
}

static fmx::geom::Mesh make_plate() {
  fmx::geom::Mesh m;
  // YZ plane, normal toward -X
  m.tris.push_back({Vec3{0, 0.5, 0.5}, Vec3{0, 0.5, -0.5}, Vec3{0, -0.5, -0.5}});
  m.tris.push_back({Vec3{0, -0.5, 0.5}, Vec3{0, 0.5, 0.5}, Vec3{0, -0.5, -0.5}});
  return m;
}

static fmx::geom::Mesh make_two_plates() {
  auto m = make_plate();
  m.tris.push_back({Vec3{0.2, 0.5, 0.5}, Vec3{0.2, 0.5, -0.5}, Vec3{0.2, -0.5, -0.5}});
  m.tris.push_back({Vec3{0.2, -0.5, 0.5}, Vec3{0.2, 0.5, 0.5}, Vec3{0.2, -0.5, -0.5}});
  return m;
}

static fmx::geom::Mesh make_cube(double a=1.0) {
  fmx::geom::Mesh m;
  double h=a*0.5;
  auto add_quad = [&](Vec3 a, Vec3 b, Vec3 c, Vec3 d){ m.tris.push_back({a,b,c}); m.tris.push_back({d,a,c}); };
  // +X
  add_quad({ h,-h,-h},{ h, h,-h},{ h, h, h},{ h,-h, h});
  // -X
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

static void write_json_result(const std::string& path, const fmx::solver::Output& out) {
  std::ofstream of(path);
  of << "{\n";
  of << "  \"F\": [" << out.F.x << ", " << out.F.y << ", " << out.F.z << "],\n";
  of << "  \"M\": [" << out.M.x << ", " << out.M.y << ", " << out.M.z << "]\n";
  of << "}\n";
}
} // namespace

int main(int argc, char** argv) {
  std::cout << "fmx CLI\n";
  std::string config_path;
  std::string validate_case;
  std::string mesh_override;
  std::string out_path;
  double theta_deg = 0.0;
  int bench_iters = 0;
  for (int i=1;i<argc;++i) {
    std::string a = argv[i];
    if (a == "--config" && i+1<argc) config_path = argv[++i];
    else if (a == "--validate" && i+1<argc) validate_case = argv[++i];
    else if (a == "--mesh" && i+1<argc) mesh_override = argv[++i];
    else if (a == "--out" && i+1<argc) out_path = argv[++i];
    else if (a == "--theta_deg" && i+1<argc) theta_deg = std::stod(argv[++i]);
    else if (a == "--bench" && i+1<argc) bench_iters = std::stoi(argv[++i]);
    else if (a == "--help") { std::cout << "Usage: fmx_cli [--config file.json] [--validate plate|two-plates|cube|torque-plate] [--mesh path] [--theta_deg deg] [--bench iters] [--out result.json]\n"; return 0; }
  }

  CliConfig cfg;
  if (!config_path.empty()) {
    try { cfg = parse_config(slurp(config_path)); }
    catch (...) { std::cerr << "Failed to parse config; using defaults\n"; }
  }

  // Construct mesh
  fmx::geom::Mesh mesh;
  if (!validate_case.empty()) {
    if (validate_case == "plate") mesh = make_plate();
    else if (validate_case == "two-plates") mesh = make_two_plates();
    else if (validate_case == "cube") mesh = make_cube();
    else if (validate_case == "torque-plate") mesh = make_plate();
    else { std::cerr << "Unknown validate case: " << validate_case << "\n"; return 1; }
  } else if (!mesh_override.empty()) {
    auto m = fmx::geom::Mesh::load(mesh_override);
    if (!m) { std::cerr << "Failed to load mesh: " << mesh_override << "\n"; return 1; }
    mesh = *m;
  } else if (!cfg.geometry.empty()) {
    auto m = fmx::geom::Mesh::load(cfg.geometry);
    if (!m) { std::cerr << "Failed to load mesh: " << cfg.geometry << "\n"; return 1; }
    mesh = *m;
  } else {
    mesh = make_plate();
  }

  auto facets = mesh.to_facets(0);
  fmx::geom::BVHOccluder occ(mesh.tris);

  // Atmosphere
  fmx::atm::AtmosphereState st{};
  if (cfg.atm_model == "Stub") {
    fmx::atm::StubAtmosphere atm;
    fmx::atm::Indices idx; idx.F10_7 = cfg.F10_7; idx.F10_7A = cfg.F10_7A; idx.Kp = cfg.Kp;
    if (cfg.has_ap7) { idx.has_ap_array = true; for (int i=0;i<7;i++) idx.Ap[i] = cfg.Ap7[i]; }
    idx.Ap_daily = cfg.Ap_daily; idx.Ap_now = cfg.Ap_now;
    st = atm.evaluate(cfg.alt_km, cfg.lat_deg, cfg.lon_deg, cfg.utc, idx);
  } else if (cfg.atm_model == "NRLMSIS2") {
    fmx::atm::NRLMSIS2Atmosphere atm;
    fmx::atm::Indices idx; idx.F10_7 = cfg.F10_7; idx.F10_7A = cfg.F10_7A; idx.Kp = cfg.Kp;
    if (cfg.has_ap7) { idx.has_ap_array = true; for (int i=0;i<7;i++) idx.Ap[i] = cfg.Ap7[i]; }
    idx.Ap_daily = cfg.Ap_daily; idx.Ap_now = cfg.Ap_now;
    st = atm.evaluate(cfg.alt_km, cfg.lat_deg, cfg.lon_deg, cfg.utc, idx);
  } else if (cfg.atm_model == "NRLMSIS2+HWM14") {
    fmx::atm::CombinedAtmosphere atm;
    fmx::atm::Indices idx; idx.F10_7 = cfg.F10_7; idx.F10_7A = cfg.F10_7A; idx.Kp = cfg.Kp;
    if (cfg.has_ap7) { idx.has_ap_array = true; for (int i=0;i<7;i++) idx.Ap[i] = cfg.Ap7[i]; }
    idx.Ap_daily = cfg.Ap_daily; idx.Ap_now = cfg.Ap_now;
    st = atm.evaluate(cfg.alt_km, cfg.lat_deg, cfg.lon_deg, cfg.utc, idx);
  } else {
    std::cerr << "Atmosphere model '" << cfg.atm_model << "' not linked; using Stub.\n";
    fmx::atm::StubAtmosphere atm;
    fmx::atm::Indices idx; idx.F10_7 = cfg.F10_7; idx.F10_7A = cfg.F10_7A; idx.Kp = cfg.Kp;
    if (cfg.has_ap7) { idx.has_ap_array = true; for (int i=0;i<7;i++) idx.Ap[i] = cfg.Ap7[i]; }
    idx.Ap_daily = cfg.Ap_daily; idx.Ap_now = cfg.Ap_now;
    st = atm.evaluate(cfg.alt_km, cfg.lat_deg, cfg.lon_deg, cfg.utc, idx);
  }

  fmx::solver::Input in;
  in.facets = facets;
  in.materials = { {cfg.alpha_n, cfg.alpha_t, cfg.alpha_E, cfg.Tw_K} };
  in.species.clear();
  for (const auto& sp : st.species) in.species.push_back({sp.rho, sp.mass});
  in.T_K = st.T_K;
  in.V_sat_ms = cfg.V_sat;
  in.wind_ms = st.wind_ms;
  // If torque-plate validation, offset CG to induce torque (about +Z)
  if (validate_case == "torque-plate") {
    in.r_CG = {0.0, 0.1, 0.0};
  } else {
    in.r_CG = cfg.cg;
  }
  in.occluder = &occ;

  // Optional mesh rotation about Z
  if (theta_deg != 0.0) {
    double th = theta_deg * M_PI/180.0, c=std::cos(th), s=std::sin(th);
    for (auto& f : in.facets) {
      auto rotv = [&](fmx::Vec3 v){ return fmx::Vec3{ v.x*c - v.y*s, v.x*s + v.y*c, v.z }; };
      f.n = rotv(f.n).normalized();
      f.r_center = rotv(f.r_center);
    }
  }

  // GSI model selection
  if (cfg.gsi_model == "CLL") in.gsi_model = fmx::solver::GsiModel::CLL;
  else in.gsi_model = fmx::solver::GsiModel::Sentman;

  // Optional CLL kernel table
  fmx::gsi::KernelSet cll_table;
  fmx::gsi::CLLQuantization qcfg;
  qcfg.theta_deg_step = cfg.rt_theta_deg_step;
  qcfg.tau_step = cfg.rt_tau_step;
  qcfg.alpha_step = cfg.rt_alpha_step;
  qcfg.Ma_anchors = cfg.rt_Ma_anchors;
  qcfg.sample_count = cfg.rt_sample_count;
  qcfg.workers = cfg.rt_workers;
  fmx::gsi::CLLRuntime cll_runtime(qcfg);
  if (in.gsi_model == fmx::solver::GsiModel::CLL && !cfg.gsi_table_path.empty()) {
    if (cll_table.load_csv(cfg.gsi_table_path)) in.cll_kernel = &cll_table;
    else std::cerr << "Failed to load CLL table at: " << cfg.gsi_table_path << " (falling back)\n";
  }
  if (in.gsi_model == fmx::solver::GsiModel::CLL) {
    in.cll_runtime = &cll_runtime;
  }

  auto solve_once = [&](const fmx::solver::Input& in_local){ return fmx::solver::solve(in_local); };
  auto out = solve_once(in);
  // Regime adapter (optional)
  fmx::solver::RegimeDiagnostics rdiag{};
  if (cfg.regime_enabled) {
    // Auto L_char as bbox diagonal if not set
    double Lchar = cfg.regime_L_char_m;
    if (Lchar <= 0.0) {
      fmx::Vec3 lo{1e300,1e300,1e300}, hi{-1e300,-1e300,-1e300};
      for (const auto& t : mesh.tris) {
        auto upd=[&](const fmx::Vec3& p){ lo.x=std::min(lo.x,p.x); lo.y=std::min(lo.y,p.y); lo.z=std::min(lo.z,p.z);
                                          hi.x=std::max(hi.x,p.x); hi.y=std::max(hi.y,p.y); hi.z=std::max(hi.z,p.z); };
        upd(t.v0); upd(t.v1); upd(t.v2);
      }
      fmx::Vec3 ext = hi - lo; Lchar = std::sqrt(ext.x*ext.x + ext.y*ext.y + ext.z*ext.z);
    }
    fmx::solver::RegimeConfig rcfg{};
    rcfg.enabled = true; rcfg.L_char_m = Lchar; rcfg.gamma = cfg.regime_gamma;
    if (cfg.regime_corr_mode == "scalar") rcfg.corr_mode = fmx::solver::RegimeConfig::CorrMode::Scalar;
    rcfg.corr_a = cfg.regime_corr_a; rcfg.corr_b = cfg.regime_corr_b;
    rdiag = fmx::solver::apply_regime_blend(st, in.species, in.T_K, rcfg, out);
    // Pass per-facet regime blending config to solver input
    in.regime = &rcfg;
    in.regime_Kn = rdiag.Kn;
    in.regime_beta = rdiag.beta;
  }
  // Diagnostics
  // Compute occluded facet count (front-facing only)
  std::size_t occluded = 0, front = 0;
  Vec3 crel = in.V_sat_ms - in.wind_ms; double cn = crel.norm(); Vec3 chat = (cn>0)?(crel/cn):Vec3{1,0,0};
  for (const auto& f : in.facets) {
    double mu = -fmx::Vec3::dot(chat, f.n);
    if (mu <= 0.0) continue; // backface or grazing
    front++;
    fmx::geom::Ray ray{f.r_center, (-chat)};
    if (occ.any_hit(ray, 1e9)) occluded++;
  }
  // Species Mach range
  double Ma_min = 1e300, Ma_max = 0.0;
  for (const auto& sp : in.species) {
    double Ma = cn / std::sqrt(fmx::units::k_B * in.T_K / sp.mass);
    Ma_min = std::min(Ma_min, Ma);
    Ma_max = std::max(Ma_max, Ma);
  }
  double tau = (in.T_K>0)? (cfg.Tw_K / in.T_K) : 1.0;

  std::cout << "F = [" << out.F.x << ", " << out.F.y << ", " << out.F.z << "] N\n";
  std::cout << "M = [" << out.M.x << ", " << out.M.y << ", " << out.M.z << "] N*m\n";
  std::cout << "facets=" << in.facets.size() << ", front=" << front << ", shadowed=" << occluded << "\n";
  std::cout << "alt_km=" << cfg.alt_km << ", T_K=" << in.T_K << ", Tw_K=" << cfg.Tw_K << ", tau=Tw/T=" << tau << "\n";
  std::cout << "Ma[min,max]=" << Ma_min << ", " << Ma_max << "\n";
  if (cfg.regime_enabled) {
    std::cout << "Kn=" << rdiag.Kn << ", beta=" << rdiag.beta << "\n";
  }
  // Bench
  if (bench_iters > 0) {
    auto t0 = std::chrono::high_resolution_clock::now();
    fmx::solver::Output tmp{};
    for (int i=0;i<bench_iters;++i) {
      tmp = solve_once(in);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1-t0).count();
    std::cout << "bench_iters=" << bench_iters << ", total_ms=" << ms << ", avg_ms=" << (ms/bench_iters) << "\n";
  }
  // UQ
  if (cfg.uq_samples > 0) {
    std::mt19937_64 rng(12345);
    std::normal_distribution<double> z(0.0, 1.0);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    auto in_base = in;
    std::vector<double> Fx, Fy, Fz, Mx, My, Mz; Fx.reserve(cfg.uq_samples);
    for (int s=0;s<cfg.uq_samples;++s) {
      auto in_s = in_base;
      // Perturb densities (log-normal)
      for (auto& sp : in_s.species) {
        double factor = std::exp(cfg.uq_sigma_rho * z(rng));
        sp.rho = std::max(0.0, sp.rho * factor);
      }
      // Perturb material
      for (auto& m : in_s.materials) {
        double da = cfg.uq_alpha_spread * u(rng);
        m.alpha_n = std::clamp(m.alpha_n + da, 0.0, 1.0);
        da = cfg.uq_alpha_spread * u(rng);
        m.alpha_t = std::clamp(m.alpha_t + da, 0.0, 1.0);
        m.Tw_K = std::max(0.0, m.Tw_K + cfg.uq_Tw_spread * u(rng));
      }
      auto out_s = solve_once(in_s);
      Fx.push_back(out_s.F.x); Fy.push_back(out_s.F.y); Fz.push_back(out_s.F.z);
      Mx.push_back(out_s.M.x); My.push_back(out_s.M.y); Mz.push_back(out_s.M.z);
    }
    auto pct = [](std::vector<double>& v, double p){ std::sort(v.begin(), v.end()); double idx = p * (v.size()-1); size_t i = (size_t)idx; double frac = idx - i; if (i+1<v.size()) return v[i]*(1-frac)+v[i+1]*frac; else return v[i]; };
    double Fx5=pct(Fx,0.05), Fx50=pct(Fx,0.50), Fx95=pct(Fx,0.95);
    double Fy5=pct(Fy,0.05), Fy50=pct(Fy,0.50), Fy95=pct(Fy,0.95);
    double Fz5=pct(Fz,0.05), Fz50=pct(Fz,0.50), Fz95=pct(Fz,0.95);
    double Mx5=pct(Mx,0.05), Mx50=pct(Mx,0.50), Mx95=pct(Mx,0.95);
    double My5=pct(My,0.05), My50=pct(My,0.50), My95=pct(My,0.95);
    double Mz5=pct(Mz,0.05), Mz50=pct(Mz,0.50), Mz95=pct(Mz,0.95);
    std::cout << "UQ F.x P5/50/95 = "<<Fx5<<", "<<Fx50<<", "<<Fx95<<"\n";
    std::cout << "UQ F.y P5/50/95 = "<<Fy5<<", "<<Fy50<<", "<<Fy95<<"\n";
    std::cout << "UQ F.z P5/50/95 = "<<Fz5<<", "<<Fz50<<", "<<Fz95<<"\n";
    std::cout << "UQ M.x P5/50/95 = "<<Mx5<<", "<<Mx50<<", "<<Mx95<<"\n";
    std::cout << "UQ M.y P5/50/95 = "<<My5<<", "<<My50<<", "<<My95<<"\n";
    std::cout << "UQ M.z P5/50/95 = "<<Mz5<<", "<<Mz50<<", "<<Mz95<<"\n";
  }
  if (!out_path.empty()) {
    std::ofstream of(out_path);
    of << "{\n";
    of << "  \"F\": [" << out.F.x << ", " << out.F.y << ", " << out.F.z << "],\n";
    of << "  \"M\": [" << out.M.x << ", " << out.M.y << ", " << out.M.z << "],\n";
    of << "  \"diagnostics\": {\n";
    of << "    \"facets\": " << in.facets.size() << ", \"front\": " << front << ", \"shadowed\": " << occluded << ",\n";
    of << "    \"alt_km\": " << cfg.alt_km << ", \"T_K\": " << in.T_K << ", \"Tw_K\": " << cfg.Tw_K << ", \"tau\": " << tau << ",\n";
    of << "    \"Ma_min\": " << Ma_min << ", \"Ma_max\": " << Ma_max << "\n";
    of << "  }\n";
    of << "}\n";
  }
  return 0;
}
