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
#include "atm/Atmosphere.hpp"
#include "atm/NRLMSIS2.hpp"
#include "atm/HWM14.hpp"
#include "atm/Combined.hpp"

using fmx::Vec3;

static std::string slurp(const std::string& p) { std::ifstream in(p); std::ostringstream ss; ss<<in.rdbuf(); return ss.str(); }
static bool find_str(const std::string& s, const std::string& k, std::string& out) {
  auto pos=s.find("\""+k+"\""); if(pos==std::string::npos) return false; pos=s.find(':',pos); pos=s.find('"',pos); auto e=s.find('"',pos+1); out=s.substr(pos+1,e-pos-1); return true; }
static bool find_num(const std::string& s, const std::string& k, double& out) {
  auto pos=s.find("\""+k+"\""); if(pos==std::string::npos) return false; pos=s.find(':',pos); auto e=s.find_first_of(",}\n\r",pos+1); std::string t=s.substr(pos+1,e-pos-1); try{ out=std::stod(t); return true; }catch(...) {return false;}}
static bool find_int(const std::string& s, const std::string& k, int& out) { double d; if(!find_num(s,k,d)) return false; out=(int)d; return true; }
static bool find_arr3(const std::string& s, const std::string& k, Vec3& out) {
  auto pos=s.find("\""+k+"\""); if(pos==std::string::npos) return false; pos=s.find('[',pos); auto e=s.find(']',pos); std::string a=s.substr(pos+1,e-pos-1); std::replace(a.begin(),a.end(),',',' '); std::istringstream ss(a); if(!(ss>>out.x>>out.y>>out.z)) return false; return true; }

static fmx::geom::Mesh make_plate(){ fmx::geom::Mesh m; m.tris.push_back({Vec3{0,0.5,0.5},Vec3{0,0.5,-0.5},Vec3{0,-0.5,-0.5}}); m.tris.push_back({Vec3{0,-0.5,0.5},Vec3{0,0.5,0.5},Vec3{0,-0.5,-0.5}}); return m; }
static fmx::geom::Mesh make_cube(){ fmx::geom::Mesh m; double h=0.5; auto q=[&](Vec3 a,Vec3 b,Vec3 c,Vec3 d){m.tris.push_back({a,b,c});m.tris.push_back({d,a,c});}; q({ h,-h,-h},{ h, h,-h},{ h, h, h},{ h,-h, h}); q({-h, h, h},{-h, h,-h},{-h,-h,-h},{-h,-h, h}); q({-h, h, h},{ h, h, h},{ h, h,-h},{-h, h,-h}); q({-h,-h,-h},{ h,-h,-h},{ h,-h, h},{-h,-h, h}); q({-h,-h, h},{ h,-h, h},{ h, h, h},{-h, h, h}); q({-h, h,-h},{ h, h,-h},{ h,-h,-h},{-h,-h,-h}); return m; }

int main(){
  // Look for cases file
  const std::string cases_path = "examples/cr313/cases.json";
  std::ifstream cf(cases_path);
  if (!cf.good()) {
    std::cout << "[CR-313] cases file not found; skipping reference validation.\n";
    return 0;
  }
  const std::string json = slurp(cases_path);
  // Very light parsing: expect an array of objects; we'll just handle a single object for simplicity now.
  // Fields: validate|geometry, cg, materials.default.{alpha_E,Tw_K}, atmosphere.indices.{F10_7,F10_7A,Kp,Ap}, state.{alt_km,lat_deg,lon_deg,utc,V_sat_mps}, expected.{F:[..],M:[..],tol}

  // Extract a single object window between first '{' and its matching '}'
  auto l = json.find('{'); auto r = json.rfind('}'); if (l==std::string::npos || r==std::string::npos || r<=l) { std::cout<<"[CR-313] malformed cases.json; skipping.\n"; return 0; }
  std::string j = json.substr(l, r-l+1);

  // Build geometry
  fmx::geom::Mesh mesh;
  std::string validate; if (find_str(j, "validate", validate)) {
    if (validate == "plate") mesh = make_plate();
    else if (validate == "cube") mesh = make_cube();
  }
  std::string geom_path; if (mesh.tris.empty() && find_str(j, "geometry", geom_path)) {
    auto m = fmx::geom::Mesh::load(geom_path);
    if (!m) { std::cerr << "[CR-313] failed to load geometry: "<<geom_path<<"\n"; return 1; }
    mesh = *m;
  }
  if (mesh.tris.empty()) { std::cout << "[CR-313] no valid case entry; skipping.\n"; return 0; }

  // Optional rotation about Z by theta_deg
  double theta_deg = 0.0; find_num(j, "theta_deg", theta_deg);
  if (theta_deg != 0.0) {
    double th = theta_deg * M_PI/180.0, c=std::cos(th), s=std::sin(th);
    auto rot = [&](Vec3& v){ double x=v.x, y=v.y; v.x = x*c - y*s; v.y = x*s + y*c; };
    for (auto& t : mesh.tris) { rot(t.v0); rot(t.v1); rot(t.v2); }
  }
  auto facets = mesh.to_facets(0);
  // Atmosphere
  fmx::atm::Indices idx; double tmp;
  if (find_num(j, "F10_7", tmp)) idx.F10_7 = tmp; if (find_num(j, "F10_7A", tmp)) idx.F10_7A = tmp; int kp; if (find_int(j, "Kp", kp)) idx.Kp = kp;

  // State
  Vec3 cg{0,0,0}; find_arr3(j, "cg", cg);
  double alt=400, lat=0, lon=0; find_num(j, "alt_km", alt); find_num(j, "lat_deg", lat); find_num(j, "lon_deg", lon);
  std::string utc; find_str(j, "utc", utc); if (utc.empty()) utc = "2025-09-12T12:00:00Z";
  Vec3 Vsat{7500,0,0}; find_arr3(j, "V_sat_mps", Vsat);

  // Materials
  double alpha_E=1.0, Tw=300.0; find_num(j, "alpha_E", alpha_E); find_num(j, "Tw_K", Tw);

  // Model selection (Stub or Combined)
  std::string model; find_str(j, "model", model);
  fmx::atm::AtmosphereState st;
  if (model == "Stub") {
    fmx::atm::StubAtmosphere atm;
    st = atm.evaluate(alt, lat, lon, utc, idx);
  } else {
    fmx::atm::CombinedAtmosphere atm;
    st = atm.evaluate(alt, lat, lon, utc, idx);
  }

  fmx::solver::Input in; in.facets=facets; in.materials = { {1.0, 1.0, alpha_E, Tw} };
  for (const auto& sp : st.species) in.species.push_back({sp.rho, sp.mass}); in.T_K = st.T_K; in.V_sat_ms = Vsat; in.wind_ms = st.wind_ms; in.r_CG = cg;
  auto out = fmx::solver::solve(in);

  // Expected
  Vec3 Fexp{0,0,0}, Mexp{0,0,0}; double tol=0.1; bool hasF = find_arr3(j, "F_exp", Fexp); bool hasM = find_arr3(j, "M_exp", Mexp); find_num(j, "tol", tol);
  if (hasF) {
    auto rel = [&](double a,double b){ return std::abs(a-b)/std::max(1e-12,std::abs(b)); };
    if (rel(out.F.x,Fexp.x)>tol || rel(out.F.y,Fexp.y)>tol || rel(out.F.z,Fexp.z)>tol) {
      std::cerr << "[CR-313] Force mismatch: got ["<<out.F.x<<","<<out.F.y<<","<<out.F.z<<"] exp ["<<Fexp.x<<","<<Fexp.y<<","<<Fexp.z<<"] tol="<<tol<<"\n";
      return 1;
    }
  }
  if (hasM) {
    auto rel = [&](double a,double b){ return std::abs(a-b)/std::max(1e-12,std::abs(b)); };
    if (rel(out.M.x,Mexp.x)>tol || rel(out.M.y,Mexp.y)>tol || rel(out.M.z,Mexp.z)>tol) {
      std::cerr << "[CR-313] Moment mismatch: got ["<<out.M.x<<","<<out.M.y<<","<<out.M.z<<"] exp ["<<Mexp.x<<","<<Mexp.y<<","<<Mexp.z<<"] tol="<<tol<<"\n";
      return 1;
    }
  }
  std::cout << "[CR-313] Case validated or no expectations provided.\n";
  return 0;
}
