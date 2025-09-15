#include "atm/NRLMSIS2.hpp"
#include "core/units.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>

extern "C" {
#if defined(FMX_WITH_NRLMSIS2_LOCAL)
void msis_eval_c(double day, double utsec, double z_km, double lat_deg, double lon_deg,
                 double f107a, double f107, double ap_daily, double ap_now,
                 double& Tn, double& dn_tot, double& dn_n2, double& dn_o2, double& dn_o, double& dn_he,
                 double& dn_h, double& dn_ar, double& dn_n, double& dn_ao, double& dn_no);
#endif
}

namespace fmx::atm {

static bool parse_iso_utc(const std::string& s, int& year, int& mon, int& mday, int& hour, int& min, double& sec) {
  // Expect formats like YYYY-MM-DDThh:mm:ssZ or with fractional seconds
  // Minimal, robust-ish parsing without locale.
  year=0; mon=0; mday=0; hour=0; min=0; sec=0.0;
  if (s.size() < 20) return false;
  if (s[4] != '-' || s[7] != 'T' || (s[10] != ':' && s[13] != ':')) {
    // Try YYYY-MM-DDThh:mm:ssZ exact positions (allow 'Z' optional)
  }
  try {
    year = std::stoi(s.substr(0,4));
    mon  = std::stoi(s.substr(5,2));
    mday = std::stoi(s.substr(8,2));
    hour = std::stoi(s.substr(11,2));
    min  = std::stoi(s.substr(14,2));
    size_t p = 17;
    // seconds may have fractional and optional trailing Z
    size_t end = s.find_first_of("Zz ", p);
    std::string ssec = s.substr(p, end==std::string::npos ? std::string::npos : (end-p));
    sec = std::stod(ssec);
    return true;
  } catch (...) { return false; }
}

static int doy_from_ymd(int y, int m, int d) {
  static const int mdays[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
  auto is_leap = [](int yy){ return (yy%4==0 && (yy%100!=0 || yy%400==0)); };
  int doy = 0;
  for (int i=1;i<m;i++) doy += (i==2 && is_leap(y)) ? 29 : mdays[i-1];
  doy += d;
  return doy;
}

static double kp_to_ap_int(int kp) {
  // Rough mapping for integer Kp
  static const int table[10] = {0,4,7,15,27,48,80,140,240,400};
  if (kp < 0) kp = 0; if (kp > 9) kp = 9;
  return static_cast<double>(table[kp]);
}

AtmosphereState NRLMSIS2Atmosphere::evaluate(double alt_km, double lat_deg, double lon_deg,
                                              const std::string& utc_iso,
                                              const Indices& idx) const {
#if defined(FMX_WITH_NRLMSIS2_LOCAL)
  AtmosphereState st{};
  // Ensure parameter file is present in CWD for MSIS
  {
    std::ifstream test("msis21.parm");
    if (!test.good()) {
      std::ifstream src("atm/models/NRLMSIS2.1/msis21.parm", std::ios::binary);
      if (src.good()) {
        std::ofstream dst("msis21.parm", std::ios::binary);
        dst << src.rdbuf();
      }
    }
  }
  // Convert UTC to day-of-year and UT seconds is left as input provided (UTC seconds of day recommended).
  // Parse utc string; on failure fallback to noon day 100
  int y, mo, md, hh, mm; double ss;
  double day = 100.0;
  double UTsec = 43200.0;
  if (parse_iso_utc(utc_iso, y, mo, md, hh, mm, ss)) {
    int doy = doy_from_ymd(y, mo, md);
    double frac = (hh*3600.0 + mm*60.0 + ss) / 86400.0;
    day = static_cast<double>(doy) + frac;
    UTsec = hh*3600.0 + mm*60.0 + ss;
  }
  double Tn, dn_tot, dn_n2, dn_o2, dn_o, dn_he, dn_h, dn_ar, dn_n, dn_ao, dn_no;
  double f107a = idx.F10_7A > 0 ? idx.F10_7A : idx.F10_7;
  double f107 = idx.F10_7;
  // Build Ap inputs
  double ap7[7] = {0,0,0,0,0,0,0};
  if (idx.has_ap_array) {
    for (int i=0;i<7;i++) ap7[i] = idx.Ap[i];
  } else if (idx.Ap_daily > 0.0 || idx.Ap_now > 0.0) {
    ap7[0] = (idx.Ap_daily > 0.0) ? idx.Ap_daily : kp_to_ap_int(idx.Kp);
    ap7[1] = (idx.Ap_now > 0.0) ? idx.Ap_now : ap7[0];
  } else {
    double ap = kp_to_ap_int(idx.Kp);
    ap7[0] = ap; ap7[1] = ap;
  }
  // Call Fortran wrapper: currently only daily and now are used
  msis_eval_c(day, UTsec, alt_km, lat_deg, lon_deg, f107a, f107, ap7[0], ap7[1],
              Tn, dn_tot, dn_n2, dn_o2, dn_o, dn_he, dn_h, dn_ar, dn_n, dn_ao, dn_no);
  st.T_K = Tn;
  st.wind_ms = {0.0, 0.0, 0.0};
  // Convert number densities [m^-3] to mass densities [kg/m^3]
  auto n2rho = [](double n, double m)->double{ return n * m; };
  st.species.clear();
  st.species.push_back({n2rho(dn_o,  fmx::units::m_O),  fmx::units::m_O});
  st.species.push_back({n2rho(dn_n2, fmx::units::m_N2), fmx::units::m_N2});
  st.species.push_back({n2rho(dn_o2, fmx::units::m_O2), fmx::units::m_O2});
  st.species.push_back({n2rho(dn_he, fmx::units::m_He), fmx::units::m_He});
  st.species.push_back({n2rho(dn_h,  fmx::units::m_H),  fmx::units::m_H});
  return st;
#else
  static bool warned = false;
  if (!warned) {
    std::cerr << "[NRLMSIS2] Library not linked; using placeholder densities.\n";
    warned = true;
  }
  AtmosphereState st{};
  // Simple parametric thermosphere temperature and densities
  double h = std::clamp(alt_km, 100.0, 800.0);
  st.T_K = 700.0 + 0.8 * (h - 100.0); // crude trend
  st.wind_ms = {0.0, 0.0, 0.0};
  double rho_O  = 5e-11 * std::exp(-(h - 200.0)/60.0);
  double rho_He = 2e-12 * std::exp(-(h - 300.0)/80.0);
  double rho_H  = 1e-12 * std::exp(-(h - 400.0)/120.0);
  st.species.push_back({rho_O,  fmx::units::m_O});
  st.species.push_back({rho_He, fmx::units::m_He});
  st.species.push_back({rho_H,  fmx::units::m_H});
  return st;
#endif
}

} // namespace fmx::atm
