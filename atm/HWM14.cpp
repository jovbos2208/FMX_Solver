#include "atm/HWM14.hpp"
#include <cmath>
#include <iostream>
#include <fstream>

extern "C" {
#if defined(FMX_WITH_HWM14_LOCAL)
void hwm14_eval_c(int doy, double utsec, double alt_km, double lat_deg, double lon_deg, double ap3hr,
                  double& w_meridional, double& w_zonal);
#endif
}

namespace fmx::atm {

static bool parse_iso_utc_hw(const std::string& s, int& year, int& mon, int& mday, int& hour, int& min, double& sec) {
  year=0;mon=0;mday=0;hour=0;min=0;sec=0.0;
  if (s.size() < 20) return false;
  try {
    year = std::stoi(s.substr(0,4));
    mon  = std::stoi(s.substr(5,2));
    mday = std::stoi(s.substr(8,2));
    hour = std::stoi(s.substr(11,2));
    min  = std::stoi(s.substr(14,2));
    size_t p = 17; size_t end = s.find_first_of("Zz ", p);
    std::string ssec = s.substr(p, end==std::string::npos ? std::string::npos : (end-p));
    sec = std::stod(ssec);
    return true;
  } catch (...) { return false; }
}

static int doy_from_ymd_hw(int y, int m, int d) {
  static const int mdays[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
  auto is_leap = [](int yy){ return (yy%4==0 && (yy%100!=0 || yy%400==0)); };
  int doy = 0; for (int i=1;i<m;i++) doy += (i==2 && is_leap(y)) ? 29 : mdays[i-1];
  return doy + d;
}

fmx::Vec3 HWM14::evaluate_wind(double alt_km, double lat_deg, double lon_deg,
                                const std::string& utc_iso, double ap3hr) {
  (void)lon_deg; (void)utc_iso;
#if defined(FMX_WITH_HWM14_LOCAL)
  // Ensure HWM14 binary climatology present in CWD
  {
    auto ensure_file = [](const char* fname){
      std::ifstream t(fname, std::ios::binary);
      if (t.good()) return;
      std::string srcp = std::string("atm/models/hwm14/data/") + fname;
      std::ifstream src(srcp, std::ios::binary);
      if (src.good()) { std::ofstream dst(fname, std::ios::binary); dst << src.rdbuf(); }
    };
    ensure_file("hwm123114.bin");
    ensure_file("dwm07b104i.dat");
    ensure_file("gd2qd.dat");
  }
  // Parse UTC to get day-of-year and UT seconds
  int y, mo, md, hh, mm; double ss;
  int doy = 100; double utsec = 43200.0;
  if (parse_iso_utc_hw(utc_iso, y, mo, md, hh, mm, ss)) {
    doy = doy_from_ymd_hw(y, mo, md);
    utsec = hh*3600.0 + mm*60.0 + ss;
  }
  double wm=0.0, wz=0.0;
  hwm14_eval_c(doy, utsec, alt_km, lat_deg, lon_deg, ap3hr, wm, wz);
  return {0.0, wz, wm};
#endif
  static bool warned = false;
  if (!warned) {
    std::cerr << "[HWM14] Library not linked; using synthesized wind field.\n";
    warned = true;
  }
  // Simple lat/alt-dependent synthetic wind: zonal winds ~ 50-150 m/s
  double phi = lat_deg * M_PI/180.0;
  double zonal = 100.0 * std::cos(phi) * std::exp(-(alt_km-100.0)/400.0);
  return {0.0, zonal, 0.0};
}

} // namespace fmx::atm
