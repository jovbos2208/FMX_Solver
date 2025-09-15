#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "gsi/Sentman.hpp"
#include "gsi/CLL.hpp"

static void usage() {
  std::cout << "Usage: gen_gsi_table --model Sentman|CLL --out table.csv\n"
               "       [--theta_deg_start 0] [--theta_deg_end 90] [--theta_deg_step 2]\n"
               "       [--Ma 0.5,1,2,4,8,12,16] [--tau 0.3,0.5,1,1.5,2]\n"
               "       [--alpha_n 0,0.25,0.5,0.75,1] [--alpha_t 0,0.25,0.5,0.75,1]\n";
}

static std::vector<double> parse_list(const std::string& s) {
  std::vector<double> v; std::string tok; for (size_t i=0,j=0; i<=s.size(); ++i) {
    if (i==s.size() || s[i]==',') { tok = s.substr(j, i-j); try{ v.push_back(std::stod(tok)); }catch(...){} j=i+1; }
  } return v;
}

int main(int argc, char** argv) {
  std::string model = "CLL";
  std::string out_path;
  double th0=0.0, th1=90.0, thstep=2.0;
  std::vector<double> ax_Ma{0.5,1,2,4,8,12,16};
  std::vector<double> ax_tau{0.3,0.5,1.0,1.5,2.0};
  std::vector<double> ax_an{0.0,0.25,0.5,0.75,1.0};
  std::vector<double> ax_at{0.0,0.25,0.5,0.75,1.0};

  for (int i=1;i<argc;++i) {
    std::string a=argv[i];
    if (a=="--model" && i+1<argc) model=argv[++i];
    else if (a=="--out" && i+1<argc) out_path=argv[++i];
    else if (a=="--theta_deg_start" && i+1<argc) th0=std::stod(argv[++i]);
    else if (a=="--theta_deg_end" && i+1<argc) th1=std::stod(argv[++i]);
    else if (a=="--theta_deg_step" && i+1<argc) thstep=std::stod(argv[++i]);
    else if (a=="--Ma" && i+1<argc) ax_Ma=parse_list(argv[++i]);
    else if (a=="--tau" && i+1<argc) ax_tau=parse_list(argv[++i]);
    else if (a=="--alpha_n" && i+1<argc) ax_an=parse_list(argv[++i]);
    else if (a=="--alpha_t" && i+1<argc) ax_at=parse_list(argv[++i]);
    else if (a=="--help") { usage(); return 0; }
  }
  if (out_path.empty()) { usage(); std::cerr << "--out is required\n"; return 1; }
  if (thstep<=0 || th1<th0) { std::cerr << "Invalid theta range/step\n"; return 1; }

  // Build theta axis
  std::vector<double> ax_th; for (double d=th0; d<=th1+1e-9; d+=thstep) ax_th.push_back(d * M_PI/180.0);
  size_t NT=ax_th.size(), NM=ax_Ma.size(), NK=ax_tau.size(), NA=ax_an.size(), NB=ax_at.size();
  std::ofstream out(out_path);
  if (!out) { std::cerr << "Failed to open output: "<<out_path<<"\n"; return 1; }

  // Header
  out << "dims:" << NT << "," << NM << "," << NK << "," << NA << "," << NB << "\n";
  out << "theta:"; for (size_t i=0;i<NT;i++) { if(i) out<<","; out<<ax_th[i]; } out<<"\n";
  out << "Ma:";    for (size_t i=0;i<NM;i++) { if(i) out<<","; out<<ax_Ma[i]; } out<<"\n";
  out << "tau:";   for (size_t i=0;i<NK;i++) { if(i) out<<","; out<<ax_tau[i]; } out<<"\n";
  out << "alpha_n:";for (size_t i=0;i<NA;i++) { if(i) out<<","; out<<ax_an[i]; } out<<"\n";
  out << "alpha_t:";for (size_t i=0;i<NB;i++) { if(i) out<<","; out<<ax_at[i]; } out<<"\n";

  // Loop grid and write CN,CT rows
  for (size_t ib=0; ib<NB; ++ib) {
    for (size_t ia=0; ia<NA; ++ia) {
      for (size_t ik=0; ik<NK; ++ik) {
        for (size_t im=0; im<NM; ++im) {
          for (size_t it=0; it<NT; ++it) {
            double theta = ax_th[it];
            double Ma    = ax_Ma[im];
            double tau   = ax_tau[ik];
            double an    = ax_an[ia];
            double at    = ax_at[ib];
            double CN=0.0, CT=0.0;
            if (model == "Sentman") {
              std::tie(CN, CT) = fmx::gsi::coefficients(theta, Ma, tau, fmx::gsi::SentmanParams{1.0});
            } else {
              std::tie(CN, CT) = fmx::gsi::coefficients(theta, Ma, tau, fmx::gsi::CLLParams{an, at});
            }
            out << CN << "," << CT << "\n";
          }
        }
      }
    }
  }
  std::cerr << "Wrote grid: "<< NT<<"x"<<NM<<"x"<<NK<<"x"<<NA<<"x"<<NB<<" to "<<out_path<<"\n";
  return 0;
}

