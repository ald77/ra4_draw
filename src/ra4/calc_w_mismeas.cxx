#include "ra4/calc_w_mismeas.hpp"

#include <random>
#include <array>
#include <iostream>
#include <string>
#include <limits>
#include <numeric>
#include <algorithm>

#include <unistd.h>
#include <getopt.h>

#include "TError.h"

#include "core/table.hpp"
#include "core/table_row.hpp"
#include "core/process.hpp"
#include "core/gamma_params.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/config_parser.hpp"

using namespace std;

namespace{
  string mismeas_scenario = "mc_as_data";
  NamedFunc reweight_cut = "mt>140. && ntruleps<=1 && mt_tru<=140 && !(type==5000 || type==13000 || type==15000 || type==16000)";

  string output_file = "txt/sys_weights.cfg";
  
  double luminosity = 35.;
  size_t num_toys = 10000000;

  bool debug = false;
}

namespace{
  mt19937_64 InitializePRNG(){
    array<int, 128> sd;
    random_device r;
    generate_n(sd.begin(), sd.size(), ref(r));
    seed_seq ss(begin(sd), end(sd));
    return mt19937_64(ss);
  }

  mt19937_64 prng = InitializePRNG();
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000;
  GetOptions(argc, argv);

  string base_folder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    base_folder = "/net/cms2"; // In laptops, you can't create a /net folder

  string mc_folder(base_folder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_met100_stdnj5/");
  string ntupletag = "met100";
  set<string> mc_files = {
    mc_folder+"*_TTJets*Lept*"+ntupletag+"*.root", mc_folder+"*_TTJets_HT*"+ntupletag+"*.root",
    mc_folder+"*_WJetsToLNu*"+ntupletag+"*.root",mc_folder+"*_ST_*"+ntupletag+"*.root",
    mc_folder+"*_TTW*"+ntupletag+"*.root",mc_folder+"*_TTZ*"+ntupletag+"*.root",
    mc_folder+"*_TTGJets*"+ntupletag+"*.root",mc_folder+"*_TTTT*"+ntupletag+"*.root",
    mc_folder+"*QCD_HT*Inf_Tune*"+ntupletag+"*.root", mc_folder+"*QCD_HT*0_Tune*"+ntupletag+"*.root",
    mc_folder+"*DYJetsToLL*"+ntupletag+"*.root",
    mc_folder+"*_ZJet*"+ntupletag+"*.root",mc_folder+"*_ttHJetTobb*"+ntupletag+"*.root",
    mc_folder+"*_WH_HToBB*"+ntupletag+"*.root",mc_folder+"*_ZH_HToBB*"+ntupletag+"*.root",
    mc_folder+"*_WWTo*"+ntupletag+"*.root",mc_folder+"*_WZ*"+ntupletag+"*.root",mc_folder+"*_ZZ_*"+ntupletag+"*.root"
  };

  string data_folder(base_folder+"/cms2r0/babymaker/babies/2016_08_10/data/unskimmed/");
  set<string> data_files = {
    data_folder+"*.root"
  };

  NamedFunc baseline = "pass && stitch && mj14>250 && nleps>=1 && st>500 && met>100 && met<=150. && njets>=6 && weight<=1";
  
  auto data = Process::MakeShared<Baby_full>("data", Process::Type::data, kBlack,
                                             data_files, baseline);
  auto mc_all = Process::MakeShared<Baby_full>("mc_all", Process::Type::data, kBlack,
                                               mc_files, baseline);
  auto mc_good = Process::MakeShared<Baby_full>("mc_good", Process::Type::background, kBlack,
                                             mc_files, baseline && !reweight_cut);
  auto mc_bad = Process::MakeShared<Baby_full>("mc_bad", Process::Type::background, kBlack,
                                               mc_files, baseline && reweight_cut);
  auto pseudo_data = UseData() ? data : mc_all;
  vector<shared_ptr<Process> > mc_procs = {mc_good, mc_bad}, data_procs = {pseudo_data};

  NamedFunc r1 = "mj14<=400. && mt<=140.";
  NamedFunc r2 = "mj14>400. && mt<=140. && nbm>=1";
  NamedFunc r3 = "mj14<=400. && mt>140.";
  NamedFunc r4 = "mj14>400. && mt>140. && nbm>=1";

  PlotMaker pm;
  NamedFunc weight = "weight";
  Table & mc_table = pm.Push<Table>("mc_yields", vector<TableRow>{
      TableRow("R1", r1, 0, 0, weight),
        TableRow("R2", r2, 0, 0, weight),
        TableRow("R3", r3, 0, 0, weight),
        TableRow("R4", r4, 0, 0, weight)
        }, mc_procs, false, false, false);
  NamedFunc mismeas_wgt = WeightWithCut()*weight;
  Table & data_table = pm.Push<Table>("data_yields", vector<TableRow>{
      TableRow("R1", r1, 0, 0, mismeas_wgt),
        TableRow("R2", r2, 0, 0, mismeas_wgt),
        TableRow("R3", r3, 0, 0, mismeas_wgt),
        TableRow("R4", r4, 0, 0, mismeas_wgt)
        }, data_procs, false, false, false);

  vector<GammaParams> data_yield, mc_good_yield, mc_bad_yield;
  if(!debug){
    pm.multithreaded_ = false;
    pm.MakePlots(luminosity);
    
    data_yield = data_table.Yield(pseudo_data.get(), UseData() ? 1. : luminosity);
    mc_good_yield = mc_table.Yield(mc_good.get(), luminosity);
    mc_bad_yield = mc_table.Yield(mc_bad.get(), luminosity);
  }else{
    data_yield.resize(4);
    data_yield.at(0).SetNEffectiveAndWeight(95785.4, 0.32526);
    data_yield.at(1).SetNEffectiveAndWeight(53560.1, 0.170445);
    data_yield.at(2).SetNEffectiveAndWeight(2158.15, 1.06504);
    data_yield.at(3).SetNEffectiveAndWeight(2724.64, 0.375241);
    mc_good_yield.resize(4);
    mc_good_yield.at(0).SetNEffectiveAndWeight(24809.3, 0.101977);
    mc_good_yield.at(1).SetNEffectiveAndWeight(13684.4, 0.0604334);
    mc_good_yield.at(2).SetNEffectiveAndWeight(12251.9, 0.0985605);
    mc_good_yield.at(3).SetNEffectiveAndWeight(6522.73, 0.0739027);
    mc_bad_yield.resize(4);
    mc_bad_yield.at(0).SetNEffectiveAndWeight(82972.8, 0.344995);
    mc_bad_yield.at(1).SetNEffectiveAndWeight(45765.6, 0.181403);
    mc_bad_yield.at(2).SetNEffectiveAndWeight(511.043, 2.1348);
    mc_bad_yield.at(3).SetNEffectiveAndWeight(838.962, 0.644067);
  }

  double omega_lo_b, omega_mid_b, omega_hi_b;
  BayesianOmega(omega_lo_b, omega_mid_b, omega_hi_b,
                data_yield.at(0), data_yield.at(1), data_yield.at(2), data_yield.at(3),
                mc_good_yield.at(0), mc_good_yield.at(1), mc_good_yield.at(2), mc_good_yield.at(3),
                mc_bad_yield.at(0), mc_bad_yield.at(1), mc_bad_yield.at(2), mc_bad_yield.at(3));
  PrintYield("Data", data_yield);
  PrintYield("MC Good", mc_good_yield);
  PrintYield("MC Bad", mc_bad_yield);

  ConfigParser cp;
  cp.SetOpt("mismeas_cut", MismeasurementCut());
  cp.SetOpt("mismeas_wgt", MismeasurementWeight());
  cp.SetOpt("reweight_cut", reweight_cut);
  cp.SetOpt("w_down", exp(omega_lo_b));
  cp.SetOpt("w_central", exp(omega_mid_b));
  cp.SetOpt("w_up", exp(omega_hi_b));
  if(debug){
    cout << "\nDEBUG MODE: NOT WRITING TO FILE:\n" << endl;
    cout << "Results: " << cp << endl;
  }else{
    cp.Save(output_file, mismeas_scenario);
    cout << "Wrote configuration " << mismeas_scenario << " to file " << output_file << endl;
  }
}

void PrintYield(const string &name, const vector<GammaParams> &yields){
  double kappa = yields.at(0).Yield()*yields.at(3).Yield()/(yields.at(1).Yield()*yields.at(2).Yield());
  cout << name << ": R1=" << yields.at(0)
       << ", R2=" << yields.at(1)
       << ", R3=" << yields.at(2)
       << ", R4=" << yields.at(3)
       << ", kappa=" << kappa << endl;
}

bool UseData(){
  return mismeas_scenario == "data";
}

NamedFunc MismeasurementCut(){
  if(mismeas_scenario == "data"){
    return 0.;
  }else if(mismeas_scenario == "mc_as_data"){
    return 0.;
  }else if(mismeas_scenario == "met_res"){
    return "(met-met_tru)/met>0.5";
  }else if(mismeas_scenario == "xsec_ttw"){
    return "type>=4000 && type<5000";
  }else if(mismeas_scenario == "xsec_ttz"){
    return "type>=5000 && type<6000";
  }else if(mismeas_scenario == "xsec_wjets"){
    return "type>=2000 && type<3000";
  }else if(mismeas_scenario == "xsec_qcd"){
    return "type>=7000 && type<8000";
  }else if(mismeas_scenario == "lep_fake_rate"){
    return "ntruleps==0";
  }else if(mismeas_scenario == "w_isr_1l"){
    return "ntruleps<=1 && type>=1000 && type<2000";
  }else if(mismeas_scenario == "isr_pt"){
    return 1.;
  }else if(mismeas_scenario == "st_model"){
    return 1.;
  }else if(mismeas_scenario == "mismeas_kappa"){
    return "mt>140. && mj14>400. && ntruleps<=1 && mt_tru<=140 && !(type==5000 || type==13000 || type==15000 || type==16000)";
  }else if(mismeas_scenario == "mismeas_mt"){
    return "mt>140. && ntruleps<=1 && mt_tru<=140 && !(type==5000 || type==13000 || type==15000 || type==16000)";
  }else{
    ERROR("Unrecognized mismeasurement scenario: "+mismeas_scenario);
    return 0.;
  }
}

NamedFunc MismeasurementWeight(){
  if(mismeas_scenario == "data"){
    return 1.;
  }else if(mismeas_scenario == "mc_as_data"){
    return 1.;
  }else if(mismeas_scenario == "met_res"){
    return 2.;
  }else if(mismeas_scenario == "xsec_ttw"){
    return 3.;
  }else if(mismeas_scenario == "xsec_ttz"){
    return 3.;
  }else if(mismeas_scenario == "xsec_wjets"){
    return 3.;
  }else if(mismeas_scenario == "xsec_qcd"){
    return 4.;
  }else if(mismeas_scenario == "lep_fake_rate"){
    return 2.;
  }else if(mismeas_scenario == "w_isr_1l"){
    return "1./w_isr";
  }else if(mismeas_scenario == "isr_pt"){
    return "(isr_tru_pt<=600.) + (isr_tru_pt>600. && isr_tru_pt<=800.)*0.5 + (isr_tru_pt>800.)*0.25";
  }else if(mismeas_scenario == "st_model"){
    return "(st<=500.) + (st>500.&&st<=650.)*0.7 + (st>650.&&st<=800.)*0.5 + (st>800.&&st<=1000.)*0.25 + (st>1000.)*0.1";
  }else if(mismeas_scenario == "mismeas_kappa"){
    return 2.;
  }else if(mismeas_scenario == "mismeas_mt"){
    return 2.;
  }else{
    ERROR("Unrecognized mismeasurement scenario: "+mismeas_scenario);
    return 1.;
  }
}

NamedFunc WeightWithCut(){
  NamedFunc cut = MismeasurementCut();
  NamedFunc wgt = MismeasurementWeight();
  string name = "w_mismeas";
  if(cut.IsScalar()){
    if(wgt.IsScalar()){
      return NamedFunc(name, [cut, wgt](const Baby &b){
          return cut.GetScalar(b) ? wgt.GetScalar(b) : 1.;
        });
    }else{
      return NamedFunc(name, [cut, wgt](const Baby &b){
          NamedFunc::VectorType wgts = wgt.GetVector(b);
          return cut.GetScalar(b) ? wgts : NamedFunc::VectorType(wgts.size(), 1.);
        });
    }
  }else{
    if(wgt.IsScalar()){
      return NamedFunc(name, [cut, wgt](const Baby &b){
          NamedFunc::VectorType cuts = cut.GetVector(b);
          NamedFunc::VectorType wgts(cuts.size());
          NamedFunc::ScalarType scalar_wgt = wgt.GetScalar(b);
          for(size_t i = 0; i < wgts.size(); ++i){
            wgts.at(i) = cuts.at(i) ? scalar_wgt : 1.;
          }
          return wgts;
        });
    }else{
      return NamedFunc(name, [cut, wgt](const Baby &b){
          NamedFunc::VectorType cuts = cut.GetVector(b);
          NamedFunc::VectorType wgts = wgt.GetVector(b);
          size_t min_size = min(cuts.size(), wgts.size());
          NamedFunc::VectorType out(min_size);
          for(size_t i = 0; i < min_size; ++i){
            out.at(i) = cuts.at(i) ? wgts.at(i) : 1.;
          }
          return out;
        });
    }
  }
}

void BayesianOmega(double &omega_lo, double &omega_mid, double &omega_hi,
                   const GammaParams &a_data, const GammaParams &b_data,
                   const GammaParams &c_data, const GammaParams &d_data,
                   const GammaParams &a_mc_good, const GammaParams &b_mc_good,
                   const GammaParams &c_mc_good, const GammaParams &d_mc_good,
                   const GammaParams &a_mc_bad, const GammaParams &b_mc_bad,
                   const GammaParams &c_mc_bad, const GammaParams &d_mc_bad){
  gamma_distribution<double> ga_data(a_data.Yield(), 1.);
  gamma_distribution<double> gb_data(b_data.Yield(), 1.);
  gamma_distribution<double> gc_data(c_data.Yield(), 1.);
  gamma_distribution<double> gd_data(d_data.Yield(), 1.);
  gamma_distribution<double> ga_mc_good(a_mc_good.NEffective(), a_mc_good.Weight());
  gamma_distribution<double> gb_mc_good(b_mc_good.NEffective(), b_mc_good.Weight());
  gamma_distribution<double> gc_mc_good(c_mc_good.NEffective(), c_mc_good.Weight());
  gamma_distribution<double> gd_mc_good(d_mc_good.NEffective(), d_mc_good.Weight());
  gamma_distribution<double> ga_mc_bad(a_mc_bad.NEffective(), a_mc_bad.Weight());
  gamma_distribution<double> gb_mc_bad(b_mc_bad.NEffective(), b_mc_bad.Weight());
  gamma_distribution<double> gc_mc_bad(c_mc_bad.NEffective(), c_mc_bad.Weight());
  gamma_distribution<double> gd_mc_bad(d_mc_bad.NEffective(), d_mc_bad.Weight());

  vector<double> omega(num_toys);

  for(size_t i = 0; i < num_toys; ++i){
    double mua_data = ga_data(prng), mub_data = gb_data(prng),
      muc_data = gc_data(prng), mud_data = gd_data(prng);
    double mua_mc_good = ga_mc_good(prng), mub_mc_good = gb_mc_good(prng),
      muc_mc_good = gc_mc_good(prng), mud_mc_good = gd_mc_good(prng);
    double mua_mc_bad = ga_mc_bad(prng), mub_mc_bad = gb_mc_bad(prng),
      muc_mc_bad = gc_mc_bad(prng), mud_mc_bad = gd_mc_bad(prng);
    omega.at(i) = GetOmega(mua_data, mub_data, muc_data, mud_data,
                           mua_mc_good, mub_mc_good, muc_mc_good, mud_mc_good,
                           mua_mc_bad, mub_mc_bad, muc_mc_bad, mud_mc_bad);
    if(fabs(omega.at(i))>1. && false){
      DBG(i << " " << omega.at(i));
      DBG("a_data: " << mua_data << " " << a_data);
      DBG("b_data: " << mub_data << " " << b_data);
      DBG("c_data: " << muc_data << " " << c_data);
      DBG("d_data: " << mud_data << " " << d_data);
      DBG("a_mc_good: " << mua_mc_good << " " << a_mc_good);
      DBG("b_mc_good: " << mub_mc_good << " " << b_mc_good);
      DBG("c_mc_good: " << muc_mc_good << " " << c_mc_good);
      DBG("d_mc_good: " << mud_mc_good << " " << d_mc_good);
      DBG("a_mc_bad: " << mua_mc_bad << " " << a_mc_bad);
      DBG("b_mc_bad: " << mub_mc_bad << " " << b_mc_bad);
      DBG("c_mc_bad: " << muc_mc_bad << " " << c_mc_bad);
      DBG("d_mc_bad: " << mud_mc_bad << " " << d_mc_bad);
    }
  }

  omega_mid = GetOmega(a_data.Yield(), b_data.Yield(), c_data.Yield(), d_data.Yield(),
		       a_mc_good.Yield(), b_mc_good.Yield(), c_mc_good.Yield(), d_mc_good.Yield(),
		       a_mc_bad.Yield(), b_mc_bad.Yield(), c_mc_bad.Yield(), d_mc_bad.Yield());

  sort(omega.begin(), omega.end());
  GetInterval(omega_lo, omega_hi, omega);
}

double GetOmega(double mua_data, double mub_data, double muc_data, double mud_data,
                double mua_mc_good, double mub_mc_good, double muc_mc_good, double mud_mc_good,
                double mua_mc_bad, double mub_mc_bad, double muc_mc_bad, double mud_mc_bad){
  double muc_x = mua_data*mud_data*(mub_mc_good+mub_mc_bad);
  double mud_x = mub_data*muc_data*(mua_mc_good+mua_mc_bad);

  double w_denom = mud_x*mud_mc_bad-muc_x*muc_mc_bad;
  if(w_denom == 0.){DBG("denom==0"); return 0.;}
  double w_numer = muc_x*muc_mc_good-mud_x*mud_mc_good;
  double w = w_numer/w_denom;

  if(w <= 0.){
    return -numeric_limits<double>::infinity();
  }else{
    return log(w);
  }
}

void GetInterval(double &x_lo, double &x_hi, const vector<double> &x_vals){
  double min_delta_x = numeric_limits<double>::infinity();
  size_t delta_i = ceil(erf(sqrt(0.5))*x_vals.size())-1;
  for(size_t i_lo = 0; i_lo + delta_i < x_vals.size(); ++i_lo){
    double this_x_lo = x_vals.at(i_lo);
    double this_x_hi = x_vals.at(i_lo+delta_i);
    double this_delta_x = this_x_hi-this_x_lo;
    if(this_delta_x < min_delta_x){
      min_delta_x = this_delta_x;
      x_lo = this_x_lo;
      x_hi = this_x_hi;
    }
  }
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"mismeas_scenario", required_argument, 0, 'm'},
      {"reweight", required_argument, 0, 'r'},
      {"file", required_argument, 0, 'f'},
      {"luminosity", required_argument, 0, 'l'},
      {"toys", required_argument, 0, 't'},
      {"debug", no_argument, 0, 0},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "m:r:f:l:t:", long_options, &option_index);

    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'm':
      mismeas_scenario = optarg;
      break;
    case 'r':
      reweight_cut = optarg;
      break;
    case 'f':
      output_file = optarg;
      break;
    case 'l':
      luminosity = atof(optarg);
      break;
    case 't':
      num_toys = atoi(optarg);
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "debug"){
        debug = true;
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
