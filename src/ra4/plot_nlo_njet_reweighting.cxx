#include <cmath>
#include <algorithm>
#include <getopt.h>

#include "TError.h"
#include "TVector2.h"
#include "TString.h"

#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist1d.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"
#include "hig/hig_functions.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace {
  string isrtype = "zcand";
  // string isrtype = "ttisr";

  bool single_thread = false;
  bool quick = false; 
  float lumi = 35.9;
}

void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  Palette colors("txt/colors.txt", "default");

  //// Processes for ISR skims
  string dir_mc_isr = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_minisrmc_"+isrtype+"/";
  string dir_mc_isr_nlo = bfolder+"/cms2r0/babymaker/babies/2017_04_01/mc/merged_minisrmc_"+isrtype+"/";
  auto ttbar = Process::MakeShared<Baby_full>("t#bar{t}", Process::Type::background, colors("tt_2l"),
    {dir_mc_isr_nlo+"*_TT_*.root"});
  auto single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {dir_mc_isr+"*_ST_*.root"});
  auto dyjets = Process::MakeShared<Baby_full>("DY+jets", Process::Type::background, kRed-3,
    {dir_mc_isr_nlo+"*DYJetsToLL_M-50*.root"});//(ptll_me<100||type!=6200)"); 
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, kTeal-8,
    {dir_mc_isr+"*_TTWJets*.root", dir_mc_isr+"*_TTZTo*.root", dir_mc_isr+"*_TTGJets*.root"});
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {dir_mc_isr+"*_ZJet*.root", dir_mc_isr+"*_WJetsToLNu*.root",
    dir_mc_isr+"*QCD_HT*0_Tune*.root", dir_mc_isr+"*QCD_HT*Inf_Tune*.root",
    dir_mc_isr+"*ggZH_HToBB*.root", dir_mc_isr+"*_ttHJetTobb*.root",
    dir_mc_isr+"*_TTTT*.root", dir_mc_isr+"*_WWTo*.root",
    dir_mc_isr+"*_WH_HToBB*.root",dir_mc_isr+"*_ZH_HToBB*.root",
    dir_mc_isr+"*_WZ*.root",dir_mc_isr+"*_ZZ_*.root"},"stitch");

  string dir_data_isr = bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_minisrdata_"+isrtype+"/*root";
  string lumi_label = RoundNumber(lumi,1).Data();
  auto data = Process::MakeShared<Baby_full>("Data "+lumi_label+" fb^{-1}", Process::Type::data, kBlack,
    {dir_data_isr}, "pass" && Higfuncs::trig_hig>0.);

  vector<shared_ptr<Process> > procs = {data, dyjets, ttbar, single_t, ttv, other};
  if (isrtype=="ttisr") procs = {data, ttbar, dyjets, single_t, ttv, other};

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  vector<PlotOpt> plot_log = {log_lumi};
  vector<PlotOpt> plot_lin = {lin_lumi};
  vector<PlotOpt> plot_vals = {log_lumi().PrintVals(true)};
  PlotMaker pm;


  NamedFunc zcand("zcand", [](const Baby &b) -> NamedFunc::ScalarType{
    if (b.nels()==2) {
      if (b.elel_m()>80 && b.elel_m()<100) return 1;
    } else if (b.nmus()==2) {
      if (b.mumu_m()>80 && b.mumu_m()<100) return 1;
    } 
    return -1;
  });


  NamedFunc baseline = "nleps==2&&leps_pt[0]>40&&((elel_m>80&&elel_m<100)||(mumu_m>80&&mumu_m<100))";
  // if (isrtype=="ttisr") baseline = "nleps==2&&leps_pt[0]>40&&nbm==2";
  if (isrtype=="ttisr") baseline = "nleps==2&&leps_pt[0]>40&&nbm==2" && zcand<0.;

  const vector<double> isr_syspt_bins = {0, 50, 100, 150, 200, 300, 400, 600, 800, 1000};
  const vector<double> nisrjet_bins = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5,7.5};

  vector<NamedFunc> weight_opts;
  // weight_opts.push_back("w_lumi");
  weight_opts.push_back(NamedFunc("nom", 
    [](const Baby &b) -> NamedFunc::ScalarType{
      if (b.w_isr()>0) return b.weight()/b.w_isr(); 
      else return b.weight();
    }) * Higfuncs::eff_higtrig);

  for (const auto &iweight: weight_opts){
    if (isrtype=="zcand") {
      pm.Push<Hist1D>(Axis(nisrjet_bins, "njets", "ISR jet multiplicity"), 
        baseline, procs, plot_vals).Weight(iweight).Tag(isrtype);
      pm.Push<Hist1D>(Axis(25,0,1250, "jetsys_pt", "ISR p_{T} [GeV]"), 
        baseline, procs, plot_log).Weight(iweight).Tag(isrtype);
    } else {
      pm.Push<Hist1D>(Axis(nisrjet_bins, "njets-2", "ISR jet multiplicity"), 
        baseline, procs, plot_vals).Weight(iweight).Tag(isrtype);
      pm.Push<Hist1D>(Axis(25,0,1250, "jetsys_nob_pt", "ISR p_{T} [GeV]"), 
        baseline, procs, plot_log).Weight(iweight).Tag(isrtype);
    }
    pm.Push<Hist1D>(Axis(20,0.,2000., "ht", "H_{T} [GeV]"), 
      baseline, procs, plot_log).Weight(iweight).Tag(isrtype);
  } // Loop over weights

  if (single_thread) pm.multithreaded_ = false;
  pm.min_print_ = true;
  pm.MakePlots(lumi);
    double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Making plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
} // main

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"type", required_argument, 0, 't'},  // Method to run on (if you just want one)
      {"quick", no_argument, 0, 0},           // Used inclusive ttbar for quick testing
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "t:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 't':
      isrtype = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "quick"){
        quick = true;
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
        exit(1);
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}