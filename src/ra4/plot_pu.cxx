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
  string dir_mc_isr = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_isrmc_zcand/";
  auto ttbar = Process::MakeShared<Baby_full>("t#bar{t}", Process::Type::background, colors("tt_1l"),
    {dir_mc_isr+"*_TTJets_SingleLeptFromT_Tune*.root", dir_mc_isr+"*_TTJets_SingleLeptFromTbar_Tune*.root",
     dir_mc_isr+"*_TTJets_DiLept_Tune*.root"}, "1");
  auto single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {dir_mc_isr+"*_ST_*.root"});
  auto dyjets = Process::MakeShared<Baby_full>("DY+jets", Process::Type::background, kOrange+1,
    {dir_mc_isr+"*DYJetsToLL_M-50_*.root"},"stitch"); // Inclusive + HT-binned DY
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, kTeal-8,
    {dir_mc_isr+"*_TTWJets*.root", dir_mc_isr+"*_TTZTo*.root", dir_mc_isr+"*_TTGJets*.root"});
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {dir_mc_isr+"*_ZJet*.root", dir_mc_isr+"*_WJetsToLNu*.root",
    dir_mc_isr+"*QCD_HT*0_Tune*.root", dir_mc_isr+"*QCD_HT*Inf_Tune*.root",
    dir_mc_isr+"*ggZH_HToBB*.root", dir_mc_isr+"*_ttHJetTobb*.root",
    dir_mc_isr+"*_TTTT*.root", dir_mc_isr+"*_WWTo*.root",
    dir_mc_isr+"*_WH_HToBB*.root",dir_mc_isr+"*_ZH_HToBB*.root",
    dir_mc_isr+"*_WZ*.root",dir_mc_isr+"*_ZZ_*.root"},"stitch");

  string dir_data_isr = bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_isrdata_zcand/";
  dir_data_isr += quick ? "mergedbaby__0*root" : "*root";
  string lumi_label = RoundNumber(lumi,1).Data();
  auto data = Process::MakeShared<Baby_full>("Data "+lumi_label+" fb^{-1}", Process::Type::data, kBlack,
    {dir_data_isr}, "pass" && Higfuncs::trig_hig>0.);

  vector<shared_ptr<Process> > procs;
  procs = {data, dyjets, ttbar, single_t, ttv, other};
  if (quick) {
    dyjets = Process::MakeShared<Baby_full>("DY+jets", Process::Type::background, kOrange+1,
                        {dir_mc_isr+"*TTJets_Tune*.root"},"1"); // Inclusive + HT-binned DY
    procs = {data, dyjets};
  }


  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .Bottom(BottomType::off)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  vector<PlotOpt> plot_types = {log_lumi, lin_lumi};
  vector<PlotOpt> plot_vals = {lin_lumi().PrintVals(true)};
  PlotMaker pm;

  NamedFunc w_pu_up("w_pu_up", [](const Baby &b) -> NamedFunc::ScalarType{
      if(b.type()<1000) return b.weight();
      else return b.weight()*b.sys_pu()->at(0);
    });
  NamedFunc w_pu_nom("w_pu_nom", [](const Baby &b) -> NamedFunc::ScalarType{
      if(b.type()<1000) return b.weight();
      else return b.weight()*b.w_pu();
    });
  NamedFunc w_pu_dn("w_pu_dn", [](const Baby &b) -> NamedFunc::ScalarType{
      if(b.type()<1000) return b.weight();
      else return b.weight()*b.sys_pu()->at(1);
    });

  pm.Push<Hist1D>(Axis(60,0,60, "npv", "NPV (no pu wgt)"), 
    "1", procs, plot_types).Weight("weight").Tag("pu");
  pm.Push<Hist1D>(Axis(60,0,60, "npv", "NPV (pu nom)"), 
    "1", procs, plot_types).Weight(w_pu_nom).Tag("pu");
  pm.Push<Hist1D>(Axis(60,0,60, "npv", "NPV (pu up)"), 
    "1", procs, plot_types).Weight(w_pu_up).Tag("pu");
  pm.Push<Hist1D>(Axis(60,0,60, "npv", "NPV (pu down)"), 
    "1", procs, plot_types).Weight(w_pu_dn).Tag("pu");

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
      {"quick", no_argument, 0, 0},           // Used inclusive ttbar for quick testing
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
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