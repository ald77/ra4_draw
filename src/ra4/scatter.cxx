#include "ra4/scatter.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TColor.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  bool single_thread = false;
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  double lumi = 35.9;

  string base_path = "";
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname,"compute-")){
    base_path = "/net/cms29";
  }
  string mc_dir = base_path+"/cms29r0/babymaker/babies/2017_01_27/mc/merged_mcbase_standard/";

  Palette colors("txt/colors.txt", "default");

  auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors.RGB(1,57,166),
    {mc_dir+"*_TTJets*Lept*.root"},
    "ntruleps<=1&&stitch_met");
  tt1l->SetMarkerStyle(23);
  tt1l->SetMarkerSize(0.8);
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors.RGB(86,160,211),
    {mc_dir+"*_TTJets*Lept*.root"},
    "ntruleps>=2&&stitch_met");
  tt2l->SetMarkerStyle(22);
  tt2l->SetMarkerSize(0.8);
  auto wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {mc_dir+"*_WJetsToLNu*.root"},"stitch");
  auto single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {mc_dir+"*_ST_*.root"});
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {mc_dir+"*_TTWJets*.root", mc_dir+"*_TTZTo*.root"});
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {mc_dir+"*DYJetsToLL*.root", mc_dir+"*_QCD_HT*.root",
        mc_dir+"*_ZJet*.root", mc_dir+"*_WWTo*.root",
        mc_dir+"*ggZH_HToBB*.root", mc_dir+"*ttHJetTobb*.root",
        mc_dir+"*_TTGJets*.root", mc_dir+"*_TTTT_*.root",
        mc_dir+"*_WH_HToBB*.root", mc_dir+"*_WZTo*.root",
        mc_dir+"*_ZH_HToBB*.root", mc_dir+"*_ZZ_*.root"});

  auto t1tttt = Process::MakeShared<Baby_full>("T1tttt(1800,100)", Process::Type::signal, colors("t1tttt"),
    {base_path+"/cms29r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/*SMS-T1tttt_mGluino-1800_mLSP-100_*.root"});
  t1tttt->SetMarkerStyle(21);
  t1tttt->SetMarkerSize(0.9);

  auto data = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {base_path+"/cms29r0/babymaker/babies/2017_02_14/data/merged_database_stdnj5/*.root"},"pass&&trig_ra4");
  data->SetMarkerStyle(20);
  data->SetMarkerSize(1.);

  vector<shared_ptr<Process> > all_procs = {data, t1tttt, tt1l, tt2l, wjets, single_t, ttv, other};
  vector<shared_ptr<Process> > tt_sig = {tt1l, tt2l, t1tttt};

  PlotOpt style("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> bkg_hist = {style().Stack(StackType::data_norm).Title(TitleType::supplementary)};
  vector<PlotOpt> bkg_pts = {style().Stack(StackType::lumi_shapes).Title(TitleType::simulation_supplementary)};

  NamedFunc baseline = "nleps==1&&st>500&&met>150&&njets>=6&&nbm>=1&&nveto==0 && met/met_calo<5.0 && pass_ra2_badmu";
  NamedFunc weight = "weight";
  vector<NamedFunc> met_bins = {"met>200","met>200&&met<=350", "met>350&&met<=500", "met>500"};
  vector<NamedFunc> nbm_bins = {"nbm==1", "nbm>=2"};

  PlotMaker pm;
  for(const auto &met_bin: met_bins){
    for(const auto &nbm_bin: nbm_bins){
      NamedFunc cut = baseline && met_bin && nbm_bin;
      pm.Push<Hist2D>(Axis(44, 0., 1100., "mj14", "M_{J} [GeV]", {250., 400.}),
                      Axis(25, 0., 700., "mt", "m_{T} [GeV]", {140.}),
                      cut, all_procs, bkg_hist);
      pm.Push<Hist2D>(Axis(48, 0., 1200., "mj14", "M_{J} [GeV]", {250., 400.}),
                      Axis(175, 0., 700., "mt", "m_{T} [GeV]", {140.}),
                      cut, tt_sig, bkg_pts);
    }
  }

  if(single_thread) pm.multithreaded_ = false;
  pm.min_print_=true;
  pm.MakePlots(lumi);
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 's'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s", long_options, &option_index);

    if( opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      single_thread = true;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(false){
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
