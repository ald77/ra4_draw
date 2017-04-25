#include "core/test.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TColor.h"
#include "TVector2.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "hig/hig_functions.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

vector<unsigned> higidx(const Baby &b);

namespace{
  bool single_thread = false;
  double lumi = 35.9;
}

int main(int argc, char *argv[]){
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info).LogMinimum(0.01).PrintVals(true);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info);
  vector<PlotOpt> plot_opts = {log_lumi_info};

  string folder_lo = "/cms2r0/babymaker/babies/2017_01_27/mc/merged_bare_dy/";
  string folder_nlo = "/cms2r0/babymaker/babies/2017_04_01/mc/merged_bare_dy/";

  Palette colors("txt/colors.txt", "default");

  vector<shared_ptr<Process> >  procs_all;
  procs_all.push_back(Process::MakeShared<Baby_full>("Inclusive", Process::Type::data, 
    kBlack, {folder_nlo+"*DYJetsToLL_M-50_Tune*.root"}, "1"));
  procs_all.push_back(Process::MakeShared<Baby_full>("Stitch p_{T}(ll)", Process::Type::background, 
    kGreen+1, {folder_nlo+"*DYJetsToLL_M-50_Tune*.root"}, "ptll_me<100"));
  procs_all.push_back(Process::MakeShared<Baby_full>("Stitch p_{T}(ll) 100-250", Process::Type::background, 
    kMagenta+2, {folder_nlo+"*DYJetsToLL_Pt-100To250*.root"}, "1"));
  procs_all.push_back(Process::MakeShared<Baby_full>("Stitch p_{T}(ll) 250-400", Process::Type::background, 
    kAzure+1, {folder_nlo+"*DYJetsToLL_Pt-250To400*.root"}, "1"));
  procs_all.push_back(Process::MakeShared<Baby_full>("Stitch p_{T}(ll) 400-650", Process::Type::background, 
    kOrange, {folder_nlo+"*DYJetsToLL_Pt-400To650*.root"}, "1"));
  procs_all.push_back(Process::MakeShared<Baby_full>("Stitch p_{T}(ll)", Process::Type::background, 
    kRed+1, {folder_nlo+"*DYJetsToLL_Pt-650ToInf*.root"}, "1"));

  vector<shared_ptr<Process> >  procs_ht;
  procs_ht.push_back(Process::MakeShared<Baby_full>("Inclusive", Process::Type::data, 
    kBlack, {folder_lo+"*DYJetsToLL_M-50_Tune*.root"}, "1"));
  
  procs_ht.push_back(Process::MakeShared<Baby_full>("HT 0-70", Process::Type::background, kOrange-3, {folder_lo+"*DYJetsToLL_M-50_Tune*.root"}, "stitch"));
  procs_ht.push_back(Process::MakeShared<Baby_full>("HT 70-100", Process::Type::background, kGreen+1, {folder_lo+"*DYJetsToLL_M-50_HT-70to100_*.root"}, "stitch"));
  procs_ht.push_back(Process::MakeShared<Baby_full>("HT 100-200", Process::Type::background, kCyan-2, {folder_lo+"*DYJetsToLL_M-50_HT-100to200_*.root"}, "stitch"));
  procs_ht.push_back(Process::MakeShared<Baby_full>("HT 200-400", Process::Type::background, kMagenta-2, {folder_lo+"*DYJetsToLL_M-50_HT-200to400_*.root"}, "stitch"));
  procs_ht.push_back(Process::MakeShared<Baby_full>("HT 400-600", Process::Type::background, kAzure+1, {folder_lo+"*DYJetsToLL_M-50_HT-400to600_*.root"}, "stitch"));
  procs_ht.push_back(Process::MakeShared<Baby_full>("HT 600-800", Process::Type::background, kOrange, {folder_lo+"*DYJetsToLL_M-50_HT-600to800_*.root"}, "stitch"));
  procs_ht.push_back(Process::MakeShared<Baby_full>("HT 800-1200", Process::Type::background, kRed+1, {folder_lo+"*DYJetsToLL_M-50_HT-800to1200_*.root"}, "stitch"));
  procs_ht.push_back(Process::MakeShared<Baby_full>("HT 1200-2500", Process::Type::background, kBlue, {folder_lo+"*DYJetsToLL_M-50_HT-1200to2500_*.root",folder_lo+"*DYJetsToLL_M-50_HT-2500toInf_*.root"}, "stitch"));

  PlotMaker pm;
  NamedFunc wlumi_dy_nlo("wlumi_dy_nlo",[&](const Baby &b){
      double wgt_ = b.w_lumi(); 
      // if (b.type()==6202) wgt_*=0.984;
      // else if (b.type()==6203) wgt_*=1.169;
      // else if (b.type()==6204) wgt_*=1.267;
      // else if (b.type()==6205) wgt_*=1.216;
      return wgt_;
    });

  NamedFunc wlumi_dy_lo("wlumi_dy_lo",[&](const Baby &b){
      double wgt_ = b.w_lumi(); 
      // if (b.type()==6100) wgt_*=0.84;
      // else if (b.type()==6101) wgt_*=1.061;
      // else if (b.type()==6102) wgt_*=0.954;
      // else if (b.type()==6103) wgt_*=0.983;
      // else if (b.type()==6104) wgt_*=0.974;
      // else if (b.type()==6105) wgt_*=0.855;
      // else if (b.type()==6106 || b.type()==6107) wgt_*=1.04;
      return wgt_;
    });

  // vector<double> ptbins = {0,100,250,400,650, 800};
  // pm.Push<Hist1D>(Axis(ptbins,"ptll_me", "ME p_{T} (ll) [GeV]", {100, 250, 400, 650}),
  //   "1", procs_all, plot_opts).Tag("stitch").Weight(wlumi_dy_nlo);
  vector<double> ptbins_fine = {0, 20, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 575, 650, 725, 875, 1000};
  pm.Push<Hist1D>(Axis(ptbins_fine,"ptll_me", "ME p_{T} (ll) [GeV]", {100, 250, 400, 650}),
    "1", procs_all, plot_opts).Tag("stitch_fine").Weight(wlumi_dy_nlo)
    .RatioTitle("Inclusive","p_{T}(ll) bins");

  // vector<double> htbins = {0,70, 100, 200, 400, 600, 800, 1200, 2500};
  // pm.Push<Hist1D>(Axis(htbins,"ht_isr_me", "Gen H_{T} [GeV]", {0,70, 100, 200, 400, 600, 800, 1200, 2500}),
  //   "1", procs_ht, plot_opts).Tag("stitch").Weight(wlumi_dy_lo);
  vector<double> htbins_fine = {0, 20, 40, 70, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 700, 800, 1000, 1200, 1500};
  pm.Push<Hist1D>(Axis(htbins_fine,"ht_isr_me", "Gen H_{T} [GeV]", {70, 100, 200, 400, 600, 800, 1200, 2500}),
    "1", procs_ht, plot_opts).Tag("stitch_fine").Weight(wlumi_dy_lo)
  .RatioTitle("Inclusive","H_{T} bins");

  if(single_thread) pm.multithreaded_ = false;
  pm.min_print_ = true;
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

vector<unsigned> higidx(const Baby &b){
  vector<unsigned> idx;
  for (unsigned i(0); i<b.mc_pt()->size(); i++){
    if (b.mc_id()->at(i)==25) idx.push_back(i);
    if (idx.size()>1) break;
  }
  return idx;
}

// vector<unsigned> bidx(const Baby &b){
//   vector<unsigned> idx;
//   for (unsigned i(0); i<b.mc_pt()->size(); i++){
//     if (b.mc_id()->at(i)==25) higidx.push_back(i);
//     if (higidx.size()>1) break;
//   }
//   return idx;
// }
