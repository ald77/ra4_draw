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
  double lumi = 1;
}

int main(int argc, char *argv[]){
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  PlotOpt log_lumi("txt/plot_styles.txt", "RatioPlot");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info).RatioMinimum(0.)
  .RatioMaximum(0.58).PrintVals(true);
  vector<PlotOpt> plot_opts = {log_lumi_info};
  Palette colors("txt/colors.txt", "default");

  string folder_dy = bfolder+"/cms2r0/babymaker/babies/2017_04_01/mc/merged_bare_dy/";
  vector<shared_ptr<Process> >  procs_dy;
  procs_dy.push_back(Process::MakeShared<Baby_full>("# negative", Process::Type::data, 
    kBlack, {folder_dy+"*DYJetsToLL_*.root"}, "(ptll_me<100||type!=6200)&&w_lumi<0"));
  procs_dy.push_back(Process::MakeShared<Baby_full>("# total", Process::Type::background, 
    kGray+2, {folder_dy+"*DYJetsToLL_*.root"}, "(ptll_me<100||type!=6200)"));

  string folder_tt = bfolder+"/cms2r0/babymaker/babies/2016_11_29/mc/unskimmed/";
  string tt_type("*TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX*.root"), tag("ttnom");
   // tt_type = "*TTJets_TuneCUETP8M1_alphaS0118_13TeV-amcatnloFXFX*.root"; tag = "ttaS0118";
  vector<shared_ptr<Process> >  procs_tt;
  procs_tt.push_back(Process::MakeShared<Baby_full>("# negative", Process::Type::data, 
    kBlack, {folder_tt+tt_type}, "w_lumi<0"));
  procs_tt.push_back(Process::MakeShared<Baby_full>("# total", Process::Type::background, 
    kGray+2, {folder_tt+tt_type}, "1"));

  tt_type = "*TTJets_TuneCUETP8M1_alphaS0118_13TeV-amcatnloFXFX*.root"; 
  vector<shared_ptr<Process> >  procs_tt_alphas;
  procs_tt_alphas.push_back(Process::MakeShared<Baby_full>("# negative", Process::Type::data, 
    kBlack, {folder_tt+tt_type}, "w_lumi<0"));
  procs_tt_alphas.push_back(Process::MakeShared<Baby_full>("# total", Process::Type::background, 
    kGray+2, {folder_tt+tt_type}, "1"));

  PlotMaker pm;

  NamedFunc wgt("wgt",[&](const Baby &b){
      return fabs(b.w_lumi()); 
    });

  // pm.Push<Hist1D>(Axis(25,0,2500,"ht_isr_me", "Gen H_{T} [GeV]", {}),
  //   "1", procs_dy, plot_opts).Tag("fneg_dy").Weight(wgt).RatioTitle("# negative","# total");
  // pm.Push<Hist1D>(Axis(20,0,1000,"ptll_me", "Gen p_{T} (ll) [GeV]", {}),
  //   "1", procs_dy, plot_opts).Tag("fneg_dy").Weight(wgt).RatioTitle("# negative","# total");
  // pm.Push<Hist1D>(Axis(12,-0.5,11.5,"njets", "Reco. N_{jets} [GeV]", {}),
  //   "1", procs_dy, plot_opts).Tag("fneg_dy").Weight(wgt).RatioTitle("# negative","# total");

  pm.Push<Hist1D>(Axis(20,0,1000,"met", "E_{T}^{miss} [GeV]", {}),
    "1", procs_tt, plot_opts).Tag("fneg_"+tag).Weight(wgt).RatioTitle("# negative","# total");
  pm.Push<Hist1D>(Axis(12,-0.5,11.5,"njets", "Reco. N_{jets} [GeV]", {}),
    "1", procs_tt, plot_opts).Tag("fneg_"+tag).Weight(wgt).RatioTitle("# negative","# total");
  pm.Push<Hist1D>(Axis(25,0,2500,"ht", "Reco. H_{T} [GeV]", {}),
    "1", procs_tt, plot_opts).Tag("fneg_"+tag).Weight(wgt).RatioTitle("# negative","# total");

  // tag = "ttaS0118";
  // pm.Push<Hist1D>(Axis(25,0,2500,"ht_isr_me", "Gen H_{T} [GeV]", {}),
  //   "1", procs_tt_alphas, plot_opts).Tag("fneg_"+tag).Weight(wgt).RatioTitle("# negative","# total");
  // pm.Push<Hist1D>(Axis(20,0,1000,"met", "E_{T}^{miss} [GeV]", {}),
  //   "1", procs_tt_alphas, plot_opts).Tag("fneg_"+tag).Weight(wgt).RatioTitle("# negative","# total");
  // pm.Push<Hist1D>(Axis(12,-0.5,11.5,"njets", "Reco. N_{jets} [GeV]", {}),
  //   "1", procs_tt_alphas, plot_opts).Tag("fneg_"+tag).Weight(wgt).RatioTitle("# negative","# total");
  // pm.Push<Hist1D>(Axis(25,0,2500,"ht", "Reco. H_{T} [GeV]", {}),
  //   "1", procs_tt_alphas, plot_opts).Tag("fneg_"+tag).Weight(wgt).RatioTitle("# negative","# total");


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
