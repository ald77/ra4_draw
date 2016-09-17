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

using namespace std;
using namespace PlotOptTypes;

namespace{
  bool single_thread = false;
}

int main(int argc, char *argv[]){
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  double lumi = 20;

  string trig_skim_mc = "/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/";
  string trig_skim_signal = "/cms2r0/babymaker/babies/2016_08_10/TChiHH/merged_higmc_higloose/";

  Palette colors("txt/colors.txt", "default");
  auto tchi = Process::MakeShared<Baby_full>("TChiHH(400,1)", Process::Type::signal, kRed,
    {trig_skim_signal+"*TChiHH_mGluino-400*.root"});
  auto tt = Process::MakeShared<Baby_full>("t#bar{t}", Process::Type::background, colors("tt_1l"),
    {trig_skim_mc+"*_TTJets*Lept*.root", trig_skim_mc+"*_TTJets_HT*.root"}, "stitch");
  auto wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {trig_skim_mc+"*_WJetsToLNu*.root"});
  //  auto znunu = Process::MakeShared<Baby_full>("Z#rightarrow#nu#nu", Process::Type::background, colors("znunu"),
  //  {trig_skim_mc+"*ZJetsToNuNu_HT*.root"}); 
  auto wt = Process::MakeShared<Baby_full>("W+top", Process::Type::background, kRed+3,
    {trig_skim_mc+"*ST_tW*.root"});
  auto qcd = Process::MakeShared<Baby_full>("QCD", Process::Type::background, colors("qcd"),
    {trig_skim_mc+"*_QCD_HT*00_Tune*.root", trig_skim_mc+"*_QCD_HT*Inf_Tune*.root"});
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, kViolet-6,
    {trig_skim_mc+"*DYJetsToLL*.root", trig_skim_mc+"*TTTT*.root", trig_skim_mc+"*_TTWJets*.root",
        trig_skim_mc+"*_WH_*.root", trig_skim_mc+"*_ttHJet*.root", trig_skim_mc+"*_ZH_*.root",
        trig_skim_mc+"*_ST_*channel*.root", trig_skim_mc+"*_TTGJets*.root", trig_skim_mc+"*_TTZTo*.root"});


  vector<shared_ptr<Process> > full_trig_skim = {tchi, tt, wjets, wt, qcd, other};

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::off)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info);
  vector<PlotOpt> all_plot_types = {log_lumi_info};

  PlotMaker pm;

  vector<NamedFunc> nbcuts;
  nbcuts.push_back(" nbt >= 2 && nbl <= 3");
  nbcuts.push_back(" nbt >= 2 && nbl >= 4");
  NamedFunc metskim("njets>=4&&njets<=5&&met>100&&nvleps==0");
  NamedFunc trkskim("njets>=4&&njets<=5&&met>250&&nvleps==0");
  NamedFunc skim("njets>=4&&njets<=5&&met>250&&nvleps==0&&ntks==0");
  NamedFunc DeltaR("hig_drmax < 2.2");
  NamedFunc AverageM("hig_am > 100 && hig_am < 140");
  NamedFunc DeltaM("hig_dm < 40");
  NamedFunc LDP("!low_dphi");

  pm.Push<Hist1D>(Axis(7,-0.5,6.5,"nbl", "N_{b-jet}^{L}", {100.}),
                  skim&&DeltaR&&AverageM&&DeltaM&&LDP, full_trig_skim, all_plot_types);
  pm.Push<Hist1D>(Axis(7,-0.5,6.5,"nbm", "N_{b-jet}^{M}", {100.}),
                  skim&&DeltaR&&AverageM&&DeltaM&&LDP, full_trig_skim, all_plot_types);
  pm.Push<Hist1D>(Axis(7,-0.5,6.5,"nbt", "N_{b-jet}^{T}", {100.}),
                  skim&&DeltaR&&AverageM&&DeltaM&&LDP, full_trig_skim, all_plot_types);

  for(auto inb: nbcuts) {
    pm.Push<Hist1D>(Axis(20,100,600,"met", "E_{T}^{miss} [GeV]", {150., 250.}),
                    inb&&metskim&&DeltaR&&AverageM&&DeltaM&&LDP, full_trig_skim, all_plot_types);
    pm.Push<Hist1D>(Axis(32,0,160,"hig_dm", "#Deltam [GeV]", {40.}),
                    inb&&skim&&DeltaR&&AverageM&&LDP, full_trig_skim, all_plot_types);
    pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m> [GeV]", {100., 140.}),
                    inb&&skim&&DeltaR&&DeltaM&&LDP, full_trig_skim, all_plot_types);
    pm.Push<Hist1D>(Axis(20,0,4,"hig_drmax", "#DeltaR_{max}", {2.2}),
                    inb&&skim&&AverageM&&DeltaM&&LDP, full_trig_skim, all_plot_types);
    pm.Push<Hist1D>(Axis(28,0,1400,"ht", "H_{T} [GeV]"),
                    inb&&skim&&DeltaR&&AverageM&&DeltaM&&LDP, full_trig_skim, all_plot_types);
    pm.Push<Hist1D>(Axis(5,2.5,7.5,"njets", "N_{jet}", {3.5, 5.5}),
                    inb&&skim&&DeltaR&&AverageM&&DeltaM&&LDP, full_trig_skim, all_plot_types);
    pm.Push<Hist1D>(Axis(5,-.5,4.5,"ntks", "N_{track}", {0.5}),
                    inb&&trkskim&&DeltaR&&AverageM&&DeltaM&&LDP, full_trig_skim, all_plot_types);
    pm.Push<Hist1D>(Axis(32,0,3.2,"dphi2", "#Delta#phi_{2}", {0.5}),
                    inb&&skim&&DeltaR&&AverageM&&DeltaM, full_trig_skim, all_plot_types);
    pm.Push<Hist1D>(Axis(32,0,3.2,"dphi3", "#Delta#phi_{3}", {0.3}),
                    inb&&skim&&DeltaR&&AverageM&&DeltaM, full_trig_skim, all_plot_types);
    pm.Push<Hist1D>(Axis(32,0,3.2,"dphi4", "#Delta#phi_{4}", {0.3}),
                    inb&&skim&&DeltaR&&AverageM&&DeltaM, full_trig_skim, all_plot_types);
    pm.Push<Hist1D>(Axis(30,0,600,"jets_pt[0]", "Jet 1 p_{T} [GeV]", {50.}),
                    inb&&skim&&DeltaR&&AverageM&&DeltaM&&LDP, full_trig_skim, all_plot_types);
    pm.Push<Hist1D>(Axis(35,0,350,"jets_pt[1]", "Jet 2 p_{T} [GeV]", {50.}),
                    inb&&skim&&DeltaR&&AverageM&&DeltaM&&LDP, full_trig_skim, all_plot_types);
    pm.Push<Hist1D>(Axis(25,0,250,"jets_pt[2]", "Jet 3 p_{T} [GeV]"),
                    inb&&skim&&DeltaR&&AverageM&&DeltaM&&LDP, full_trig_skim, all_plot_types);
    pm.Push<Hist1D>(Axis(20,0,200,"jets_pt[3]", "Jet 4 p_{T} [GeV]"),
                    inb&&skim&&DeltaR&&AverageM&&DeltaM&&LDP, full_trig_skim, all_plot_types);
  }

  pm.Push<Table>("cutflow", vector<TableRow>{
      TableRow("$MET > 100$, $\\text{2M b-tags}$, $\\text{4 or 5 jets}$, $0\\ell$", "1"),
        TableRow("$\\Delta\\phi_{\\text{min}}$", LDP),
        TableRow("$\\Delta m < 40$", DeltaM&&LDP,1,0),
        TableRow("$\\left< m \\right> \\in (100,140)$", AverageM&&LDP),
        TableRow("$\\Delta R_{\\text{max}} < 2.2$", DeltaR&&LDP)
        },full_trig_skim,0);

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
