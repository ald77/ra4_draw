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
#include "core/event_scan.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"

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

  auto t1tttt_nc = Process::MakeShared<Baby_full>("T1tttt(1800,100)", Process::Type::signal, colors("t1tttt"),
    {base_path+"/cms29r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/*SMS-T1tttt_mGluino-1800_mLSP-100_*.root"});
  t1tttt_nc->SetMarkerStyle(21);
  t1tttt_nc->SetMarkerSize(0.9);

  auto t1tttt_c = Process::MakeShared<Baby_full>("T1tttt(1800,100)", Process::Type::signal, colors("t1tttt"),
    {base_path+"/cms29r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/*SMS-T1tttt_mGluino-1800_mLSP-100_*.root"});
  t1tttt_c->SetLineStyle(2);
  t1tttt_c->SetMarkerStyle(21);
  t1tttt_c->SetMarkerSize(0.9);

  auto data = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {base_path+"/cms29r0/babymaker/babies/2017_02_14/data/merged_database_stdnj5/*.root"},"pass&&trig_ra4");
  data->SetMarkerStyle(20);
  data->SetMarkerSize(1.);

  vector<shared_ptr<Process> > full_trig_skim = {data, t1tttt_nc, t1tttt_c, tt1l, tt2l, wjets, single_t, ttv, other};

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::ratio)
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
  vector<PlotOpt> all_plot_types = {log_lumi, lin_lumi, log_shapes, lin_shapes,
                                    log_lumi_info, lin_lumi_info, log_shapes_info, lin_shapes_info};
  PlotOpt style2D("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> bkg_hist = {style2D().Stack(StackType::data_norm).Title(TitleType::preliminary)};
  vector<PlotOpt> bkg_pts = {style2D().Stack(StackType::lumi_shapes).Title(TitleType::info)};

  PlotMaker pm;
  pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nleps", "Num. Leptons", {0.5, 1.5}),
                  "st>500&&met>200&&njets>=6&&nbm>=1", full_trig_skim, all_plot_types);
  pm.Push<Hist1D>(Axis(40, 0, 2000., "st", "S_{T} [GeV]", {500.}),
                  "nleps==1&&met>200&&njets>=6&&nbm>=1", full_trig_skim, all_plot_types);
  pm.Push<Hist1D>(Axis(40, 0, 1000., "met", "MET [GeV]", {200., 400.}),
                  "nleps==1&&st>500&&njets>=6&&nbm>=1", full_trig_skim, all_plot_types);
  pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "Num. AK4 Jets", {5.5, 8.5}),
                  "nleps==1&&st>500&&met>200&&nbm>=1", full_trig_skim, all_plot_types);
  pm.Push<Hist1D>(Axis(11, -0.5, 10.5, "nbm", "Num. b-Tagged Jets", {0.5, 1.5, 2.5}),
                  "nleps==1&&st>500&&met>200&&njets>=6", full_trig_skim, all_plot_types);
  pm.Push<Hist1D>(Axis(24, 0., 1200., "mj14", "M_{J} [GeV]", {250., 400.}),
                  "nleps==1&&st>500&&met>200&&njets>=6&&nbm>=1", full_trig_skim, all_plot_types);
  pm.Push<Hist1D>(Axis(12, 0., 420., "mt", "m_{T} [GeV]", {140.}),
                  "nleps==1&&st>500&&met>200&&njets>=6&&nbm>=1", full_trig_skim, all_plot_types);
  pm.Push<Hist1D>(Axis(15, 0., 1500., "mj14", "M_{J} [GeV]", {400.}),
                  "nleps==1&&st>500&&met>200", full_trig_skim, all_plot_types)
    .Tag("changing_tags_and_weights").Weight("1.2345*weight").RatioTitle("Numerator","Denominator");
  Table & cutflow = pm.Push<Table>("cutflow", vector<TableRow>{
      TableRow("Baseline"),
        TableRow("No Selection", "1"),
        TableRow("$1\\ell$, $H_{T}>500$, $E_{\\text{T}}^{\\text{miss}}>200$", "nleps==1&&st>500&&met>200"),
        TableRow("$N_{\\text{jets}}\\geq6$", "nleps==1&&st>500&&met>200&&njets>=6"),
        TableRow("$N_{b}\\geq1$", "nleps==1&&st>500&&met>200&&njets>=6&&nbm>=1"),
        TableRow("$M_{J}>250$", "nleps==1&&st>500&&met>200&&njets>=6&&nbm>=1&&mj14>250", 1, 0),
        TableRow("ABCD Signal Region"),
        TableRow("$m_{T}>140$", "nleps==1&&st>500&&met>200&&njets>=6&&nbm>=1&&mj14>250&&mt>140"),
        TableRow("$M_{J}>400$", "nleps==1&&st>500&&met>200&&njets>=6&&nbm>=1&&mj14>400&&mt>140"),
        TableRow("Binning"),
        TableRow("$E_{\\text{T}}^{\\text{miss}}>500$", "nleps==1&&st>500&&met>500&&njets>=6&&nbm>=1&&mj14>400&&mt>140"),
        TableRow("$N_{\\text{jets}}\\geq9$", "nleps==1&&st>500&&met>500&&njets>=9&&nbm>=1&&mj14>400&&mt>140"),
        TableRow("$N_{b}\\geq3$", "nleps==1&&st>500&&met>500&&njets>=9&&nbm>=3&&mj14>400&&mt>140")
        }, full_trig_skim);

  NamedFunc mm_wgt = Functions::MismeasurementWeight("txt/sys_weights.cfg", "off");
  NamedFunc mm_crt = Functions::MismeasurementCorrection("txt/sys_weights.cfg", "off",
                                                         Functions::Variation::central);
  pm.Push<EventScan>("scan", true, vector<NamedFunc>{"weight", "met", mm_wgt, mm_crt},
                     vector<shared_ptr<Process> >{tt1l});
  pm.Push<Hist2D>(Axis(48, 0., 1200., "mj14", "M_{J} [GeV]", {250., 400.}),
                  Axis(25, 0., 700., "mt", "m_{T} [GeV]", {140.}),
                  "nleps==1&&st>500&&met>200&&njets>=6&&nbm>=1",
                  full_trig_skim, bkg_hist);
  pm.Push<Hist2D>(Axis(48, 0., 1200., "mj14", "M_{J} [GeV]", {250., 400.}),
                  Axis(25, 0., 700., "mt", "m_{T} [GeV]", {140.}),
                  "nleps==1&&st>500&&met>200&&njets>=6&&nbm>=1",
                  vector<shared_ptr<Process> >{tt1l, tt2l, t1tttt_nc}, bkg_pts);

  if(single_thread) pm.multithreaded_ = false;
  pm.MakePlots(lumi);

  vector<GammaParams> yields = cutflow.BackgroundYield(lumi);
  for(const auto &yield: yields){
    cout << yield << endl;
  }
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
