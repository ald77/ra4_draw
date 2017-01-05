
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

int main(){
  gErrorIgnoreLevel = 6000;
  

  double lumi = 36202;

  string base_path = "";
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname,"compute-")){
    base_path = "/net/cms29";
  }
  string mc_dir = base_path+"/cms29r0/babymaker/babies/2016_08_10/mc/unskimmed/";
  string sig_dir = base_path+"/cms29r0/babymaker/babies/2016_08_10/T1tttt/unskimmed/";
  Palette colors("txt/colors.txt", "default");

  auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {mc_dir+"*_TTJets*Lept*.root", mc_dir+"*_TTJets_HT*.root"},
    "ntruleps<=1&&stitch");
  tt1l->SetMarkerStyle(23);
  tt1l->SetMarkerSize(0.8);
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {mc_dir+"*_TTJets*Lept*.root", mc_dir+"*_TTJets_HT*.root"},
    "ntruleps>=2&&stitch");
  tt1l->SetMarkerStyle(22);
  tt1l->SetMarkerSize(0.8);
  auto wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {mc_dir+"*_WJetsToLNu*.root"});
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
    {sig_dir+"*SMS-T1tttt_mGluino-1800_mLSP-100_Tu*.root"});
  t1tttt_nc->SetMarkerStyle(21);
  t1tttt_nc->SetMarkerSize(0.9);
  auto t1tttt_c = Process::MakeShared<Baby_full>("T1tttt(1400,1000)", Process::Type::signal, colors("t1tttt"),
    {sig_dir+"*SMS-T1tttt_mGluino-1400_mLSP-1000_Tu*.root"});

  t1tttt_c->SetLineStyle(2);
  t1tttt_c->SetMarkerStyle(21);
  t1tttt_c->SetMarkerSize(0.9);

  vector<shared_ptr<Process> > full_trig_skim = {t1tttt_nc, t1tttt_c, tt1l, tt2l};

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
 
  Table & cutflow = pm.Push<Table>("cutflow", vector<TableRow>{
      TableRow("Baseline"),
     
      TableRow("No Selection, low PU", "ntrupv<20"),
        TableRow("Baseline", "nleps==1&&st>500&&met>200&&njets>=6&&nbm>=1&&nveto==0&&ntrupv<20"),
        TableRow("Baseline and high MJ14", "nleps==1&&st>500&&met>200&&njets>=6&&nbm>=1&&nveto==0&&mj14>400&&ntrupv<20"),
	TableRow("Baseline and high MJ12", "nleps==1&&st>500&&met>200&&njets>=6&&nbm>=1&&nveto==0&&mj12>400&&ntrupv<20"),


	TableRow("No Selection, high PU", "ntrupv>35"),
	TableRow("Baseline", "nleps==1&&st>500&&met>200&&njets>=6&&nbm>=1&&nveto==0&&ntrupv>35"),
	TableRow("Baseline and high MJ14", "nleps==1&&st>500&&met>200&&njets>=6&&nbm>=1&&nveto==0&&mj14>400&&ntrupv>35"),
	TableRow("Baseline and high MJ12", "nleps==1&&st>500&&met>200&&njets>=6&&nbm>=1&&nveto==0&&mj12>400&&ntrupv>35")
        }, full_trig_skim);


 

  if(single_thread) pm.multithreaded_ = false;
  pm.MakePlots(lumi);

  vector<GammaParams> yields = cutflow.BackgroundYield(lumi);
  for(const auto &yield: yields){
    cout << yield << endl;
  }
}
