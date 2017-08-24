#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "TError.h"
#include "TColor.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist1d.hpp"

using namespace std;
using namespace PlotOptTypes;

int main(){
  gErrorIgnoreLevel = 6000;
  string samp("2017");
  string data_path(""),data_cuts("");
  double lumi(0);
//If running on cms#, need to add /net/cms2/ to data_path names
  if(samp == "2017") { 
	lumi = 8.3;
	data_path = "/cms2r0/babymaker/babies/2017_08_14/data/merged_database_standard/*.root";
	data_cuts = "(met/met_calo<5.0)&&pass_ra2_badmu&&(trig[3]||trig[7]||trig[9]||trig[15]||trig[20]||trig[21]||trig[23]||trig[24])&&(pass_hbhe&&pass_hbheiso&&pass_goodv&&pass_ecaldeadcell)";
	}
  else if(samp == "2016") {
	lumi = 36.2;
	data_path = "/cms2r0/babymaker/babies/2017_02_14/data/skim_standard/*.root";
	data_cuts = "(met/met_calo<5.0)&&pass&&trig_ra4";
	}
  else if(samp == "2017-5ifb") { //First 5.8 ifb of 2017 data
        lumi = 5.829;
        data_path = "/cms2r0/babymaker/babies/2017_08_14/data/merged_database_standard/*100.root";
        data_cuts = "(met/met_calo<5.0)&&pass_ra2_badmu&&(trig[3]||trig[7]||trig[9]||trig[15]||trig[20]||trig[21]||trig[23]||trig[24])&&(pass_hbhe&&pass_hbheiso&&pass_goodv&&pass_ecaldeadcell)";
        }
  else if(samp == "2017-extra") { //Additional 2.5 ifb
        lumi = 2.488;
        data_path = "/cms2r0/babymaker/babies/2017_08_14/data/merged_database_standard/*53.root";
        data_cuts = "(met/met_calo<5.0)&&pass_ra2_badmu&&(trig[3]||trig[7]||trig[9]||trig[15]||trig[20]||trig[21]||trig[23]||trig[24])&&(pass_hbhe&&pass_hbheiso&&pass_goodv&&pass_ecaldeadcell)";
        }
  else if(samp == "2017-Halo") { //2017 data with cschalo filter applied
        lumi = 8.3;
        data_path = "/cms2r0/babymaker/babies/2017_08_14/data/merged_database_standard/*.root";
        data_cuts = "(met/met_calo<5.0)&&pass_ra2_badmu&&(trig[3]||trig[7]||trig[9]||trig[15]||trig[20]||trig[21]||trig[23]||trig[24])&&(pass_hbhe&&pass_hbheiso&&pass_goodv&&pass_ecaldeadcell&&pass_cschalo)";
        }
  string trig_skim1l_mc = "/cms2r0/babymaker/babies/2017_01_27/mc/merged_mcbase_standard/";
  string trig_skim2l_mc = "/cms2r0/babymaker/babies/2017_01_27/mc/merged_mcbase_standard/"; 

  Palette colors("txt/colors.txt", "default");
  //
  //1-Lepton Plots
  //
  auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {trig_skim1l_mc+"*_TTJets*SingleLept*.root"}, "ntruleps<=1&&stitch_met");
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {trig_skim1l_mc+"*_TTJets*DiLept*.root"}, "ntruleps>=2&&stitch_met");
  auto wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {trig_skim1l_mc+"*_WJetsToLNu*.root"},"stitch_met");
  auto single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {trig_skim1l_mc+"*_ST_*.root"});
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {trig_skim1l_mc+"*_TTWJets*.root", trig_skim1l_mc+"*_TTZ*.root"});
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {trig_skim1l_mc+"*DYJetsToLL*.root", trig_skim1l_mc+"*QCD_HT*0_Tune*.root", trig_skim1l_mc+"*QCD_HT*Inf_Tune*.root",
        trig_skim1l_mc+"*_ZJet*.root", trig_skim1l_mc+"*_ttHTobb_M125_*.root",
        trig_skim1l_mc+"*_TTGJets*.root", trig_skim1l_mc+"*_TTTT_*.root",
        trig_skim1l_mc+"*_WH_HToBB*.root", trig_skim1l_mc+"*_ZH_HToBB*.root", 
        trig_skim1l_mc+"*_WWTo*.root", trig_skim1l_mc+"*_WZ*.root",
        trig_skim1l_mc+"_ZZ_*.root"}, "stitch_met");
 
  auto data_1l = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {data_path},data_cuts);
  vector<shared_ptr<Process> > full_trig_skim_1l = {data_1l, tt1l, tt2l, wjets, single_t, ttv, other};

  //
  // 2-lepton Plots
  //
  auto tt1l_2l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {trig_skim2l_mc+"*_TTJets*Lept*.root"},
    "ntruleps<=1");
  auto tt2l_2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {trig_skim2l_mc+"*_TTJets*Lept*.root"},
    "ntruleps>=2");
  auto dy_2l = Process::MakeShared<Baby_full>("DY", Process::Type::background, colors("wjets"),
    {trig_skim2l_mc+"*_DYJetsToLL*HT*.root"});
  auto other_2l = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("ttv"),
    {trig_skim2l_mc+"*_TTWJet*root",  trig_skim2l_mc+"*_TTZTo*.root",
        trig_skim2l_mc+"*_ZJet*.root", trig_skim2l_mc+"*_WWTo*.root",
        trig_skim2l_mc+"*ggZH_HToBB*.root", trig_skim2l_mc+"*ttHJetTobb*.root",
        trig_skim2l_mc+"*_TTGJets*.root", trig_skim2l_mc+"*_TTTT_*.root",
        trig_skim2l_mc+"*_WH_HToBB*.root", trig_skim2l_mc+"*_WZTo*.root",
        trig_skim2l_mc+"*_ST_*.root", trig_skim2l_mc+"*_WJetsToLNu*.root",
        trig_skim2l_mc+"*_ZH_HToBB*.root", trig_skim2l_mc+"_ZZ_*.root",
        trig_skim2l_mc+"*_QCD_HT700to1000*.root",trig_skim2l_mc+"*_QCD_HT1000to1500*.root",
        trig_skim2l_mc+"*_QCD_HT1500to2000*.root",trig_skim2l_mc+"*_QCD_HT2000toInf*.root"});

  auto data_2l = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {data_path},data_cuts); // 2016 
  vector<shared_ptr<Process> > full_trig_skim_2l = {data_2l, tt1l_2l, tt2l_2l, dy_2l, other_2l};

  //
  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .Bottom(BottomType::off)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info);
  vector<PlotOpt> all_plot_types = {log_lumi_info, lin_lumi_info};

  PlotMaker pm;


//Good Plots for ttbar
  //1 Lep
	//Total
  pm.Push<Hist1D>(Axis(15, 500, 2000., "st", "S_{T} [GeV]", {500.}),
                  "nleps==1&&st>500&&met>200&&mt<140&&njets>=6&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
                  "nleps==1&&st>500&&met>200&&mt<140&&njets>=6&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}", {5.5, 8.5}),
                  "nleps==1&&st>500&&met>200&&mt<140&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm_moriond", "N_{b}", {0.5, 1.5, 2.5}),
                  "nleps==1&&st>500&&met>200&&mt<140&&njets>=6&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(15, 0., 2000., "mj14", "M_{J} [GeV]", {400.}),
                  "nleps==1&&st>500&&met>200&&mt<140&&njets>=6&&nbm_moriond>=1&&nveto==0",
                   full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(14, 0., 280., "mt", "m_{T} [GeV]", {140.}),
                  "nleps==1&&st>500&&met>200&&mj14<400&&njets>=6&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(14, 0., 280., "mt_nohf", "m_{T} [GeV]", {140.}),
                  "nleps==1&&st>500&&met>200&&mj14<400&&njets>=6&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
        //Muons
  pm.Push<Hist1D>(Axis(15, 500, 2000., "st", "S_{T} [GeV]", {500.}),
                  "nmus==1&&nleps==1&&st>500&&met>200&&mt<140&&njets>=6&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
                  "nmus==1&&nleps==1&&st>500&&met>200&&mt<140&&njets>=6&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}", {5.5, 8.5}),
                  "nmus==1&&nleps==1&&st>500&&met>200&&mt<140&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm_moriond", "N_{b}", {0.5, 1.5, 2.5}),
                  "nmus==1&&nleps==1&&st>500&&met>200&&mt<140&&njets>=6&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(15, 0., 2000., "mj14", "M_{J} [GeV]", {400.}),
                  "nmus==1&&nleps==1&&st>500&&met>200&&mt<140&&njets>=6&&nbm_moriond>=1&&nveto==0",
                   full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(14, 0., 280., "mt", "m_{T} [GeV]", {140.}),
                  "nmus==1&&nleps==1&&st>500&&met>200&&mj14<400&&njets>=6&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(14, 0., 280., "mt_nohf", "m_{T} [GeV]", {140.}),
                  "nmus==1&&nleps==1&&st>500&&met>200&&mj14<400&&njets>=6&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
        //Electron
  pm.Push<Hist1D>(Axis(15, 500, 2000., "st", "S_{T} [GeV]", {500.}),
                  "nels==1&&nleps==1&&st>500&&met>200&&mt<140&&njets>=6&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
                  "nels==1&&nleps==1&&st>500&&met>200&&mt<140&&njets>=6&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}", {5.5, 8.5}),
                  "nels==1&&nleps==1&&st>500&&met>200&&mt<140&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm_moriond", "N_{b}", {0.5, 1.5, 2.5}),
                  "nels==1&&nleps==1&&st>500&&met>200&&mt<140&&njets>=6&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(15, 0., 2000., "mj14", "M_{J} [GeV]", {400.}),
                  "nels==1&&nleps==1&&st>500&&met>200&&mt<140&&njets>=6&&nbm_moriond>=1&&nveto==0",
                   full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(14, 0., 280., "mt", "m_{T} [GeV]", {140.}),
                  "nels==1&&nleps==1&&st>500&&met>200&&mj14<400&&njets>=6&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(14, 0., 280., "mt_nohf", "m_{T} [GeV]", {140.}),
                  "nels==1&&nleps==1&&st>500&&met>200&&mj14<400&&njets>=6&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  //2 Lep
  pm.Push<Hist1D>(Axis(15, 500, 2000., "st", "S_{T} [GeV]", {500.}),
          "st>500&&nleps==2&&njets>=5&&nbm_moriond<=2&&met>200&&met<500",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
          "st>500&&nleps==2&&njets>=5&&nbm_moriond<=2&&met>200",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}", {4.5, 7.5}),
          "st>500&&nleps==2&&nbm_moriond<=2&&met>200&&met<500",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm_moriond", "N_{b}", {0.5, 1.5, 2.5}),
          "st>500&&nleps==2&&njets>=5&&met>200&&met<500",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(14, 0., 280., "mt", "m_{T} [GeV]", {140.}),
          "st>500&&nleps==2&&njets>=5&&nbm_moriond<=2&&met>200&&met<500",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(14, 0., 280., "mt_nohf", "m_{T} [GeV]", {140.}),
          "st>500&&nleps==2&&njets>=5&&nbm_moriond<=2&&met>200&&met<500",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(15, 0., 2000., "mj14", "M_{J} [GeV]", {400.}),
                  "st>500&&nleps==2&&njets>=5&&nbm_moriond<=2&&met>200&&met<500",
                  full_trig_skim_1l, all_plot_types);
//Looser Selection on njets(>=4) for both 1/2 Lep regions
  //1 Lep
        //Combined
  pm.Push<Hist1D>(Axis(15, 500, 2000., "st", "S_{T} [GeV]", {500.}),
                  "nleps==1&&st>500&&met>200&&mt<140&&njets>=4&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
                  "nleps==1&&st>500&&met>200&&mt<140&&njets>=4&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}", {5.5, 8.5}),
                  "nleps==1&&st>500&&met>200&&mt<140&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm_moriond", "N_{b}", {0.5, 1.5, 2.5}),
                  "nleps==1&&st>500&&met>200&&mt<140&&njets>=4&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(15, 0., 2000., "mj14", "M_{J} [GeV]", {400.}),
                  "nleps==1&&st>500&&met>200&&mt<140&&njets>=4&&nbm_moriond>=1&&nveto==0",
                   full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(14, 0., 280., "mt", "m_{T} [GeV]", {140.}),
                  "nleps==1&&st>500&&met>200&&mj14<400&&njets>=4&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(14, 0., 280., "mt_nohf", "m_{T} [GeV]", {140.}),
                  "nleps==1&&st>500&&met>200&&mj14<400&&njets>=4&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
        //Muons
  pm.Push<Hist1D>(Axis(15, 500, 2000., "st", "S_{T} [GeV]", {500.}),
                  "nmus==1&&nleps==1&&st>500&&met>200&&mt<140&&njets>=4&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
                  "nmus==1&&nleps==1&&st>500&&met>200&&mt<140&&njets>=4&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}", {5.5, 8.5}),
                  "nmus==1&&nleps==1&&st>500&&met>200&&mt<140&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm_moriond", "N_{b}", {0.5, 1.5, 2.5}),
                  "nmus==1&&nleps==1&&st>500&&met>200&&mt<140&&njets>=4&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(15, 0., 2000., "mj14", "M_{J} [GeV]", {400.}),
                  "nmus==1&&nleps==1&&st>500&&met>200&&mt<140&&njets>=4&&nbm_moriond>=1&&nveto==0",
                   full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(14, 0., 280., "mt", "m_{T} [GeV]", {140.}),
                  "nmus==1&&nleps==1&&st>500&&met>200&&mj14<400&&njets>=4&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(14, 0., 280., "mt_nohf", "m_{T} [GeV]", {140.}),
                  "nmus==1&&nleps==1&&st>500&&met>200&&mj14<400&&njets>=4&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
        //Electron
  pm.Push<Hist1D>(Axis(15, 500, 2000., "st", "S_{T} [GeV]", {500.}),
                  "nels==1&&nleps==1&&st>500&&met>200&&mt<140&&njets>=4&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
                  "nels==1&&nleps==1&&st>500&&met>200&&mt<140&&njets>=4&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}", {5.5, 8.5}),
                  "nels==1&&nleps==1&&st>500&&met>200&&mt<140&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm_moriond", "N_{b}", {0.5, 1.5, 2.5}),
                  "nels==1&&nleps==1&&st>500&&met>200&&mt<140&&njets>=4&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(15, 0., 2000., "mj14", "M_{J} [GeV]", {400.}),
                  "nels==1&&nleps==1&&st>500&&met>200&&mt<140&&njets>=4&&nbm_moriond>=1&&nveto==0",
                   full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(14, 0., 280., "mt", "m_{T} [GeV]", {140.}),
                  "nels==1&&nleps==1&&st>500&&met>200&&mj14<400&&njets>=4&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(14, 0., 280., "mt_nohf", "m_{T} [GeV]", {140.}),
                  "nels==1&&nleps==1&&st>500&&met>200&&mj14<400&&njets>=4&&nbm_moriond>=1&&nveto==0",
                  full_trig_skim_1l, all_plot_types);
  //2 Lep
  pm.Push<Hist1D>(Axis(15, 500, 2000., "st", "S_{T} [GeV]", {500.}),
          "st>500&&nleps==2&&njets>=4&&nbm_moriond<=2&&met>200&&met<500",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
          "st>500&&nleps==2&&njets>=4&&nbm_moriond<=2&&met>200",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}", {4.5, 7.5}),
          "st>500&&nleps==2&&nbm_moriond<=2&&met>200&&met<500",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm_moriond", "N_{b}", {0.5, 1.5, 2.5}),
          "st>500&&nleps==2&&njets>=4&&met>200&&met<500",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(14, 0., 280., "mt", "m_{T} [GeV]", {140.}),
          "st>500&&nleps==2&&njets>=4&&nbm_moriond<=2&&met>200&&met<500",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(14, 0., 280., "mt_nohf", "m_{T} [GeV]", {140.}),
          "st>500&&nleps==2&&njets>=4&&nbm_moriond<=2&&met>200&&met<500",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(15, 0., 2000., "mj14", "M_{J} [GeV]", {400.}),
                  "st>500&&nleps==2&&njets>=4&&nbm_moriond<=2&&met>200&&met<500",
                  full_trig_skim_1l, all_plot_types);
  pm.MakePlots(lumi);

}
