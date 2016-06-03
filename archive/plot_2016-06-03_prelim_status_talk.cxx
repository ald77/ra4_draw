#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "TError.h"
#include "TColor.h"

#include "baby.hpp"
#include "process.hpp"
#include "named_func.hpp"
#include "plot_maker.hpp"
#include "plot_opt.hpp"
#include "palette.hpp"

using namespace std;
using namespace PlotOptTypes;

template<typename T>
shared_ptr<Process> Proc(const string process_name, Process::Type type,
                         int color, const set<string> &files, const string &cut = "1"){
  return make_shared<Process>(process_name, type, color,
                              unique_ptr<Baby>(new T(files)),
                              cut);
}

int main(){
  gErrorIgnoreLevel = 6000;

  double lumi = 0.589;
  
  // 74x
  string trig_skim_mc = "/net/cms27/cms27r0/babymaker/2016_04_29/mc/merged_1lht500met200/";
  string trig_mc = "/net/cms27/cms27r0/babymaker/2016_04_29/mc/unskimmed/";

  Palette colors("txt/colors.txt", "default");

  auto tt1l = Proc<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {trig_skim_mc+"*_TTJets*Lept*.root", trig_skim_mc+"*_TTJets_HT*.root"},
    "ntruleps<=1&&stitch");
  auto tt2l = Proc<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {trig_skim_mc+"*_TTJets*Lept*.root", trig_skim_mc+"*_TTJets_HT*.root"},
    "ntruleps>=2&&stitch");
  auto wjets = Proc<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {trig_skim_mc+"*_WJetsToLNu*.root"});
  auto single_t = Proc<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {trig_skim_mc+"*_ST_*.root"});
  auto ttv = Proc<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {trig_skim_mc+"*_TTWJets*.root", trig_skim_mc+"*_TTZTo*.root"});
  auto other = Proc<Baby_full>("Other", Process::Type::background, colors("other"),
    {trig_skim_mc+"*DYJetsToLL*.root", trig_skim_mc+"*_QCD_HT*.root",
        trig_skim_mc+"*_ZJet*.root", trig_skim_mc+"*_WWTo*.root",
        trig_skim_mc+"*ggZH_HToBB*.root", trig_skim_mc+"*ttHJetTobb*.root",
        trig_skim_mc+"*_TTGJets*.root", trig_skim_mc+"*_TTTT_*.root",
        trig_skim_mc+"*_WH_HToBB*.root", trig_skim_mc+"*_WZTo*.root",
        trig_skim_mc+"*_ZH_HToBB*.root", trig_skim_mc+"_ZZ_*.root"});

  auto t1tttt_nc = Proc<Baby_full>("T1tttt(1500,100)", Process::Type::signal, colors("t1tttt"),
    {trig_mc+"*SMS-T1tttt_mGluino-1500_mLSP-100*.root"});
  auto t1tttt_c = Proc<Baby_full>("T1tttt(1200,800)", Process::Type::signal, colors("t1tttt"),
    {trig_mc+"*SMS-T1tttt_mGluino-1200_mLSP-800*.root"});
  t1tttt_c->SetLineStyle(2);

  auto data_1l = Proc<Baby_full>("Data", Process::Type::data, kBlack,
    {"/net/cms2/cms2r0/babymaker/babies/2016_05_31/data/met/*.root"},"pass&&trig[28]"); // 0.589 fb-1
  vector<shared_ptr<Process> > full_trig_skim_1l = {data_1l, t1tttt_nc, t1tttt_c, tt1l, tt2l, wjets, single_t, ttv, other}; 
  
  // 
  // 2l  / 0l
  //  
  auto tt1l_unskim = Proc<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {trig_mc+"*_TTJets*Lept*.root"},
    "ntruleps<=1");
  auto tt2l_unskim = Proc<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {trig_mc+"*_TTJets*Lept*.root"},
    "ntruleps>=2");
  auto tt_unskim = Proc<Baby_full>("t#bar{t}", Process::Type::background, colors("tt_1l"),
    {trig_mc+"*_TTJets_TuneCUETP8M1_13TeV*.root"});
  auto dy_unskim = Proc<Baby_full>("DY", Process::Type::background, colors("wjets"),
    {trig_mc+"*DYJetsToLL*HT*.root"});
  auto qcd_unskim = Proc<Baby_full>("QCD", Process::Type::background, colors("qcd"),
    {trig_mc+"*_QCD_HT700to1000*.root",trig_mc+"*_QCD_HT1000to1500*.root",
	trig_mc+"*_QCD_HT1500to2000*.root",trig_mc+"*_QCD_HT2000toInf*.root"});
  auto other_unskim = Proc<Baby_full>("Other", Process::Type::background, colors("ttv"),
    {trig_mc+"*_TTWJet*root",  trig_mc+"*_TTZTo*.root",
        trig_mc+"*_ZJet*.root", trig_mc+"*_WWTo*.root",
        trig_mc+"*ggZH_HToBB*.root", trig_mc+"*ttHJetTobb*.root",
        trig_mc+"*_TTGJets*.root", trig_mc+"*_TTTT_*.root",
        trig_mc+"*_WH_HToBB*.root", trig_mc+"*_WZTo*.root",
        trig_skim_mc+"*_ST_*.root", trig_skim_mc+"*_WJetsToLNu*.root",
        trig_mc+"*_ZH_HToBB*.root", trig_mc+"_ZZ_*.root"});

  auto data_2l = Proc<Baby_full>("Data", Process::Type::data, kBlack,
    {"/net/cms2/cms2r0/babymaker/babies/2016_05_31/data/singlelep/combined/*.root"},"pass&&(trig[20]||trig[22])"); // 0.589 fb-1
  vector<shared_ptr<Process> > full_trig_skim_2l = {data_2l, t1tttt_nc, t1tttt_c, tt1l_unskim, tt2l_unskim, dy_unskim, other_unskim};
  
  auto data_0l = Proc<Baby_full>("Data", Process::Type::data, kBlack,
    {"/net/cms2/cms2r0/babymaker/babies/2016_05_31/data/jetht/*.root"},"pass&&trig[12]"); // FIXME: trigger
  vector<shared_ptr<Process> > full_trig_skim_0l = {data_0l, t1tttt_nc, t1tttt_c, qcd_unskim, tt_unskim, other_unskim};

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


  string commoncut = "mt<140&&ht>500&&met>200";
  PlotMaker pm;  


  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////// ttbar 
  //////////////////////////////////////////////////////////////////////////////

  for(int ilep=2; ilep<3; ilep++) 
    {
      string lepcut;
      if(ilep==0)  lepcut=commoncut+"&&nels==1&&nmus==0"; 
      if(ilep==1)  lepcut=commoncut+"&&nels==0&&nmus==1"; 
      if(ilep==2)  lepcut=commoncut+"&&nleps==1"; 
 
      //
      // event-level variables except ones related to fatjet
      // 
      pm.AddPlot(HistoDef(15, 500, 2000., "ht", "H_{T} [GeV]",
			  lepcut+"&&njets>=5&&nbm>=1", "weight", {500.}),
		 full_trig_skim_1l, all_plot_types);
      pm.AddPlot(HistoDef(10, 200, 500., "met", "E_{T}^{miss} [GeV]",
			  lepcut+"&&njets>=5&&nbm>=1", "weight", {200., 350., 500.}),
		 full_trig_skim_1l, all_plot_types);
      pm.AddPlot(HistoDef(16, -0.5, 15.5, "njets", "N_{jets}",
			  lepcut+"&&nbm>=1", "weight", {5.5, 8.5}),
		 full_trig_skim_1l, all_plot_types);

      pm.AddPlot(HistoDef(16, 0., 280., "mt", "m_{T} [GeV]",
			  lepcut+"&&njets>=5&&nbm>=1", "weight", {140.}),
		 full_trig_skim_1l, all_plot_types); 

      //
      // fatjets
      //
      pm.AddPlot(HistoDef(10, 0, 500., "fjets14_m[0]", "m_{J1} [GeV]",
			  lepcut+"&&njets>=5&&nbm>=1", "weight", {999.}),
		 full_trig_skim_1l, all_plot_types); 
      pm.AddPlot(HistoDef(20, 0, 400., "fjets14_m[0]", "m_{J1} [GeV]",
			  lepcut+"&&njets<=5&&nbm>=1", "weight", {999.}),
		 full_trig_skim_1l, all_plot_types); 
  
      // event-level variables related to fatjet
      pm.AddPlot(HistoDef(10, 0., 1500., "mj14", "M_{J} [GeV]",
			  lepcut+"&&njets>=5&&nbm>=1", "weight", {400.}),
		 full_trig_skim_1l, all_plot_types);
  
    }  


  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////// Dilepton: Drelll-Yan and ttbar 
  //////////////////////////////////////////////////////////////////////////////

  string mll="(mumu_m*(mumu_m>0&&mumu_pt1>30)+elel_m*(elel_m>0&&elel_pt1>30))";
  string mllcut="(mumu_m*(mumu_m>0&&mumu_pt1>30)+elel_m*(elel_m>0&&elel_pt1>30))>80&&(mumu_m*(mumu_m>0&&mumu_pt1>30)+elel_m*(elel_m>0&&elel_pt1>30))<100";

  pm.AddPlot(HistoDef(10, 0., 1000., "mj14", "M_{J} [GeV]",
                      "nels==1&&nmus==1&&ht>300&&met>100&&njets>=2&&met<400&&nbm<3&&nbm>=1&&((leps_pt[0]>30&&leps_eta[0]<2.1&&leps_eta[0]>-2.1)||(leps_pt[1]>30&&leps_eta[1]<2.1&&leps_eta[1]>-2.1))", "weight", {400.}),
	     full_trig_skim_2l, all_plot_types);

  pm.AddPlot(HistoDef(10, 0., 1000., "mj14", "M_{J} [GeV]",
                      "(nmus>=2||nels>=2)&&ht>350&&njets>=4&&"+mllcut, "weight", {400.}),
	     full_trig_skim_2l, all_plot_types);

  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////// QCD
  //////////////////////////////////////////////////////////////////////////////

  pm.AddPlot(HistoDef(30, 0., 1500., "mj14", "M_{J} [GeV]",
                      "nleps==0&&met<50&&ht>1000&&njets>=9", "weight", {400.}),
	     full_trig_skim_0l, all_plot_types);

  //
  //
  //
  pm.MakePlots(lumi); 
}
