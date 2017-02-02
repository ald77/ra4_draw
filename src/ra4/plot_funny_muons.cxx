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
#include "core/functions.hpp"

using namespace std;
using namespace Functions;
using namespace PlotOptTypes;
NamedFunc::ScalarType DilepMass(const Baby &b);
NamedFunc::ScalarType DilepMassWindow(const Baby &b);
int main(){
  gErrorIgnoreLevel = 6000;
  double lumi =1.0;
  cout<<"Luminosity is "<<lumi<<endl;

  //Only produces plots that correspond to predefined lumis. Options are:
  //  vector<double> lumis = {1.0  ,  20.2        ,  20.2        , 7.72        ,  16.6};

  Palette colors("txt/colors.txt", "default");

  string rereco_dir = "/net/cms2/cms2r0/babymaker/babies/2017_01_27/data/merged_database_stdnj5/";
  string prompt_dir = "/net/cms2/cms2r0/babymaker/babies/2016_08_10/data/merged_database_standard/"; 
  string prompt_dir2 = "/net/cms2/cms2r0/babymaker/babies/2016_10_26/data/merged_database_standard/";
 
  string rereco_2l_dir = "/net/cms2/cms2r0/babymaker/babies/2017_01_27/data/merged_database_zcand/";
  string prompt_2l_dir = "/net/cms2/cms2r0/babymaker/babies/2016_08_10/data/merged_database_zcand/";
  string prompt_2l_dir2 = "/net/cms2/cms2r0/babymaker/babies/2016_10_26/data/merged_database_zcand/";

  string trig_mc        = "/net/cms2/cms2r0/babymaker/babies/2016_08_10/mc/unskimmed/";
  string trig_skim1l_mc = "/net/cms2/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_standard/";
  //string trig_skim0l_mc = "/net/cms2/cms2r0/babymaker/babies/2016_08_10/mc/merged_qcd/";
  string trig_skim2l_mc = "/net/cms2/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_zcand/";


  ///1l Data

  auto rereco = Process::MakeShared<Baby_full>("ReReco", Process::Type::data, kBlack,
    { rereco_dir+"*.root"},
    "trig_ra4");

 
  //Runs B through F, 20.236 ifb

  auto rerecoBF = Process::MakeShared<Baby_full>("ReReco, Runs B-F", Process::Type::data, kBlack,
    { rereco_dir+"*.root"},
    "trig_ra4 && run<=278808");


  //RunH, 8.496 ifb

  auto rerecoGH = Process::MakeShared<Baby_full>("ReReco, Run G-H", Process::Type::data, kBlack,
    { rereco_dir+"*RunG*.root",rereco_dir+"*RunH*.root"},
    "trig_ra4 && run>=278820"); 


  auto promptreco = Process::MakeShared<Baby_full>("PromptReco", Process::Type::background, kAzure+2,
    { prompt_dir+"*.root", prompt_dir2+"*.root"},
    "trig_ra4");

  auto promptrecoBF = Process::MakeShared<Baby_full>("PromptReco, Runs B-F", Process::Type::data, kBlack,
    { prompt_dir+"*.root", prompt_dir2+"*.root"},
    "trig_ra4 && run<=278808");

  auto promptrecoG = Process::MakeShared<Baby_full>("PromptReco, Run G", Process::Type::data, kBlack,
    { prompt_dir+"*.root", prompt_dir2+"*.root"},
    "trig_ra4 && run>=278820&&run<=280385");

  

  ///2l Data


  auto rereco_2l = Process::MakeShared<Baby_full>("ReReco", Process::Type::data, kBlack,
    { rereco_2l_dir+"*.root"},
    "trig_ra4");

 auto rereco_2l_no_bad_any = Process::MakeShared<Baby_full>("ReReco, no bad, duplicate, or bad track muons", Process::Type::data, kBlack,
    { rereco_2l_dir+"*.root"},
    "trig_ra4"&&n_mus_bad==0.&&n_mus_bad_dupl==0.&&n_mus_bad_trkmu==0.);

  auto rereco_2l_no_bad = Process::MakeShared<Baby_full>("ReReco, no bad or duplicate muons", Process::Type::data, kBlack,
    { rereco_2l_dir+"*.root"},
    "trig_ra4"&&n_mus_bad==0.&&n_mus_bad_dupl==0.);

  auto rereco_2l_no_bad_dupl = Process::MakeShared<Baby_full>("ReReco, no duplicate muons", Process::Type::data, kBlack,
    { rereco_2l_dir+"*.root"},
    "trig_ra4"&&n_mus_bad_dupl==0.);

  auto rereco_2l_no_bad_trk = Process::MakeShared<Baby_full>("ReReco, no bad track muons", Process::Type::data, kBlack,
    { rereco_2l_dir+"*.root"},
    "trig_ra4"&&n_mus_bad_trkmu==0.);



  auto rereco_2lBF = Process::MakeShared<Baby_full>("ReReco, Runs B-F", Process::Type::data, kBlack,
    { rereco_2l_dir+"*RunB*.root", rereco_2l_dir+"*RunC*.root", rereco_2l_dir+"*RunD*.root", rereco_2l_dir+"*RunE*.root", rereco_2l_dir+"*RunF*.root"},
    "trig_ra4 && run<=278808");

  auto rereco_2lGH = Process::MakeShared<Baby_full>("ReReco, Run G-H", Process::Type::data, kBlack,
    { rereco_2l_dir+"*RunG*.root", rereco_2l_dir+"*RunH*.root"},
    "trig_ra4 && run>=278820");


  auto promptreco_2l = Process::MakeShared<Baby_full>("PromptReco", Process::Type::background, kAzure+2,
    { prompt_2l_dir+"*.root", prompt_2l_dir2+"*.root"},
    "trig_ra4");

  auto promptreco_2lBF = Process::MakeShared<Baby_full>("PromptReco, Runs B-F", Process::Type::data, kBlack,
    { prompt_2l_dir+"*.root", prompt_2l_dir2+"*.root"},
    "trig_ra4 && run<=278808");

  auto promptreco_2lG = Process::MakeShared<Baby_full>("PromptReco, Run G", Process::Type::data, kBlack,
    { prompt_2l_dir+"*.root", prompt_2l_dir2+"*.root"},
    "trig_ra4 && run>=278820&&run<=280385");

  

  auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {trig_skim1l_mc+"*_TTJets*SingleLept*.root", trig_skim1l_mc+"*_TTJets_HT*.root"}, "ntruleps<=1&&stitch");
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {trig_skim1l_mc+"*_TTJets*DiLept*.root", trig_skim1l_mc+"*_TTJets_HT*.root"}, "ntruleps>=2&&stitch");
  auto wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {trig_skim1l_mc+"*_WJetsToLNu*.root"},"stitch");
  auto single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {trig_skim1l_mc+"*_ST_*.root"});
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {trig_skim1l_mc+"*_TTWJets*.root", trig_skim1l_mc+"*_TTZTo*.root"});
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {trig_skim1l_mc+"*DYJetsToLL*.root", trig_skim1l_mc+"*_QCD_HT*.root",
        trig_skim1l_mc+"*_ZJet*.root", trig_skim1l_mc+"*_WWTo*.root",
        trig_skim1l_mc+"*ggZH_HToBB*.root", trig_skim1l_mc+"*ttHJetTobb*.root",
        trig_skim1l_mc+"*_TTGJets*.root", trig_skim1l_mc+"*_TTTT_*.root",
        trig_skim1l_mc+"*_WH_HToBB*.root", trig_skim1l_mc+"*_WZTo*.root",
        trig_skim1l_mc+"*_ZH_HToBB*.root", trig_skim1l_mc+"_ZZ_*.root"});


  auto t1tttt_nc = Process::MakeShared<Baby_full>("T1tttt(1500,100)", Process::Type::signal, colors("t1tttt"),
    {trig_mc+"*SMS-T1tttt_mGluino-1500_mLSP-100*.root"});
  auto t1tttt_c = Process::MakeShared<Baby_full>("T1tttt(1200,800)", Process::Type::signal, colors("t1tttt"),
    {trig_mc+"*SMS-T1tttt_mGluino-1200_mLSP-800*.root"});
  t1tttt_c->SetLineStyle(2);



 //
  // 2-lepton plots
  //
  auto tt1l_2l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {trig_skim2l_mc+"*_TTJets*SingleLept*.root"},
    "ntruleps<=1");
  auto tt2l_2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {trig_skim2l_mc+"*_TTJets*DiLept*.root"},
    "ntruleps>=2");
  auto dy_2l = Process::MakeShared<Baby_full>("DY", Process::Type::background, colors("wjets"),
    {trig_skim2l_mc+"*DYJetsToLL_M-50_Tune*.root"});
  auto other_2l = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("ttv"),
    {trig_skim2l_mc+"*_TTWJet*root",  trig_skim2l_mc+"*_TTZTo*.root",
        trig_skim2l_mc+"*_ZJet*.root", trig_skim2l_mc+"*_WWTo*.root",
        trig_skim2l_mc+"*ggZH_HToBB*.root", trig_skim2l_mc+"*ttHJetTobb*.root",
        trig_skim2l_mc+"*_TTGJets*.root", trig_skim2l_mc+"*_TTTT_*.root",
        trig_skim2l_mc+"*_WH_HToBB*.root", trig_skim2l_mc+"*_WZTo*.root",
        trig_skim2l_mc+"*_ST_*.root", trig_skim2l_mc+"*_WJetsToLNu_Tune*.root",
        trig_skim2l_mc+"*_ZH_HToBB*.root", trig_skim2l_mc+"_ZZ_*.root",
        trig_skim2l_mc+"*_QCD_HT700to1000*.root",trig_skim2l_mc+"*_QCD_HT1000to1500*.root",
        trig_skim2l_mc+"*_QCD_HT1500to2000*.root",trig_skim2l_mc+"*_QCD_HT2000toInf*.root"});

  

  vector<shared_ptr<Process> > rr_mc_2l_BF = {rereco_2lBF, /*t1tttt_nc, t1tttt_c,*/ tt1l_2l, tt2l_2l, dy_2l, other_2l};
  vector<shared_ptr<Process> > rr_mc_2l_GH = {rereco_2lGH, /*t1tttt_nc, t1tttt_c,*/ tt1l_2l, tt2l_2l, dy_2l, other_2l};

  vector<shared_ptr<Process> > pr_mc_2l_BF = {promptreco_2lBF, /*t1tttt_nc, t1tttt_c,*/ tt1l_2l, tt2l_2l, dy_2l, other_2l};
  vector<shared_ptr<Process> > pr_mc_2l_G = {promptreco_2lG, /*t1tttt_nc, t1tttt_c,*/ tt1l_2l, tt2l_2l, dy_2l, other_2l};

  vector<shared_ptr<Process> > rr_mc_BF = {rerecoBF, /*t1tttt_nc, t1tttt_c,*/ tt1l, tt2l, wjets, single_t, ttv, other};
  vector<shared_ptr<Process> > rr_mc_GH = {rerecoGH, /*t1tttt_nc, t1tttt_c,*/ tt1l, tt2l, wjets, single_t, ttv, other};

  vector<shared_ptr<Process> > pr_mc_BF = {promptrecoBF, /*t1tttt_nc, t1tttt_c,*/ tt1l, tt2l, wjets, single_t, ttv, other};

  vector<shared_ptr<Process> > pr_mc_G = {promptrecoG, /*t1tttt_nc, t1tttt_c,*/ tt1l, tt2l, wjets, single_t, ttv, other};     


  vector<shared_ptr<Process> > rr_pr_2l = {rereco_2l, promptreco_2l};
  vector<shared_ptr<Process> > rr_pr_2l_no_bad_any = {rereco_2l_no_bad_any, promptreco_2l};
  vector<shared_ptr<Process> > rr_pr_2l_no_bad = {rereco_2l_no_bad, promptreco_2l};
  vector<shared_ptr<Process> > rr_pr_2l_no_bad_dupl = {rereco_2l_no_bad_dupl, promptreco_2l};
  vector<shared_ptr<Process> > rr_pr_2l_no_bad_trk = {rereco_2l_no_bad_trk, promptreco_2l};

  vector<shared_ptr<Process> > rr_pr = {rereco, promptreco};


  
  vector<double> lumis = { 1.0  ,  1.0  ,  1.0  ,  1.0  ,  1.0  ,  20.2        ,  20.2        , 7.72        ,  16.6};
  vector<string> tags =  {"data", "data_no_bad_any", "data_no_bad", "data_no_bad_dupl", "data_no_bad_trk", "prompt_mc_BF","rereco_mc_BF","prompt_mc_G","rereco_mc_GH"};
  vector<string> numers = {"ReReco","ReReco","ReReco","ReReco","ReReco","Data","Data","Data","Data"};
  vector<string> denoms = {"Prompt","Prompt","Prompt","Prompt","Prompt","MC","MC","MC","MC"};

  vector< vector<shared_ptr<Process> > > DY = {rr_pr_2l,rr_pr_2l_no_bad_any,rr_pr_2l_no_bad,rr_pr_2l_no_bad_dupl,rr_pr_2l_no_bad_trk, pr_mc_2l_BF,rr_mc_2l_BF,pr_mc_2l_G,rr_mc_2l_GH};
  // vector< vector<shared_ptr<Process> > > DY_mc_H = {pr_mc_2l_H,rr_mc_2l_H};

  vector< vector<shared_ptr<Process> > > single = {rr_pr,pr_mc_BF,rr_mc_BF,pr_mc_G,rr_mc_GH};
  // vector< vector<shared_ptr<Process> > > single_mc_H = {pr_mc_H,rr_mc_H};
  
  // vector< vector<shared_ptr<Process> > > DY = {rr_mc_2l,rr_pr_2l};
  // vector< vector<shared_ptr<Process> > > single = {rr_mc,rr_pr};
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
  PlotOpt log_lumi_shapes = log_lumi_info().Stack(StackType::lumi_shapes);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info);
  vector<PlotOpt> all_plot_types = {log_lumi_info, lin_lumi_info, log_lumi_shapes};

  //NamedFunc lowmtcut = "ht>500&&met>200&&mt<140";
  //  NamedFunc lowmjcut = "ht>500&&met>200&&mj14<400";
  PlotMaker pm;
  pm.min_print_ = true;

  // vector<TString> runs({"run>=272007&&run<=275376", "run>=275657&&run<=276283", "run>=276315&&run<=276811", 
  // 	"run>=276831&&run<=277420 ", "run>=277772&&run<=278808 ", "run>=278820&&run<=280385",   
  // 	"run>=280919&&run<=2830590"});
  // vector<TString> runNames({"RunB", "RunC", "RunD", "RunE", "RunF", "RunG", "RunH"});


  vector<NamedFunc> runrange = {"1","run<=278808","run>=278820&&run<=280385"};

  vector<NamedFunc> sels = {"nleps==1&&nveto==0&&st>500",
			  "nels==1&&nleps==1&&nveto==0&&st>500",
			  "nmus==1&&nleps==1&&nveto==0&&st>500",
			  "nleps==1&&nveto==0&&st>500&&mt>140",
			  "nels==1&&nleps==1&&nveto==0&&st>500&&mt>140",
			  "nmus==1&&nleps==1&&nveto==0&&st>500&&mt>140"};

  NamedFunc ll_m("ll_m", DilepMass);
  NamedFunc ll_window("ll_window", DilepMassWindow);
  vector<NamedFunc> DY_sels = {"nleps==2&&leps_pt[0]>40"&&ll_window,"nleps==2&&leps_pt[0]>40&&nmus==2&&mumu_m>80&&mumu_m<100","nleps==2&&leps_pt[0]>40&&nels==2&&elel_m>80&&elel_m<100"};
  vector<NamedFunc> njets_sels = {"1","njets==0","njets==1","njets>=2"};


  vector<NamedFunc> filters = {"pass","pass&&pass_ra2_badmu&&(met/met_calo)<5"};


  for(unsigned iproc=0;iproc<DY.size();iproc++){
    if(lumis[iproc]!=lumi) continue;
    for(unsigned int irun=0; irun < runrange.size();irun++) {
      //  if(iproc>0&&irun>0) continue;
      for(auto baseline: DY_sels){
	for(auto njets : njets_sels){
	  for(auto filter : filters){

	    continue;

	    pm.Push<Hist1D>(Axis(8, 70, 110., ll_m, "Dilepton invariant mass [GeV]",{80,100}),
			    runrange[irun]&&njets&&filter,
			    DY[iproc], all_plot_types).Tag(tags[iproc]).RatioTitle(numers[iproc],denoms[iproc]);

	    pm.Push<Hist1D>(Axis(17, 0, 850., "met", "E_{T}^{miss} [GeV]"),
			    baseline&&runrange[irun]&&njets&&filter,
			    DY[iproc], all_plot_types).Tag(tags[iproc]).RatioTitle(numers[iproc],denoms[iproc]);

	     pm.Push<Hist1D>(Axis(20, 0, 400., "met", "E_{T}^{miss} [GeV]"),
			    baseline&&runrange[irun]&&njets&&filter,
			    DY[iproc], all_plot_types).Tag(tags[iproc]+"_zoom").RatioTitle(numers[iproc],denoms[iproc]);

	    pm.Push<Hist1D>(Axis(17, 0, 850., "leps_pt", "Lepton p_{T} [GeV]"),
			    baseline&&runrange[irun]&&njets&&filter,
			    DY[iproc], all_plot_types).Tag(tags[iproc]).RatioTitle(numers[iproc],denoms[iproc]);
	    

	    pm.Push<Hist1D>(Axis(16, 0, 8., "met/met_calo", "E_{T}^{miss} / Calo MET"),
			    baseline&&runrange[irun]&&njets&&filter,
			    DY[iproc], all_plot_types).Tag(tags[iproc]).RatioTitle(numers[iproc],denoms[iproc]);

	    pm.Push<Hist1D>(Axis(16, 0, 8., "met/met_calo", "E_{T}^{miss} / Calo MET"),
			    baseline&&runrange[irun]&&njets&&filter&&"met>150",
			    DY[iproc], all_plot_types).Tag(tags[iproc]).RatioTitle(numers[iproc],denoms[iproc]);

	    pm.Push<Hist1D>(Axis(16, 0, 8., "met/met_calo", "E_{T}^{miss} / Calo MET"),
			    baseline&&runrange[irun]&&njets&&filter&&"met>100",
			    DY[iproc], all_plot_types).Tag(tags[iproc]).RatioTitle(numers[iproc],denoms[iproc]);

	  }
	}
      }
    }
    // pm.MakePlots(lumis[iproc]);

  };


  vector<NamedFunc> filtered = {/*"1","pass","pass&&!pass_ra2_badmu",*/"pass","pass&&pass_ra2_badmu&&(met/met_calo)<5"};
  for(unsigned iproc=0;iproc<single.size();iproc++){
    if(lumis[iproc]!=lumi) continue;
    for(unsigned int irun=0; irun < runrange.size();irun++) {
      //if(iproc>0&&irun>0) continue;
      for(auto baseline: sels){
	for(auto filter : filtered){

	  

	  	pm.Push<Hist1D>(Axis(16, 500, 2100., "st", "S_{T} [GeV]", {500.}),
			runrange[irun]&&baseline&&filter&&"met>200&&njets>=6&&nbm_moriond>=1",
			single[iproc], all_plot_types).Tag(tags[iproc]).RatioTitle(numers[iproc],denoms[iproc]);
	

	pm.Push<Hist1D>(Axis(20, 0, 1000., "leps_pt", "pT_{lep} [GeV]", {-999.}),
			runrange[irun]&&baseline&&filter&&"met>200&&njets>=6&&nbm_moriond>=1",
			single[iproc], all_plot_types).Tag(tags[iproc]).RatioTitle(numers[iproc],denoms[iproc]);
	  
	pm.Push<Hist1D>(Axis(16, 0, 8., "met/met_calo", "E_{T}^{miss} / Calo MET"),
			runrange[irun]&&baseline&&filter&&"met>200&&njets>=6&&nbm_moriond>=1",
			single[iproc], all_plot_types).Tag(tags[iproc]).RatioTitle(numers[iproc],denoms[iproc]);
	// pm.Push<Hist1D>(Axis(15, 0, 600., "fjets14_m", "m_{J} [GeV]", {-999}),
	// 		runrange[irun]&&baseline&&filter&&"met>200&&njets>=6&&nbm_moriond>=1",
	// 		single[iproc], all_plot_types).Tag(proctag);
	
	pm.Push<Hist1D>(Axis(24, 0, 1200., "mj14", "M_{J} [GeV]", {250.,400}),
			runrange[irun]&&baseline&&filter&&"met>200&&njets>=6&&nbm_moriond>=1",
			single[iproc], all_plot_types).Tag(tags[iproc]).RatioTitle(numers[iproc],denoms[iproc]);
	
	pm.Push<Hist1D>(Axis(15, 100, 850., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
			runrange[irun]&&baseline&&filter&&"met>100&&njets>=6&&nbm_moriond>=1",
			single[iproc], all_plot_types).Tag(tags[iproc]).RatioTitle(numers[iproc],denoms[iproc]);
	
	pm.Push<Hist1D>(Axis(35, 0, 700., "mt", "m_{T} [GeV]", {-999}),
			runrange[irun]&&baseline&&filter&&"met>200&&njets>=6&&nbm_moriond>=1",
			single[iproc], all_plot_types).Tag(tags[iproc]).RatioTitle(numers[iproc],denoms[iproc]);
	pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}", {5.5, 8.5}),
			runrange[irun]&&baseline&&filter&&"met>200&&nbm_moriond>=1",
			single[iproc], all_plot_types).Tag(tags[iproc]).RatioTitle(numers[iproc],denoms[iproc]);
	pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm_moriond", "N_{b}", {-999}),
			runrange[irun]&&baseline&&filter&&"met>200&&njets>=6",
			single[iproc], all_plot_types).Tag(tags[iproc]).RatioTitle(numers[iproc],denoms[iproc]);
	
	}
      }
    }
    // pm.MakePlots(lumis[iproc]);
  }

  pm.MakePlots(lumi);
 

}


NamedFunc::ScalarType DilepMass(const Baby &b){
  if (b.nmus()>=2) return b.mumu_m();
  if (b.nels()>=2) return b.elel_m();
  else return 0;

}

NamedFunc::ScalarType DilepMassWindow(const Baby &b){
  if (b.nmus()>=2&&b.mumu_m()>80&&b.mumu_m()<100) return true;
  if (b.nels()>=2&&b.elel_m()>80&&b.elel_m()<100) return true;
  else return false;

}
