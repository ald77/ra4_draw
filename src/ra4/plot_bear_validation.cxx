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
  double lumi = 35.9;

  Palette colors("txt/colors.txt", "default");

  string bear_dir = "/net/cms2/cms2r0/babymaker/babies/2017_01_27/data/merged_database_stdnj5/";
  string capy_dir = "/net/cms2/cms2r0/babymaker/babies/2016_11_08/data/merged_database_standard/"; 

  bool noskim=false;
  if(noskim) {
     bear_dir = "/net/cms2/cms2r0/babymaker/babies/2017_01_27/data/merged_database_met100/";
     capy_dir = "/net/cms2/cms2r0/babymaker/babies/2016_11_08/data/skim_met100/"; 
  }
  string bear_mc = "/net/cms2/cms2r0/babymaker/babies/2017_01_27/mc/merged_mcbase_stdnj5/";



  auto bear = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    { bear_dir+"*.root"},
    "pass && trig_ra4");

  auto capy = Process::MakeShared<Baby_full>("Data with ICHEP JEC", Process::Type::background, kAzure+2,
    { capy_dir+"*.root"},
    "pass && trig_ra4");

  auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {bear_mc+"*_TTJets*SingleLept*.root"}, "ntruleps<=1&&stitch_met");
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {bear_mc+"*_TTJets*DiLept*.root"}, "ntruleps>=2&&stitch_met");
  auto wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {bear_mc+"*_WJetsToLNu*.root"},"stitch");
  auto single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {bear_mc+"*_ST_*.root"});
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {bear_mc+"*_TTWJets*.root", bear_mc+"*_TTZTo*.root"});
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {bear_mc+"*DYJetsToLL*.root", bear_mc+"*_QCD_HT*.root",
        bear_mc+"*_ZJet*.root", bear_mc+"*_WWTo*.root",
        bear_mc+"*ggZH_HToBB*.root", bear_mc+"*ttHJetTobb*.root",
        bear_mc+"*_TTGJets*.root", bear_mc+"*_TTTT_*.root",
        bear_mc+"*_WH_HToBB*.root", bear_mc+"*_WZTo*.root",
        bear_mc+"*_ZH_HToBB*.root", bear_mc+"_ZZ_*.root"});


  vector<shared_ptr<Process> > procs;
  procs = {bear, capy};

  vector< vector<shared_ptr<Process> > > proc_pairs = {{bear,tt1l, tt2l, wjets, single_t, ttv, other}/*,{bear,capy}*/};
  vector<string> proc_tags = {"bear_data_mc","bear_capy_data"};

  vector<string> numers = {"Data","Moriond JEC"};
  vector<string> denoms = {"MC","ICHEP JEC"};

 //procs = {rereco_singlelept, prompt_singlelept};

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

 // if ("Run2016F1" in taskname): config.Data.runRange = '277772-278801'
 // 	    elif ("Run2016F2" in taskname): config.Data.runRange = '278802-278808'
// vector<TString> runs({"run>=272007&&run<=275376", "run>=275657&&run<=276283", "run>=276315&&run<=276811", 
  // 	"run>=276831&&run<=277420 ", "run>=277772&&run<=278808 ", "run>=278820&&run<=280385",   
  // 	"run>=280919&&run<=2830590"});
  // vector<TString> runNames({"RunB", "RunC", "RunD", "RunE", "RunF", "RunG", "RunH"});
  

  vector<TString> runs = {""/*,"run<=276811&&","run>=276831&&run<=278801&&","run>=278802&&run<=280385&&","run>=280919&&run<284044&&"*/};
  vector<TString> runNames={"All","Runs B-D" ,"Runs E-F1"               ,"Runs F2-G"             , "Run H"                 };

 
  vector<TString> sels = {"nleps==1&&nveto==0&&st>500&&"/*,"nleps==1&&nveto==0&&st>500&&mt>140&&","nleps==1&&nveto==0&&st>500&&mj14>400&&",*/
			  ,"nels==1&&nleps==1&&nveto==0&&st>500&&"/*,"nels==1&&nleps==1&&nveto==0&&st>500&&mt>140&&","nels==1&&nleps==1&&nveto==0&&st>500&&mj14>400&&",*/
			  ,"nmus==1&&nleps==1&&nveto==0&&st>500&&"/*,"nmus==1&&nleps==1&&nveto==0&&st>500&&mt>140&&","nmus==1&&nleps==1&&nveto==0&&st>500&&mj14>400&&"*/};

  vector<TString> filters = {"pass&&pass_ra2_badmu&&(met/met_calo)<5"};
  for(unsigned int i=0; i<proc_pairs.size();i++){
    for(auto irun : runs) {
      for(auto baseline: sels){
	for(auto filter : filters){

	  if(!noskim){
	    pm.Push<Hist1D>(Axis(16, 400, 2000., "ht", "H_{T} [GeV]", {500.}),
			  irun+baseline+filter+"&&met>200&&njets>=6&&nbm_moriond>=1",
			proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);
	pm.Push<Hist1D>(Axis(16, 500, 2100., "st", "S_{T} [GeV]", {500.}),
			irun+baseline+filter+"&&met>200&&njets>=6&&nbm_moriond>=1",
			proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);

	pm.Push<Hist1D>(Axis(34, 0, 1020., "jets_pt", "pT_{jets} [GeV]", {-999}),
			irun+baseline+filter+"&&met>200&&njets>=6&&nbm_moriond>=1&&!jets_islep",
			proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);

        pm.Push<Hist1D>(Axis(33, 30, 1020., "jets_pt", "pT_{jets} [GeV]", {-999}),
                        irun+baseline+filter+"&&met>200&&njets>=6&&nbm_moriond>=1&&!jets_islep&&jets_pt>30",
                        proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);

	pm.Push<Hist1D>(Axis(20, 0, 1000., "leps_pt", "pT_{lep} [GeV]", {-999.}),
			irun+baseline+filter+"&&met>200&&njets>=6&&nbm_moriond>=1",
			proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);


	pm.Push<Hist1D>(Axis(15, 0, 600., "fjets14_m", "m_{J} [GeV]", {-999}),
			irun+baseline+filter+"&&met>200&&njets>=6&&nbm_moriond>=1",
			proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);
	// pm.Push<Hist1D>(Axis(15, 0, 600., "fjets14_m[0]", "m_{J1} [GeV]", {-999}),
	// 		irun+baseline+filter+"&&met>200&&njets>=6&&nbm_moriond>=1",
	// 		proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);
	pm.Push<Hist1D>(Axis(24, 0, 1200., "mj14", "M_{J} [GeV]", {250.,400}),
			irun+baseline+filter+"&&met>200&&njets>=6&&nbm_moriond>=1",
			proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);

	pm.Push<Hist1D>(Axis(15, 100, 850., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
			irun+baseline+filter+"&&met>100&&njets>=6&&nbm_moriond>=1",
			proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);
	
	pm.Push<Hist1D>(Axis(20, 0, 700., "mt", "m_{T} [GeV]", {-999}),
			irun+baseline+filter+"&&met>200&&njets>=6&&nbm_moriond>=1",
			proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);
	pm.Push<Hist1D>(Axis(11, 4.5, 15.5, "njets", "N_{jets}", {5.5, 8.5}),
			irun+baseline+filter+"&&met>200&&nbm_moriond>=1&&njets>=5",
			proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);
	pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm_moriond", "N_{b}", {-999}),
			irun+baseline+filter+"&&met>200&&njets>=6",
			proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);

	pm.Push<Hist1D>(Axis(20, 0.54, 1, "jets_csv", "CSV", {0.8484,0.9535}),
			irun+baseline+filter+"&&met>200&&njets>=6&&jets_csv>0.54&&jets_pt>30&&!jets_islep",
			proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);
	pm.Push<Hist1D>(Axis(40, 0, 1, "jets_csv", "CSV",{0.5426,0.8484,0.9535}),
			irun+baseline+filter+"&&met>200&&njets>=6&&jets_pt>30&&!jets_islep",
			proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);

	pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nels", "N_{e}", {0.5}),
			irun+filter+"&&pass&&nleps>=1&&st>500&&met>200&&njets>=6&&nbm_moriond>=1",
			proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);
	pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nmus", "N_{mu}", {0.5}),
			irun+filter+"&&pass&&nleps>=1&&st>500&&met>200&&njets>=6&&nbm_moriond>=1",
			proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);
	pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nleps", "N_{leps}", {0.5}),
			irun+filter+"&&pass&&nleps>=1&&st>500&&met>200&&njets>=6&&nbm_moriond>=1",
			proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);
	pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nveto", "N_{veto}", {0.5}),
			irun+filter+"&&pass&&nleps==1&&st>500&&met>200&&njets>=6&&nbm_moriond>=1",
			proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);
	}
	
	else{
	  pm.Push<Hist1D>(Axis(20, 0, 2000., "ht", "H_{T} [GeV]", {500.}),
			irun+filter+"&&nleps==1&&nveto==0&&met>200&&njets>=6&&nbm_moriond>=1",
			proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);

	  pm.Push<Hist1D>(Axis(17, 0, 850., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
			  irun+baseline+filter+"&&njets>=6&&nbm_moriond>=1",
			  proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);




	  pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}", {5.5, 8.5}),
			  irun+baseline+filter+"&&met>200&&nbm_moriond>=1",
			  proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);
	  pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}", {5.5, 8.5}),
			  irun+baseline+filter+"&&met>100&&nbm_moriond>=1",
			  proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);


	  pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nleps", "N_{leps}", {0.5}),
			  irun+filter+"&&pass&&st>500&&met>200&&njets>=6&&nbm_moriond>=1",
			  proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);


	  pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nels", "N_{e}", {0.5}),
			  irun+filter+"&&pass&&st>500&&met>200&&njets>=6&&nbm_moriond>=1",
			  proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);

	  pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nmus", "N_{mu}", {0.5}),
			  irun+filter+"&&pass&&st>500&&met>200&&njets>=6&&nbm_moriond>=1",
			proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]);



	}

      }
    }
    if(noskim){
	 
      pm.Push<Hist1D>(Axis(24, 0, 1200., "mj14", "M_{J} [GeV]", {250.,400}),
		      irun+"1",
		      proc_pairs[i], all_plot_types).RatioTitle("Moriond JEC","ICHEP JEC").Tag(proc_tags[i]);
      pm.Push<Hist1D>(Axis(20, 0, 700., "mt", "m_{T} [GeV]", {-999}),
		      irun+"1",
		      proc_pairs[i], all_plot_types).RatioTitle("Moriond JEC","ICHEP JEC").Tag(proc_tags[i]);
      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nmus", "N_{mu}", {0.5}),
		      irun+"1",
		      proc_pairs[i], all_plot_types).RatioTitle("Moriond JEC","ICHEP JEC").Tag(proc_tags[i]);	 
      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nels", "N_{e}", {0.5}),
		      irun+"1",
		      proc_pairs[i], all_plot_types).RatioTitle("Moriond JEC","ICHEP JEC").Tag(proc_tags[i]);
      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nleps", "N_{leps}", {0.5}),
		      irun+"1",
		      proc_pairs[i], all_plot_types).RatioTitle("Moriond JEC","ICHEP JEC").Tag(proc_tags[i]);
      pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}", {5.5, 8.5}),
		      irun+"1",
		      proc_pairs[i], all_plot_types).RatioTitle("Moriond JEC","ICHEP JEC").Tag(proc_tags[i]);
      pm.Push<Hist1D>(Axis(17, 0, 850., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
		      irun+"1",
		      proc_pairs[i], all_plot_types).RatioTitle("Moriond JEC","ICHEP JEC").Tag(proc_tags[i]);
      pm.Push<Hist1D>(Axis(20, 0, 2000., "ht", "H_{T} [GeV]", {500.}),
		      irun+"1",
		      proc_pairs[i], all_plot_types).RatioTitle("Moriond JEC","ICHEP JEC").Tag(proc_tags[i]);
      pm.Push<Hist1D>(Axis(20, 0, 2000., "st", "S_{T} [GeV]", {500.}),
		      irun+"1",
		      proc_pairs[i], all_plot_types).RatioTitle("Moriond JEC","ICHEP JEC").Tag(proc_tags[i]);

      pm.Push<Hist1D>(Axis(20, 0, 1000., "leps_pt", "pT_{lep} [GeV]", {-999.}),
		      irun+"1",
		      proc_pairs[i], all_plot_types).RatioTitle("Moriond JEC","ICHEP JEC").Tag(proc_tags[i]);
      pm.Push<Hist1D>(Axis(33, 30, 1020., "jets_pt", "pT_{jets} [GeV]", {-999}),
		      irun+"1",
		      proc_pairs[i], all_plot_types).RatioTitle("Moriond JEC","ICHEP JEC").Tag(proc_tags[i]);
    }

  }



 
  }
  pm.min_print_ = true;
  pm.MakePlots(lumi);

}
