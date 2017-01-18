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
  bool teddybear=true;
  double lumi = 1.0;

  Palette colors("txt/colors.txt", "default");

  //Can't use skims
  string rereco_dir = "/net/cms2/cms2r0/babymaker/babies/2016_11_08/data/skim_standard/";
  string prompt_dir = "/net/cms2/cms2r0/babymaker/babies/2016_08_10/data/skim_standard/"; 
  string prompt_dir2 = "/net/cms2/cms2r0/babymaker/babies/2016_10_26/data/skim_standard/";
 


  auto rereco = Process::MakeShared<Baby_full>("ReReco", Process::Type::data, kBlack,
    { rereco_dir+"*.root"},
    "pass && trig_ra4");

  auto promptreco = Process::MakeShared<Baby_full>("PromptReco", Process::Type::background, kAzure+2,
    { prompt_dir+"*.root", prompt_dir2+"*.root"},
    "pass && trig_ra4");

  auto rereco_singlelept = Process::MakeShared<Baby_full>("Moriond MC", Process::Type::data, kBlack,
    { "/net/cms29/cms29r0/heller/teddybear/unskimmed/*.root"},
    "pass");
  
  auto prompt_singlelept = Process::MakeShared<Baby_full>("Spring16 MC", Process::Type::background, kAzure+2,
    { "/net/cms2/cms2r0/babymaker/babies/2016_08_10/mc/unskimmed/*SingleLeptFromT_*.root"},
    "pass");
  
  vector<shared_ptr<Process> > procs;
  if(!teddybear) procs = {rereco, promptreco};
  else procs = {rereco_singlelept, prompt_singlelept};

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

  //NamedFunc lowmtcut = "ht>500&&met>200&&mt<140";
  //  NamedFunc lowmjcut = "ht>500&&met>200&&mj14<400";
  PlotMaker pm;

  
  vector<TString> runrange = {"mt<140&&run<=278808&&","mj14<400&&run<=278808&&","json12p9&&"};
  if(teddybear) runrange = {"1&&"};
  vector<TString> sels = {"nleps==1&&nveto==0&&st>500&&","nleps==1&&nveto==0&&st>500&&mt>140&&","nleps==1&&nveto==0&&st>500&&mj14>400&&",
			  "nels==1&&nleps==1&&nveto==0&&st>500&&","nels==1&&nleps==1&&nveto==0&&st>500&&mt>140&&","nels==1&&nleps==1&&nveto==0&&st>500&&mj14>400&&",
			  "nmus==1&&nleps==1&&nveto==0&&st>500&&","nmus==1&&nleps==1&&nveto==0&&st>500&&mt>140&&","nmus==1&&nleps==1&&nveto==0&&st>500&&mj14>400&&"};

  vector<TString> filters = {"1","pass_ra2_badmu","pass_ra2_badmu&&(met/met_calo)<5"};
  if(teddybear) filters = {"1"};
  for(auto irun : runrange) {
    //string baseline = "nleps==1&&nveto==0&&st>500";
    for(auto baseline: sels){
      if(baseline.Contains("mt>140") && irun.Contains("mt<140")) continue;
      for(auto filter : filters){
	pm.Push<Hist1D>(Axis(16, 400, 2000., "ht", "H_{T} [GeV]", {500.}),
			irun+baseline+filter+"&&met>200&&njets>=6&&nbm>=1",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(16, 500, 2100., "st", "S_{T} [GeV]", {500.}),
			irun+baseline+filter+"&&met>200&&njets>=6&&nbm>=1",
			procs, all_plot_types);

	pm.Push<Hist1D>(Axis(34, 0, 1020., "jets_pt", "pT_{jets} [GeV]", {-999}),
			irun+baseline+filter+"&&met>200&&njets>=6&&nbm>=1&&!jets_islep",
			procs, all_plot_types);

        pm.Push<Hist1D>(Axis(33, 30, 1020., "jets_pt", "pT_{jets} [GeV]", {-999}),
                        irun+baseline+filter+"&&met>200&&njets>=6&&nbm>=1&&!jets_islep&&jets_pt>30",
                        procs, all_plot_types);

	pm.Push<Hist1D>(Axis(20, 0, 1000., "leps_pt", "pT_{lep} [GeV]", {-999.}),
			irun+baseline+filter+"&&met>200&&njets>=6&&nbm>=1",
			procs, all_plot_types);


	pm.Push<Hist1D>(Axis(15, 0, 600., "fjets14_m", "m_{J} [GeV]", {-999}),
			irun+baseline+filter+"&&met>200&&njets>=6&&nbm>=1",
			procs, all_plot_types);
	// pm.Push<Hist1D>(Axis(15, 0, 600., "fjets14_m[0]", "m_{J1} [GeV]", {-999}),
	// 		irun+baseline+filter+"&&met>200&&njets>=6&&nbm>=1",
	// 		procs, all_plot_types);
	pm.Push<Hist1D>(Axis(24, 0, 1200., "mj14", "M_{J} [GeV]", {250.,400}),
			irun+baseline+filter+"&&met>200&&njets>=6&&nbm>=1",
			procs, all_plot_types);
	if(!teddybear){
	  pm.Push<Hist1D>(Axis(14, 150, 850., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
			  irun+baseline+filter+"&&met>150&&njets>=6&&nbm>=1",
			  procs, all_plot_types);}
	
	else{
	  pm.Push<Hist1D>(Axis(17, 0, 850., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
			irun+baseline+filter+"&&njets>=6&&nbm>=1",
			procs, all_plot_types);

	}
	pm.Push<Hist1D>(Axis(35, 0, 700., "mt", "m_{T} [GeV]", {-999}),
			irun+baseline+filter+"&&met>200&&njets>=6&&nbm>=1",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}", {5.5, 8.5}),
			irun+baseline+filter+"&&met>200&&nbm>=1",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm", "N_{b}", {-999}),
			irun+baseline+filter+"&&met>200&&njets>=6",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbl", "N_{l}", {-999}),
			irun+baseline+filter+"&&met>200&&njets>=6",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbt", "N_{t}", {-999}),
			irun+baseline+filter+"&&met>200&&njets>=6",
			procs, all_plot_types);

	pm.Push<Hist1D>(Axis(20, 0.46, 1, "jets_csv", "CSV", {0.8,0.935}),
			irun+baseline+filter+"&&met>200&&njets>=6&&jets_csv>0.46&&jets_pt>30&&!jets_islep",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(40, 0, 1, "jets_csv", "CSV", {0.46,0.8,0.935}),
			irun+baseline+filter+"&&met>200&&njets>=6&&jets_pt>30&&!jets_islep",
			procs, all_plot_types);

	pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nels", "N_{e}", {0.5}),
			irun+filter+"&&pass&&nleps>=1&&st>500&&met>200&&(nleps>1||mt<140)&&njets>=6&&nbm>=1",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nmus", "N_{mu}", {0.5}),
			irun+filter+"&&pass&&nleps>=1&&st>500&&met>200&&(nleps>1||mt<140)&&njets>=6&&nbm>=1",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nleps", "N_{leps}", {0.5}),
			irun+filter+"&&pass&&nleps>=1&&st>500&&met>200&&(nleps>1||mt<140)&&njets>=6&&nbm>=1",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nveto", "N_{veto}", {0.5}),
			irun+filter+"&&pass&&nleps==1&&st>500&&met>200&&mt<140&&njets>=6&&nbm>=1",
			procs, all_plot_types);
      }
    }
  }

  if(teddybear){

    	pm.Push<Hist1D>(Axis(40, 0, 2000., "ht", "H_{T} [GeV]", {500.}),
			"1",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(40, 0, 2000., "st", "S_{T} [GeV]", {500.}),
			"1",
			procs, all_plot_types);

	pm.Push<Hist1D>(Axis(34, 0, 1020., "jets_pt", "pT_{jets} [GeV]", {-999}),
			"!jets_islep",
			procs, all_plot_types);

	pm.Push<Hist1D>(Axis(33, 30, 1020., "jets_pt", "pT_{jets} [GeV]", {-999}),
                        "!jets_islep&&jets_pt>30",
                        procs, all_plot_types);


	pm.Push<Hist1D>(Axis(50, 0, 1000., "leps_pt", "pT_{lep} [GeV]", {-999.}),
			"1",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(50, 0, 1000., "els_pt", "pT_{lep} [GeV]", {-999.}),
			"1",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(50, 0, 1000., "mus_pt", "pT_{lep} [GeV]", {-999.}),
			"1",
			procs, all_plot_types);


	pm.Push<Hist1D>(Axis(30, 0, 600., "fjets14_m", "m_{J} [GeV]", {-999}),
			"1",
			procs, all_plot_types);
	//	pm.Push<Hist1D>(Axis(30, 0, 600., "fjets14_m[0]", "m_{J1} [GeV]", {-999}),
	///		"1",
	//		procs, all_plot_types);
	pm.Push<Hist1D>(Axis(24, 0, 1200., "mj14", "M_{J} [GeV]", {250.,400}),
			"1",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(32, 0, 800., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
			"1",
			procs, all_plot_types);

	pm.Push<Hist1D>(Axis(32, 0, 800., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
			"mt>140",
			procs, all_plot_types);

	pm.Push<Hist1D>(Axis(32, 0, 800., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
			"nleps&&1&&nmus==1&&mt>140",
			procs, all_plot_types);

	pm.Push<Hist1D>(Axis(32, 0, 800., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
			"nleps&&1&&nels==1&&mt>140",
			procs, all_plot_types);
	
	pm.Push<Hist1D>(Axis(35, 0, 700., "mt", "m_{T} [GeV]", {-999}),
			"1",
			procs, all_plot_types);

	pm.Push<Hist1D>(Axis(35, 0, 700., "mt", "m_{T} [GeV]", {-999}),
			"nleps==1&&nels==1",
			procs, all_plot_types);

	pm.Push<Hist1D>(Axis(35, 0, 700., "mt", "m_{T} [GeV]", {-999}),
			"nleps==1&&nmus==1",
			procs, all_plot_types);

	pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}", {5.5, 8.5}),
			"1",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm", "N_{b}", {-999}),
			"1",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbl", "N_{l}", {-999}),
			"1",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbt", "N_{t}", {-999}),
			"1",
			procs, all_plot_types);

	pm.Push<Hist1D>(Axis(20, 0.46, 1, "jets_csv", "CSV", {0.8,0.935}),
			"jets_csv>0.46&&jets_pt>30&&!jets_islep",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(40, 0, 1, "jets_csv", "CSV", {0.46,0.8,0.935}),
			"jets_pt>30&&!jets_islep",
			procs, all_plot_types);

	pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nels", "N_{e}", {0.5}),
			"1",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nmus", "N_{mu}", {0.5}),
			"1",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nleps", "N_{leps}", {0.5}),
			"1",
			procs, all_plot_types);
	pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nveto", "N_{veto}", {0.5}),
			"1",
			procs, all_plot_types);







  }


  /*for(int ilep=0; ilep<3; ilep++){
    NamedFunc lepcut(true);
    if(ilep==0)  lepcut="nels==1&&nmus==0";
    if(ilep==1)  lepcut="nels==0&&nmus==1";
    if(ilep==2)  lepcut="nleps==1";

    for(int inveto=0; inveto<2; inveto++){
    NamedFunc nvetocut(true);
    if(inveto==0)  nvetocut="nveto==0";
    if(inveto==1)  nvetocut="1";
    if(inveto==1)  continue;

    for(int injets=0; injets<2; injets++){
    NamedFunc njetscut(true);
        if(injets==0)  njetscut="njets>=5";
        if(injets==1)  njetscut="njets>=6";
        if(injets==0)  continue;

        for(int inb=0; inb<2; inb++){
          NamedFunc nbcut(true);
          if(inb==0)  nbcut="nbm>=1";
          if(inb==1)  nbcut="nbm>=2";
          if(inb==1)  continue;

          //
          // event-level variables except ones related to fatjet
          //
          pm.Push<Hist1D>(Axis(15, 500, 2000., "ht", "H_{T} [GeV]", {500.}),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
                          lepcut&&lowmjcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}", {5.5, 8.5}),
                          lepcut&&lowmtcut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm", "N_{b}", {0.5, 1.5, 2.5}),
                          lepcut&&lowmtcut&&njetscut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);


          if(ilep==0){
            // only barrel or endcap electron
            pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
                            lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"leps_eta[0]<1.479&&leps_eta[0]>-1.479",
                            full_trig_skim_1l, all_plot_types);
            pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
                            lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"leps_eta[0]>1.479||leps_eta[0]<-1.479",
                            full_trig_skim_1l, all_plot_types);
          }

          //
          // leptons
          //
          if(ilep==0){// electron
            pm.Push<Hist1D>(Axis(10, 0, 200., "leps_pt[0]", "p_{T}(electron) [GeV]"),
                            lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                            full_trig_skim_1l, all_plot_types);
            pm.Push<Hist1D>(Axis(10, -2.5, 2.5, "leps_eta[0]", "#eta(electron) [GeV]"),
                            lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                            full_trig_skim_1l, all_plot_types);
          }else if(ilep==1){// muon
            pm.Push<Hist1D>(Axis(10, 0, 200., "leps_pt[0]", "p_{T}(muon) [GeV]"),
                            lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                            full_trig_skim_1l, all_plot_types);
            pm.Push<Hist1D>(Axis(10, -2.5, 2.5, "leps_eta[0]", "#eta(muon) [GeV]"),
                            lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                            full_trig_skim_1l, all_plot_types);
          }else if(ilep==2){  // electron or muon
            pm.Push<Hist1D>(Axis(10, 0, 200., "leps_pt[0]", "p_{T}(lepton) [GeV]"),
                            lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                            full_trig_skim_1l, all_plot_types);
            pm.Push<Hist1D>(Axis(10, -2.5, 2.5, "leps_eta[0]", "#eta(lepton) [GeV]"),
                            lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                            full_trig_skim_1l, all_plot_types);
          }

          //
          // jets
          //
          pm.Push<Hist1D>(Axis(10, 30, 630., "jets_pt[0]", "p_{T}(leading jet) [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"!jets_islep[0]",
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, -2.5, 2.5, "jets_eta[0]", "#eta(leading jet) [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"!jets_islep[0]",
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 0, 200., "jets_m[0]", "mass(leading jet) [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"!jets_islep[0]",
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 30, 630., "jets_pt", "p_{T}(jet) [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"!jets_islep",
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, -2.5, 2.5, "jets_eta", "#eta(jet) [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"!jets_islep",
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 0, 200., "jets_m", "m(jet) [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"!jets_islep",
                          full_trig_skim_1l, all_plot_types);

          //
          // fatjets
          //
          pm.Push<Hist1D>(Axis(10, 0, 1000., "fjets14_pt[0]", "p_{T}(J1) [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 0, 500., "fjets14_m[0]", "m_{J1} [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 0, 500., "fjets14_m[0]", "m_{J1} [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"leps_pt[0]>160",
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 0, 500., "fjets14_m[0]", "m_{J1} [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"leps_pt[0]<=160",
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, -2.5, 2.5, "fjets14_eta[0]", "#eta(J1)"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(7, 0.5, 7.5, "fjets14_nconst[0]", "N_{constituents}(J1)"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 0, 1000., "fjets14_pt", "p_{T}(J) [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 0, 500., "fjets14_m", "m_{J} [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, -2.5, 2.5, "fjets14_eta", "#eta(J)"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(7, 0.5, 7.5, "fjets14_nconst", "N_{constituents}(J)"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);

          pm.Push<Hist1D>(Axis(20, 0, 400., "fjets14_m[0]", "m_{J1} [GeV]"),
                          lepcut&&lowmtcut&&"njets<=5&&nbm>=1",
                          full_trig_skim_1l, all_plot_types);

          // event-level variables related to fatjet
          pm.Push<Hist1D>(Axis(15, 0., 1500., "mj14", "M_{J} [GeV]", {400.}),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(7, 0.5, 7.5, "nfjets14", "N_{J} [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);

        } //for(int inb=0; inb<2; inb++)
      } //for(int injets=0; injets<2; injets++)
    } //for(int inveto=0; inveto<2; inveto++)
  }
  */
  
  pm.min_print_ = true;
  pm.MakePlots(lumi);

}
