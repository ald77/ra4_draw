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
  double lumi = 1.0;

  Palette colors("txt/colors.txt", "default");

 

  string bear_mc_dir = "/net/cms2/cms2r0/babymaker/babies/2017_01_27/mc/merged_mcbase_met100/";
  string capy_mc_dir = "/net/cms2/cms2r0/babymaker/babies/2016_08_10/mc/skim_met100/";

 bool noskim=false;
 if(noskim) {
   bear_mc_dir = "/net/cms2/cms2r0/babymaker/babies/2017_01_27/mc/unskimmed/";
   capy_mc_dir = "/net/cms2/cms2r0/babymaker/babies/2016_11_08/mc/unskimmed/"; 
 }



  auto bear_met_ttbar = Process::MakeShared<Baby_full>("Bear MC, MET-binned ttbar", Process::Type::data, kBlack,
    {bear_mc_dir+"*_TTJets*Lept*.root"},
    "stitch_met&&pass");



 // auto bear_inclusive_ttbar = Process::MakeShared<Baby_full>("Bear MC, inclusive ttbar", Process::Type::data, kBlack,
 //     {bear_mc_dir+"*_TTJets_Tune*"},
 //    "ntruleps>=1");

  auto bear_HT_ttbar_dat = Process::MakeShared<Baby_full>("Bear MC, HT-binned ttbar", Process::Type::data, kBlack,
    {bear_mc_dir+"*_TTJets*DiLept_Tune*.root",bear_mc_dir+"*_TTJets*LeptFromT_Tune*.root",bear_mc_dir+"*_TTJets*LeptFromTbar_Tune*.root",bear_mc_dir+"*_TTJets*HT*"},
    "stitch&&pass&&ntruleps>=1");

  auto bear_HT_ttbar = Process::MakeShared<Baby_full>("Bear MC, HT-binned ttbar", Process::Type::background, kAzure+2,
    {bear_mc_dir+"*_TTJets*DiLept_Tune*.root",bear_mc_dir+"*_TTJets*LeptFromT_Tune*.root",bear_mc_dir+"*_TTJets*LeptFromTbar_Tune*.root",bear_mc_dir+"*_TTJets*HT*"},
    "stitch&&pass&&ntruleps>=1");

  auto capy_HT_ttbar = Process::MakeShared<Baby_full>("Capybara MC, HT-binned ttbar", Process::Type::background, kAzure+2,
    {capy_mc_dir+"*_TTJets*DiLept_Tune*.root",capy_mc_dir+"*_TTJets*LeptFromT_Tune*.root",capy_mc_dir+"*_TTJets*LeptFromTbar_Tune*.root",capy_mc_dir+"*_TTJets*HT*"},
    "stitch&&pass&&ntruleps>=1&&w_btag>0.");


  auto bear_other = Process::MakeShared<Baby_full>("Bear MC, non-ttbar", Process::Type::data, kBlack,
    {	bear_mc_dir+"*_WJetsToLNu*.root",bear_mc_dir+"*_ST_*.root",
	bear_mc_dir+"*_TTWJets*.root", bear_mc_dir+"*_TTZTo*.root",
	bear_mc_dir+"*DYJetsToLL*.root", bear_mc_dir+"*_QCD_HT*.root",
        bear_mc_dir+"*_ZJet*.root", bear_mc_dir+"*_WWTo*.root",
        bear_mc_dir+"*ggZH_HToBB*.root", bear_mc_dir+"*ttHJetTobb*.root",
        bear_mc_dir+"*_TTGJets*.root", bear_mc_dir+"*_TTTT_*.root",
        bear_mc_dir+"*_WH_HToBB*.root", bear_mc_dir+"*_WZTo*.root",
	bear_mc_dir+"*_ZH_HToBB*.root", bear_mc_dir+"_ZZ_*.root"},
    "stitch&&pass");

  auto capy_other = Process::MakeShared<Baby_full>("Capybara MC, non-ttbar", Process::Type::background, kAzure+2,
    {	capy_mc_dir+"*_WJetsToLNu*.root",capy_mc_dir+"*_ST_*.root",
	capy_mc_dir+"*_TTWJets*.root", capy_mc_dir+"*_TTZTo*.root",
	capy_mc_dir+"*DYJetsToLL*.root", capy_mc_dir+"*_QCD_HT*.root",
        capy_mc_dir+"*_ZJet*.root", capy_mc_dir+"*_WWTo*.root",
        capy_mc_dir+"*ggZH_HToBB*.root", capy_mc_dir+"*ttHJetTobb*.root",
        capy_mc_dir+"*_TTGJets*.root", capy_mc_dir+"*_TTTT_*.root",
        capy_mc_dir+"*_WH_HToBB*.root", capy_mc_dir+"*_WZTo*.root",
	capy_mc_dir+"*_ZH_HToBB*.root", capy_mc_dir+"_ZZ_*.root"},
    "stitch&&pass");


  // vector< vector<shared_ptr<Process> > > proc_pairs = {{bear_met_ttbar,bear_HT_ttbar},{bear_inclusive_ttbar,bear_HT_ttbar}};
  // vector<string> proc_tags = {"noskim_bear_bear_ttbar","noskim_inclu_bear_bear_ttbar"};
  // vector<string> numers = {"MET-binned","Inclusive"};
  // vector<string> denoms = {"HT-binned","HT-binned"};
 

  vector< vector<shared_ptr<Process> > > proc_pairs = {{bear_met_ttbar,capy_HT_ttbar},{bear_met_ttbar,bear_HT_ttbar},{bear_HT_ttbar_dat,capy_HT_ttbar},{bear_other,capy_other}};
  vector<string> proc_tags = {"bear_capy_ttbar","bear_bear_ttbar","bear_capy_HT_ttbar","bear_capy_other"};

  vector<string> numers = {"Bear","MET-binned","Bear HT-binned","Bear"};
  vector<string> denoms = {"Capybara","HT-binned","Capy HT-binned","Capybara"};
  
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

  string baseline = "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6&&nbm_moriond>=1&&met/met_calo<5.";
  vector<string> weights = {"weight","weight/(w_btag)"};
  for(unsigned int i=0; i<proc_pairs.size();i++){
    for(unsigned int iw=0; iw<weights.size();iw++){
      pm.Push<Hist1D>(Axis(20, 0, 2000., "st", "S_{T} [GeV]", {500.}),
		      "nleps==1&&nveto==0&&met>200&&njets>=6&&nbm_moriond>=1&&met/met_calo<5.",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);
      pm.Push<Hist1D>(Axis(10, 0, 2000., "st", "S_{T} [GeV]", {500.}),
		      "1",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);

      pm.Push<Hist1D>(Axis(33, 30, 1020., "jets_pt", "pT_{jets} [GeV]", {-999}),
		      baseline+"&&!jets_islep&&jets_pt>30",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);
      pm.Push<Hist1D>(Axis(33, 30, 1020., "jets_pt", "pT_{jets} [GeV]", {-999}),
		      "!jets_islep&&jets_pt>30",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);

      pm.Push<Hist1D>(Axis(20, 0, 1000., "leps_pt", "pT_{lep} [GeV]", {-999.}),
		      baseline,
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);
      pm.Push<Hist1D>(Axis(20, 0, 1000., "leps_pt", "pT_{lep} [GeV]", {-999.}),
		      "1",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);

      pm.Push<Hist1D>(Axis(15, 0, 600., "fjets14_m", "m_{J} [GeV]", {-999}),
		      baseline,
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);
      pm.Push<Hist1D>(Axis(15, 0, 600., "fjets14_m", "m_{J} [GeV]", {-999}),
		      "1",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);

      pm.Push<Hist1D>(Axis(24, 0, 1200., "mj14", "M_{J} [GeV]", {250.,400}),
		      baseline,
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);
      pm.Push<Hist1D>(Axis(24, 0, 1200., "mj14", "M_{J} [GeV]", {250.,400}),
		      "1",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);

      pm.Push<Hist1D>(Axis(15, 100, 850., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
		      "nleps==1&&nveto==0&&st>500&&met>100&&njets>=6&&nbm_moriond>=1&&met/met_calo<5.",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);
      pm.Push<Hist1D>(Axis(15, 100, 850., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
		      "1",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);
	
      pm.Push<Hist1D>(Axis(28, 0, 700., "met_tru", "True E_{T}^{miss} [GeV]", {150.}),
		      "1",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);

      pm.Push<Hist1D>(Axis(28, 0, 700., "met_tru", "True E_{T}^{miss} [GeV]", {150.}),
		      baseline,
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);

      pm.Push<Hist1D>(Axis(30, 0, 1500., "ht_isr_me", "True ISR H_{T} [GeV]", {600.}),
		      "1",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);

      pm.Push<Hist1D>(Axis(30, 0, 1500., "ht_isr_me", "True ISR H_{T} [GeV]", {600.}),
		      "ntruleps>=2",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);

      pm.Push<Hist1D>(Axis(30, 0, 1500., "ht_isr_me", "True ISR H_{T} [GeV]", {600.}),
		      "ntruleps==1",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);
 

      pm.Push<Hist1D>(Axis(20, 0, 700., "mt", "m_{T} [GeV]", {-999}),
		      baseline,
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);

      pm.Push<Hist1D>(Axis(20, 0, 700., "mt", "m_{T} [GeV]", {-999}),
		      "1",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);

      pm.Push<Hist1D>(Axis(12, -0.5, 11.5, "njets", "N_{jets}", {5.5, 8.5}),
		      "nleps==1&&nveto==0&&st>500&&met>200&&nbm_moriond>=1&&met/met_calo<5.",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);
      pm.Push<Hist1D>(Axis(12, -0.5, 11.5, "njets", "N_{jets}", {5.5, 8.5}),
		      "1",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);


      pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm_moriond", "N_{b}", {-999}),
		      "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6&&met/met_calo<5.",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);
      pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm_moriond", "N_{b}", {-999}),
		      "1",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);

      pm.Push<Hist1D>(Axis(20, 0.54, 1, "jets_csv", "CSV", {0.8484,0.9535}),
		      baseline+"&&jets_pt>30&&!jets_islep",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);
      pm.Push<Hist1D>(Axis(20, 0.54, 1, "jets_csv", "CSV", {0.8484,0.9535}),
		      "jets_pt>30&&!jets_islep",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);

      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nels", "N_{e}", {0.5}),
		      baseline,
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);
      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nels", "N_{e}", {0.5}),
		      "1",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);

      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nmus", "N_{mu}", {0.5}),
		      baseline,
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);
      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nmus", "N_{mu}", {0.5}),
		      "1",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);

      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nleps", "N_{leps}", {0.5}),
		      "st>500&&met>200&&njets>=6&&nbm_moriond>=1&&met/met_calo<5.",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);
      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nleps", "N_{leps}", {0.5}),
		      "1",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);

      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nveto", "N_{veto}", {0.5}),
		      "nleps==1&&st>500&&met>200&&njets>=6&&nbm_moriond>=1&&met/met_calo<5.",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);
      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nveto", "N_{veto}", {0.5}),
		      "1",
		      proc_pairs[i], all_plot_types).RatioTitle(numers[i],denoms[i]).Tag(proc_tags[i]).Weight(weights[iw]);
    }
  }
   
 
  
  pm.min_print_ = true;
  pm.MakePlots(lumi);

}
