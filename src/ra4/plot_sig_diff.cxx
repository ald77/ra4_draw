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

 

  string bear_T1tttt_dir = "/net/cms29/cms29r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/";
  string capy_T1tttt_dir = "/net/cms29/cms29r0/babymaker/babies/2016_08_10/T1tttt/unskimmed/";




  auto bear_NC = Process::MakeShared<Baby_full>("New T1tttt(1800,100)", Process::Type::signal, colors("t1tttt"),
    {bear_T1tttt_dir+"*T1tttt_mGluino-1800_mLSP-100_*.root"},
    "pass&&st<100000");

  auto capy_NC = Process::MakeShared<Baby_full>("Old T1tttt(1800,100)", Process::Type::signal, colors("t1tttt"),
    {capy_T1tttt_dir+"*T1tttt_mGluino-1800_mLSP-100_*.root"},
    "pass");
  capy_NC->SetLineStyle(2);


  auto bear_C = Process::MakeShared<Baby_full>("New T1tttt(1400,1000)", Process::Type::signal, kAzure+2,
    {bear_T1tttt_dir+"*T1tttt_mGluino-1400_mLSP-1000_*.root"},
    "pass&&st<100000");

  auto capy_C = Process::MakeShared<Baby_full>("Old T1tttt(1400,1000)", Process::Type::signal, kAzure+2,
    {capy_T1tttt_dir+"*T1tttt_mGluino-1400_mLSP-1000_*.root"},
    "pass");
  capy_C->SetLineStyle(2);


  vector<shared_ptr<Process> > procs = {bear_NC,capy_NC,bear_C,capy_C};
  
  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    //    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::lumi_shapes);
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
  vector<string> weights = {"weight"/*,"weight/(w_btag)"*/};
  for(unsigned iw=0; iw<weights.size(); iw++){
    pm.Push<Hist1D>(Axis(20, 0, 2000., "st", "S_{T} [GeV]", {500.}),
		    "nleps==1&&nveto==0&&met>200&&njets>=6&&nbm_moriond>=1&&met/met_calo<5.",
		   procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");
    pm.Push<Hist1D>(Axis(20, 0, 2000., "st", "S_{T} [GeV]", {500.}),
		    "1",
		   procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");

    pm.Push<Hist1D>(Axis(33, 30, 1020., "jets_pt", "pT_{jets} [GeV]", {-999}),
		    baseline+"&&!jets_islep&&jets_pt>30",
		    procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");
    pm.Push<Hist1D>(Axis(33, 30, 1020., "jets_pt", "pT_{jets} [GeV]", {-999}),
		    "!jets_islep&&jets_pt>30",
		 procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");
    
      pm.Push<Hist1D>(Axis(20, 0, 1000., "leps_pt", "pT_{lep} [GeV]", {-999.}),
		      baseline,
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");
      pm.Push<Hist1D>(Axis(20, 0, 1000., "leps_pt", "pT_{lep} [GeV]", {-999.}),
		      "1",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");

      pm.Push<Hist1D>(Axis(15, 0, 600., "fjets14_m", "m_{J} [GeV]", {-999}),
		      baseline,
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");
      pm.Push<Hist1D>(Axis(15, 0, 600., "fjets14_m", "m_{J} [GeV]", {-999}),
		      "1",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");

      pm.Push<Hist1D>(Axis(24, 0, 1200., "mj14", "M_{J} [GeV]", {250.,400}),
		      baseline,
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");
      pm.Push<Hist1D>(Axis(24, 0, 1200., "mj14", "M_{J} [GeV]", {250.,400}),
		      "1",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");

      pm.Push<Hist1D>(Axis(15, 100, 850., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
		      "nleps==1&&nveto==0&&st>500&&met>100&&njets>=6&&nbm_moriond>=1&&met/met_calo<5.",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");
      pm.Push<Hist1D>(Axis(15, 100, 850., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
		      "1",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");
	
      pm.Push<Hist1D>(Axis(28, 0, 700., "met_tru", "True E_{T}^{miss} [GeV]", {150.}),
		      "1",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");

      pm.Push<Hist1D>(Axis(28, 0, 700., "met_tru", "True E_{T}^{miss} [GeV]", {150.}),
		      baseline,
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");

      pm.Push<Hist1D>(Axis(30, 0, 1500., "ht_isr_me", "True ISR H_{T} [GeV]", {600.}),
		      "1",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");

      pm.Push<Hist1D>(Axis(30, 0, 1500., "ht_isr_me", "True ISR H_{T} [GeV]", {600.}),
		      "ntruleps>=2",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");

      pm.Push<Hist1D>(Axis(30, 0, 1500., "ht_isr_me", "True ISR H_{T} [GeV]", {600.}),
		      "ntruleps==1",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");
 

      pm.Push<Hist1D>(Axis(20, 0, 700., "mt", "m_{T} [GeV]", {-999}),
		      baseline,
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");

      pm.Push<Hist1D>(Axis(20, 0, 700., "mt", "m_{T} [GeV]", {-999}),
		      "1",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");

      pm.Push<Hist1D>(Axis(12, -0.5, 11.5, "njets", "N_{jets}", {5.5, 8.5}),
		      "nleps==1&&nveto==0&&st>500&&met>200&&nbm_moriond>=1&&met/met_calo<5.",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");
      pm.Push<Hist1D>(Axis(12, -0.5, 11.5, "njets", "N_{jets}", {5.5, 8.5}),
		      "1",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");


      pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm_moriond", "N_{b}", {-999}),
		      "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6&&met/met_calo<5.",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");
      pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm_moriond", "N_{b}", {-999}),
		      "1",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");

      pm.Push<Hist1D>(Axis(20, 0.54, 1, "jets_csv", "CSV", {0.8484,0.9535}),
		      baseline+"&&jets_pt>30&&!jets_islep",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");
      pm.Push<Hist1D>(Axis(20, 0.54, 1, "jets_csv", "CSV", {0.8484,0.9535}),
		      "jets_pt>30&&!jets_islep",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");

      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nels", "N_{e}", {0.5}),
		      baseline,
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");
      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nels", "N_{e}", {0.5}),
		      "1",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");

      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nmus", "N_{mu}", {0.5}),
		      baseline,
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");
      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nmus", "N_{mu}", {0.5}),
		      "1",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");

      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nleps", "N_{leps}", {0.5}),
		      "st>500&&met>200&&njets>=6&&nbm_moriond>=1&&met/met_calo<5.",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");
      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nleps", "N_{leps}", {0.5}),
		      "1",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");

      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nveto", "N_{veto}", {0.5}),
		      "nleps==1&&st>500&&met>200&&njets>=6&&nbm_moriond>=1&&met/met_calo<5.",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");
      pm.Push<Hist1D>(Axis(5, -0.5, 4.5, "nveto", "N_{veto}", {0.5}),
		      "1",
		     procs, all_plot_types).Weight(weights[iw]).Tag("sig_diff_v2");
    
  }
   
 
  
  pm.min_print_ = true;
  pm.MakePlots(lumi);

}
