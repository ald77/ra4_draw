#include <cmath>

#include "TError.h"
#include "TVector2.h"

#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist1d.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  NamedFunc isr_weight("isr_weight", [](const Baby &b) -> NamedFunc::ScalarType{
      if(b.ntrupv()<0) return b.weight();
      return b.SampleType() == 20
	? Functions::njets_weights_ttisr.GetScalar(b)
	: b.weight()/(b.eff_trig()*b.w_toppt());
    });

  NamedFunc full_weight("full_weight", [](const Baby &b) -> NamedFunc::ScalarType{
      if(b.ntrupv()<0) return b.weight();
      double orig = isr_weight.GetScalar(b);
      double mj = b.mj14();
      if(mj <= 300.) return orig * (1.3 - 0.4*mj/300.);
      else return 0.9 * orig;
    });
}

int main(){
  gErrorIgnoreLevel = 6000;

  double lumi = 2.6;

  Palette colors("txt/colors.txt", "default");

  string folder_mc = "/net/cms29/cms29r0/babymaker/babies/mismeasured_v2/2016_06_14/mc/unskimmed/";
  set<string> file_names = {folder_mc+"*_TTJets*Lept*.root", folder_mc+"*_TTJets_HT*.root",
			    folder_mc+"*_WJetsToLNu*.root", folder_mc+"*_ST_*.root",
			    folder_mc+"*_TTWJets*.root", folder_mc+"*_TTZTo*.root",
			    folder_mc+"*DYJetsToLL*.root", folder_mc+"*_QCD_HT*.root",
			    folder_mc+"*ggZH_HToBB*.root", folder_mc+"*ttHJetTobb*.root",
			    folder_mc+"*_TTGJets*.root", folder_mc+"*_TTTT_*.root",
			    folder_mc+"*_WH_HToBB*.root", folder_mc+"*_WZTo*.root",
			    folder_mc+"*_ZH_HToBB*.root", folder_mc+"_ZZ_*.root"};
  file_names = {folder_mc+"*_TTJets*Lept*.root", folder_mc+"*_TTJets_HT*.root"};

  vector<int> scenarios = {0,2,5};
  PlotMaker pm;
  vector<string> tags = {"lownj", "highnj", "allnj"};
  vector<float> lownj = {5.5, 8.5, 5.5};
  vector<float> highnj = {8.5, 9999., 9999.};
  for(size_t i = 0; i < tags.size(); ++i){
    for(const auto &scenario: scenarios){
      string s = "["+to_string(scenario)+"]";
      NamedFunc baseline("pass&&stitch&&mm_nleps"+s+">=1&&mm_met"+s+">200&&mm_ht"+s+">500");
      NamedFunc lowmt_cut("mm_nleps"+s+"==1&&mm_mt"+s+"<=140&&mm_njets"+s+">="+to_string(lownj.at(i))+"&&mm_njets"+s+"<"+to_string(highnj.at(i))+"&&mm_nbm"+s+">=1");
      NamedFunc highmt_cut("mm_nleps"+s+"==1&&mm_mt"+s+">140&&mm_njets"+s+">="+to_string(lownj.at(i))+"&&mm_njets"+s+"<"+to_string(highnj.at(i))+"&&mm_nbm"+s+">=1&&ntruleps==2");
      NamedFunc dilep_cut_lep("mm_nleps"+s+"==2&&mm_njets"+s+">="+to_string(lownj.at(i)-1)+"&&mm_njets"+s+"<"+to_string(highnj.at(i)-1)+"&&mm_nbm"+s+">=0");
      NamedFunc dilep_cut_nolep("mm_nleps"+s+"==2&&mm_njets"+s+">="+to_string(lownj.at(i))+"&&mm_njets"+s+"<"+to_string(highnj.at(i))+"&&mm_nbm"+s+">=0");

      auto lowmt = Process::MakeShared<Baby_full>("m_{T}<=140. GeV", Process::Type::background, colors.RGB(0, 255, 0),
						  file_names, baseline&&lowmt_cut);
      auto highmt = Process::MakeShared<Baby_full>("m_{T}>140 GeV, 2 True Leps", Process::Type::background, colors.RGB(0, 0, 255),
						   file_names, baseline&&highmt_cut);
      auto dilep_lep = Process::MakeShared<Baby_full>("Dilepton", Process::Type::background, colors.RGB(255, 0, 0),
						      file_names, baseline&&dilep_cut_lep);//njets5
      auto dilep_nolep = Process::MakeShared<Baby_full>("Dilepton", Process::Type::background, colors.RGB(255, 0, 0),
							file_names, baseline&&dilep_cut_nolep);//njets6

      vector<shared_ptr<Process> > procs_lep = {lowmt, highmt, dilep_lep};
      vector<shared_ptr<Process> > procs_nolep = {lowmt, highmt, dilep_nolep};

      string tag = "rohanplot_"+tags.at(i);

      PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
      log_lumi.Title(TitleType::info)
	.Bottom(BottomType::ratio)
	.YAxis(YAxisType::log)
	.Stack(StackType::shapes);
      PlotOpt noinfo = log_lumi().Title(TitleType::preliminary);
      vector<PlotOpt> plot_types = {log_lumi, noinfo};

      NamedFunc mjlep("mm_mj14_lep"+s+">250.");
      NamedFunc mjnolep("mm_mj14_nolep"+s+">250.");
      NamedFunc low_met("mm_met"+s+"<=350");
      NamedFunc med_met("mm_met"+s+">350&&mm_met"+s+"<=500");
      NamedFunc high_met("mm_met"+s+">500");

      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]", {250., 400.}),
                      true, procs_lep, plot_types).Tag(tag);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]", {250., 400.}),
                      true, procs_nolep, plot_types).Tag(tag);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]", {250., 400.}),
                      low_met, procs_lep, plot_types).Tag(tag);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]", {250., 400.}),
                      low_met, procs_nolep, plot_types).Tag(tag);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]", {250., 400.}),
                      med_met, procs_lep, plot_types).Tag(tag);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]", {250., 400.}),
                      med_met, procs_nolep, plot_types).Tag(tag);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]", {250., 400.}),
                      high_met, procs_lep, plot_types).Tag(tag);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]", {250., 400.}),
                      high_met, procs_nolep, plot_types).Tag(tag);

      NamedFunc w_no_toppt("weight/(eff_trig*w_toppt)");
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]", {250., 400.}),
                      true, procs_lep, plot_types).Tag(tag).Weight(w_no_toppt);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]", {250., 400.}),
                      true, procs_nolep, plot_types).Tag(tag).Weight(w_no_toppt);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]", {250., 400.}),
                      low_met, procs_lep, plot_types).Tag(tag).Weight(w_no_toppt);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]", {250., 400.}),
                      low_met, procs_nolep, plot_types).Tag(tag).Weight(w_no_toppt);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]", {250., 400.}),
                      med_met, procs_lep, plot_types).Tag(tag).Weight(w_no_toppt);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]", {250., 400.}),
                      med_met, procs_nolep, plot_types).Tag(tag).Weight(w_no_toppt);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]", {250., 400.}),
                      high_met, procs_lep, plot_types).Tag(tag).Weight(w_no_toppt);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]", {250., 400.}),
                      high_met, procs_nolep, plot_types).Tag(tag).Weight(w_no_toppt);

      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]", {250., 400.}),
                      true, procs_lep, plot_types).Tag(tag).Weight(isr_weight);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]", {250., 400.}),
                      true, procs_nolep, plot_types).Tag(tag).Weight(isr_weight);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]", {250., 400.}),
                      low_met, procs_lep, plot_types).Tag(tag).Weight(isr_weight);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]", {250., 400.}),
                      low_met, procs_nolep, plot_types).Tag(tag).Weight(isr_weight);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]", {250., 400.}),
                      med_met, procs_lep, plot_types).Tag(tag).Weight(isr_weight);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]", {250., 400.}),
                      med_met, procs_nolep, plot_types).Tag(tag).Weight(isr_weight);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]", {250., 400.}),
                      high_met, procs_lep, plot_types).Tag(tag).Weight(isr_weight);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]", {250., 400.}),
                      high_met, procs_nolep, plot_types).Tag(tag).Weight(isr_weight);

      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]", {250., 400.}),
                      true, procs_lep, plot_types).Tag(tag).Weight(full_weight);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]", {250., 400.}),
                      true, procs_nolep, plot_types).Tag(tag).Weight(full_weight);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]", {250., 400.}),
                      low_met, procs_lep, plot_types).Tag(tag).Weight(full_weight);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]", {250., 400.}),
                      low_met, procs_nolep, plot_types).Tag(tag).Weight(full_weight);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]", {250., 400.}),
                      med_met, procs_lep, plot_types).Tag(tag).Weight(full_weight);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]", {250., 400.}),
                      med_met, procs_nolep, plot_types).Tag(tag).Weight(full_weight);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]", {250., 400.}),
                      high_met, procs_lep, plot_types).Tag(tag).Weight(full_weight);
      pm.Push<Hist1D>(Axis(30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]", {250., 400.}),
                      high_met, procs_nolep, plot_types).Tag(tag).Weight(full_weight);

      NamedFunc mc_cut = "mm_mj14_nolep"+s+">250.&&mm_mj14_lep"+s+"<=400&&mm_met"+s+"<=500&&mm_nbm"+s+"<=2&&mm_nleps"+s+"==1&&mm_ht"+s+">500&&mm_met"+s+">200&&mm_njets"+s+">=6&&mm_nbm"+s+">=1";
      auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
	{folder_mc+"*_TTJets*Lept*.root", folder_mc+"*_TTJets_HT*.root"}, mc_cut&&"ntruleps<=1&&stitch");
      auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
	{folder_mc+"*_TTJets*Lept*.root", folder_mc+"*_TTJets_HT*.root"}, mc_cut&&"ntruleps>=2&&stitch");
      auto wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
	{folder_mc+"*_WJetsToLNu*.root"}, mc_cut);
      auto single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
	{folder_mc+"*_ST_*.root"}, mc_cut);
      auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
	{folder_mc+"*_TTWJets*.root", folder_mc+"*_TTZTo*.root"}, mc_cut);
      auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
	{folder_mc+"*DYJetsToLL*.root", folder_mc+"*_QCD_HT*.root",
	    folder_mc+"*_ZJet*.root", folder_mc+"*_WWTo*.root",
	    folder_mc+"*ggZH_HToBB*.root", folder_mc+"*ttHJetTobb*.root",
	    folder_mc+"*_TTGJets*.root", folder_mc+"*_TTTT_*.root",
	    folder_mc+"*_WH_HToBB*.root", folder_mc+"*_WZTo*.root",
	    folder_mc+"*_ZH_HToBB*.root", folder_mc+"_ZZ_*.root"}, mc_cut);
      auto data_2016 = Process::MakeShared<Baby_full>("2016 Data", Process::Type::data, kBlack,
	{"/net/cms2/cms2r0/babymaker/babies/2016_06_26/data/merged_standard/*.root"},
	"json2p6&&pass&&(trig[4]||trig[8]||trig[13]||trig[33])&&mj14>250.&&mj14<=400&&met<=500&&nbm<=2&&nleps==1&&ht>500&&met>200&&njets>=6&&nbm>=1");
      vector<shared_ptr<Process> > procs = {tt1l, tt2l, wjets, single_t, ttv, other, data_2016};
      NamedFunc mt("mt",[scenario](const Baby &b){
	  bool is_data = false;
	  for(const auto &name: b.FileNames()){
	    if(Contains(name, "data")){
	      is_data = true;
	      break;
	    }
	  }
	  if(is_data) return b.mt();
	  else return b.mm_mt()->at(scenario);
	});
      pm.Push<Hist1D>(Axis(30, 0., 1500., mt, "m_{T} [GeV]", {140.}), true, procs,
                      vector<PlotOpt>{log_lumi().Stack(StackType::data_norm),
                          noinfo().Stack(StackType::data_norm)}).Tag("mm"+to_string(scenario));
    }
  }
  pm.multithreaded_ = false;
  pm.MakePlots(lumi);
}
