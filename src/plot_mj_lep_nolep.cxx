#include <cmath>

#include "TError.h"
#include "TVector2.h"

#include "plot_maker.hpp"
#include "plot_opt.hpp"
#include "palette.hpp"
#include "histo_stack.hpp"
#include "event_scan.hpp"

using namespace std;
using namespace PlotOptTypes;

template<typename T>
shared_ptr<Process> Proc(const string process_name, Process::Type type,
                         int color, const set<string> &files, const NamedFunc &cut = true){
  return make_shared<Process>(process_name, type, color,
                              unique_ptr<Baby>(new T(files)),
                              cut);
}

int main(){
  gErrorIgnoreLevel = 6000;

  double lumi = 2.6;

  Palette colors("txt/colors.txt", "default");

  string folder_mc = "/net/cms26/cms26r0/babymaker/babies/mismeasured_v2/2016_06_14/mc/merged_mm_std_nj5mj250/";
  set<string> file_names = {folder_mc+"*_TTJets*Lept*.root", folder_mc+"*_TTJets_HT*.root",
			    folder_mc+"*_WJetsToLNu*.root", folder_mc+"*_ST_*.root",
			    folder_mc+"*_TTWJets*.root", folder_mc+"*_TTZTo*.root",
			    folder_mc+"*DYJetsToLL*.root", folder_mc+"*_QCD_HT*.root",
			    folder_mc+"*ggZH_HToBB*.root", folder_mc+"*ttHJetTobb*.root",
			    folder_mc+"*_TTGJets*.root", folder_mc+"*_TTTT_*.root",
			    folder_mc+"*_WH_HToBB*.root", folder_mc+"*_WZTo*.root",
			    folder_mc+"*_ZH_HToBB*.root", folder_mc+"_ZZ_*.root"};
  file_names = {folder_mc+"*_TTJets*Lept*.root", folder_mc+"*_TTJets_HT*.root"};

  vector<int> scenarios = {0, 2};
  PlotMaker pm;
  for(const auto &scenario: scenarios){
    string s = "["+to_string(scenario)+"]";
    NamedFunc baseline("pass&&stitch&&mm_nleps"+s+">=1&&mm_met"+s+">200&&mm_ht"+s+">500");
    NamedFunc lowmt_cut("mm_nleps"+s+"==1&&mm_mt"+s+"<=140&&mm_njets"+s+">=6&&mm_njets"+s+"<9&&mm_nbm"+s+">=1");
    NamedFunc highmt_cut("mm_nleps"+s+"==1&&mm_mt"+s+">140&&mm_njets"+s+">=6&&mm_njets"+s+"<9&&mm_nbm"+s+">=1&&ntruleps==2");
    NamedFunc dilep_cut_lep("mm_nleps"+s+"==2&&mm_njets"+s+">=5&&mm_njets"+s+"<8&&mm_nbm"+s+">=0");
    NamedFunc dilep_cut_nolep("mm_nleps"+s+"==2&&mm_njets"+s+">=6&&mm_njets"+s+"<9&&mm_nbm"+s+">=0");

    auto lowmt = Proc<Baby_full>("m_{T}<=140. GeV", Process::Type::background, colors.RGB(0, 255, 0),
				      file_names, baseline&&lowmt_cut);
    auto highmt = Proc<Baby_full>("m_{T}>140 GeV, 2 True Leps", Process::Type::background, colors.RGB(0, 0, 255),
				       file_names, baseline&&highmt_cut);
    auto dilep_lep = Proc<Baby_full>("Dilepton", Process::Type::background, colors.RGB(255, 0, 0),
				     file_names, baseline&&dilep_cut_lep);//njets5
    auto dilep_nolep = Proc<Baby_full>("Dilepton", Process::Type::background, colors.RGB(255, 0, 0),
				       file_names, baseline&&dilep_cut_nolep);//njets6

    vector<shared_ptr<Process> > procs_lep = {lowmt, highmt, dilep_lep};
    vector<shared_ptr<Process> > procs_nolep = {lowmt, highmt, dilep_nolep};

    string tag = "rohanplot_lownj";

    PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
    log_lumi.Title(TitleType::info)
      .Bottom(BottomType::ratio)
      .YAxis(YAxisType::log)
      .Stack(StackType::shapes);
    PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
    vector<PlotOpt> plot_types = {log_lumi, lin_lumi};

    NamedFunc mjlep("mm_mj14_lep"+s+">250.");
    NamedFunc mjnolep("mm_mj14_nolep"+s+">250.");
    NamedFunc low_met("mm_met"+s+"<=350");
    NamedFunc med_met("mm_met"+s+">350&&mm_met"+s+"<=500");
    NamedFunc high_met("mm_met"+s+">500");

    pm.Push<HistoStack>(HistoDef(tag, 30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]",
				 mjlep, "weight", {250., 400.}), procs_lep, plot_types);
    pm.Push<HistoStack>(HistoDef(tag, 30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]",
				 mjnolep, "weight", {250., 400.}), procs_nolep, plot_types);
    pm.Push<HistoStack>(HistoDef(tag, 30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]",
				 mjlep&&low_met, "weight", {250., 400.}), procs_lep, plot_types);
    pm.Push<HistoStack>(HistoDef(tag, 30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]",
				 mjnolep&&low_met, "weight", {250., 400.}), procs_nolep, plot_types);
    pm.Push<HistoStack>(HistoDef(tag, 30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]",
				 mjlep&&med_met, "weight", {250., 400.}), procs_lep, plot_types);
    pm.Push<HistoStack>(HistoDef(tag, 30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]",
				 mjnolep&&med_met, "weight", {250., 400.}), procs_nolep, plot_types);
    pm.Push<HistoStack>(HistoDef(tag, 30, 0., 1500., "mm_mj14_lep"+s, "M_{J} (with lep) [GeV]",
				 mjlep&&high_met, "weight", {250., 400.}), procs_lep, plot_types);
    pm.Push<HistoStack>(HistoDef(tag, 30, 0., 1500., "mm_mj14_nolep"+s, "M_{J} (no lep) [GeV]",
				 mjnolep&&high_met, "weight", {250., 400.}), procs_nolep, plot_types);
  }
  pm.MakePlots(lumi);
}
