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

  double lumi = 2.3;

  string trig_skim_mc = "/net/cms27/cms27r0/babymaker/2016_04_29/mc/merged_1lht500met200/";

  auto tt1l = Proc<Baby_full>("t#bar{t} (1l)", Process::Type::background, TColor::GetColor(86, 160, 211),
    {trig_skim_mc+"*_TTJets_*.root"}, "nleps<=1");
  auto tt2l = Proc<Baby_full>("t#bar{t} (2l)", Process::Type::background, TColor::GetColor(1, 57, 166),
    {trig_skim_mc+"*_TTJets_*.root"}, "nleps>=2");
  auto wjets = Proc<Baby_full>("W+jets", Process::Type::background, TColor::GetColor(0, 79, 39),
    {trig_skim_mc+"*_WJetsTo*.root"});
  auto single_t = Proc<Baby_full>("Single t", Process::Type::background, TColor::GetColor(52, 42, 123),
    {trig_skim_mc+"*_ST_*.root"});
  auto ttv = Proc<Baby_full>("t#bar{t}V", Process::Type::background, TColor::GetColor(149, 0, 26),
    {trig_skim_mc+"*_TTWJetsTo*.root", trig_skim_mc+"*_TTZTo*.root"});
  auto other = Proc<Baby_full>("Other", Process::Type::background, TColor::GetColor(255, 200, 47),
    {trig_skim_mc+"*_DYJets*.root", trig_skim_mc+"*_QCD_HT*.root",
        trig_skim_mc+"*_TTGJets*.root", trig_skim_mc+"*_TTTT_*.root",
        trig_skim_mc+"*_WH_HToBB*.root", trig_skim_mc+"*_WWTo*.root",
        trig_skim_mc+"*_WZTo*.root", trig_skim_mc+"*_ZH_HToBB*.root",
        trig_skim_mc+"*_ZZ_*.root", trig_skim_mc+"*_ttHJetTobb*.root"});

  auto t1tttt_nc = Proc<Baby_full>("T1tttt(1500,100)", Process::Type::signal, kRed,
    {trig_skim_mc+"*SMS-T1tttt_mGluino-1500_mLSP-100*.root"});
  auto t1tttt_c = Proc<Baby_full>("T1tttt(1200,800)", Process::Type::signal, kRed,
    {trig_skim_mc+"*SMS-T1tttt_mGluino-1200_mLSP-800*.root"});
  t1tttt_c->SetLineStyle(2);

  auto data = Proc<Baby_full>("Data", Process::Type::data, kBlack,
    {"/net/cms27/cms27r0/babymaker/2016_04_29/data/merged_1lht500met200/*.root"});

  vector<shared_ptr<Process> > full_trig_skim = {other, ttv, single_t, wjets, tt2l, tt1l, t1tttt_c, t1tttt_nc, data};

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::signal_overlay);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .Bottom(BottomType::off)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info);
  vector<PlotOpt> all_plot_types = {log_lumi, lin_lumi, log_shapes, lin_shapes,
                                    log_lumi_info, lin_lumi_info, log_shapes_info, lin_shapes_info};

  PlotMaker pm;
  pm.AddPlot(HistoDef(7, -0.5, 6.5, "nleps", "Num. Leptons",
                      "ht>500&&met>200&&njets>=6&&nbm>=1", "weight", {0.5, 1.5}),
             full_trig_skim, all_plot_types);
  pm.AddPlot(HistoDef(40, 0, 2000., "ht", "H_{T} [GeV]",
                      "nleps==1&&met>200&&njets>=6&&nbm>=1", "weight", {500.}),
             full_trig_skim, all_plot_types);
  pm.AddPlot(HistoDef(40, 0, 1000., "met", "MET [GeV]",
                      "nleps==1&&ht>500&&njets>=6&&nbm>=1", "weight", {200., 400.}),
             full_trig_skim, all_plot_types);
  pm.AddPlot(HistoDef(16, -0.5, 15.5, "njets", "Num. AK4 Jets",
                      "nleps==1&&ht>500&&met>200&&nbm>=1", "weight", {5.5, 8.5}),
             full_trig_skim, all_plot_types);
  pm.AddPlot(HistoDef(11, -0.5, 10.5, "nbm", "Num. b-Tagged Jets",
                      "nleps==1&&ht>500&&met>200&&njets>=6", "weight", {0.5, 1.5, 2.5}),
             full_trig_skim, all_plot_types);
  pm.AddPlot(HistoDef(24, 0., 1200., "mj", "M_{J} [GeV]",
                      "nleps==1&&ht>500&&met>200&&njets>=6&&nbm>=1", "weight", {250., 400.}),
             full_trig_skim, all_plot_types);
  pm.AddPlot(HistoDef(12, 0., 420., "mt", "m_{T} [GeV]",
                      "nleps==1&&ht>500&&met>200&&njets>=6&&nbm>=1", "weight", {140.}),
             full_trig_skim, all_plot_types);
  pm.MakePlots(lumi);
}
