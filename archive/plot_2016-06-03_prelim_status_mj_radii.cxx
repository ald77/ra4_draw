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

  double lumi = 10.0;

  string trig_skim_mc = "/net/cms27/cms27r0/babymaker/2016_04_29/mc/merged_baseline/";
  string trig_skim_data = "/net/cms27/cms27r0/babymaker/2016_04_29/mc/T1tttt/skim_baseline/";

  Palette colors("txt/colors.txt", "default");

  auto tt = Proc<Baby_full>("t#bar{t}", Process::Type::background, colors("tt_1l"),
    {trig_skim_mc+"*_TTJets*Lept*.root", trig_skim_mc+"*_TTJets_HT*.root"},
    "stitch");

  auto t1tttt_nc = Proc<Baby_full>("T1tttt(1800,200)", Process::Type::signal, colors("t1tttt"),
    {trig_skim_data+"*SMS-T1tttt_mGluino-1800_mLSP-200_*.root"});
  auto t1tttt_c = Proc<Baby_full>("T1tttt(1400,1000)", Process::Type::signal, colors("t1tttt"),
    {trig_skim_data+"*SMS-T1tttt_mGluino-1400_mLSP-1000_*.root"});
  t1tttt_c->SetLineStyle(2);

  vector<shared_ptr<Process> > full_trig_skim = {t1tttt_nc, t1tttt_c, tt};

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm)
    .FileExtensions({"pdf","root"});
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
  pm.AddPlot(HistoDef(30, 0., 1500., "mj", "M_{J}^{R=1.2} [GeV]",
                      "nleps==1&&ht>500&&met>200&&njets>=6&&nbm>=1", "weight", {250., 400.}),
             full_trig_skim, all_plot_types);
  pm.AddPlot(HistoDef(30, 0., 1500., "mj14", "M_{J}^{R=1.4} [GeV]",
                      "nleps==1&&ht>500&&met>200&&njets>=6&&nbm>=1", "weight", {250., 400.}),
             full_trig_skim, all_plot_types);

  pm.MakePlots(lumi);
}
