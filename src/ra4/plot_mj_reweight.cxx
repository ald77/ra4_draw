#include <string>
#include <memory>

#include "TError.h"

#include "core/baby_full.hpp"
#include "core/utilities.hpp"
#include "core/process.hpp"
#include "core/plot_maker.hpp"
#include "core/histo_stack.hpp"
#include "core/palette.hpp"
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
      else return orig*0.9;
    });
}

int main(){
  gErrorIgnoreLevel = 6000;

  Palette colors("txt/colors.txt", "default");

  string folder_mc="/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/skim_nleps2/";

  NamedFunc mll = "mumuv_m*(mumuv_m>0&&mumu_pt1>25)+elelv_m*(elelv_m>0&&elel_pt1>30)";
  NamedFunc baseline = mll>80. && mll<100. && "(nels>=2||nmus>=2)&&njets>=4&&ht>500&&pass";

  auto zjets = Process::MakeShared<Baby_full>("Z+jets", Process::Type::background, kGreen+1,
    {folder_mc+"*_DYJetsToLL*.root"},
    baseline);
  auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {folder_mc+"*_TTJets*Lept*.root", folder_mc+"*_TTJets_HT*.root"},
    baseline&&"ntruleps<=1&&stitch");
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {folder_mc+"*_TTJets*Lept*.root", folder_mc+"*_TTJets_HT*.root"},
    baseline&&"ntruleps>=2&&stitch");
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {folder_mc+"*_TTWJetsTo*.root", folder_mc+"*_TTZTo*.root"},
    baseline);
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {folder_mc+"*_WJetsToLNu_HT-*.root", folder_mc+"*_QCD_HT*.root",
        folder_mc+"*_TTGJets*.root", folder_mc+"*_TTTT*.root",
        folder_mc+"*_WH_HToBB*.root", folder_mc+"*_WWTo*.root",
        folder_mc+"*_WZTo*.root", folder_mc+"*_ZH*.root",
        folder_mc+"*_ZZ*.root", folder_mc+"*_ttHJetTobb*.root",
        folder_mc+"*_ST*channel*.root"},
    baseline);
  auto data = Process::MakeShared<Baby_full>("2016 Data", Process::Type::data, kBlack,
    {"/net/cms2/cms2r0/babymaker/babies/2016_06_26/data/skim_nleps2/*.root"},
    baseline&&"json2p6&&pass&&(trig[4]||trig[8]||trig[13]||trig[33])");
  vector<shared_ptr<Process> > procs = {data, zjets, tt1l, tt2l, ttv, other};
  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  vector<PlotOpt> plot_types = {log_lumi};

  PlotMaker pm;
  pm.Push<HistoStack>(HistoDef("orig", 24, 0., 600, "mj14", "M_{J} (with lep) [GeV]",
                               true, "weight", {250., 400.}), procs, plot_types);
  pm.Push<HistoStack>(HistoDef("orig", 24, 0., 600, "mj14", "M_{J} (with lep) [GeV]",
                               true, "weight/(eff_trig*w_toppt)", {250., 400.}), procs, plot_types);
  pm.Push<HistoStack>(HistoDef("orig", 24, 0., 600, "mj14", "M_{J} (with lep) [GeV]",
                               true, isr_weight, {250., 400.}), procs, plot_types);
  pm.Push<HistoStack>(HistoDef("orig", 24, 0., 600, "mj14", "M_{J} (with lep) [GeV]",
                               true, full_weight, {250., 400.}), procs, plot_types);
  pm.MakePlots(2.6, "reweight");
}
