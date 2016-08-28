#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "TError.h"
#include "TColor.h"

#include "core/baby.hpp"
#include "core/histo_stack.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

int main(){
  gErrorIgnoreLevel = 6000;

  double lumi = 2.6;

  string mc_amcatnlo_lo = "/homes/cawest/links/TTJets_TuneCUETP8M1_13TeV-madgraphMLM/";
  string mc_powheg = "/homes/cawest/links/TT_TuneCUETP8M1_13TeV-powheg-pythia8/";
  string mc_amcatnlo_nlo = "/homes/cawest/links/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/";

  Palette colors("txt/colors.txt", "default");

  //tt amcanlo lo
  auto tt_amcatnlo_lo = Process::MakeShared<Baby_full>("t#bar{t}, aMC@NLO (LO, MLM)", Process::Type::background, kBlack,
    {mc_amcatnlo_lo+"*TTJets_TuneCUETP8M1*madgraphMLM*.root"});
  //tt Powheg
  auto tt_powheg = Process::MakeShared<Baby_full>("t#bar{t}, Powheg", Process::Type::background, kRed,
    {mc_powheg+"*_TT_*powheg*.root"});
  tt_powheg->SetLineStyle(2);
  //tt amcanlo nlo
  auto tt_amcatnlo_nlo = Process::MakeShared<Baby_full>("t#bar{t}, aMC@NLO (NLO, FxFx)", Process::Type::background, kBlue,
    {mc_amcatnlo_nlo+"*TTJets_TuneCUETP8M1*amcatnlo*.root"});
  tt_amcatnlo_nlo->SetLineStyle(3);

  vector<shared_ptr<Process> > tt_sams = {tt_amcatnlo_lo, tt_powheg, tt_amcatnlo_nlo};

  //
  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm)
    .FileExtensions({"pdf"});

  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .Bottom(BottomType::ratio)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info).Bottom(BottomType::ratio).Stack(StackType::data_norm);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::ratio).Stack(StackType::data_norm);
  vector<PlotOpt> plot_types = {log_shapes};

  PlotMaker pm;

  // 0-lepton selection
  pm.Push<HistoStack>(HistoDef(6, 0, 6, "nbm", "N_{b}",
                               "nleps==0","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);
  pm.Push<HistoStack>(HistoDef(25, 0, 2500, "mj", "M_{J} [GeV]",
                               "nleps==0","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);
  pm.Push<HistoStack>(HistoDef(20, 0, 20, "njets", "N_{jets}",
                               "nleps==0","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);
  pm.Push<HistoStack>(HistoDef(40, 0, 4000, "ht", "H_{T} [GeV]",
                               "nleps==0","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);

  pm.Push<HistoStack>(HistoDef(6, 0, 6, "nbm", "N_{b}",
                               "nleps==0&&ht>1500&&njets>=4&&mj>500&&nbm>=1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);
  pm.Push<HistoStack>(HistoDef(25, 0, 2500, "mj", "M_{J} [GeV]",
                               "nleps==0&&ht>1500&&njets>=4&&mj>500&&nbm>=1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);
  pm.Push<HistoStack>(HistoDef(20, 0, 20, "njets", "N_{jets}",
                               "nleps==0&&ht>1500&&njets>=4&&mj>500&&nbm>=1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);
  pm.Push<HistoStack>(HistoDef(40, 0, 4000, "ht", "H_{T} [GeV]",
                               "nleps==0&&ht>1500&&njets>=4&&mj>500&&nbm>=1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);

  pm.Push<HistoStack>(HistoDef(6, 0, 6, "nbm", "N_{b}",
                               "nleps==0&&ht>1500&&njets>=8&&mj>500&&nbm>=1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);
  pm.Push<HistoStack>(HistoDef(25, 0, 2500, "mj", "M_{J} [GeV]",
                               "nleps==0&&ht>1500&&njets>=8&&mj>500&&nbm>=1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);
  pm.Push<HistoStack>(HistoDef(20, 0, 20, "njets", "N_{jets}",
                               "nleps==0&&ht>1500&&njets>=8&&mj>500&&nbm>=1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);
  pm.Push<HistoStack>(HistoDef(40, 0, 4000, "ht", "H_{T} [GeV]",
                               "nleps==0&&ht>1500&&njets>=8&&mj>500&&nbm>=1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);

  // 1-lepton selectio
  pm.Push<HistoStack>(HistoDef(6, 0, 6, "nbm", "N_{b}",
                               "nleps==1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);
  pm.Push<HistoStack>(HistoDef(25, 0, 2500, "mj", "M_{J} [GeV]",
                               "nleps==1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);
  pm.Push<HistoStack>(HistoDef(20, 0, 20, "njets", "N_{jets}",
                               "nleps==1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);
  pm.Push<HistoStack>(HistoDef(40, 0, 4000, "ht", "H_{T} [GeV]",
                               "nleps==1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);

  pm.Push<HistoStack>(HistoDef(6, 0, 6, "nbm", "N_{b}",
                               "nleps==1&&ht>1200&&njets>=4&&mj>500&&nbm>=1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);
  pm.Push<HistoStack>(HistoDef(25, 0, 2500, "mj", "M_{J} [GeV]",
                               "nleps==1&&ht>1200&&njets>=4&&mj>500&&nbm>=1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);
  pm.Push<HistoStack>(HistoDef(20, 0, 20, "njets", "N_{jets}",
                               "nleps==1&&ht>1200&&njets>=4&&mj>500&&nbm>=1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);
  pm.Push<HistoStack>(HistoDef(40, 0, 4000, "ht", "H_{T} [GeV]",
                               "nleps==1&&ht>1200&&njets>=4&&mj>500&&nbm>=1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);

  pm.Push<HistoStack>(HistoDef(6, 0, 6, "nbm", "N_{b}",
                               "nleps==1&&ht>1200&&njets>=8&&mj>500&&nbm>=1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);
  pm.Push<HistoStack>(HistoDef(25, 0, 2500, "mj", "M_{J} [GeV]",
                               "nleps==1&&ht>1200&&njets>=8&&mj>500&&nbm>=1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);
  pm.Push<HistoStack>(HistoDef(20, 0, 20, "njets", "N_{jets}",
                               "nleps==1&&ht>1200&&njets>=8&&mj>500&&nbm>=1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);
  pm.Push<HistoStack>(HistoDef(40, 0, 4000, "ht", "H_{T} [GeV]",
                               "nleps==1&&ht>1200&&njets>=8&&mj>500&&nbm>=1","weight*w_pu_rpv/eff_trig", {}), tt_sams, plot_types);

  pm.MakePlots(lumi);
}
