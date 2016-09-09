#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "TError.h"
#include "TColor.h"
#include "TVector2.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/histo_stack.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace {
  double lumi = 12.9;
}

int main(){
  gErrorIgnoreLevel = 6000;
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  Palette colors("txt/colors.txt", "default");

  string ntupletag = "*.root";
  
  string fmc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_standard/";
  string fdata = bfolder+"/cms2r0/babymaker/babies/2016_08_10/data/merged_database_standard/";

  auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {fmc+"*_TTJets*Lept"+ntupletag, fmc+"*_TTJets_HT"+ntupletag}, "ntruleps<=1&&stitch");
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {fmc+"*_TTJets*Lept"+ntupletag, fmc+"*_TTJets_HT"+ntupletag}, "ntruleps>=2&&stitch");
  auto wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {fmc+"*_WJetsToLNu"+ntupletag},"stitch");
  auto single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {fmc+"*_ST_"+ntupletag});
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {fmc+"*_TTWJets"+ntupletag, fmc+"*_TTZTo"+ntupletag});
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {fmc+"*DYJetsToLL"+ntupletag, 
      fmc+"*_QCD_HT*00_Tune"+ntupletag, fmc+"*_QCD_HT*Inf_Tune"+ntupletag,
      fmc+"*_ZJet"+ntupletag, fmc+"*_WWTo"+ntupletag,
      fmc+"*ggZH_HToBB"+ntupletag, fmc+"*ttHJetTobb"+ntupletag,
      fmc+"*_TTGJets"+ntupletag, fmc+"*_TTTT_"+ntupletag,
      fmc+"*_WH_HToBB"+ntupletag, fmc+"*_WZTo"+ntupletag,
      fmc+"*_ZH_HToBB"+ntupletag, fmc+"_ZZ_"+ntupletag});

  auto data_highmt = Process::MakeShared<Baby_full>("Data 1l, m_{T} > 140", Process::Type::data, kBlack,
    {fdata+ntupletag},"pass && json12p9 && trig_ra4 && st>500 && mt>140 && nleps==1 && nveto==0 && njets>=6 && nbm>=1");
  auto data_lowmt = Process::MakeShared<Baby_full>("Data 1l, m_{T} #leq 140", Process::Type::background, kBlack,
    {fdata+ntupletag},"pass && json12p9 && trig_ra4 && st>500 && mt<=140 && nleps==1 && nveto==0 && njets>=6 && nbm>=1");
  data_lowmt->SetFillColor(kWhite);
  data_lowmt->SetLineColor(kBlue-7);
  data_lowmt->SetLineWidth(2);

  auto data2lveto = Process::MakeShared<Baby_full>("Data 2l or l+trk", Process::Type::data, kBlue+2,
    {fdata+ntupletag},
    "pass && json12p9 && trig_ra4 && st>500 && ((nleps==2 && njets>=5 && nbm<=2) || (nleps==1 && nveto==1 && njets>=6 && nbm>=1 && mt>140))");

  auto data2l = Process::MakeShared<Baby_full>("Data 2l", Process::Type::data, kMagenta+3,
    {fdata+ntupletag},
    "pass && json12p9 && trig_ra4 && st>500 && (nleps==2 && njets>=5 && nbm<=2)");


  auto data = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {fdata+ntupletag},"pass && json12p9 && trig_ra4");

  vector<shared_ptr<Process> > all_procs;
  all_procs = {data, tt1l, tt2l, wjets, single_t, ttv, other};

  vector<shared_ptr<Process> > data1l_procs = {data_highmt, data_lowmt};
  vector<shared_ptr<Process> > data2lveto_procs = {data2lveto, data_lowmt};
  vector<shared_ptr<Process> > data2l_procs = {data2l, data_lowmt};

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
  vector<PlotOpt> log = {log_lumi_info};
  vector<PlotOpt> lin = {lin_lumi_info};

  NamedFunc dphi_met_lep("dphi_met_lep",[](const Baby &b) -> NamedFunc::ScalarType{
    return fabs(TVector2::Phi_mpi_pi(b.met_phi()-b.leps_phi()->at(0)));
  });

  string baseline = "met>200 && nleps==1 && nveto==0 && st>500 && njets>=6 && nbm>=1";
  PlotMaker pm;

  //MJ modeling
  pm.Push<HistoStack>(HistoDef(20, 500., 2500., "st", "S_{T} [GeV]",
                               baseline, "weight"), all_procs, log);
  pm.Push<HistoStack>(HistoDef(20, 0., 1000., "mj14", "M_{J} [GeV]",
                               baseline +"&&mt<=140", "weight", {250.,400.}), all_procs, lin);
  pm.Push<HistoStack>(HistoDef(16, -0.5, 15.5, "njets", "N_{jets}",
                               CopyReplaceAll(baseline, "&& njets>=6", ""), "weight"), all_procs, lin);
  pm.Push<HistoStack>(HistoDef(20, 0., 700., "fjets14_m[0]", "leading fat jet mass [GeV]",
                               baseline, "weight", {5.5, 8.5}), all_procs, lin);

  //mT modeling
  pm.Push<HistoStack>(HistoDef(11, 150., 700., "met", "MET [GeV]",
                               CopyReplaceAll(baseline, "met>200", "met>150"), "weight", {200.,350.,500.}), all_procs, log);
  pm.Push<HistoStack>(HistoDef(20, 20., 520., "leps_pt[0]", "lepton p_{T} [GeV]",
                               baseline, "weight"), all_procs, log);
  pm.Push<HistoStack>(HistoDef(15, 0., 525., "mt", "m_{T} [GeV]",
                               baseline, "weight", {140.}), all_procs, log);
  pm.Push<HistoStack>(HistoDef(15, 2, 3.2, Functions::max_dphi_met_jet, "max #Delta#phi(j,MET)",
                               baseline, "weight"), all_procs, log);
  pm.Push<HistoStack>(HistoDef(15, 0., 1.5, Functions::min_dphi_met_jet, "min #Delta#phi(j,MET)",
                               baseline, "weight"), all_procs, log);
  pm.Push<HistoStack>(HistoDef(16, 0., 3.2, dphi_met_lep, "#Delta#phi(l,MET)",
                               baseline, "weight"), all_procs, log);

  //other
  pm.Push<HistoStack>(HistoDef(7, -0.5, 6.5, "nbm", "N_{b}",
                               CopyReplaceAll(baseline, "&& nbm>=1", ""), "weight", {0.5, 1.5, 2.5}), all_procs, lin);

  //data-to-data
  vector<string> metbins = {"met>150 && met<=500", "met>150 && met<=200", "met>200 && met<=350", "met>350 && met<=500", "met>200 && met<=500","met>500","met>200"};
  for (auto &imet: metbins){
    pm.Push<HistoStack>(HistoDef("data1l1b", 20,0.,1000., "mj14", "M_{J} [GeV]",
                                 imet + "&&nbm==1", "weight", {250.,400.}), data1l_procs, lin);
    pm.Push<HistoStack>(HistoDef("data1l2b", 20,0.,1000., "mj14", "M_{J} [GeV]",
                                 imet + "&&nbm>=2", "weight", {250.,400.}), data1l_procs, lin);
    pm.Push<HistoStack>(HistoDef("data2l", 20,0.,1000., "mj14", "M_{J} [GeV]",
                                 imet, "weight", {250.,400.}), data2l_procs, lin);
    pm.Push<HistoStack>(HistoDef("data2lveto", 20,0.,1000., "mj14", "M_{J} [GeV]",
                                 imet, "weight", {250.,400.}), data2lveto_procs, lin);
  } 

  pm.Push<HistoStack>(HistoDef(20, 0., 1000., "mj14", "M_{J} [GeV]",
                               baseline +"&&mt>140", "weight", {250.,400.}), all_procs, lin);
  pm.Push<HistoStack>(HistoDef(20, 0., 1000., "mj14", "M_{J} [GeV]",
                               baseline +"&&mt>140 && nbm==1", "weight", {250.,400.}), all_procs, lin);
  pm.Push<HistoStack>(HistoDef(20, 0., 1000., "mj14", "M_{J} [GeV]",
                               baseline +"&&mt>140 && nbm>=2", "weight", {250.,400.}), all_procs, lin);

  pm.min_print_ = true;
  pm.MakePlots(lumi);

}
