#include "test.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TColor.h"

#include "baby.hpp"
#include "process.hpp"
#include "named_func.hpp"
#include "plot_maker.hpp"
#include "plot_opt.hpp"
#include "palette.hpp"
#include "table.hpp"
#include "histo_stack.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  bool single_thread = false;
}

template<typename T>
shared_ptr<Process> Proc(const string process_name, Process::Type type,
                         int color, const set<string> &files, const string &cut = "1"){
  return make_shared<Process>(process_name, type, color,
                              unique_ptr<Baby>(new T(files)),
                              cut);
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  double lumi = 2.3;

  string trig_skim_mc = "/net/cms27/cms27r0/babymaker/2016_04_29/mc/merged_1lht500met200/";

  Palette colors("txt/colors.txt", "default");

  auto tt1l = Proc<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {trig_skim_mc+"*_TTJets*Lept*.root", trig_skim_mc+"*_TTJets_HT*.root"},
    "ntruleps<=1&&stitch");
  auto tt2l = Proc<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {trig_skim_mc+"*_TTJets*Lept*.root", trig_skim_mc+"*_TTJets_HT*.root"},
    "ntruleps>=2&&stitch");
  auto wjets = Proc<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {trig_skim_mc+"*_WJetsToLNu*.root"});
  auto single_t = Proc<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {trig_skim_mc+"*_ST_*.root"});
  auto ttv = Proc<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {trig_skim_mc+"*_TTWJets*.root", trig_skim_mc+"*_TTZTo*.root"});
  auto other = Proc<Baby_full>("Other", Process::Type::background, colors("other"),
    {trig_skim_mc+"*DYJetsToLL*.root", trig_skim_mc+"*_QCD_HT*.root",
        trig_skim_mc+"*_ZJet*.root", trig_skim_mc+"*_WWTo*.root",
        trig_skim_mc+"*ggZH_HToBB*.root", trig_skim_mc+"*ttHJetTobb*.root",
        trig_skim_mc+"*_TTGJets*.root", trig_skim_mc+"*_TTTT_*.root",
        trig_skim_mc+"*_WH_HToBB*.root", trig_skim_mc+"*_WZTo*.root",
        trig_skim_mc+"*_ZH_HToBB*.root", trig_skim_mc+"_ZZ_*.root"});

  auto t1tttt_nc = Proc<Baby_full>("T1tttt(1500,100)", Process::Type::signal, colors("t1tttt"),
    {trig_skim_mc+"*SMS-T1tttt_mGluino-1500_mLSP-100*.root"});
  auto t1tttt_c = Proc<Baby_full>("T1tttt(1200,800)", Process::Type::signal, colors("t1tttt"),
    {trig_skim_mc+"*SMS-T1tttt_mGluino-1200_mLSP-800*.root"});
  t1tttt_c->SetLineStyle(2);

  auto data = Proc<Baby_full>("Data", Process::Type::data, kBlack,
    {"/net/cms27/cms27r0/babymaker/2016_04_29/data/merged_1lht500met200/*.root"},"pass&&(trig[4]||trig[8])");

  vector<shared_ptr<Process> > full_trig_skim = {data, t1tttt_nc, t1tttt_c, tt1l, tt2l, wjets, single_t, ttv, other};

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
  vector<PlotOpt> all_plot_types = {log_lumi, lin_lumi, log_shapes, lin_shapes,
                                    log_lumi_info, lin_lumi_info, log_shapes_info, lin_shapes_info};

  PlotMaker pm;
  pm.Push<HistoStack>(HistoDef(7, -0.5, 6.5, "nleps", "Num. Leptons",
                               "ht>500&&met>200&&njets>=6&&nbm>=1", "weight", {0.5, 1.5}),
                      full_trig_skim, all_plot_types);
  pm.Push<HistoStack>(HistoDef(40, 0, 2000., "ht", "H_{T} [GeV]",
                               "nleps==1&&met>200&&njets>=6&&nbm>=1", "weight", {500.}),
                      full_trig_skim, all_plot_types);
  pm.Push<HistoStack>(HistoDef(40, 0, 1000., "met", "MET [GeV]",
                               "nleps==1&&ht>500&&njets>=6&&nbm>=1", "weight", {200., 400.}),
                      full_trig_skim, all_plot_types);
  pm.Push<HistoStack>(HistoDef(16, -0.5, 15.5, "njets", "Num. AK4 Jets",
                               "nleps==1&&ht>500&&met>200&&nbm>=1", "weight", {5.5, 8.5}),
                      full_trig_skim, all_plot_types);
  pm.Push<HistoStack>(HistoDef(11, -0.5, 10.5, "nbm", "Num. b-Tagged Jets",
                               "nleps==1&&ht>500&&met>200&&njets>=6", "weight", {0.5, 1.5, 2.5}),
                      full_trig_skim, all_plot_types);
  pm.Push<HistoStack>(HistoDef(24, 0., 1200., "mj", "M_{J} [GeV]",
                               "nleps==1&&ht>500&&met>200&&njets>=6&&nbm>=1", "weight", {250., 400.}),
                      full_trig_skim, all_plot_types);
  pm.Push<HistoStack>(HistoDef(12, 0., 420., "mt", "m_{T} [GeV]",
                               "nleps==1&&ht>500&&met>200&&njets>=6&&nbm>=1", "weight", {140.}),
                      full_trig_skim, all_plot_types);
  pm.Push<HistoStack>(HistoDef(15, 0., 1500., "mj08", "M_{J}^{0.8} [GeV]",
                               "nleps==1&&ht>500&&met>200", "weight", {400.}),
                      full_trig_skim, all_plot_types);
  vector<TableRow> rows = {
    TableRow("Baseline"),
    TableRow("No Selection", "1"),
    TableRow("$1\\ell$, $H_{T}>500$, $\\slashed{E}_{T}>200$", "nleps==1&&ht>500&&met>200"),
    TableRow("$N_{\\text{jets}}\\geq6$", "nleps==1&&ht>500&&met>200&&njets>=6"),
    TableRow("$N_{b}\\geq1$", "nleps==1&&ht>500&&met>200&&njets>=6&&nbm>=1"),
    TableRow("$M_{J}>250$", "nleps==1&&ht>500&&met>200&&njets>=6&&nbm>=1&&mj14>250", 0, 1),
    TableRow("ABCD Signal Region"),
    TableRow("$m_{T}>140$", "nleps==1&&ht>500&&met>200&&njets>=6&&nbm>=1&&mj14>250&&mt>140"),
    TableRow("$M_{J}>400$", "nleps==1&&ht>500&&met>200&&njets>=6&&nbm>=1&&mj14>400&&mt>140"),
    TableRow("Binning"),
    TableRow("$\\slashed{E}_{T}>500$", "nleps==1&&ht>500&&met>500&&njets>=6&&nbm>=1&&mj14>400&&mt>140"),
    TableRow("$N_{\\text{jets}}\\geq9$", "nleps==1&&ht>500&&met>500&&njets>=9&&nbm>=1&&mj14>400&&mt>140"),
    TableRow("$N_{b}\\geq3$", "nleps==1&&ht>500&&met>500&&njets>=9&&nbm>=3&&mj14>400&&mt>140"),
  };
  pm.Push<Table>("cutflow", rows, full_trig_skim);

  if(single_thread) pm.multithreaded_ = false;
  pm.MakePlots(lumi);
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 's'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s", long_options, &option_index);
    
    if( opt == -1) break;
    
    string optname;
    switch(opt){
    case 's':
      single_thread = true;
      break;
    case '0':
      optname = long_options[option_index].name;
      if(false){
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
