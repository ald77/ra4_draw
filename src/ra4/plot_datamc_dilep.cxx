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
#include "core/utilities.hpp"
#include "core/histo_stack.hpp"

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

  vector<double> lumi = {2.3, 2.6}; // 2015 and 2016 luminosity

  // 80X (ttbar, qcd, dy[ht=100-600], wjets) + 74X (rest)

  string mc_2015 = "/net/cms2/cms2r0/babymaker/babies/2016_04_29/mc/merged_1lht500met200/";
  string data_2015 = "/net/cms2/cms2r0/babymaker/babies/2016_04_29/data/merged_1lht500met200/";
  string mc_2016 = "/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/merged_standard/";
  string data_2016 = "/net/cms2/cms2r0/babymaker/babies/2016_06_26/data/merged_standard/";

  Palette colors("txt/colors.txt", "default");

  //2015 Samples
  auto tt1l_2015 = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {mc_2015+"*_TTJets*Lept*.root", mc_2015+"*_TTJets_HT*.root"}, "ntruleps<=1&&stitch");
  auto tt2l_2015 = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {mc_2015+"*_TTJets*Lept*.root", mc_2015+"*_TTJets_HT*.root"}, "ntruleps>=2&&stitch");
  auto wjets_2015 = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {mc_2015+"*_WJetsToLNu*.root"});
  auto single_t_2015 = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {mc_2015+"*_ST_*.root"});
  auto ttv_2015 = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {mc_2015+"*_TTWJets*.root", mc_2015+"*_TTZTo*.root"});
  auto other_2015 = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {mc_2015+"*DYJetsToLL*.root", mc_2015+"*_QCD_HT*.root",
        mc_2015+"*_ZJet*.root", mc_2015+"*_WWTo*.root",
        mc_2015+"*ggZH_HToBB*.root", mc_2015+"*ttHJetTobb*.root",
        mc_2015+"*_TTGJets*.root", mc_2015+"*_TTTT_*.root",
        mc_2015+"*_WH_HToBB*.root", mc_2015+"*_WZTo*.root",
        mc_2015+"*_ZH_HToBB*.root", mc_2015+"_ZZ_*.root"});

  auto data2015 = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {data_2015+"*.root"},"pass&&(trig[4]||trig[8])");

  //2016 Samples
  auto tt1l_2016 = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {mc_2016+"*_TTJets*Lept*.root", mc_2016+"*_TTJets_HT*.root"}, "ntruleps<=1&&stitch");
  auto tt2l_2016 = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {mc_2016+"*_TTJets*Lept*.root", mc_2016+"*_TTJets_HT*.root"}, "ntruleps>=2&&stitch");
  auto wjets_2016 = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {mc_2016+"*_WJetsToLNu*.root"});
  auto single_t_2016 = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {mc_2016+"*_ST_*.root"});
  auto ttv_2016 = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {mc_2016+"*_TTWJets*.root", mc_2016+"*_TTZTo*.root"});
  auto other_2016 = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {mc_2016+"*DYJetsToLL*.root", mc_2016+"*_QCD_HT*.root",
        mc_2016+"*_ZJet*.root", mc_2016+"*_WWTo*.root",
        mc_2016+"*ggZH_HToBB*.root", mc_2016+"*ttHJetTobb*.root",
        mc_2016+"*_TTGJets*.root", mc_2016+"*_TTTT_*.root",
        mc_2016+"*_WH_HToBB*.root", mc_2016+"*_WZTo*.root",
        mc_2016+"*_ZH_HToBB*.root", mc_2016+"_ZZ_*.root"});

  auto data2016 = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {data_2016+"*.root"},"pass&&(trig[4]||trig[8]||trig[13]||trig[33])&&json2p6");

  vector<shared_ptr<Process> > sam2015 = {data2015, tt1l_2015, tt2l_2015, wjets_2015, single_t_2015, ttv_2015, other_2015};
  vector<shared_ptr<Process> > sam2016 = {data2016, tt1l_2016, tt2l_2016, wjets_2016, single_t_2016, ttv_2016, other_2016};

  //
  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm)
    .FileExtensions({"pdf"});

  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .Bottom(BottomType::off)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info);
  //  vector<PlotOpt> all_plot_types = {log_lumi_info, lin_lumi_info};
  vector<PlotOpt> all_plot_types = {log_lumi_info};
  vector<PlotOpt> mc_plot_types = {log_shapes_info, lin_shapes_info};

  vector<PlotMaker>  pm(2);
  vector<vector<shared_ptr<Process> > > sams = {sam2015, sam2016};  vector<string> sam_tag = {"sam2015","sam2016"};
  vector<string> leps = {"nleps==2","nels==2&&nmus==0","nels==0&&nmus==2","nels==1&&nmus==1","nleps==1","nels==1&&nmus==0","nels==0&&nmus==1"};
  vector<string> lep_tag = {"2leps","2els","2mus","1el1mu","1lep","1el","1mu"};

  // Loop over 2015 and 2016 samples
  for(unsigned int iyr=0; iyr<2; iyr++){
    // Loop over lepton flavors
    for(unsigned int ilep=0; ilep<leps.size(); ilep++){

      string cut2l_base = "ht>500&&met>200&&"+leps[ilep]+"&&njets>=5&&nbm<=2";
      if(ilep>3) cut2l_base = "ht>500&&met>200&&"+leps[ilep]+"&&njets>=6&&nbm>=1&&mt<=140";

      //MJ
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 1000., "mj14", "M_{J}^{1.4} [GeV]",
                                        cut2l_base+"&&met<=350","weight", {250, 400}),sams[iyr], all_plot_types);
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 1000., "mj14", "M_{J}^{1.4} [GeV]",
                                        cut2l_base+"&&met>350&&met<=500","weight", {250, 400}),sams[iyr], all_plot_types);

      //MJ low mt di-lepton
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 1000., "mj14", "M_{J}^{1.4} [GeV]",
                                        cut2l_base+"&&met<=350&&mt<=140","weight", {250, 400}),sams[iyr], all_plot_types);
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 1000., "mj14", "M_{J}^{1.4} [GeV]",
                                        cut2l_base+"&&met>350&&met<=500&&mt<=140","weight", {250, 400}),sams[iyr], all_plot_types);

      //MJ low mt di-lepton
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 1000., "mj14", "M_{J}^{1.4} [GeV]",
                                        cut2l_base+"&&met<=350&&mt<=140","weight", {250, 400}),sams[iyr], all_plot_types);
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 1000., "mj14", "M_{J}^{1.4} [GeV]",
                                        cut2l_base+"&&met>350&&met<=500&&mt<=140","weight", {250, 400}),sams[iyr], all_plot_types);

      //MJ high mt di-lepton
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 1000., "mj14", "M_{J}^{1.4} [GeV]",
                                        cut2l_base+"&&met<=350&&mt>140","weight", {250, 400}),sams[iyr], all_plot_types);
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 1000., "mj14", "M_{J}^{1.4} [GeV]",
                                        cut2l_base+"&&met>350&&met<=500&&mt>140","weight", {250, 400}),sams[iyr], all_plot_types);

      //D3 and D4
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 700., "mt", "m_{T} [GeV]",
                                        "mj14>250&&"+cut2l_base+"&&met<=350","weight", {140.}),sams[iyr], all_plot_types);
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 700., "mt", "m_{T} [GeV]",
                                        "mj14>250&&"+cut2l_base+"&&met>350&&met<=500","weight", {140.}),sams[iyr], all_plot_types);
      //D3
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 700., "mt", "m_{T} [GeV]",
                                        "mj14>250&&"+cut2l_base+"&&mj14<=400&&met<=350","weight", {140.}),sams[iyr], all_plot_types);
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 700., "mt", "m_{T} [GeV]",
                                        "mj14>250&&"+cut2l_base+"&&mj14<=400&&met>350&&met<=500","weight", {140.}),sams[iyr], all_plot_types);
      //D4
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 700., "mt", "m_{T} [GeV]",
                                        "mj14>250&&"+cut2l_base+"&&mj14>400&&met<=350","weight", {140.}),sams[iyr], all_plot_types);
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 700., "mt", "m_{T} [GeV]",
                                        "mj14>250&&"+cut2l_base+"&&mj14>400&&met>350&&met<=500","weight", {140.}),sams[iyr], all_plot_types);

      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 12, 200., 500., "met", "MET [GeV]",
                                        "mj14>250&&"+cut2l_base+"&&met<=500","weight", {200}),sams[iyr], all_plot_types);

      //MJ12 Plots
      ReplaceAll(cut2l_base,"mj14","mj");

      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 1000., "mj", "M_{J}^{1.2} [GeV]",
                                        cut2l_base+"&&met<=350","weight", {250, 400.}),sams[iyr], all_plot_types);
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 1000., "mj", "M_{J}^{1.2} [GeV]",
                                        cut2l_base+"&&met>350&&met<=500","weight", {250, 400.}),sams[iyr], all_plot_types);
      //MJ low mt di-lepton
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 1000., "mj", "M_{J}^{1.2} [GeV]",
                                        cut2l_base+"&&met<=350&&mt<=140","weight", {250, 400}),sams[iyr], all_plot_types);
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 1000., "mj", "M_{J}^{1.2} [GeV]",
                                        cut2l_base+"&&met>350&&met<=500&&mt<=140","weight", {250, 400}),sams[iyr], all_plot_types);
      //MJ high mt di-lepton
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 1000., "mj", "M_{J}^{1.2} [GeV]",
                                        cut2l_base+"&&met<=350&&mt>140","weight", {250, 400}),sams[iyr], all_plot_types);
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 1000., "mj", "M_{J}^{1.2} [GeV]",
                                        cut2l_base+"&&met>350&&met<=500&&mt>140","weight", {250, 400}),sams[iyr], all_plot_types);

      //D3 and D4
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 700., "mt", "m_{T} [GeV]",
                                        "mj>250&&"+cut2l_base+"&&met<=350","weight", {140.}),sams[iyr], all_plot_types);
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 700., "mt", "m_{T} [GeV]",
                                        "mj>250&&"+cut2l_base+"&&met>350&&met<=500","weight", {140.}),sams[iyr], all_plot_types);
      //D3
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 700., "mt", "m_{T} [GeV]",
                                        "mj>250&&"+cut2l_base+"&&mj<=400&&met<=350","weight", {140.}),sams[iyr], all_plot_types);
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 700., "mt", "m_{T} [GeV]",
                                        "mj>250&&"+cut2l_base+"&&mj<=400&&met>350&&met<=500","weight", {140.}),sams[iyr], all_plot_types);
      //D4
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 700., "mt", "m_{T} [GeV]",
                                        "mj>250&&"+cut2l_base+"&&mj>400&&met<=350","weight", {140.}),sams[iyr], all_plot_types);
      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 20, 0., 700., "mt", "m_{T} [GeV]",
                                        "mj>250&&"+cut2l_base+"&&mj>400&&met>350&&met<=500","weight", {140.}),sams[iyr], all_plot_types);

      pm[iyr].Push<HistoStack>(HistoDef(sam_tag[iyr]+"_"+lep_tag[ilep], 12, 200., 500., "met", "MET [GeV]",
                                        "mj>250&&"+cut2l_base+"&&met<=500","weight", {200}),sams[iyr], all_plot_types);
    }
    //    pm[iyr].MakePlots(lumi[iyr]);
  }

  //MJ MC Comparison
  //2015 Samples
  auto allmc_lowmt_2015 = Process::MakeShared<Baby_full>("2015, low m_{T}", Process::Type::background, colors("tt_2l"),
    {mc_2015+"*_TTJets*Lept*.root", mc_2015+"*_TTJets_HT*.root",
        mc_2015+"*_TTJets*Lept*.root", mc_2015+"*_TTJets_HT*.root",
        mc_2015+"*_WJetsToLNu*.root",
        mc_2015+"*_ST_*.root",
        mc_2015+"*_TTWJets*.root", mc_2015+"*_TTZTo*.root",
        mc_2015+"*DYJetsToLL*.root", mc_2015+"*_QCD_HT*.root",
        mc_2015+"*_ZJet*.root", mc_2015+"*_WWTo*.root",
        mc_2015+"*ggZH_HToBB*.root", mc_2015+"*ttHJetTobb*.root",
        mc_2015+"*_TTGJets*.root", mc_2015+"*_TTTT_*.root",
        mc_2015+"*_WH_HToBB*.root", mc_2015+"*_WZTo*.root",
        mc_2015+"*_ZH_HToBB*.root", mc_2015+"_ZZ_*.root"},
    "stitch&&mt<=140");

  auto allmc_highmt_2015 = Process::MakeShared<Baby_full>("2015, low m_{T}", Process::Type::data, kBlack,
    {mc_2015+"*_TTJets*Lept*.root", mc_2015+"*_TTJets_HT*.root",
        mc_2015+"*_TTJets*Lept*.root", mc_2015+"*_TTJets_HT*.root",
        mc_2015+"*_WJetsToLNu*.root",
        mc_2015+"*_ST_*.root",
        mc_2015+"*_TTWJets*.root", mc_2015+"*_TTZTo*.root",
        mc_2015+"*DYJetsToLL*.root", mc_2015+"*_QCD_HT*.root",
        mc_2015+"*_ZJet*.root", mc_2015+"*_WWTo*.root",
        mc_2015+"*ggZH_HToBB*.root", mc_2015+"*ttHJetTobb*.root",
        mc_2015+"*_TTGJets*.root", mc_2015+"*_TTTT_*.root",
        mc_2015+"*_WH_HToBB*.root", mc_2015+"*_WZTo*.root",
        mc_2015+"*_ZH_HToBB*.root", mc_2015+"_ZZ_*.root"},
    "stitch&&mt>=140");

  //2016 Samples
  auto allmc_lowmt_2016 = Process::MakeShared<Baby_full>("2016, low m_{T}", Process::Type::background, colors("tt_2l"),
    {mc_2016+"*_TTJets*Lept*.root", mc_2016+"*_TTJets_HT*.root",
        mc_2016+"*_TTJets*Lept*.root", mc_2016+"*_TTJets_HT*.root",
        mc_2016+"*_WJetsToLNu*.root",
        mc_2016+"*_ST_*.root",
        mc_2016+"*_TTWJets*.root", mc_2016+"*_TTZTo*.root",
        mc_2016+"*DYJetsToLL*.root", mc_2016+"*_QCD_HT*.root",
        mc_2016+"*_ZJet*.root", mc_2016+"*_WWTo*.root",
        mc_2016+"*ggZH_HToBB*.root", mc_2016+"*ttHJetTobb*.root",
        mc_2016+"*_TTGJets*.root", mc_2016+"*_TTTT_*.root",
        mc_2016+"*_WH_HToBB*.root", mc_2016+"*_WZTo*.root",
        mc_2016+"*_ZH_HToBB*.root", mc_2016+"_ZZ_*.root"},
    "stitch&&mt<=140");

  auto allmc_highmt_2016 = Process::MakeShared<Baby_full>("2016, high m_{T}", Process::Type::data, kBlack,
    {mc_2016+"*_TTJets*Lept*.root", mc_2016+"*_TTJets_HT*.root",
        mc_2016+"*_TTJets*Lept*.root", mc_2016+"*_TTJets_HT*.root",
        mc_2016+"*_WJetsToLNu*.root",
        mc_2016+"*_ST_*.root",
        mc_2016+"*_TTWJets*.root", mc_2016+"*_TTZTo*.root",
        mc_2016+"*DYJetsToLL*.root", mc_2016+"*_QCD_HT*.root",
        mc_2016+"*_ZJet*.root", mc_2016+"*_WWTo*.root",
        mc_2016+"*ggZH_HToBB*.root", mc_2016+"*ttHJetTobb*.root",
        mc_2016+"*_TTGJets*.root", mc_2016+"*_TTTT_*.root",
        mc_2016+"*_WH_HToBB*.root", mc_2016+"*_WZTo*.root",
        mc_2016+"*_ZH_HToBB*.root", mc_2016+"_ZZ_*.root"},
    "stitch&&mt>=140");

  vector<shared_ptr<Process> > mc2015 = {allmc_lowmt_2015, allmc_highmt_2015};
  vector<shared_ptr<Process> > mc2016 = {allmc_lowmt_2016, allmc_highmt_2016};

  vector<PlotMaker>  pm_mc(2);
  vector<vector<shared_ptr<Process> > > sams_mc = {mc2015, mc2016};  vector<string> sam_mc_tag = {"mc2015","mc2016"};
  vector<string> leps_mc = {"nleps==2","nels==2&&nmus==0","nels==0&&nmus==2","nels==1&&nmus==1"};
  vector<string> lep_mc_tag = {"2leps","2els","2mus","1el1mu"};

  // Loop over 2015 and 2016 samples
  for(unsigned int iyr=0; iyr<2; iyr++){
    // Loop over lepton flavors
    for(unsigned int ilep=0; ilep<leps_mc.size(); ilep++){

      string cut2l_base = "ht>500&&met>200&&"+leps_mc[ilep]+"&&njets>=5&&nbm<=2";

      //MJ
      pm_mc[iyr].Push<HistoStack>(HistoDef(sam_mc_tag[iyr]+"_"+lep_mc_tag[ilep], 20, 0., 1000., "mj14", "M_{J}^{1.4} [GeV]",
                                           cut2l_base+"&&met<=350","weight", {250, 400}),sams_mc[iyr], mc_plot_types);
      pm_mc[iyr].Push<HistoStack>(HistoDef(sam_mc_tag[iyr]+"_"+lep_mc_tag[ilep], 20, 0., 1000., "mj14", "M_{J}^{1.4} [GeV]",
                                           cut2l_base+"&&met>350&&met<=500","weight", {250, 400}),sams_mc[iyr], mc_plot_types);

      ReplaceAll(cut2l_base,"mj14","mj");

      pm_mc[iyr].Push<HistoStack>(HistoDef(sam_mc_tag[iyr]+"_"+lep_mc_tag[ilep], 20, 0., 1000., "mj", "M_{J}^{1.2} [GeV]",
                                           cut2l_base+"&&met<=350","weight", {250, 400.}),sams_mc[iyr], mc_plot_types);
      pm_mc[iyr].Push<HistoStack>(HistoDef(sam_mc_tag[iyr]+"_"+lep_mc_tag[ilep], 20, 0., 1000., "mj", "M_{J}^{1.2} [GeV]",
                                           cut2l_base+"&&met>350&&met<=500","weight", {250, 400.}),sams_mc[iyr], mc_plot_types);
    }
    pm_mc[iyr].MakePlots(lumi[iyr]);
  }
}
