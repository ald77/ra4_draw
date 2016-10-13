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

double getWeight(const Baby &b, int bin);
NamedFunc::ScalarType nEls25(const Baby &b);
NamedFunc::ScalarType nElsTight25(const Baby &b);
NamedFunc::ScalarType nMus25(const Baby &b);

int main(){
  gErrorIgnoreLevel = 6000;

  double lumi = 12.9;

  string mc_folder = "/net/cms29/cms29r0/babymaker/babies/2016_08_10/mc/skim_nleps1__metG250__htG250/";
  string data_folder = "/net/cms29/cms29r0/babymaker/babies/2016_08_10/data/skim_nleps1__metG250__htG250/";

  Palette colors("txt/colors.txt", "default");

  auto data = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {data_folder+"*.root"},"pass&&(trig[13]||trig[33])&&json12p9");

  auto tt_1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {mc_folder+"*_TTJets*Lept*.root"},"pass&&ntruleps<=1");
  auto tt_2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {mc_folder+"*_TTJets*Lept*.root"},"pass&&ntruleps>=2");
  auto wjets = Process::MakeShared<Baby_full>("W+Jets", Process::Type::background, colors("wjets"),
    {mc_folder+"*_WJets*HT*.root"},"pass");
  auto st = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {mc_folder+"*_ST_s-channel*.root", mc_folder+"*_ST_t-channel*.root", mc_folder+"_ST_tW_*NoFullyHadronicDecays*.root"},"pass");
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {mc_folder+"*_TTWJets*.root", mc_folder+"*_TTZTo*.root"}, "pass");
  auto vv = Process::MakeShared<Baby_full>("VV", Process::Type::background, kRed,
    {mc_folder+"_ZZ_*.root", mc_folder+"*_WWTo*.root",mc_folder+"*_WZTo*.root"},"pass");
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {mc_folder+"*DYJetsToLL*.root",    mc_folder+"*_QCD_HT*.root", 
        mc_folder+"*_ZJet*.root",      mc_folder+"*ggZH_HToBB*.root", 
	mc_folder+"*ttHJetTobb*.root", mc_folder+"*_TTGJets*.root",
	mc_folder+"*_TTTT_*.root",     mc_folder+"*_WH_HToBB*.root", 
	mc_folder+"*_ZH_HToBB*.root"},"pass");


  //vector<shared_ptr<Process> > full_sam = {data, tt_1l, tt_2l, wjets, st};
  vector<shared_ptr<Process> > full_sam = {data, tt_1l, tt_2l, wjets, st, ttv, vv, other};
  
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
  vector<PlotOpt> all_plot_types = {log_lumi_info};

  PlotMaker pm;

  NamedFunc nels25("nels",nEls25);
  NamedFunc nmus25("nmus",nMus25);
  //  NamedFunc nelstight25("n_{e}^{t}",nElsTight25);

  NamedFunc baseline = "ht>250 && met>250 && njets>=2";

  vector<NamedFunc> lepcuts;
  lepcuts.push_back(nels25+nmus25==1);
  //  lepcuts.push_back(nels25==1);
  //  lepcuts.push_back(nmus25==1);
  
  //  lepcuts.push_back(nelstight25+nmus25==1);
  //  lepcuts.push_back(nelstight25==1);


  for(const auto lepcut: lepcuts){      

    vector<NamedFunc> ctrlregs;
    ctrlregs.push_back("nbm==0");
    ctrlregs.push_back("nbm>=1");
    ctrlregs.push_back("nbm>=2");

    for(const auto ctrlreg: ctrlregs){      

      vector<NamedFunc> weights;
      //      weights.push_back(NamedFunc("No SF",              [](const Baby &b) -> NamedFunc::ScalarType{return getWeight(b,0);}));
      //      weights.push_back(NamedFunc("w_lep",              [](const Baby &b) -> NamedFunc::ScalarType{return getWeight(b,1);}));
      //      weights.push_back(NamedFunc("w_btag",             [](const Baby &b) -> NamedFunc::ScalarType{return getWeight(b,2);}));
      //      weights.push_back(NamedFunc("w_isr",              [](const Baby &b) -> NamedFunc::ScalarType{return getWeight(b,3);}));
      weights.push_back(NamedFunc("w_lep*w_btag",       [](const Baby &b) -> NamedFunc::ScalarType{return getWeight(b,4);}));
      //      weights.push_back(NamedFunc("w_lep*w_isr",        [](const Baby &b) -> NamedFunc::ScalarType{return getWeight(b,5);}));
      //weights.push_back(NamedFunc("w_btag*w_isr",       [](const Baby &b) -> NamedFunc::ScalarType{return getWeight(b,6);}));
      weights.push_back(NamedFunc("w_lep*w_btag*w_isr", [](const Baby &b) -> NamedFunc::ScalarType{return getWeight(b,7);}));
      
      for(const auto &iwght: weights){      
	pm.Push<Hist1D>(Axis(25, 250., 1500., "ht", "H_{T} [GeV]"), baseline && ctrlreg && lepcut, full_sam, all_plot_types).Weight(iwght);
	pm.Push<Hist1D>(Axis(20, 250., 1250., "met", "MET [GeV]"), baseline && ctrlreg && lepcut, full_sam, all_plot_types).Weight(iwght);
	pm.Push<Hist1D>(Axis(10, 2., 12., "njets", "N_{jets}"), baseline && ctrlreg && lepcut, full_sam, all_plot_types).Weight(iwght);
	pm.Push<Hist1D>(Axis(5, 0., 5., "nbm", "N_{b}^{med}"), baseline && ctrlreg && lepcut, full_sam, all_plot_types).Weight(iwght);
	pm.Push<Hist1D>(Axis(20, 0, 1000, "mj14", "M_{J} [GeV]"), baseline && ctrlreg && lepcut, full_sam, all_plot_types).Weight(iwght);
	pm.Push<Hist1D>(Axis(12, 0, 420, "mt", "m_{T} [GeV]"), baseline && ctrlreg && lepcut, full_sam, all_plot_types).Weight(iwght);
      }
    }
  }
  pm.MakePlots(lumi);
}


double getWeight(const Baby &b, int bin){

  double weight = b.w_lumi()*b.w_pu();

  if(bin==1) 
    weight *= b.w_lep();
  else if(bin==2) 
    weight *= b.w_btag();
  else if(bin==3){ 
    if(b.ntrupv()>=0) weight *= b.w_isr(); //Only apply to non-data samples
  }
  else if(bin==4) 
    weight *= b.w_lep()*b.w_btag();
  else if(bin==5){
    weight *= b.w_lep();
    if(b.ntrupv()>=0) weight *= b.w_isr();
  }
  else if(bin==6){
    weight *= b.w_btag();
    if(b.ntrupv()>=0) weight *= b.w_isr();
  }
  else if(bin==7){
    weight *= b.w_lep()*b.w_btag();
    if(b.ntrupv()>=0) weight *= b.w_isr();
  }

  return weight;
}

NamedFunc::ScalarType nEls25(const Baby &b) {
  int nels = 0;
  for (size_t iel(0); iel<b.els_pt()->size(); iel++){
    if (b.els_pt()->at(iel)>25 && b.els_sig()->at(iel) && b.els_miniso()->at(iel)<0.1 && fabs(b.els_eta()->at(iel))<=2.4)
      nels++;
  }
  return nels;
}

NamedFunc::ScalarType nElsTight25(const Baby &b) {
  int nels = 0;
  for (size_t iel(0); iel<b.els_pt()->size(); iel++){
    if (b.els_pt()->at(iel)>25 && b.els_tight()->at(iel) && b.els_miniso()->at(iel)<0.1 && fabs(b.els_eta()->at(iel))<=2.4)
      nels++;
  }
  return nels;
}

NamedFunc::ScalarType nMus25(const Baby &b) {
  int nmus = 0;
  for (size_t imu(0); imu<b.mus_pt()->size(); imu++){
    if (b.mus_pt()->at(imu)>25 && b.mus_sig()->at(imu) && b.mus_miniso()->at(imu)<0.2 && fabs(b.mus_eta()->at(imu))<=2.4)
       nmus++;
  }
  return nmus;
}
