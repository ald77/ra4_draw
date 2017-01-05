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
#include "core/hist2d.hpp"

using namespace std;
using namespace PlotOptTypes;
NamedFunc::ScalarType LeptonicWMass(const Baby &b);
NamedFunc::ScalarType LeptonicWTransverseMass(const Baby &b);

int main(){
  gErrorIgnoreLevel = 6000;

  double lumi = 36.2;
  Palette colors("txt/colors.txt", "default");
  
  string standard_mc = "/net/cms2/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_standard/";
  string powheg = "/net/cms2/cms2r0/babymaker/babies/2016_11_29/mc/merged_mcbase_standard/";

  string standard_mc_unskim = "/net/cms2/cms2r0/babymaker/babies/2016_08_10/mc/unskimmed/";
  string powheg_unskim = "/net/cms2/cms2r0/babymaker/babies/2016_11_29/mc/unskimmed/";
  
  auto tt1l_mg = Process::MakeShared<Baby_full>("Madgraph LO t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {standard_mc+"*_TTJets*SingleLept*.root", standard_mc+"*_TTJets_HT*.root"}, "ntruleps<=1&&stitch");
  auto tt1l_powheg = Process::MakeShared<Baby_full>("Powheg t#bar{t} (1l)", Process::Type::background, colors("ttv"),
    {powheg+"*TTToSemiLeptonic_TuneCUETP8M1_alphaS01273_13TeV-powheg*.root"}, "ntruleps<=1&&stitch");
  auto wjets = Process::MakeShared<Baby_full>("Madgraph LO W+jets", Process::Type::background, colors("wjets"),
    {standard_mc+"*_WJetsToLNu*.root"},"stitch");


   auto tt1l_mg_unskim = Process::MakeShared<Baby_full>("Madgraph LO t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {standard_mc_unskim+"*_TTJets*SingleLept*.root", standard_mc+"*_TTJets_HT*.root"}, "ntruleps<=1&&stitch");
  auto tt1l_powheg_unskim = Process::MakeShared<Baby_full>("Powheg t#bar{t} (1l)", Process::Type::background, colors("ttv"),
    {powheg_unskim+"*TTToSemiLeptonic_TuneCUETP8M1_alphaS01273_13TeV-powheg*.root"}, "ntruleps<=1&&stitch");
  auto wjets_unskim = Process::MakeShared<Baby_full>("Madgraph LO W+jets", Process::Type::background, colors("wjets"),
    {standard_mc_unskim+"*_WJetsToLNu*.root"},"stitch");

  vector<shared_ptr<Process> > ttbar_wjets = {wjets, tt1l_powheg, tt1l_mg};
  vector<shared_ptr<Process> > powheg_wjets = {wjets, tt1l_powheg};
  vector<shared_ptr<Process> > ttbar_wjets_unskim = {wjets_unskim, tt1l_powheg_unskim, tt1l_mg_unskim};
  //
  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary).YAxis(YAxisType::log).Stack(StackType::lumi_shapes).ShowBackgroundError(false);
    //.Bottom(BottomType::ratio)
     
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .Bottom(BottomType::off)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info);
  vector<PlotOpt> all_plot_types = {lin_lumi_info,log_lumi_info, log_shapes_info};

  PlotMaker pm;
  
  NamedFunc abcd = "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6&&nbm>=1&&mj14>250&&(mc_id*mc_id)==576";
  NamedFunc highmet= "met>500";
  NamedFunc highnjets= "njets>=9";
  NamedFunc highmass= "mc_mass>125";

  NamedFunc mW_lep("mW_lep", LeptonicWMass);
  NamedFunc mTW_tru("mTW_tru", LeptonicWTransverseMass);

  /* pm.Push<Hist1D>(Axis(20, 0, 280., mW_lep, "m_{W} leptonic [GeV]", {140.}),
		  "1",
		  ttbar_wjets_unskim, all_plot_types);

  pm.Push<Hist1D>(Axis(20, 0, 280., "mt_tru_nuw", "True W m_{T} [GeV]", {140.}),
		  "1",
		  ttbar_wjets_unskim, all_plot_types);

  pm.Push<Hist1D>(Axis(20, 0, 280., mTW_tru , "True W m_{T} [GeV]", {140.}),
		  "1",
		  ttbar_wjets_unskim, all_plot_types);

 // pm.Push<Hist2D>(Axis(40, 0, 280., mW_lep, "m_{W} leptonic [GeV]", {140.}), Axis(40, 0, 280., "mt_tru_nuw", "True W m_{T} [GeV]", {140.}),
 // 		  "1",
 // 		  ttbar_wjets_unskim, all_plot_types);


 //  pm.Push<Hist2D>(Axis(40, 0, 280., mW_lep, "m_{W} leptonic [GeV]", {140.}), Axis(40, 0, 280., mTW_tru, "True W m_{T} [GeV]", {140.}),
 // 		  "1",
 // 		  ttbar_wjets_unskim, all_plot_types);

  pm.Push<Hist1D>(Axis(20, 0, 280., mW_lep, "m_{W} leptonic [GeV]", {140.}),
		  "mt>140",
		  ttbar_wjets_unskim, all_plot_types);

  pm.Push<Hist1D>(Axis(20, 0, 280., "mt", "m_{T} [GeV]", {140.}),
		  mW_lep>140,
		  ttbar_wjets_unskim, all_plot_types);
*/
  pm.Push<Hist1D>(Axis(20, 0, 280., "mt", "m_{T} [GeV]", {140.}),
		  "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6&&nbm>=1&&mj14>250"&&mW_lep>140,
		  ttbar_wjets, all_plot_types);


  
  pm.Push<Hist1D>(Axis(20, 0, 280., mW_lep, "m_{W} leptonic [GeV]", {140.}),
		  "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6&&nbm>=1&&mj14>250",
		  ttbar_wjets, all_plot_types);

  
  pm.Push<Hist1D>(Axis(20, 0, 280., mW_lep, "m_{W} leptonic [GeV]", {140.}),
		  "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6&&nbm>=1&&mj14>250&&mt>140",
		  ttbar_wjets, all_plot_types);

  pm.Push<Hist1D>(Axis(20, 0, 280., mW_lep, "m_{W} leptonic [GeV]", {140.}),
		  "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6&&nbm>=1&&mj14>250&&met>500",
		  ttbar_wjets, all_plot_types);

  
  pm.Push<Hist1D>(Axis(20, 0, 280., mW_lep, "m_{W} leptonic [GeV]", {140.}),
		  "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6&&nbm>=1&&mj14>250&&met>500&&mt>140",
		  ttbar_wjets, all_plot_types);


  pm.Push<Hist1D>(Axis(20, 0, 280.,mTW_tru, "True W m_{T} [GeV]", {140.}),
		  "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6&&nbm>=1&&mj14>250",
		  ttbar_wjets, all_plot_types);

  pm.Push<Hist1D>(Axis(20, 0, 280.,mTW_tru, "True W m_{T} [GeV]", {140.}),
		  "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6&&nbm>=1&&mj14>250&&met>500",
		  ttbar_wjets, all_plot_types);


  pm.Push<Hist1D>(Axis(20, 0, 280., "mt_tru_nuw", "True W m_{T} [GeV]", {140.}),
		  "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6&&nbm>=1&&mj14>250",
		  ttbar_wjets, all_plot_types);

  pm.Push<Hist1D>(Axis(20, 0, 280., "mt_tru_nuw", "True W m_{T} [GeV]", {140.}),
		  "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6&&nbm>=1&&mj14>250&&mt>140",
		  ttbar_wjets, all_plot_types);
  

    pm.Push<Hist1D>(Axis(20, 0, 280., "mt_tru", "True m_{T} [GeV]", {140.}),
		  "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6&&nbm>=1&&mj14>250",
		  ttbar_wjets, all_plot_types);
  

    pm.Push<Hist2D>(Axis(40, 0, 280., mW_lep, "m_{W} leptonic [GeV]", {140.}), Axis(40, 0, 280., "mt_tru_nuw", "True W m_{T}", {140.}),
		  "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6&&nbm>=1&&mj14>250",
		  ttbar_wjets, all_plot_types);

    pm.Push<Hist2D>(Axis(40, 0, 280., mW_lep, "m_{W} leptonic [GeV]", {140.}), Axis(40, 0, 280.,mTW_tru , "True W m_{T} (by hand)", {140.}),
		  "nleps==1&&nveto==0&&st>500&&met>200&&njets>=6&&nbm>=1&&mj14>250",
		  ttbar_wjets, all_plot_types);

    
    /* pm.Push<Hist1D>(Axis(20, 0, 280., "mc_mass", "m_{W} [GeV]", {140.}),
		  abcd,
		  ttbar_wjets, all_plot_types);

  pm.Push<Hist1D>(Axis(20, 0, 280., "mc_mass", "m_{W} [GeV]", {140.}),
		  abcd&&highmet,
		  ttbar_wjets, all_plot_types);

  pm.Push<Hist1D>(Axis(20, 0, 280., "mc_mass", "m_{W} [GeV]", {140.}),
		  abcd&&highnjets,
		  ttbar_wjets, all_plot_types);

  pm.Push<Hist1D>(Axis(20, 125, 225., "mc_mass", "m_{W} [GeV]", {140.}),
		  abcd&&highmass,
		  powheg_wjets, all_plot_types).Tag("zoom");

  pm.Push<Hist1D>(Axis(20, 125, 225., "mc_mass", "m_{W} [GeV]", {140.}),
		  abcd&&highmet&&highmass,
		  powheg_wjets, all_plot_types).Tag("zoom");

  pm.Push<Hist1D>(Axis(20, 125, 225., "mc_mass", "m_{W} [GeV]", {140.}),
		  abcd&&highnjets&&highmass,
		  powheg_wjets, all_plot_types).Tag("zoom");
    */
  pm.MakePlots(lumi);

}


//Find mass of leptonic W in event
NamedFunc::ScalarType LeptonicWMass(const Baby &b){
  float mass=0;
  for(size_t imc = 0; imc < b.mc_pt()->size(); ++imc){
    if(abs(b.mc_mom()->at(imc))==24 && (abs(b.mc_id()->at(imc)) == 11 || abs(b.mc_id()->at(imc)) == 13 || abs(b.mc_id()->at(imc)) == 15)){
      mass = b.mc_mass()->at(b.mc_momidx()->at(imc));
      break;
    }
  }
  return mass;
}


NamedFunc::ScalarType LeptonicWTransverseMass(const Baby &b){
  float pt_neu(0),phi_neu(0),pt_lep(0),phi_lep(0);
  for(size_t imc = 0; imc < b.mc_pt()->size(); ++imc){
    if(abs(b.mc_mom()->at(imc))==24 && (abs(b.mc_id()->at(imc)) == 11 || abs(b.mc_id()->at(imc)) == 13 || abs(b.mc_id()->at(imc)) == 15)){
      pt_lep=b.mc_pt()->at(imc);
      phi_lep=b.mc_phi()->at(imc);
    }
    if(abs(b.mc_mom()->at(imc))==24 && (abs(b.mc_id()->at(imc)) == 12 || abs(b.mc_id()->at(imc)) == 14 || abs(b.mc_id()->at(imc)) == 16)){
      pt_neu=b.mc_pt()->at(imc);
      phi_neu=b.mc_phi()->at(imc);
    }
  }

  return sqrt(2.*pt_neu*pt_lep*(1.-cos(phi_neu-phi_lep)));
}

