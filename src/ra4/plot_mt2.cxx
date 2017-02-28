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
NamedFunc::VectorType tk_mt2(const Baby &b, int flavor=0);
NamedFunc::VectorType tk_truth(const Baby &b, int flavor=0);
NamedFunc::ScalarType n_tk_fake(const Baby &b, int flavor=0);

int main(){
  gErrorIgnoreLevel = 6000;

  double lumi = 35.9;


  string trig_skim1l_mc = "/net/cms2/cms2r0/babymaker/babies/2016_08_10/mc/skim_nlepsGE1__stG500__metG200__njetsGE5/";
  
  string unskimmed = "/net/cms29/cms29r0/babymaker/babies/2017_01_27/mc/unskimmed/";

  Palette colors("txt/colors.txt", "default");

  auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {trig_skim1l_mc+"*_TTJets*Lept*.root", trig_skim1l_mc+"*_TTJets_HT*.root"}, "ntruleps<=1&&stitch");
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {trig_skim1l_mc+"*_TTJets*Lept*.root", trig_skim1l_mc+"*_TTJets_HT*.root"}, "ntruleps>=2&&stitch");
  auto wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {trig_skim1l_mc+"*_WJetsToLNu*.root"},"stitch");
  auto single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {trig_skim1l_mc+"*_ST_*.root"});
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {trig_skim1l_mc+"*_TTWJets*.root", trig_skim1l_mc+"*_TTZTo*.root"});
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {trig_skim1l_mc+"*DYJetsToLL*.root", trig_skim1l_mc+"*_QCD_HT*.root",
        trig_skim1l_mc+"*_ZJet*.root", trig_skim1l_mc+"*_WWTo*.root",
        trig_skim1l_mc+"*ggZH_HToBB*.root", trig_skim1l_mc+"*ttHJetTobb*.root",
        trig_skim1l_mc+"*_TTGJets*.root", trig_skim1l_mc+"*_TTTT_*.root",
        trig_skim1l_mc+"*_WH_HToBB*.root", trig_skim1l_mc+"*_WZTo*.root",
        trig_skim1l_mc+"*_ZH_HToBB*.root", trig_skim1l_mc+"_ZZ_*.root"});

  // auto t1tttt_c = Process::MakeShared<Baby_full>("T1tttt(1200,800)", Process::Type::signal, colors("t1tttt"),
  //   {trig_mc+"*SMS-T1tttt_mGluino-1200_mLSP-800*.root"});
  // t1tttt_c->SetLineStyle(2);

  auto data_1l = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {"/net/cms2/cms2r0/babymaker/babies/2017_01_27/data/skim_stdnj5/*.root"},"pass&&trig_ra4");
  vector<shared_ptr<Process> > full_trig_skim_1l = {data_1l, /*t1tttt_nc, t1tttt_c,*/ tt1l, tt2l, wjets, single_t, ttv, other};

  NamedFunc ntrk_fake("ntrk_fake",[&](const Baby &b){
      return n_tk_fake(b);
     
    });
  NamedFunc ntrk_lep_fake("ntrk_lep_fake",[&](const Baby &b){
      int n = n_tk_fake(b);
      if(b.nleps()>b.nleps_tm()) n+=b.nleps()-b.nleps_tm();
      return n;
    });

  auto tt2l_bear = Process::MakeShared<Baby_full>("t#bar{t} (2l), all", Process::Type::background, colors("tt_2l"),
    {unskimmed+"*_TTJets*DiLept*.root"}, "stitch_met");
  auto tt2l_bear_no_fake = Process::MakeShared<Baby_full>("t#bar{t} (2l), no fake track or lepton", Process::Type::background, kGreen-3,
    {unskimmed+"*_TTJets*DiLept*.root"}, "stitch_met"&&ntrk_lep_fake==0.); // && (met-met_tru_nuw)*(met-met_tru_nuw)<=225)

  auto tt2l_bear_fake_tk= Process::MakeShared<Baby_full>("t#bar{t} (2l), has fake track or lepton", Process::Type::background, colors("other"),
    {unskimmed+"*_TTJets*DiLept*.root"}, "stitch_met"&&ntrk_lep_fake>0.);

  //auto tt2l_bear_fake_met = Process::MakeShared<Baby_full>("t#bar{t} (2l), fake met but no fake tracks", Process::Type::background, kRed,
  //{unskimmed+"*_TTJets*DiLept*.root"}, "stitch_met && (met-met_tru)*(met-met_tru)>225 "&&ntrk_fake==0.);//|| (met-met_tru_nuw)*(met-met_tru_nuw)>225)

  //auto tt2l_bear_no_b_met = Process::MakeShared<Baby_full>("t#bar{t} (2l), no non-prompt met", Process::Type::background, kGreen-3,
  // {unskimmed+"*_TTJets*DiLept*.root"}, "stitch_met && (met-met_tru_nuw)*(met-met_tru_nuw)<400 ");


  //  auto t1tttt_nc = Process::MakeShared<Baby_full>("T1tttt(1800,100)", Process::Type::signal, colors("t1tttt"),
  //{"/net/cms29/cms29r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/*SMS-T1tttt_mGluino-1800_mLSP-100_*.root"});

  vector<shared_ptr<Process> > tt2l_variations = {tt2l_bear,tt2l_bear_fake_tk,tt2l_bear_no_fake};

  NamedFunc mt2("mt2",[&](const Baby &b){
      return tk_mt2(b);
    });

  NamedFunc mt2_lep("mt2_lep",[&](const Baby &b){
      return tk_mt2(b,1);
    });

  NamedFunc mt2_el("mt2_el",[&](const Baby &b){
      return tk_mt2(b,11);
    });

  NamedFunc mt2_mu("mt2_mu",[&](const Baby &b){
      return tk_mt2(b,13);
    });

  NamedFunc mt2_had("mt2_had",[&](const Baby &b){
      return tk_mt2(b,211);
    });


  NamedFunc isTM_lep("isTM_lep",[&](const Baby &b){
      return tk_truth(b,1);
    });

  NamedFunc isTM_el("isTM_el",[&](const Baby &b){
      return tk_truth(b,11);
    });

  NamedFunc isTM_mu("isTM_mu",[&](const Baby &b){
      return tk_truth(b,13);
    });

  NamedFunc isTM_had("isTM_had",[&](const Baby &b){
      return tk_truth(b,211);
    });

  NamedFunc ntrk("ntrk",[&](const Baby &b){
      return tk_mt2(b).size();
    });
  NamedFunc ntrk_lep("ntrk_lep",[&](const Baby &b){
      return tk_mt2(b,1).size();
    });

  NamedFunc ntrk_had("ntrk_had",[&](const Baby &b){
      return tk_mt2(b,211).size();
    });



  /*  NamedFunc ntrk_lep_fake("ntrk_lep_fake",[&](const Baby &b){
      int nfake = 0;
      vector<double> tksTM = tk_truth(b,1);
      for(unsigned int i=0; i<tksTM.size();i++){
	if(tksTM==0.) nfake++;
      }
      return nfake;
    });
  NamedFunc ntrk_had_fake("ntrk_had_fake",[&](const Baby &b){
      int nfake = 0;
      vector<double> tksTM = tk_truth(b,211);
      for(unsigned int i=0; i<tksTM.size();i++){
	if(tksTM==0.) nfake++;
      }
      return nfake;
      });*/

  //
  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    //.Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::lumi_shapes);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .Bottom(BottomType::off)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info);
  vector<PlotOpt> all_plot_types = {log_lumi_info, lin_lumi_info,log_shapes_info, lin_shapes_info};
  vector<PlotOpt> shapes = {log_shapes_info, lin_shapes_info};
  
  PlotMaker pm;
  vector<NamedFunc> baselines = {"nleps==1","nleps==1&&met>80&&njets>=2&&nbm>=1","nleps==1&&met<80","nleps==1&&met>200","nleps==1&&st>500&&met>200&&nbm>=1&&njets>=6&&pass&&met/met_calo<5.","mt>140&&mj14>250&&nleps==1&&st>500&&met>200&&nbm>=1&&njets>=6&&pass&&met/met_calo<5."};
    //{"nleps==1&&st>500&&met>200&&nbm_moriond>=1&&njets>=6&&pass&&met/met_calo<5.","nleps==1&&mt>140&&st>500&&met>200&&nbm_moriond>=1&&njets>=6&&pass&&met/met_calo<5."};
  
  for(auto baseline: baselines){
    // pm.Push<Hist1D>(Axis(8, 0 ,160, mt2, "m_{T2} [GeV]",{60.,80}),
    //		    tt2l,
    //		    baseline_variations, all_plot_types).Tag("tracks/");
    pm.Push<Hist1D>(Axis(30, 0 ,300, mt2_lep, "Leptonic track m_{T2} [GeV]",{80.}),
		    baseline,
		    tt2l_variations, all_plot_types).Tag("tracks");

    pm.Push<Hist1D>(Axis(17, 0 ,340, mt2_lep, "Leptonic track m_{T2} [GeV]",{80.}),
		    baseline,
		    tt2l_variations, all_plot_types).Tag("tracks_sus17001");

    pm.Push<Hist1D>(Axis(7, 0 ,140, mt2_lep, "Leptonic track m_{T2} [GeV]",{80.}),
		    baseline,
		    tt2l_variations, all_plot_types).Tag("tracks_an_binning");

 pm.Push<Hist1D>(Axis(30, 0 ,300, mt2_mu, "Muon track m_{T2} [GeV]",{80.}),
		    baseline,
		    tt2l_variations, all_plot_types).Tag("tracks");

 pm.Push<Hist1D>(Axis(30, 0 ,300, mt2_el, "Electron track m_{T2} [Gev]",{80.}),
		    baseline,
		    tt2l_variations, all_plot_types).Tag("tracks");


    pm.Push<Hist1D>(Axis(30, 0 ,300, mt2_had, "Hadronic track m_{T2} [GeV]",{60.}),
		    baseline,
		    tt2l_variations, all_plot_types).Tag("tracks");
   pm.Push<Hist1D>(Axis(7, 0 ,140, mt2_had, "Hadronic track m_{T2} [GeV]",{60.}),
		    baseline,
		    tt2l_variations, all_plot_types).Tag("tracks_an_binning");



    /* pm.Push<Hist1D>(Axis(3, -0.5, 2.5, ntrk, "N_{trks}"),
		     baseline,
		     full_trig_skim_1l, all_plot_types).Tag("tracks/");

     pm.Push<Hist1D>(Axis(3, -0.5, 2.5, ntrk_lep, "N_{trks, lep}"),
		     baseline,
		     full_trig_skim_1l, all_plot_types).Tag("tracks/");

     pm.Push<Hist1D>(Axis(3, -0.5, 2.5, ntrk_had, "N_{trks, had}"),
		     baseline,
		     full_trig_skim_1l, all_plot_types).Tag("tracks/");
     
     pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
		     baseline,
		     full_trig_skim_1l, all_plot_types).Tag("tracks/");
	
       pm.Push<Hist1D>(Axis(12, 0, 420., "mt", "m_{T} [GeV]", {140}),
		       baseline,
		       full_trig_skim_1l, all_plot_types).Tag("tracks/");*/


   }
  pm.min_print_ = true;
  pm.MakePlots(lumi);




}
 

NamedFunc::VectorType tk_mt2(const Baby &b, int flavor){
  vector<double> mt2;
  for(unsigned int itk=0;itk< b.tks_pt()->size();itk++){

    if(flavor==0 ||flavor==211){
      if (fabs(b.tks_pdg()->at(itk))==211  && b.tks_pt()->at(itk)>15. && b.tks_miniso()->at(itk)<0.1 &&  b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05) mt2.push_back(b.tks_mt2()->at(itk));						   }
   
    if(flavor==13 || flavor==0 || flavor==1){
      if (fabs(b.tks_pdg()->at(itk))==13 && b.tks_pt()->at(itk)>10. && b.tks_miniso()->at(itk)<0.2 && b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05) mt2.push_back(b.tks_mt2()->at(itk));	
    }
	
    if(flavor==11 || flavor==0 || flavor==1){

      if (fabs(b.tks_pdg()->at(itk))==11 && b.tks_pt()->at(itk)>10. && b.tks_miniso()->at(itk)<0.2 && b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05) mt2.push_back(b.tks_mt2()->at(itk));	
    }

  }

  return mt2;
}

NamedFunc::VectorType tk_truth(const Baby &b, int flavor){
  vector<double> isTM;
  for(unsigned int itk=0;itk< b.tks_pt()->size();itk++){

    if(flavor==0 ||flavor==211){
      if (fabs(b.tks_pdg()->at(itk))==211  && b.tks_pt()->at(itk)>15. && b.tks_miniso()->at(itk)<0.1 &&  b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05) isTM.push_back(b.tks_tm()->at(itk));						   }
   
    if(flavor==13 || flavor==0 || flavor==1){
      if (fabs(b.tks_pdg()->at(itk))==13 && b.tks_pt()->at(itk)>10. && b.tks_miniso()->at(itk)<0.2 && b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05) isTM.push_back(b.tks_tm()->at(itk));	
    }
	
    if(flavor==11 || flavor==0 || flavor==1){

      if (fabs(b.tks_pdg()->at(itk))==11 && b.tks_pt()->at(itk)>10. && b.tks_miniso()->at(itk)<0.2 && b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05) isTM.push_back(b.tks_tm()->at(itk));	
    }

  }

  return isTM;
}


NamedFunc::ScalarType n_tk_fake(const Baby &b, int flavor){
  int n_fake=0;
  for(unsigned int itk=0;itk< b.tks_pt()->size();itk++){

    if(flavor==0 ||flavor==211){
      if (fabs(b.tks_pdg()->at(itk))==211  && b.tks_pt()->at(itk)>15. && b.tks_miniso()->at(itk)<0.1 &&  b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05 && !b.tks_tm()->at(itk)) n_fake++; 						   }
   
    if(flavor==13 || flavor==0 || flavor==1){
      if (fabs(b.tks_pdg()->at(itk))==13 && b.tks_pt()->at(itk)>10. && b.tks_miniso()->at(itk)<0.2 && b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05 && !b.tks_tm()->at(itk)) n_fake++;	
    }
	
    if(flavor==11 || flavor==0 || flavor==1){

      if (fabs(b.tks_pdg()->at(itk))==11 && b.tks_pt()->at(itk)>10. && b.tks_miniso()->at(itk)<0.2 && b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05 && !b.tks_tm()->at(itk)) n_fake++;	
    }

  }

  return n_fake;
}




