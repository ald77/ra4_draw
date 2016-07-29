#include "test.hpp"
#include "plot_dy_mismeasurement_study.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TColor.h"
#include "TVector2.h"

#include "baby.hpp"
#include "process.hpp"
#include "named_func.hpp"
#include "plot_maker.hpp"
#include "plot_opt.hpp"
#include "palette.hpp"
#include "table.hpp"
#include "histo_stack.hpp"
#include "event_scan.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  bool single_thread = false;
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  double lumi = 2.7;

  string trig_skim_mc = "/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/skim_nleps2/";

  Palette colors("txt/colors.txt", "default");

  auto dy = Process::MakeShared<Baby_full>("Drell Yan", Process::Type::background, colors("wjets"),
    {trig_skim_mc+"*DYJetsToLL*.root"},"stitch");
  auto diboson = Process::MakeShared<Baby_full>("WW,ZZ,WZ", Process::Type::background, colors("single_t"),
    {trig_skim_mc+"*_WWTo*.root",trig_skim_mc+"*_WZTo*.root", trig_skim_mc+"_ZZ_*.root"});
  auto tt = Process::MakeShared<Baby_full>("t#bar{t}", Process::Type::background, colors("tt_2l"),
    {trig_skim_mc+"*_TTJets*Lept*.root", trig_skim_mc+"*_TTJets_HT*.root"},
    "stitch");
  auto wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("ttv"),
    {trig_skim_mc+"*_WJetsToLNu*.root"});

  //auto wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
  //  {trig_skim_mc+"*_WJetsToLNu*.root"});
  //auto single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
  // {trig_skim_mc+"*_ST_*.root"});
  /* auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
     {trig_skim_mc+"*_TTWJets*.root", trig_skim_mc+"*_TTZTo*.root"});*/
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {/*trig_skim_mc+"*_WJetsToLNu*.root",*/ trig_skim_mc+"*_QCD_HT*.root",
        trig_skim_mc+"*_ZJet*.root",/* trig_skim_mc+"*_WWTo*.root",*/
        trig_skim_mc+"*ggZH_HToBB*.root", trig_skim_mc+"*ttHJetTobb*.root",
        trig_skim_mc+"*_TTGJets*.root", trig_skim_mc+"*_TTTT_*.root",
        trig_skim_mc+"*_WH_HToBB*.root"/*, trig_skim_mc+"*_WZTo*.root"*/,
        trig_skim_mc+"*_ZH_HToBB*.root",trig_skim_mc+"*_ST_*.root"/*, trig_skim_mc+"_ZZ_*.root"*/
        ,trig_skim_mc+"*_TTWJets*.root", trig_skim_mc+"*_TTZTo*.root"});

  auto t1tttt_nc = Process::MakeShared<Baby_full>("T1tttt(1500,100)", Process::Type::signal, colors("t1tttt"),
    {trig_skim_mc+"*SMS-T1tttt_mGluino-1500_mLSP-100*.root"});
  auto t1tttt_c = Process::MakeShared<Baby_full>("T1tttt(1200,800)", Process::Type::signal, colors("t1tttt"),
    {trig_skim_mc+"*SMS-T1tttt_mGluino-1200_mLSP-800*.root"});
  t1tttt_c->SetLineStyle(2);

  auto data = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {"/net/cms2/cms2r0/babymaker/babies/2016_06_26/data/skim_nleps2/*.root"},"json2p6&&pass&&(trig[20]||trig[21]||trig[23]||trig[29]||trig[24]||trig[4]||trig[8]||trig[13]||trig[33])");

  auto data_ee = Process::MakeShared<Baby_full>("Data, ee", Process::Type::data, kBlack,
    {"/net/cms2/cms2r0/babymaker/babies/2016_06_26/data/skim_nleps2/*.root"},"nels==2&&elel_m>60&&json2p6&&pass&&(trig[20]||trig[21]||trig[23]||trig[29]||trig[24]||trig[4]||trig[8]||trig[13]||trig[33])");

  auto data_mumu = Process::MakeShared<Baby_full>("Data, #mu#mu", Process::Type::background, 30,
    {"/net/cms2/cms2r0/babymaker/babies/2016_06_26/data/skim_nleps2/*.root"},"nmus==2&&mumu_m>60&&json2p6&&pass&&(trig[20]||trig[21]||trig[23]||trig[29]||trig[24]||trig[4]||trig[8]||trig[13]||trig[33])");

  vector<shared_ptr<Process> > full_trig_skim = {data, t1tttt_nc, t1tttt_c, dy, diboson, tt, wjets, other};
  vector<shared_ptr<Process> > ee_vs_mumu = {data_ee,data_mumu};
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
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::ratio);
  vector<PlotOpt> all_plot_types = {log_lumi, lin_lumi, log_shapes, lin_shapes,
                                    log_lumi_info, lin_lumi_info, log_shapes_info, lin_shapes_info};
  vector<PlotOpt> norms = {log_lumi_info, lin_lumi_info};
  vector<PlotOpt> shapes = {lin_shapes_info,log_shapes_info};

  PlotMaker pm;
  double pi = acos(-1.);

  NamedFunc dphi_ll("dphi_ll", [](const Baby &b) -> NamedFunc::ScalarType{
      double phi1, eta1, phi2, eta2;
      GetAngles(b, phi1, eta1, phi2, eta2);
      if(phi1 == -999 || phi2 == -999) return -1.;
      return fabs(TVector2::Phi_mpi_pi(phi2-phi1));
    });
  NamedFunc min_dphi_lmet("min_dphi_lmet", MinDeltaPhiLepMet);
  NamedFunc max_dphi_lmet("max_dphi_lmet", MaxDeltaPhiLepMet);
  NamedFunc min_dphi_lj("min_dphi_lj", MinDeltaPhiLepJet);
  NamedFunc max_dphi_lj("max_dphi_lj", MaxDeltaPhiLepJet);
  NamedFunc min_dphi_metj("min_dphi_metj", MinDeltaPhiMetJet);
  NamedFunc max_dphi_metj("max_dphi_metj", MaxDeltaPhiMetJet);

  NamedFunc dr_ll("dr_ll", [](const Baby &b) -> NamedFunc::ScalarType{
      double phi1, eta1, phi2, eta2;
      GetAngles(b, phi1, eta1, phi2, eta2);
      if(phi1 == -999 || phi2 == -999) return -1.;
      return hypot(TVector2::Phi_mpi_pi(phi2-phi1), eta2-eta1);
    });
  NamedFunc min_dr_lj("min_dr_lj", MinDeltaRLepJet);
  NamedFunc max_dr_lj("max_dr_lj", MaxDeltaRLepJet);

  string selections[] = {"nleps==2&&ht<200&&njets<=2&&nbl==0&&met<500&&leps_pt[0]>120",
                         "nleps==2&&ht>200&&njets<=2&&nbl==0&&met<500&&leps_pt[0]>120",
                         "nleps==2&&ht<200&&njets<=2&&nbl==0&&met<500&&leps_pt[0]>120&&met>200",
                         "nleps==2&&ht>200&&njets<=2&&nbl==0&&met<500&&leps_pt[0]>120&&met>200",
                         /*"nleps==2&&ht>200&&njets<=2&&nbl==0&&met<500&&(elel_m+mumu_m+999)>200&&leps_pt[0]>120"*/};

  string leps[] = {"&&nels==1&&nmus==1","&&nels==2&&elel_m>60","&&nmus==2&&mumu_m>60"/*,"&&nmus==2&&mumu_m>200","&&nels==2&&elel_m>200","&&nels==1&&nmus==1"*/};

  pm.Push<HistoStack>(HistoDef(50,0,500, "elel_m", "m_{el el} [GeV]", "nels==2&&ht>200&&njets<=2&&nbl==0&&met<500&&leps_pt[0]>120", "weight", {60}), full_trig_skim, norms);
  pm.Push<HistoStack>(HistoDef(20,0,1000, "mumu_m", "m_{mu mu} [GeV]", "nmus==2&&ht>200&&njets<=2&&nbl==0&&met<500&&leps_pt[0]>120", "weight", {60}), full_trig_skim, norms);

  for(unsigned int ilep=0;ilep<1;ilep++){
    pm.Push<HistoStack>(HistoDef(30,0,600, "leps_pt[0]", "Leading lepton p_{T}","nleps==2&&ht>200&&nbl==0&&njets<=2&&met<500"+leps[ilep], "weight", {120}), full_trig_skim, norms);
    pm.Push<HistoStack>(HistoDef(20,0,1000, "ht", "H_{T}","nleps==2&&nbl==0&&njets<=2&&met<500&&leps_pt[0]>120"+leps[ilep], "weight", {200}), full_trig_skim, norms);
    pm.Push<HistoStack>(HistoDef(5,-0.5,4.5, "nbl", "N_{b, loose}","nleps==2&&ht>200&&njets<=2&&met<500&&leps_pt[0]>120"+leps[ilep], "weight", {0.5}), full_trig_skim, norms);
    pm.Push<HistoStack>(HistoDef(9,-0.5,8.5, "njets", "N_{jets}","nleps==2&&ht>200&&nbl==0&&met<500&&leps_pt[0]>120"+leps[ilep], "weight", {2.5}), full_trig_skim, norms);
  }

  pm.Push<HistoStack>(HistoDef(50,0,500, "elel_m", "m_{el el} [GeV]", "nels==2&&ht<200&&njets<=2&&nbl==0&&met<500&&leps_pt[0]>120", "weight", {60}), full_trig_skim, norms);
  pm.Push<HistoStack>(HistoDef(20,0,1000, "mumu_m", "m_{mu mu} [GeV]", "nmus==2&&ht<200&&njets<=2&&nbl==0&&met<500&&leps_pt[0]>120", "weight", {60}), full_trig_skim, norms);

  for(unsigned int ilep=0;ilep<2;ilep++){
    pm.Push<HistoStack>(HistoDef(30,0,600, "leps_pt[0]", "Leading lepton p_{T}","nleps==2&&ht<200&&nbl==0&&njets<=2&&met<500"+leps[ilep], "weight", {120}), full_trig_skim, norms);
    pm.Push<HistoStack>(HistoDef(20,0,1000, "ht", "H_{T}","nleps==2&&nbl==0&&njets<=2&&met<500&&leps_pt[0]>120"+leps[ilep], "weight", {200}), full_trig_skim, norms);
    pm.Push<HistoStack>(HistoDef(5,-0.5,4.5, "nbl", "N_{b, loose}","nleps==2&&ht<200&&njets<=2&&met<500&&leps_pt[0]>120"+leps[ilep], "weight", {0.5}), full_trig_skim, norms);
    pm.Push<HistoStack>(HistoDef(9,-0.5,8.5, "njets", "N_{jets}","nleps==2&&ht<200&&nbl==0&&met<500&&leps_pt[0]>120"+leps[ilep], "weight", {2.5}), full_trig_skim, norms);
  }

  for(unsigned int isel=0;isel<2;isel++){
    for(unsigned int ilep=0;ilep<2;ilep++){
      pm.Push<HistoStack>(HistoDef(8,0,800, "mj14", "M_{J} [GeV]", selections[isel]+leps[ilep], "weight", {250,400}), full_trig_skim, norms);
      pm.Push<HistoStack>(HistoDef(10,0,500, "met", "MET [GeV]", selections[isel]+leps[ilep], "weight", {350}), full_trig_skim, norms);
      pm.Push<HistoStack>(HistoDef(20,0,1000, "leps_pt[0]", "Leading lepton p_{T} [GeV]", selections[isel]+leps[ilep], "weight", {20}), full_trig_skim, norms);
      pm.Push<HistoStack>(HistoDef(25,0,500, "mt", "m_{T} [GeV]", selections[isel]+leps[ilep], "weight", {140}), full_trig_skim, norms);
      pm.Push<HistoStack>(HistoDef(20, 0., pi, min_dphi_lmet, "Min. #Delta#phi(l,MET)", selections[isel]+leps[ilep]),  full_trig_skim, norms);
      pm.Push<HistoStack>(HistoDef(20, 0., pi, max_dphi_lmet, "Max. #Delta#phi(l,MET)", selections[isel]+leps[ilep]),  full_trig_skim, norms);
      //pm.Push<HistoStack>(HistoDef(20,0,40, "npv", "N_{PV}", selections[isel]+leps[ilep], "weight", {-1}), full_trig_skim, norms);
      if(ilep==1) pm.Push<HistoStack>(HistoDef(20,0,1000, "elel_m", "m_{el el} [GeV]", selections[isel]+leps[ilep], "weight", {90}), full_trig_skim, norms);
      if(ilep==0) pm.Push<HistoStack>(HistoDef(20,0,1000, "elmu_m", "m_{el mu} [GeV]", selections[isel]+leps[ilep], "weight", {90}), full_trig_skim, norms);
      if(ilep==2) pm.Push<HistoStack>(HistoDef(20,0,1000, "mumu_m", "m_{mu mu} [GeV]", selections[isel]+leps[ilep], "weight", {90}), full_trig_skim, norms);
    }
  }

  PlotMaker pm_data;
  for(unsigned int isel=0;isel<2;isel++){
    pm_data.Push<HistoStack>(HistoDef(10,0,500, "mj14", "MJ 1.4 with leptons [GeV]", selections[isel], "weight", {350}), ee_vs_mumu, norms);
    pm_data.Push<HistoStack>(HistoDef(10,0,500, "met", "MET [GeV]", selections[isel], "weight", {350}), ee_vs_mumu, norms);
    pm_data.Push<HistoStack>(HistoDef(20,0,1000, "leps_pt[0]", "Leading lepton p_{T} [GeV]", selections[isel], "weight", {20}), ee_vs_mumu, norms);
    pm_data.Push<HistoStack>(HistoDef(25,0,500, "mt", "m_{T} [GeV]", selections[isel], "weight", {140}), ee_vs_mumu, norms);
    pm_data.Push<HistoStack>(HistoDef(20, 0., pi, min_dphi_lmet, "Min. #Delta#phi(l,MET)", selections[isel]), ee_vs_mumu, norms);
    pm_data.Push<HistoStack>(HistoDef(20, 0., pi, max_dphi_lmet, "Max. #Delta#phi(l,MET)", selections[isel]), ee_vs_mumu, norms);
    // if(ilep==0)
    pm_data.Push<HistoStack>(HistoDef(20,0,1000, "elel_m+mumu_m+999", "m_{lep lep} [GeV]", selections[isel], "weight", {90}), ee_vs_mumu, norms);
    // if(ilep==1) pm.Push<HistoStack>(HistoDef(20,0,1000, "elmu_m", "m_{el mu} [GeV]", selections[isel]+leps[ilep], "weight", {90}), ee_vs_mumu, shapes);
    // if(ilep==2) pm.Push<HistoStack>(HistoDef(20,0,1000, "mumu_m", "m_{mu mu} [GeV]", selections[isel]+leps[ilep], "weight", {90}), ee_vs_mumu, shapes);

  }

  if(single_thread) pm.multithreaded_ = false;
  pm.MakePlots(lumi);
  cout<<lumi<<endl;
  if(single_thread) pm_data.multithreaded_ = false;
  pm_data.MakePlots(1.0);
}

bool IsGoodMuon(const Baby &b, size_t imu){
  return imu<b.mus_pt()->size()
    && b.mus_pt()->at(imu)>20.
    && fabs(b.mus_eta()->at(imu))<2.4
                                  && b.mus_sigid()->at(imu)
    && b.mus_miniso()->at(imu) >= 0.
    && b.mus_miniso()->at(imu) < 0.2;
}

bool IsGoodElectron(const Baby &b, size_t iel){
  return iel<b.els_pt()->size()
    && b.els_pt()->at(iel)>20.
    && fabs(b.els_sceta()->at(iel))<2.5
                                    && b.els_sigid()->at(iel)
    && b.els_miniso()->at(iel) >= 0.
    && b.els_miniso()->at(iel) < 0.1;
}

bool IsGoodTrack(const Baby &b, size_t itk){
  if(itk >= b.tks_pt()->size()) return false;

  if(fabs(b.tks_pdg()->at(itk))==211  && b.tks_pt()->at(itk)>15. && b.tks_miniso()->at(itk)<0.1 && b.tks_mt2()->at(itk)<60 && b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05 ){
    return true;
  }else if(fabs(b.tks_pdg()->at(itk))==13 && b.tks_pt()->at(itk)>10. && b.tks_miniso()->at(itk)<0.2 && b.tks_mt2()->at(itk)<80 && b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05){
    return true;
  }else if(fabs(b.tks_pdg()->at(itk))==11 && b.tks_pt()->at(itk)>10. && b.tks_miniso()->at(itk)<0.2 && b.tks_mt2()->at(itk)<80 && b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05){
    return true;
  }else{
    return false;
  }
}

bool IsGoodJet(const Baby &b, size_t ijet){
  if(ijet >= b.jets_pt()->size()) return false;

  return b.jets_pt()->at(ijet) > 30.
    && b.jets_eta()->at(ijet) < 2.
                                && !b.jets_islep()->at(ijet);
}

void GetAngles(const Baby &b,
               double &phi1, double &eta1,
               double &phi2, double &eta2){
  phi1 = -999.; eta1 = -999.;
  phi2 = -999.; eta2 = -999.;
  bool h1=false, h2=false;
  if(b.nels()==2 && b.nmus()==0){
    for(size_t iel = 0; iel < b.els_pt()->size() && !(h1&&h2); ++iel){
      if(IsGoodElectron(b, iel)){
        if(!h1){
          phi1 = b.els_phi()->at(iel);
          eta1 = b.els_sceta()->at(iel);
          h1 = true;
        }else if(!h2){
          phi2 = b.els_phi()->at(iel);
          eta2 = b.els_sceta()->at(iel);
          h2 = true;
        }
      }
    }
  }else if(b.nels()==1 && b.nmus()==1){
    for(size_t iel = 0; iel < b.els_pt()->size() && !h1; ++iel){
      if(IsGoodElectron(b, iel)){
        if(!h1){
          phi1 = b.els_phi()->at(iel);
          eta1 = b.els_sceta()->at(iel);
          h1 = true;
        }else if(!h2){
          phi2 = b.els_phi()->at(iel);
          eta2 = b.els_sceta()->at(iel);
          h2 = true;
        }
      }
    }
    for(size_t imu = 0; imu < b.mus_pt()->size() && h1 && !h2; ++imu){
      if(IsGoodMuon(b, imu)){
        if(!h1){
          phi1 = b.mus_phi()->at(imu);
          eta1 = b.mus_eta()->at(imu);
          h1 = true;
        }else if(!h2){
          phi2 = b.mus_phi()->at(imu);
          eta2 = b.mus_eta()->at(imu);
          h2 = true;
        }
      }
    }
  }else if(b.nels()==0 && b.nmus()==2){
    for(size_t imu = 0; imu < b.mus_pt()->size() && !(h1&&h2); ++imu){
      if(IsGoodMuon(b, imu)){
        if(!h1){
          phi1 = b.mus_phi()->at(imu);
          eta1 = b.mus_eta()->at(imu);
          h1 = true;
        }else if(!h2){
          phi2 = b.mus_phi()->at(imu);
          eta2 = b.mus_eta()->at(imu);
          h2 = true;
        }
      }
    }
  }else if(b.nels()==1 && b.nmus()==0 && b.nveto()==1){
    for(size_t iel = 0; iel < b.els_pt()->size() && !h1; ++iel){
      if(IsGoodElectron(b, iel)){
        if(!h1){
          phi1 = b.els_phi()->at(iel);
          eta1 = b.els_sceta()->at(iel);
          h1 = true;
        }else if(!h2){
          phi2 = b.els_phi()->at(iel);
          eta2 = b.els_sceta()->at(iel);
          h2 = true;
        }
      }
    }
    for(size_t itk = 0; itk < b.tks_pt()->size() && h1 && !h2; ++itk){
      if(IsGoodTrack(b, itk)){
        if(!h1){
          phi1 = b.tks_phi()->at(itk);
          eta1 = b.tks_eta()->at(itk);
          h1 = true;
        }else if(!h2){
          phi2 = b.tks_phi()->at(itk);
          eta2 = b.tks_eta()->at(itk);
          h2 = true;
        }
      }
    }
  }else if(b.nels()==0 && b.nmus()==1 && b.nveto()==1){
    for(size_t imu = 0; imu < b.mus_pt()->size() && !h1; ++imu){
      if(IsGoodMuon(b, imu)){
        if(!h1){
          phi1 = b.mus_phi()->at(imu);
          eta1 = b.mus_eta()->at(imu);
          h1 = true;
        }else if(!h2){
          phi2 = b.mus_phi()->at(imu);
          eta2 = b.mus_eta()->at(imu);
          h2 = true;
        }
      }
    }
    for(size_t itk = 0; itk < b.tks_pt()->size() && h1 && !h2; ++itk){
      if(IsGoodTrack(b, itk)){
        if(!h1){
          phi1 = b.tks_phi()->at(itk);
          eta1 = b.tks_eta()->at(itk);
          h1 = true;
        }else if(!h2){
          phi2 = b.tks_phi()->at(itk);
          eta2 = b.tks_eta()->at(itk);
          h2 = true;
        }
      }
    }
  }
}

NamedFunc::ScalarType MinDeltaPhiLepMet(const Baby &b){
  double phi1, eta1, phi2, eta2;
  GetAngles(b, phi1, eta1, phi2, eta2);
  double dphi1 = fabs(TVector2::Phi_mpi_pi(phi1-b.met_phi()));
  double dphi2 = fabs(TVector2::Phi_mpi_pi(phi2-b.met_phi()));
  if(phi1 != -999 && phi2 != -999){
    return std::min(dphi1,dphi2);
  }else if(phi1 != -999 && phi2 == -999){
    return dphi1;
  }else if(phi1 == -999 && phi2 != -999){
    return dphi2;
  }
  return -1;
}

NamedFunc::ScalarType MaxDeltaPhiLepMet(const Baby &b){
  double phi1, eta1, phi2, eta2;
  GetAngles(b, phi1, eta1, phi2, eta2);
  double dphi1 = fabs(TVector2::Phi_mpi_pi(phi1-b.met_phi()));
  double dphi2 = fabs(TVector2::Phi_mpi_pi(phi2-b.met_phi()));
  if(phi1 != -999 && phi2 != -999){
    return std::max(dphi1,dphi2);
  }else if(phi1 != -999 && phi2 == -999){
    return dphi1;
  }else if(phi1 == -999 && phi2 != -999){
    return dphi2;
  }
  return -1;
}

NamedFunc::ScalarType MinDeltaPhiLepJet(const Baby &b){
  double phi1, eta1, phi2, eta2;
  GetAngles(b, phi1, eta1, phi2, eta2);
  double minphi = -1.;
  for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
    if(!IsGoodJet(b,ijet)) continue;
    double dphi1 = fabs(TVector2::Phi_mpi_pi(phi1-b.jets_phi()->at(ijet)));
    double dphi2 = fabs(TVector2::Phi_mpi_pi(phi2-b.jets_phi()->at(ijet)));
    double thisdphi = -1;
    if(phi1 != -999 && phi2 != -999){
      thisdphi = std::min(dphi1, dphi2);
    }else if(phi1 != -999 && phi2 == -999){
      thisdphi = dphi1;
    }else if(phi1 == -999 && phi2 != -999){
      thisdphi = dphi2;
    }
    if(minphi < 0. || thisdphi < minphi){
      minphi = thisdphi;
    }
  }
  return minphi;
}

NamedFunc::ScalarType MaxDeltaPhiLepJet(const Baby &b){
  double phi1, eta1, phi2, eta2;
  GetAngles(b, phi1, eta1, phi2, eta2);
  double maxphi = -1.;
  for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
    if(!IsGoodJet(b,ijet)) continue;
    double dphi1 = fabs(TVector2::Phi_mpi_pi(phi1-b.jets_phi()->at(ijet)));
    double dphi2 = fabs(TVector2::Phi_mpi_pi(phi2-b.jets_phi()->at(ijet)));
    double thisdphi = -1;
    if(phi1 != -999 && phi2 != -999){
      thisdphi = std::max(dphi1, dphi2);
    }else if(phi1 != -999 && phi2 == -999){
      thisdphi = dphi1;
    }else if(phi1 == -999 && phi2 != -999){
      thisdphi = dphi2;
    }
    if(maxphi < 0. || thisdphi > maxphi){
      maxphi = thisdphi;
    }
  }
  return maxphi;
}

NamedFunc::ScalarType MinDeltaPhiMetJet(const Baby &b){
  double minphi = -1.;
  for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
    if(!IsGoodJet(b,ijet)) continue;
    double thisdphi = fabs(TVector2::Phi_mpi_pi(b.met_phi()-b.jets_phi()->at(ijet)));
    if(minphi < 0. || thisdphi < minphi){
      minphi = thisdphi;
    }
  }
  return minphi;
}

NamedFunc::ScalarType MaxDeltaPhiMetJet(const Baby &b){
  double maxphi = -1.;
  for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
    if(!IsGoodJet(b,ijet)) continue;
    double thisdphi = fabs(TVector2::Phi_mpi_pi(b.met_phi()-b.jets_phi()->at(ijet)));
    if(maxphi < 0. || thisdphi > maxphi){
      maxphi = thisdphi;
    }
  }
  return maxphi;
}

NamedFunc::ScalarType MinDeltaRLepJet(const Baby &b){
  double phi1, eta1, phi2, eta2;
  GetAngles(b, phi1, eta1, phi2, eta2);
  double minr = -1.;
  for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
    if(!IsGoodJet(b,ijet)) continue;
    double dr1 = hypot(TVector2::Phi_mpi_pi(phi1-b.jets_phi()->at(ijet)), eta2-eta1);
    double dr2 = hypot(TVector2::Phi_mpi_pi(phi2-b.jets_phi()->at(ijet)), eta2-eta1);
    double thisdr = -1;
    if(phi1 != -999 && phi2 != -999){
      thisdr = std::min(dr1, dr2);
    }else if(phi1 != -999 && phi2 == -999){
      thisdr = dr1;
    }else if(phi1 == -999 && phi2 != -999){
      thisdr = dr2;
    }
    if(minr < 0. || thisdr < minr){
      minr = thisdr;
    }
  }
  return minr;
}

NamedFunc::ScalarType MaxDeltaRLepJet(const Baby &b){
  double phi1, eta1, phi2, eta2;
  GetAngles(b, phi1, eta1, phi2, eta2);
  double maxr = -1.;
  for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
    if(!IsGoodJet(b,ijet)) continue;
    double dr1 = hypot(TVector2::Phi_mpi_pi(phi1-b.jets_phi()->at(ijet)), eta2-eta1);
    double dr2 = hypot(TVector2::Phi_mpi_pi(phi2-b.jets_phi()->at(ijet)), eta2-eta1);
    double thisdr = -1;
    if(phi1 != -999 && phi2 != -999){
      thisdr = std::max(dr1, dr2);
    }else if(phi1 != -999 && phi2 == -999){
      thisdr = dr1;
    }else if(phi1 == -999 && phi2 != -999){
      thisdr = dr2;
    }
    if(maxr < 0. || thisdr > maxr){
      maxr = thisdr;
    }
  }
  return maxr;
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
    case 0:
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
