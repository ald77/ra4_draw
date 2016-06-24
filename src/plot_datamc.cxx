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

  double lumi = 2.6;
  
  // 80X (ttbar, qcd, dy[ht=100-600], wjets) + 74X (rest)
  string trig_mc        = "/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/unskimmed/";
  string trig_skim1l_mc = "/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/merged_standard/";
  string trig_skim0l_mc = "/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/merged_qcd/";     
  string trig_skim2l_mc = "/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/merged_dy_ht300/"; 

  Palette colors("txt/colors.txt", "default");

  auto tt1l = Proc<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {trig_skim1l_mc+"*_TTJets*Lept*.root", trig_skim1l_mc+"*_TTJets_HT*.root"}, "ntruleps<=1&&stitch");
  auto tt2l = Proc<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {trig_skim1l_mc+"*_TTJets*Lept*.root", trig_skim1l_mc+"*_TTJets_HT*.root"}, "ntruleps>=2&&stitch");
  auto wjets = Proc<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {trig_skim1l_mc+"*_WJetsToLNu*.root"});
  auto single_t = Proc<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {trig_skim1l_mc+"*_ST_*.root"});
  auto ttv = Proc<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {trig_skim1l_mc+"*_TTWJets*.root", trig_skim1l_mc+"*_TTZTo*.root"});
  auto other = Proc<Baby_full>("Other", Process::Type::background, colors("other"),
    {trig_skim1l_mc+"*DYJetsToLL*.root", trig_skim1l_mc+"*_QCD_HT*.root",
     trig_skim1l_mc+"*_ZJet*.root", trig_skim1l_mc+"*_WWTo*.root",
     trig_skim1l_mc+"*ggZH_HToBB*.root", trig_skim1l_mc+"*ttHJetTobb*.root",
     trig_skim1l_mc+"*_TTGJets*.root", trig_skim1l_mc+"*_TTTT_*.root",
     trig_skim1l_mc+"*_WH_HToBB*.root", trig_skim1l_mc+"*_WZTo*.root",
     trig_skim1l_mc+"*_ZH_HToBB*.root", trig_skim1l_mc+"_ZZ_*.root"});

  auto t1tttt_nc = Proc<Baby_full>("T1tttt(1500,100)", Process::Type::signal, colors("t1tttt"),
    {trig_mc+"*SMS-T1tttt_mGluino-1500_mLSP-100*.root"});
  auto t1tttt_c = Proc<Baby_full>("T1tttt(1200,800)", Process::Type::signal, colors("t1tttt"),
    {trig_mc+"*SMS-T1tttt_mGluino-1200_mLSP-800*.root"});
  t1tttt_c->SetLineStyle(2);

  auto data_1l = Proc<Baby_full>("Data", Process::Type::data, kBlack,
    {"/net/cms2/cms2r0/babymaker/babies/2016_06_21/data/skim_standard/*.root"},"pass&&(trig[4]||trig[8]||trig[13]||trig[33])");
    vector<shared_ptr<Process> > full_trig_skim_1l = {data_1l, /*t1tttt_nc, t1tttt_c,*/ tt1l, tt2l, wjets, single_t, ttv, other}; 
  
  // 
  // 0-lepton plots
  //  
  auto tt_0l = Proc<Baby_full>("t#bar{t}", Process::Type::background, colors("tt_1l"),
    {trig_skim0l_mc+"*_TTJets_TuneCUETP8M1_13TeV*.root"});
  auto qcd_0l = Proc<Baby_full>("QCD", Process::Type::background, colors("qcd"),
    {trig_skim0l_mc+"*_QCD_HT700to1000*.root",trig_skim0l_mc+"*_QCD_HT1000to1500*.root",
     trig_skim0l_mc+"*_QCD_HT1500to2000*.root",trig_skim0l_mc+"*_QCD_HT2000toInf*.root"});
  auto other_0l = Proc<Baby_full>("Other", Process::Type::background, colors("ttv"),
    {trig_skim0l_mc+"*_TTWJet*root",  trig_skim0l_mc+"*_TTZTo*.root",
     trig_skim0l_mc+"*_ZJet*.root", trig_skim0l_mc+"*_WWTo*.root",
     trig_skim0l_mc+"*ggZH_HToBB*.root", trig_skim0l_mc+"*ttHJetTobb*.root",
     trig_skim0l_mc+"*_TTGJets*.root", trig_skim0l_mc+"*_TTTT_*.root",
     trig_skim0l_mc+"*_WH_HToBB*.root", trig_skim0l_mc+"*_WZTo*.root",
     trig_skim0l_mc+"*_ST_*.root", trig_skim0l_mc+"*_WJetsToLNu*.root",
     trig_skim0l_mc+"*_ZH_HToBB*.root", trig_skim0l_mc+"_ZZ_*.root", 
     trig_skim0l_mc+"*_DYJetsToLL*HT*.root"});

  auto data_0l = Proc<Baby_full>("Data", Process::Type::data, kBlack,
    {"/net/cms2/cms2r0/babymaker/babies/2016_06_21/data/skim_ht900/*.root"},"pass&&trig[12]");
  vector<shared_ptr<Process> > full_trig_skim_0l = {data_0l, /*t1tttt_nc, t1tttt_c,*/ qcd_0l, tt_0l, other_0l};

  //
  // 2-lepton plots
  //
  auto tt1l_2l = Proc<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {trig_skim2l_mc+"*_TTJets*Lept*.root"},
    "ntruleps<=1");
  auto tt2l_2l = Proc<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {trig_skim2l_mc+"*_TTJets*Lept*.root"},
    "ntruleps>=2");
  auto dy_2l = Proc<Baby_full>("DY", Process::Type::background, colors("wjets"),
    {trig_skim2l_mc+"*_DYJetsToLL*HT*.root"});
  auto other_2l = Proc<Baby_full>("Other", Process::Type::background, colors("ttv"),
    {trig_skim2l_mc+"*_TTWJet*root",  trig_skim2l_mc+"*_TTZTo*.root",
     trig_skim2l_mc+"*_ZJet*.root", trig_skim2l_mc+"*_WWTo*.root",
     trig_skim2l_mc+"*ggZH_HToBB*.root", trig_skim2l_mc+"*ttHJetTobb*.root",
     trig_skim2l_mc+"*_TTGJets*.root", trig_skim2l_mc+"*_TTTT_*.root",
     trig_skim2l_mc+"*_WH_HToBB*.root", trig_skim2l_mc+"*_WZTo*.root",
     trig_skim2l_mc+"*_ST_*.root", trig_skim2l_mc+"*_WJetsToLNu*.root",
     trig_skim2l_mc+"*_ZH_HToBB*.root", trig_skim2l_mc+"_ZZ_*.root", 
     trig_skim2l_mc+"*_QCD_HT700to1000*.root",trig_skim2l_mc+"*_QCD_HT1000to1500*.root",
     trig_skim2l_mc+"*_QCD_HT1500to2000*.root",trig_skim2l_mc+"*_QCD_HT2000toInf*.root"});

  auto data_2l = Proc<Baby_full>("Data", Process::Type::data, kBlack,
    {"/net/cms2/cms2r0/babymaker/babies/2016_06_21/data/skim_nm1nj2/*.root"},"pass&&(trig[20]||trig[22])");
  vector<shared_ptr<Process> > full_trig_skim_2l = {data_2l, /*t1tttt_nc, t1tttt_c,*/ tt1l_2l, tt2l_2l, dy_2l, other_2l};

  //
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
  vector<PlotOpt> all_plot_types = {log_lumi_info, lin_lumi_info};


  string lowmtcut = "ht>500&&met>200&&mt<140";
  string lowmjcut = "ht>500&&met>200&&mj14<400";
  PlotMaker pm;  

  for(int ilep=0; ilep<3; ilep++) 
  {
      string lepcut;
      if(ilep==0)  lepcut="nels==1&&nmus==0"; 
      if(ilep==1)  lepcut="nels==0&&nmus==1"; 
      if(ilep==2)  lepcut="nleps==1"; 
  
      for(int inveto=0; inveto<2; inveto++)  
      { 
          string nvetocut;
          if(inveto==0)  nvetocut="nveto==0"; 
          if(inveto==1)  nvetocut="1"; 
          if(inveto==1)  continue;

          for(int injets=0; injets<2; injets++)  
          { 
              string njetscut;
              if(injets==0)  njetscut="njets>=5"; 
              if(injets==1)  njetscut="njets>=6"; 
              if(injets==0)  continue;
              
              for(int inb=0; inb<2; inb++)  
              { 
                  string nbcut;
                  if(inb==0)  nbcut="nbm>=1"; 
                  if(inb==1)  nbcut="nbm>=2"; 
                  if(inb==1)  continue;

                  //
                  // event-level variables except ones related to fatjet
                  // 
                  pm.Push<HistoStack>(HistoDef(15, 500, 2000., "ht", "H_{T} [GeV]",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {500.}),
                          full_trig_skim_1l, all_plot_types);
                  pm.Push<HistoStack>(HistoDef(10, 200, 700., "met", "E_{T}^{miss} [GeV]",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {200., 350., 500.}),
                          full_trig_skim_1l, all_plot_types);
                  pm.Push<HistoStack>(HistoDef(10, 200, 700., "met", "E_{T}^{miss} [GeV]",
                              lepcut+"&&"+lowmjcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {200., 350., 500.}),
                          full_trig_skim_1l, all_plot_types); 
                  pm.Push<HistoStack>(HistoDef(16, -0.5, 15.5, "njets", "N_{jets}",
                              lepcut+"&&"+lowmtcut+"&&"+nbcut+"&&"+nvetocut, "weight", {5.5, 8.5}),
                          full_trig_skim_1l, all_plot_types);
                  pm.Push<HistoStack>(HistoDef(7, -0.5, 6.5, "nbm", "N_{b}",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nvetocut, "weight", {0.5, 1.5, 2.5}),
                          full_trig_skim_1l, all_plot_types);
                  pm.Push<HistoStack>(HistoDef(14, 0., 280., "mt", "m_{T} [GeV]",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {140.}),
                          full_trig_skim_1l, all_plot_types); 
                  pm.Push<HistoStack>(HistoDef(14, 0., 280., "mt", "m_{T} [GeV]",
                              lepcut+"&&"+lowmjcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {140.}),
                          full_trig_skim_1l, all_plot_types); 
                  pm.Push<HistoStack>(HistoDef(14, 0., 280., "mt", "m_{T} [GeV]",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&nbm==0"+"&&"+nvetocut, "weight", {140.}),
                          full_trig_skim_1l, all_plot_types);


                  if(ilep==0) 
                  {
                      // only barrel or endcap electron 
                      pm.Push<HistoStack>(HistoDef(10, 200, 700., "met", "E_{T}^{miss} [GeV]",
                                  lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut+"&&(leps_eta[0]<1.479&&leps_eta[0]>-1.479)", "weight", {200., 350., 500.}),
                              full_trig_skim_1l, all_plot_types);
                      pm.Push<HistoStack>(HistoDef(10, 200, 700., "met", "E_{T}^{miss} [GeV]",
                                  lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut+"&&(leps_eta[0]>1.479||leps_eta[0]<-1.479)", "weight", {200., 350., 500.}),
                              full_trig_skim_1l, all_plot_types);
                  }

                  //
                  // leptons
                  //
                  if(ilep==0) // electron 
                  {
                      pm.Push<HistoStack>(HistoDef(10, 0, 200., "leps_pt[0]", "p_{T}(electron) [GeV]",
                                  lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {999.}),
                              full_trig_skim_1l, all_plot_types);
                      pm.Push<HistoStack>(HistoDef(10, -2.5, 2.5, "leps_eta[0]", "#eta(electron) [GeV]",
                                  lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {999.}),
                              full_trig_skim_1l, all_plot_types);
                  } 
                  else if(ilep==1)  // muon
                  { 
                      pm.Push<HistoStack>(HistoDef(10, 0, 200., "leps_pt[0]", "p_{T}(muon) [GeV]",
                                  lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {999.}),
                              full_trig_skim_1l, all_plot_types);
                      pm.Push<HistoStack>(HistoDef(10, -2.5, 2.5, "leps_eta[0]", "#eta(muon) [GeV]",
                                  lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {999.}),
                              full_trig_skim_1l, all_plot_types);
                  } 
                  else if(ilep==2)  // electron or muon
                  { 
                      pm.Push<HistoStack>(HistoDef(10, 0, 200., "leps_pt[0]", "p_{T}(lepton) [GeV]",
                                  lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {999.}),
                              full_trig_skim_1l, all_plot_types);
                      pm.Push<HistoStack>(HistoDef(10, -2.5, 2.5, "leps_eta[0]", "#eta(lepton) [GeV]",
                                  lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {999.}),
                              full_trig_skim_1l, all_plot_types);
                  } 

                  //
                  // jets
                  //
                  pm.Push<HistoStack>(HistoDef(10, 30, 630., "jets_pt[0]", "p_{T}(leading jet) [GeV]",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut+"&&(!jets_islep[0])", "weight", {999.}),
                          full_trig_skim_1l, all_plot_types);
                  pm.Push<HistoStack>(HistoDef(10, -2.5, 2.5, "jets_eta[0]", "#eta(leading jet) [GeV]",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut+"&&(!jets_islep[0])", "weight", {999.}),
                          full_trig_skim_1l, all_plot_types);
                  pm.Push<HistoStack>(HistoDef(10, 0, 200., "jets_m[0]", "mass(leading jet) [GeV]",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut+"&&(!jets_islep[0])", "weight", {999.}),
                          full_trig_skim_1l, all_plot_types);
                  pm.Push<HistoStack>(HistoDef(10, 30, 630., "jets_pt", "p_{T}(jet) [GeV]",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut+"&&(!jets_islep)", "weight", {999.}),
                          full_trig_skim_1l, all_plot_types);
                  pm.Push<HistoStack>(HistoDef(10, -2.5, 2.5, "jets_eta", "#eta(jet) [GeV]",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut+"&&(!jets_islep)", "weight", {999.}),
                          full_trig_skim_1l, all_plot_types);
                  pm.Push<HistoStack>(HistoDef(10, 0, 200., "jets_m", "m(jet) [GeV]",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut+"&&(!jets_islep)", "weight", {999.}),
                          full_trig_skim_1l, all_plot_types);

                  //
                  // fatjets
                  //
                  pm.Push<HistoStack>(HistoDef(10, 0, 1000., "fjets14_pt[0]", "p_{T}(J1) [GeV]",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {999.}),
                          full_trig_skim_1l, all_plot_types); 
                  pm.Push<HistoStack>(HistoDef(10, 0, 500., "fjets14_m[0]", "m_{J1} [GeV]",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {999.}),
                          full_trig_skim_1l, all_plot_types); 
                  pm.Push<HistoStack>(HistoDef(10, 0, 500., "fjets14_m[0]", "m_{J1} [GeV]",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut+"&&leps_pt[0]>160", "weight", {999.}),
                          full_trig_skim_1l, all_plot_types); 
                  pm.Push<HistoStack>(HistoDef(10, 0, 500., "fjets14_m[0]", "m_{J1} [GeV]",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut+"&&leps_pt[0]<=160", "weight", {999.}),
                          full_trig_skim_1l, all_plot_types); 
                  pm.Push<HistoStack>(HistoDef(10, -2.5, 2.5, "fjets14_eta[0]", "#eta(J1)",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {999.}),
                          full_trig_skim_1l, all_plot_types);
                  pm.Push<HistoStack>(HistoDef(7, 0.5, 7.5, "fjets14_nconst[0]", "N_{constituents}(J1)",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {999.}),
                          full_trig_skim_1l, all_plot_types);
                  pm.Push<HistoStack>(HistoDef(10, 0, 1000., "fjets14_pt", "p_{T}(J) [GeV]",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {999.}),
                          full_trig_skim_1l, all_plot_types); 
                  pm.Push<HistoStack>(HistoDef(10, 0, 500., "fjets14_m", "m_{J} [GeV]",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {999.}),
                          full_trig_skim_1l, all_plot_types);
                  pm.Push<HistoStack>(HistoDef(10, -2.5, 2.5, "fjets14_eta", "#eta(J)",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {999.}),
                          full_trig_skim_1l, all_plot_types);
                  pm.Push<HistoStack>(HistoDef(7, 0.5, 7.5, "fjets14_nconst", "N_{constituents}(J)",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {999.}),
                          full_trig_skim_1l, all_plot_types);

                  pm.Push<HistoStack>(HistoDef(20, 0, 400., "fjets14_m[0]", "m_{J1} [GeV]",
                              lepcut+"&&"+lowmtcut+"&&njets<=5&&nbm>=1", "weight", {999.}),
                          full_trig_skim_1l, all_plot_types); 

                  // event-level variables related to fatjet
                  pm.Push<HistoStack>(HistoDef(15, 0., 1500., "mj14", "M_{J} [GeV]",
  //
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {400.}),
                          full_trig_skim_1l, all_plot_types);
                  pm.Push<HistoStack>(HistoDef(7, 0.5, 7.5, "nfjets14", "N_{J} [GeV]",
                              lepcut+"&&"+lowmtcut+"&&"+njetscut+"&&"+nbcut+"&&"+nvetocut, "weight", {999.}),
                          full_trig_skim_1l, all_plot_types);

              } //for(int inb=0; inb<2; inb++)  
          } //for(int injets=0; injets<2; injets++)  
     } //for(int inveto=0; inveto<2; inveto++)  
  }  


  //
  // 0-lepton events
  //
  pm.Push<HistoStack>(HistoDef(30, 0., 1500., "mj14", "M_{J} [GeV]",
                      "nleps==0&&met<50&&ht>1000&&njets>=4&&njets<=5", "weight", {400.}),
                       full_trig_skim_0l, all_plot_types);
  pm.Push<HistoStack>(HistoDef(30, 0., 500., "fjets14_m", "m_{J} [GeV] [GeV]",
                      "nleps==0&&met<50&&ht>1000&&njets>=4&&njets<=5", "weight", {400.}),
                       full_trig_skim_0l, all_plot_types); 
  
  pm.Push<HistoStack>(HistoDef(30, 0., 1500., "mj14", "M_{J} [GeV]",
                      "nleps==0&&met<50&&ht>1000&&njets>=6&&njets<=8", "weight", {400.}),
                       full_trig_skim_0l, all_plot_types); 
  pm.Push<HistoStack>(HistoDef(30, 0., 500., "fjets14_m", "m_{J} [GeV]",
                      "nleps==0&&met<50&&ht>1000&&njets>=6&&njets<=8", "weight", {400.}),
                       full_trig_skim_0l, all_plot_types); 

  pm.Push<HistoStack>(HistoDef(30, 0., 1500., "mj14", "M_{J} [GeV]",
                      "nleps==0&&met<50&&ht>1000&&njets>=9", "weight", {400.}),
                       full_trig_skim_0l, all_plot_types);
  pm.Push<HistoStack>(HistoDef(30, 0., 500., "fjets14_m", "m_{J} [GeV]",
                      "nleps==0&&met<50&&ht>1000&&njets>=9", "weight", {400.}),
                       full_trig_skim_0l, all_plot_types); 


  //
  // di-lepton events
  //
  string mll="(mumu_m*(mumu_m>0&&mumu_pt1>30)+elel_m*(elel_m>0&&elel_pt1>30))";
  string mllcut="(mumu_m*(mumu_m>0&&mumu_pt1>30)+elel_m*(elel_m>0&&elel_pt1>30))>80&&(mumu_m*(mumu_m>0&&mumu_pt1>30)+elel_m*(elel_m>0&&elel_pt1>30))<100";
  // ttbar
  pm.Push<HistoStack>(HistoDef(15, 0., 1500., "mj14", "M_{J} [GeV]",
                      "nels==1&&nmus==1&&ht>300&&met>100&&njets>=2&&njets<=3&&met<400&&nbm<3&&nbm>=1&&((leps_pt[0]>30&&leps_eta[0]<2.1&&leps_eta[0]>-2.1)||(leps_pt[1]>30&&leps_eta[1]<2.1&&leps_eta[1]>-2.1))", "weight", {400.}),
                       full_trig_skim_2l, all_plot_types);
  pm.Push<HistoStack>(HistoDef(15, 0., 500., "fjets14_m", "m_{J} [GeV]",
                      "nels==1&&nmus==1&&ht>300&&met>100&&njets>=2&&njets<=3&&met<400&&nbm<3&&nbm>=1&&((leps_pt[0]>30&&leps_eta[0]<2.1&&leps_eta[0]>-2.1)||(leps_pt[1]>30&&leps_eta[1]<2.1&&leps_eta[1]>-2.1))", "weight", {400.}),
                       full_trig_skim_2l, all_plot_types);
  pm.Push<HistoStack>(HistoDef(15, 0., 1500., "mj14", "M_{J} [GeV]",
                      "nels==1&&nmus==1&&ht>300&&met>100&&njets>=4&&met<400&&nbm<3&&nbm>=1&&((leps_pt[0]>30&&leps_eta[0]<2.1&&leps_eta[0]>-2.1)||(leps_pt[1]>30&&leps_eta[1]<2.1&&leps_eta[1]>-2.1))", "weight", {400.}),
                       full_trig_skim_2l, all_plot_types);
  pm.Push<HistoStack>(HistoDef(15, 0., 500., "fjets14_m", "m_{J} [GeV]",
                      "nels==1&&nmus==1&&ht>300&&met>100&&njets>=4&&met<400&&nbm<3&&nbm>=1&&((leps_pt[0]>30&&leps_eta[0]<2.1&&leps_eta[0]>-2.1)||(leps_pt[1]>30&&leps_eta[1]<2.1&&leps_eta[1]>-2.1))", "weight", {400.}),
                       full_trig_skim_2l, all_plot_types);
  // DY
  pm.Push<HistoStack>(HistoDef(15, 0., 1500., "mj14", "M_{J} [GeV]",
                      "(nmus>=2||nels>=2)&&ht>350&&njets>=2&&njets<=3&&"+mllcut, "weight", {400.}),
                       full_trig_skim_2l, all_plot_types);
  pm.Push<HistoStack>(HistoDef(15, 0., 500., "fjets14_m", "m_{J} [GeV]",
                      "(nmus>=2||nels>=2)&&ht>350&&njets>=2&&njets<=3&&"+mllcut, "weight", {400.}),
                       full_trig_skim_2l, all_plot_types);
  pm.Push<HistoStack>(HistoDef(15, 0., 1500., "mj14", "M_{J} [GeV]",
                      "(nmus>=2||nels>=2)&&ht>350&&njets>=4&&"+mllcut, "weight", {400.}),
                       full_trig_skim_2l, all_plot_types);
  pm.Push<HistoStack>(HistoDef(15, 0., 500., "fjets14_m", "m_{J} [GeV]",
                      "(nmus>=2||nels>=2)&&ht>350&&njets>=4&&"+mllcut, "weight", {400.}),
                       full_trig_skim_2l, all_plot_types);
  
  pm.MakePlots(lumi); 

}
