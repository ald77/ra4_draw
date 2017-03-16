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

int main(){
  gErrorIgnoreLevel = 6000;

  double lumi = 35.9;

  string trig_mc        = "/net/cms29/cms29r0/babymaker/babies/2017_01_27/mc/unskimmed/";
  string trig_skim1l_mc = "/net/cms29/cms29r0/babymaker/babies/2017_01_27/mc/merged_mcbase_stdnj5/";
  //string trig_skim1l_mc = "/net/cms29/cms29r0/babymaker/babies/2017_01_27/mc/skim_met100/"; // for njets plots
  string trig_skim0l_mc = "/net/cms29/cms29r0/babymaker/babies/2017_01_27/mc/merged_qcd/"; // this does not exist
  string trig_skim2l_mc = "/net/cms29/cms29r0/babymaker/babies/2017_01_27/mc/merged_dy_ht300/"; // this does not exist 

  Palette colors("txt/colors.txt", "default");

  auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {trig_skim1l_mc+"*_TTJets*SingleLept*.root"}, "ntruleps<=1&&stitch_met");
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {trig_skim1l_mc+"*_TTJets*DiLept*.root"}, "ntruleps>=2&&stitch_met");
  auto wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {trig_skim1l_mc+"*_WJetsToLNu*.root"},"stitch_met");
  auto single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {trig_skim1l_mc+"*_ST_*.root"});
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {trig_skim1l_mc+"*_TTWJets*.root", trig_skim1l_mc+"*_TTZ*.root"});
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {trig_skim1l_mc+"*DYJetsToLL*.root", trig_skim1l_mc+"*QCD_HT*0_Tune*.root", trig_skim1l_mc+"*QCD_HT*Inf_Tune*.root",
        trig_skim1l_mc+"*_ZJet*.root", trig_skim1l_mc+"*_ttHTobb_M125_*.root",
        trig_skim1l_mc+"*_TTGJets*.root", trig_skim1l_mc+"*_TTTT_*.root",
        trig_skim1l_mc+"*_WH_HToBB*.root", trig_skim1l_mc+"*_ZH_HToBB*.root", 
        trig_skim1l_mc+"*_WWTo*.root", trig_skim1l_mc+"*_WZ*.root",
        trig_skim1l_mc+"_ZZ_*.root"}, "stitch_met");

  auto t1tttt_nc = Process::MakeShared<Baby_full>("T1tttt(1800,100)", Process::Type::signal, colors("t1tttt"),
    {"/net/cms29/cms29r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/*SMS-T1tttt_mGluino-1800_mLSP-100_*.root"});
  auto t1tttt_c = Process::MakeShared<Baby_full>("T1tttt(1400,1000)", Process::Type::signal, colors("t1tttt"),
    {"/net/cms29/cms29r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/*SMS-T1tttt_mGluino-1400_mLSP-1000*.root"});
  t1tttt_c->SetLineStyle(2);

  auto data_1l = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    //{"/net/cms29/cms29r0/babymaker/babies/2017_02_14/data/skim_met100/*.root"},"pass&&(met/met_calo<5.0)&&pass_ra2_badmu&&trig_ra4"); // for njets plots
    {"/net/cms29/cms29r0/babymaker/babies/2017_02_14/data/merged_database_stdnj5/*.root"},"pass&&(met/met_calo<5.0)&&pass_ra2_badmu&&trig_ra4");
  vector<shared_ptr<Process> > full_trig_skim_1l = {data_1l, t1tttt_nc, t1tttt_c, tt1l, tt2l, wjets, single_t, ttv, other};

  //
  // 0-lepton plots
  //
  auto tt_0l = Process::MakeShared<Baby_full>("t#bar{t}", Process::Type::background, colors("tt_1l"),
    {trig_skim0l_mc+"*_TTJets_TuneCUETP8M1_13TeV*.root"});
  auto qcd_0l = Process::MakeShared<Baby_full>("QCD", Process::Type::background, colors("qcd"),
    {trig_skim0l_mc+"*_QCD_HT700to1000*.root",trig_skim0l_mc+"*_QCD_HT1000to1500*.root",
        trig_skim0l_mc+"*_QCD_HT1500to2000*.root",trig_skim0l_mc+"*_QCD_HT2000toInf*.root"});
  auto other_0l = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("ttv"),
    {trig_skim0l_mc+"*_TTWJet*root",  trig_skim0l_mc+"*_TTZTo*.root",
        trig_skim0l_mc+"*_ZJet*.root", trig_skim0l_mc+"*_WWTo*.root",
        trig_skim0l_mc+"*ggZH_HToBB*.root", trig_skim0l_mc+"*ttHJetTobb*.root",
        trig_skim0l_mc+"*_TTGJets*.root", trig_skim0l_mc+"*_TTTT_*.root",
        trig_skim0l_mc+"*_WH_HToBB*.root", trig_skim0l_mc+"*_WZTo*.root",
        trig_skim0l_mc+"*_ST_*.root", trig_skim0l_mc+"*_WJetsToLNu*.root",
        trig_skim0l_mc+"*_ZH_HToBB*.root", trig_skim0l_mc+"_ZZ_*.root",
        trig_skim0l_mc+"*_DYJetsToLL*HT*.root"});

  auto data_0l = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {"/net/cms2/cms2r0/babymaker/babies/2016_06_21/data/skim_ht900/*.root"},"pass&&trig[12]");
  vector<shared_ptr<Process> > full_trig_skim_0l = {data_0l, /*t1tttt_nc, t1tttt_c,*/ qcd_0l, tt_0l, other_0l};

  //
  // 2-lepton plots
  //
  auto tt1l_2l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {trig_skim2l_mc+"*_TTJets*Lept*.root"},
    "ntruleps<=1");
  auto tt2l_2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {trig_skim2l_mc+"*_TTJets*Lept*.root"},
    "ntruleps>=2");
  auto dy_2l = Process::MakeShared<Baby_full>("DY", Process::Type::background, colors("wjets"),
    {trig_skim2l_mc+"*_DYJetsToLL*HT*.root"});
  auto other_2l = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("ttv"),
    {trig_skim2l_mc+"*_TTWJet*root",  trig_skim2l_mc+"*_TTZTo*.root",
        trig_skim2l_mc+"*_ZJet*.root", trig_skim2l_mc+"*_WWTo*.root",
        trig_skim2l_mc+"*ggZH_HToBB*.root", trig_skim2l_mc+"*ttHJetTobb*.root",
        trig_skim2l_mc+"*_TTGJets*.root", trig_skim2l_mc+"*_TTTT_*.root",
        trig_skim2l_mc+"*_WH_HToBB*.root", trig_skim2l_mc+"*_WZTo*.root",
        trig_skim2l_mc+"*_ST_*.root", trig_skim2l_mc+"*_WJetsToLNu*.root",
        trig_skim2l_mc+"*_ZH_HToBB*.root", trig_skim2l_mc+"_ZZ_*.root",
        trig_skim2l_mc+"*_QCD_HT700to1000*.root",trig_skim2l_mc+"*_QCD_HT1000to1500*.root",
        trig_skim2l_mc+"*_QCD_HT1500to2000*.root",trig_skim2l_mc+"*_QCD_HT2000toInf*.root"});

  auto data_2l = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
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

  // No need to be blind anymore
  //NamedFunc lowmtcut = "st>500&&met>200&&mt<140"; 
  //NamedFunc lowmjcut = "st>500&&met>200&&mj14<400";
  NamedFunc lowmtcut = "st>500&&met>200";
  NamedFunc lowmjcut = "st>500&&met>200";
  PlotMaker pm;

  for(int ilep=0; ilep<3; ilep++){
    NamedFunc lepcut(true);
    if(ilep==0)  lepcut="nels==1&&nmus==0";
    if(ilep==1)  lepcut="nels==0&&nmus==1";
    if(ilep==2)  lepcut="nleps==1";
    if(ilep!=2) continue;

    for(int inveto=0; inveto<2; inveto++){
      NamedFunc nvetocut(true);
      if(inveto==0)  nvetocut="nveto==0";
      if(inveto==1)  nvetocut="1";
      if(inveto==1)  continue;

      for(int injets=0; injets<2; injets++){
        NamedFunc njetscut(true);
        if(injets==0)  njetscut="njets>=5";
        if(injets==1)  njetscut="njets>=6";
        if(injets==0)  continue;

        for(int inb=0; inb<2; inb++){
          NamedFunc nbcut(true);
          if(inb==0)  nbcut="nbm>=1";
          if(inb==1)  nbcut="nbm>=2";
          if(inb==1)  continue;

// 12.9 fb-1 nbm plot with rereco
// - need to change lumi and add json12p9 to data
// - note that this is fully unblinded, i.e., no mT<140 cut
// 
//          pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm", "N_{b}", {0.5, 1.5, 2.5}),
//                          lepcut&&njetscut&&nvetocut&&"st>500&&met>200",
//                          full_trig_skim_1l, all_plot_types);

          //
          // event-level variables except ones related to fatjet
          //
          pm.Push<Hist1D>(Axis(15, 500, 2000., "st", "S_{T} [GeV]", {500.}),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
                          lepcut&&lowmjcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}", {5.5, 8.5}),
                          lepcut&&lowmtcut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm", "N_{b}", {0.5, 1.5, 2.5}),
                          lepcut&&lowmtcut&&njetscut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(14, 0., 280., "mt", "m_{T} [GeV]", {140.}),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(14, 0., 280., "mt", "m_{T} [GeV]", {140.}),
                          lepcut&&lowmjcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(14, 0., 280., "mt", "m_{T} [GeV]", {140.}),
                          lepcut&&lowmtcut&&njetscut&&nvetocut&&"nbm==0",
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(14, 0., 280., "mt", "m_{T} [GeV]", {140.}),
                          lepcut&&lowmjcut&&njetscut&&nvetocut&&"nbm==0",
                          full_trig_skim_1l, all_plot_types);

          if(ilep==0){
            // only barrel or endcap electron
            pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
                            lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"leps_eta[0]<1.479&&leps_eta[0]>-1.479",
                            full_trig_skim_1l, all_plot_types);
            pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
                            lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"leps_eta[0]>1.479||leps_eta[0]<-1.479",
                            full_trig_skim_1l, all_plot_types);
          }

          //
          // leptons
          //
          if(ilep==0){// electron
            pm.Push<Hist1D>(Axis(10, 0, 200., "leps_pt[0]", "p_{T}(electron) [GeV]"),
                            lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                            full_trig_skim_1l, all_plot_types);
            pm.Push<Hist1D>(Axis(10, -2.5, 2.5, "leps_eta[0]", "#eta(electron) [GeV]"),
                            lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                            full_trig_skim_1l, all_plot_types);
          }else if(ilep==1){// muon
            pm.Push<Hist1D>(Axis(10, 0, 200., "leps_pt[0]", "p_{T}(muon) [GeV]"),
                            lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                            full_trig_skim_1l, all_plot_types);
            pm.Push<Hist1D>(Axis(10, -2.5, 2.5, "leps_eta[0]", "#eta(muon) [GeV]"),
                            lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                            full_trig_skim_1l, all_plot_types);
          }else if(ilep==2){  // electron or muon
            pm.Push<Hist1D>(Axis(10, 0, 200., "leps_pt[0]", "p_{T}(lepton) [GeV]"),
                            lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                            full_trig_skim_1l, all_plot_types);
            pm.Push<Hist1D>(Axis(10, -2.5, 2.5, "leps_eta[0]", "#eta(lepton) [GeV]"),
                            lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                            full_trig_skim_1l, all_plot_types);
          }

          //
          // jets
          //
          pm.Push<Hist1D>(Axis(10, 30, 630., "jets_pt[0]", "p_{T}(leading jet) [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"!jets_islep[0]",
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, -2.5, 2.5, "jets_eta[0]", "#eta(leading jet) [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"!jets_islep[0]",
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 0, 200., "jets_m[0]", "mass(leading jet) [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"!jets_islep[0]",
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 30, 630., "jets_pt", "p_{T}(jet) [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"!jets_islep",
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, -2.5, 2.5, "jets_eta", "#eta(jet) [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"!jets_islep",
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 0, 200., "jets_m", "m(jet) [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"!jets_islep",
                          full_trig_skim_1l, all_plot_types);

          //
          // fatjets
          //
          pm.Push<Hist1D>(Axis(10, 0, 1000., "fjets14_pt[0]", "p_{T}(J1) [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 0, 500., "fjets14_m[0]", "m_{J1} [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 0, 500., "fjets14_m[0]", "m_{J1} [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"leps_pt[0]>160",
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 0, 500., "fjets14_m[0]", "m_{J1} [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut&&"leps_pt[0]<=160",
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, -2.5, 2.5, "fjets14_eta[0]", "#eta(J1)"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(7, 0.5, 7.5, "fjets14_nconst[0]", "N_{constituents}(J1)"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 0, 1000., "fjets14_pt", "p_{T}(J) [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, 0, 500., "fjets14_m", "m_{J} [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(10, -2.5, 2.5, "fjets14_eta", "#eta(J)"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(7, 0.5, 7.5, "fjets14_nconst", "N_{constituents}(J)"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);

          pm.Push<Hist1D>(Axis(40, 0, 400., "fjets14_m[0]", "m_{J1} [GeV]"),
                          lepcut&&lowmtcut&&"njets<=5&&nbm>=1",
                          full_trig_skim_1l, all_plot_types);

          // event-level variables related to fatjet
          pm.Push<Hist1D>(Axis(15, 0., 1500., "mj14", "M_{J} [GeV]", {400.}),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
          pm.Push<Hist1D>(Axis(7, 0.5, 7.5, "nfjets14", "N_{J} [GeV]"),
                          lepcut&&lowmtcut&&njetscut&&nbcut&&nvetocut,
                          full_trig_skim_1l, all_plot_types);
        } //for(int inb=0; inb<2; inb++)
      } //for(int injets=0; injets<2; injets++)
    } //for(int inveto=0; inveto<2; inveto++)
  }
/*
  //
  // 0-lepton events
  //
  pm.Push<Hist1D>(Axis(30, 0., 1500., "mj14", "M_{J} [GeV]", {400.}),
                  "nleps==0&&met<50&&ht>1000&&njets>=4&&njets<=5",
                  full_trig_skim_0l, all_plot_types);
  pm.Push<Hist1D>(Axis(30, 0., 500., "fjets14_m", "m_{J} [GeV] [GeV]", {400.}),
                  "nleps==0&&met<50&&ht>1000&&njets>=4&&njets<=5",
                  full_trig_skim_0l, all_plot_types);

  pm.Push<Hist1D>(Axis(30, 0., 1500., "mj14", "M_{J} [GeV]", {400.}),
                  "nleps==0&&met<50&&ht>1000&&njets>=6&&njets<=8",
                  full_trig_skim_0l, all_plot_types);
  pm.Push<Hist1D>(Axis(30, 0., 500., "fjets14_m", "m_{J} [GeV]", {400.}),
                  "nleps==0&&met<50&&ht>1000&&njets>=6&&njets<=8",
                  full_trig_skim_0l, all_plot_types);

  pm.Push<Hist1D>(Axis(30, 0., 1500., "mj14", "M_{J} [GeV]", {400.}),
                  "nleps==0&&met<50&&ht>1000&&njets>=9",
                  full_trig_skim_0l, all_plot_types);
  pm.Push<Hist1D>(Axis(30, 0., 500., "fjets14_m", "m_{J} [GeV]", {400.}),
                  "nleps==0&&met<50&&ht>1000&&njets>=9",
                  full_trig_skim_0l, all_plot_types);

  //
  // di-lepton events
  //
  string mll="(mumu_m*(mumu_m>0&&mumu_pt1>30)+elel_m*(elel_m>0&&elel_pt1>30))";
  string mllcut="(mumu_m*(mumu_m>0&&mumu_pt1>30)+elel_m*(elel_m>0&&elel_pt1>30))>80&&(mumu_m*(mumu_m>0&&mumu_pt1>30)+elel_m*(elel_m>0&&elel_pt1>30))<100";
  // ttbar
  pm.Push<Hist1D>(Axis(15, 0., 1500., "mj14", "M_{J} [GeV]", {400.}),
                  "nels==1&&nmus==1&&ht>300&&met>100&&njets>=2&&njets<=3&&met<400&&nbm<3&&nbm>=1&&((leps_pt[0]>30&&leps_eta[0]<2.1&&leps_eta[0]>-2.1)||(leps_pt[1]>30&&leps_eta[1]<2.1&&leps_eta[1]>-2.1))",
                  full_trig_skim_2l, all_plot_types);
  pm.Push<Hist1D>(Axis(15, 0., 500., "fjets14_m", "m_{J} [GeV]", {400.}),
                  "nels==1&&nmus==1&&ht>300&&met>100&&njets>=2&&njets<=3&&met<400&&nbm<3&&nbm>=1&&((leps_pt[0]>30&&leps_eta[0]<2.1&&leps_eta[0]>-2.1)||(leps_pt[1]>30&&leps_eta[1]<2.1&&leps_eta[1]>-2.1))",
                  full_trig_skim_2l, all_plot_types);
  pm.Push<Hist1D>(Axis(15, 0., 1500., "mj14", "M_{J} [GeV]", {400.}),
                  "nels==1&&nmus==1&&ht>300&&met>100&&njets>=4&&met<400&&nbm<3&&nbm>=1&&((leps_pt[0]>30&&leps_eta[0]<2.1&&leps_eta[0]>-2.1)||(leps_pt[1]>30&&leps_eta[1]<2.1&&leps_eta[1]>-2.1))",
                  full_trig_skim_2l, all_plot_types);
  pm.Push<Hist1D>(Axis(15, 0., 500., "fjets14_m", "m_{J} [GeV]", {400.}),
                  "nels==1&&nmus==1&&ht>300&&met>100&&njets>=4&&met<400&&nbm<3&&nbm>=1&&((leps_pt[0]>30&&leps_eta[0]<2.1&&leps_eta[0]>-2.1)||(leps_pt[1]>30&&leps_eta[1]<2.1&&leps_eta[1]>-2.1))",
                  full_trig_skim_2l, all_plot_types);
  // DY
  pm.Push<Hist1D>(Axis(15, 0., 1500., "mj14", "M_{J} [GeV]", {400.}),
                  "(nmus>=2||nels>=2)&&ht>350&&njets>=2&&njets<=3"&&mllcut,
                  full_trig_skim_2l, all_plot_types);
  pm.Push<Hist1D>(Axis(15, 0., 500., "fjets14_m", "m_{J} [GeV]", {400.}),
                  "(nmus>=2||nels>=2)&&ht>350&&njets>=2&&njets<=3"&&mllcut,
                  full_trig_skim_2l, all_plot_types);
  pm.Push<Hist1D>(Axis(15, 0., 1500., "mj14", "M_{J} [GeV]", {400.}),
                  "(nmus>=2||nels>=2)&&ht>350&&njets>=4"&&mllcut,
                  full_trig_skim_2l, all_plot_types);
  pm.Push<Hist1D>(Axis(15, 0., 500., "fjets14_m", "m_{J} [GeV]", {400.}),
                  "(nmus>=2||nels>=2)&&ht>350&&njets>=4"&&mllcut,
                  full_trig_skim_2l, all_plot_types);
*/ 

 // nleps==1 && nveto==1
  pm.Push<Hist1D>(Axis(10, 0, 1000., "fjets14_pt[0]", "p_{T}(J1) [GeV]"),
                  "st>500&&nleps==1&&nveto==1&&njets>=6&&mt>140&&nbm>=1&&met>200", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, 0, 500., "fjets14_m[0]", "m_{J1} [GeV]"),
                  "st>500&&nleps==1&&nveto==1&&njets>=6&&mt>140&&nbm>=1&&met>200", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, -2.5, 2.5, "fjets14_eta[0]", "#eta(J1)"),
                  "st>500&&nleps==1&&nveto==1&&njets>=6&&mt>140&&nbm>=1&&met>200", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(7, 0.5, 7.5, "fjets14_nconst[0]", "N_{constituents}(J1)"),
                  "st>500&&nleps==1&&nveto==1&&njets>=6&&mt>140&&nbm>=1&&met>200", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, 0, 1000., "fjets14_pt", "p_{T}(J) [GeV]"),
                  "st>500&&nleps==1&&nveto==1&&njets>=6&&mt>140&&nbm>=1&&met>200", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, 0, 500., "fjets14_m", "m_{J} [GeV]"),
                  "st>500&&nleps==1&&nveto==1&&njets>=6&&mt>140&&nbm>=1&&met>200", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, -2.5, 2.5, "fjets14_eta", "#eta(J)"),
                  "st>500&&nleps==1&&nveto==1&&njets>=6&&mt>140&&nbm>=1&&met>200", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(7, 0.5, 7.5, "fjets14_nconst", "N_{constituents}(J)"),
                  "st>500&&nleps==1&&nveto==1&&njets>=6&&mt>140&&nbm>=1&&met>200", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(15, 0., 2000., "mj14", "M_{J} [GeV]", {400.}),
                  "st>500&&nleps==1&&nveto==1&&njets>=6&&mt>140&&nbm>=1&&met>200",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(7, 0.5, 7.5, "nfjets14", "N_{J} [GeV]"),
                  "st>500&&nleps==1&&nveto==1&&njets>=6&&mt>140&&nbm>=1&&met>200", 
                  full_trig_skim_1l, all_plot_types); 
  pm.Push<Hist1D>(Axis(15, 500, 2000., "st", "S_{T} [GeV]", {500.}),
          "st>500&&nleps==1&&nveto==1&&njets>=6&&mt>140&&nbm>=1&&met>200", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
          "st>500&&nleps==1&&nveto==1&&njets>=6&&mt>140&&nbm>=1&&met>200", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}", {5.5, 8.5}),
          "st>500&&nleps==1&&nveto==1&&mt>140&&nbm>=1&&met>200", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm", "N_{b}", {0.5, 1.5, 2.5}),
          "st>500&&nleps==1&&nveto==1&&nbm>=0&&njets>=6&&mt>140&&met>200", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(14, 0., 280., "mt", "m_{T} [GeV]", {140.}),
          "st>500&&nleps==1&&nveto==1&&njets>=6&&mt>140&&nbm>=1&&met>200",
                  full_trig_skim_1l, all_plot_types);

 // nleps==2
  pm.Push<Hist1D>(Axis(10, 0, 1000., "fjets14_pt[0]", "p_{T}(J1) [GeV]"),
                  "st>500&&nleps==2&&njets>=5&&nbm<=2&&met>200&&met<500", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, 0, 500., "fjets14_m[0]", "m_{J1} [GeV]"),
                  "st>500&&nleps==2&&njets>=5&&nbm<=2&&met>200&&met<500", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, -2.5, 2.5, "fjets14_eta[0]", "#eta(J1)"),
                  "st>500&&nleps==2&&njets>=5&&nbm<=2&&met>200&&met<500", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(7, 0.5, 7.5, "fjets14_nconst[0]", "N_{constituents}(J1)"),
                  "st>500&&nleps==2&&njets>=5&&nbm<=2&&met>200&&met<500", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, 0, 1000., "fjets14_pt", "p_{T}(J) [GeV]"),
                  "st>500&&nleps==2&&njets>=5&&nbm<=2&&met>200&&met<500", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, 0, 500., "fjets14_m", "m_{J} [GeV]"),
                  "st>500&&nleps==2&&njets>=5&&nbm<=2&&met>200&&met<500", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, -2.5, 2.5, "fjets14_eta", "#eta(J)"),
                  "st>500&&nleps==2&&njets>=5&&nbm<=2&&met>200&&met<500", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(7, 0.5, 7.5, "fjets14_nconst", "N_{constituents}(J)"),
                  "st>500&&nleps==2&&njets>=5&&nbm<=2&&met>200&&met<500", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(15, 0., 2000., "mj14", "M_{J} [GeV]", {400.}),
                  "st>500&&nleps==2&&njets>=5&&nbm<=2&&met>200&&met<500",
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(7, 0.5, 7.5, "nfjets14", "N_{J} [GeV]"),
                  "st>500&&nleps==2&&njets>=5&&nbm<=2&&met>200&&met<500", 
                  full_trig_skim_1l, all_plot_types); 
  pm.Push<Hist1D>(Axis(15, 500, 2000., "st", "S_{T} [GeV]", {500.}),
          "st>500&&nleps==2&&njets>=5&&nbm<=2&&met>200&&met<500", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(10, 200, 700., "met", "E_{T}^{miss} [GeV]", {200., 350., 500.}),
          "st>500&&nleps==2&&njets>=5&&nbm<=2&&met>200&&met<500", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}", {4.5, 7.5}),
          "st>500&&nleps==2&&nbm<=2&&met>200&&met<500", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(7, -0.5, 6.5, "nbm", "N_{b}", {0.5, 1.5, 2.5}),
          "st>500&&nleps==2&&njets>=5&&nbm<=2&&met>200&&met<500", 
                  full_trig_skim_1l, all_plot_types);
  pm.Push<Hist1D>(Axis(14, 0., 280., "mt", "m_{T} [GeV]", {140.}),
          "st>500&&nleps==2&&njets>=5&&nbm<=2&&met>200&&met<500", 
                  full_trig_skim_1l, all_plot_types);

  pm.MakePlots(lumi);

}
