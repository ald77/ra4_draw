#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TColor.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  bool single_thread = true;
  TString method = "";
}

void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  double lumi = 2.6;

  string bfolder("");
  string hostname = execute("HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder
  string lowmet_fndata(bfolder+"/cms2r0/babymaker/babies/2016_06_26/data/merged_met150/");
  string himet_fndata(bfolder+"/cms2r0/babymaker/babies/2016_06_26/data/merged_standard/");
  string himet_fodata(bfolder+"/cms2r0/babymaker/babies/2016_04_29/data/merged_1lht500met200/");
  string lowmet_fmc(bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_met150/");
  string himet_fmc(bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_standard/");

  string baseline("nleps>=1 && ht>500 && met>150 && met<500 && pass && njets>=5");

  // vector<TString> abcdcuts_veto = {"mt<=140 && mj14<=400 && nleps==1 && nveto==0 && nbm>=1  &&  njets>=6",
  //                                  "mt<=140 && mj14>400  && nleps==1 && nveto==0 && nbm>=1  &&  njets>=6 && njets<=8",
  //                                  "mt>140  && mj14<=400 && nleps==1 && nveto==1 && nbm>=1 && nbm<=2  &&  njets>=6",
  //                                  "mt>140  && mj14>400  && nleps==1 && nveto==1 && nbm>=1 && nbm<=2  &&  njets>=6 && njets<=8"};

  vector<string> abcdcuts_2l = {"mt<=140 && mj14<=400 && nleps==1 && nveto==0 && nbm>=1 &&  njets>=6",
                                "mt<=140 && mj14>400  && nleps==1 && nveto==0 && nbm>=1 &&  njets>=6",
                                "           mj14<=400 && nels==1 && nmus==1 && nbm<=2 &&  njets>=5",
                                "           mj14>400  && nels==1 && nmus==1 && nbm<=2  &&  njets>=5"};

  Palette colors("txt/colors.txt", "default");

  // auto r1data = Process::MakeShared<Baby_full>("R1 data", Process::Type::background, kAzure+1, {himet_fndata, lowmet_fndata}, baseline+"&&"+abcdcuts_2l[0]);

  // auto r2data = Process::MakeShared<Baby_full>("R2 data", Process::Type::background, kOrange, {himet_fndata, lowmet_fndata}, baseline+"&&"+abcdcuts_2l[1]);

  // auto d3data = Process::MakeShared<Baby_full>("D3 data", Process::Type::background, kGreen+2, {himet_fndata, lowmet_fndata}, baseline+"&&"+abcdcuts_2l[2]);

  // auto d4data = Process::MakeShared<Baby_full>("D4 data", Process::Type::background, kRed+1, {himet_fndata, lowmet_fndata}, baseline+"&&"+abcdcuts_2l[3]);

  // vector<shared_ptr<Process> > all_regs = {r1data, d3data, r2data, d4data};

  auto ndata = Process::MakeShared<Baby_full>("Data (2016, 2.6 fb^{-1})", Process::Type::data, kBlack, {lowmet_fndata+"/*.root", himet_fndata+"/*.root"},
                                              baseline+"&&json2p6 && (trig[4]||trig[8]||trig[13]||trig[33])");
  auto odata = Process::MakeShared<Baby_full>("Data (2015, 2.3 fb^{-1})", Process::Type::data, kRed+1, {himet_fodata+"/*.root"},
                                              baseline+" && (trig[4]||trig[8]||trig[28]||trig[14])");

  auto tt1l = Process::MakeShared<Baby_full>("tt 1lep", Process::Type::background, colors("tt_1l"),
    {lowmet_fmc+"*_TTJets*Lept*.root", lowmet_fmc+"*_TTJets_HT*.root",
        himet_fmc+"*_TTJets*Lept*.root", himet_fmc+"*_TTJets_HT*.root"},
    baseline+" && stitch && ntruleps==1");
  auto tt2l = Process::MakeShared<Baby_full>("tt 2lep", Process::Type::background, colors("tt_2l"),
    {lowmet_fmc+"*_TTJets*Lept*.root", lowmet_fmc+"*_TTJets_HT*.root",
        himet_fmc+"*_TTJets*Lept*.root", himet_fmc+"*_TTJets_HT*.root"},
    baseline+" && stitch && ntruleps==2");
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {lowmet_fmc+"*_WJetsToLNu*.root",lowmet_fmc+"*_ST_*.root",
        lowmet_fmc+"*_TTW*.root",lowmet_fmc+"*_TTZ*.root",
        lowmet_fmc+"*DYJetsToLL*.root",lowmet_fmc+"*QCD_HT*.root",
        lowmet_fmc+"*_ZJet*.root",lowmet_fmc+"*_ttHJetTobb*.root",
        lowmet_fmc+"*_TTGJets*.root",lowmet_fmc+"*_TTTT*.root",
        lowmet_fmc+"*_WH_HToBB*.root",lowmet_fmc+"*_ZH_HToBB*.root",
        lowmet_fmc+"*_WWTo*.root",lowmet_fmc+"*_WZ*.root",lowmet_fmc+"*_ZZ_*.root",
        himet_fmc+"*_WJetsToLNu*.root",himet_fmc+"*_ST_*.root",
        himet_fmc+"*_TTW*.root",himet_fmc+"*_TTZ*.root",
        himet_fmc+"*DYJetsToLL*.root",himet_fmc+"*QCD_HT*.root",
        himet_fmc+"*_ZJet*.root",himet_fmc+"*_ttHJetTobb*.root",
        himet_fmc+"*_TTGJets*.root",himet_fmc+"*_TTTT*.root",
        himet_fmc+"*_WH_HToBB*.root",himet_fmc+"*_ZH_HToBB*.root",
        himet_fmc+"*_WWTo*.root",himet_fmc+"*_WZ*.root",himet_fmc+"*_ZZ_*.root"},
    baseline+" && stitch");

  vector<shared_ptr<Process> > lowmet_procs = {ndata, tt1l, tt2l, other};
  vector<shared_ptr<Process> > himet_procs = {ndata, odata, tt1l, tt2l, other};
  vector<shared_ptr<Process> > allmet_procs = {ndata, odata, tt1l, tt2l, other};

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
  // vector<PlotOpt> all_plot_types = {log_lumi, lin_lumi, log_shapes, lin_shapes,
  //                                   log_lumi_info, lin_lumi_info, log_shapes_info, lin_shapes_info};

  vector<PlotOpt> plot_types = {lin_lumi_info};

  PlotMaker pm;
  // pm.Push<Hist1D>(Axis(25, 0., 600., "mus_pt[0]", "Muon pT"), "nmus==1", all_regs, plot_types).Weight(1./2.6);
  // pm.Push<Hist1D>(Axis(20, 0., 0.1, "mus_pterr[0]/mus_pt[0]", "Muon relative pT error"), "nmus==1", all_regs, plot_types).Weight(1./2.6);
  // pm.Push<Hist1D>(Axis(20, 0., 20, "mus_trk_algo[0]", "Muon reco algorithm"), "nmus==1", all_regs, plot_types).Weight(1./2.6);
  // pm.Push<Hist1D>(Axis(3, 0., 3, "mus_trk_nholes_in[0]", "Holes in inner track"), "nmus==1", all_regs, plot_types).Weight(1./2.6);
  // pm.Push<Hist1D>(Axis(5, 0., 5, "mus_trk_nholes_out[0]", "Holes in outer track"), "nmus==1", all_regs, plot_types).Weight(1./2.6);
  // pm.Push<Hist1D>(Axis(2, 0., 2, "mus_trk_quality[0]", "Muon is high purity"), "nmus==1", all_regs, plot_types).Weight(1./2.6);
  // pm.Push<Hist1D>(Axis(2, 0., 2, "mus_tight[0]", "Muon is tight"), "nmus==1", all_regs, plot_types).Weight(1./2.6);
  // pm.Push<Hist1D>(Axis(3, -1., 2, "mus_charge[0]", "Muon charge"), "nmus==1", all_regs, plot_types).Weight(1./2.6);
  // pm.Push<Hist1D>(Axis(30, 0., 5., "mus_em_e[0]", "Muon E deposited in ECAL"), "nmus==1", all_regs, plot_types).Weight(1./2.6);
  // pm.Push<Hist1D>(Axis(30, 0., 20., "mus_had_e[0]", "Muon E deposited in HCAL"), "nmus==1", all_regs, plot_types).Weight(1./2.6);

  // pm.Push<Hist1D>(Axis(15, 0., 600., "elmu_m", "elmu_m"), true, all_regs, plot_types).Weight(1./2.6);
  // pm.Push<Hist1D>(Axis(15, 0., 600., "elmu_pt", "elmu_pt"), true, all_regs, plot_types).Weight(1./2.6);
  // pm.Push<Hist1D>(Axis(16, 0., 3.2, "dphi_lmet", "dphi_lmet"), true, all_regs, plot_types).Weight(1./2.6);
  // pm.Push<Hist1D>(Axis(10, 0., 1000, "jetsys_pt", "jetsys_pt"), true, all_regs, plot_types).Weight(1./2.6);
  // pm.Push<Hist1D>(Axis(10, 0., 1000, "jetsys_nob_pt", "jetsys_nob_pt"), true, all_regs, plot_types).Weight(1./2.6);
  // pm.Push<Hist1D>(Axis(20, 0., 1, "mus_miniso[0]", "mus_miniso"), "nmus==1", all_regs, plot_types).Weight(1./2.6);
  // pm.Push<Hist1D>(Axis(20, 0., 1, "mus_reliso[0]", "mus_reliso"), "nmus==1", all_regs, plot_types).Weight(1./2.6);
  // pm.Push<Hist1D>(Axis(20, 0., 200, "mus_reliso[0]*mus_pt[0]", "mus_absiso"), "nmus==1", all_regs, plot_types).Weight(1./2.6);
  // pm.Push<Hist1D>(Axis(21, 0, 21, "run-274400", "Run number"), "nmus==1", all_regs, plot_types).Weight(1./2.6);
  // pm.Push<Hist1D>(Axis(30, 274100, 274421, "run", "Run number"), true, all_regs, plot_types).Weight(1./2.6);

  // ------------------- MET 150-200 ----------------------
  pm.Push<Hist1D>(Axis(13, 25, 1000, "mj14", "M_{J} [GeV]", {250.,400.}), "met>150&&met<200 && nels==1 && nmus==1 && nbm<=2 && njets>=5", lowmet_procs, plot_types);
  // pm.Push<Hist1D>(Axis(13, 25, 1000, "mj14", "M_{J} [GeV]", {250.,400.}), "met>150&&met<200 && nels==2 && nbm<=2 && njets>=5", lowmet_procs, plot_types);
  // pm.Push<Hist1D>(Axis(13, 25, 1000, "mj14", "M_{J} [GeV]", {250.,400.}), "met>150&&met<200 && nmus==2 && nbm<=2 && njets>=5", lowmet_procs, plot_types);
  // pm.Push<Hist1D>(Axis(13, 25, 1000, "mj14", "M_{J} [GeV]", {250.,400.}), "met>150&&met<200 && nleps==2 && nbm<=2 && njets>=5", lowmet_procs, plot_types);

  pm.Push<Hist1D>(Axis(10,0.,280., "mt", "m_{T} [GeV]", {140.}), "met>150&&met<200 && nels==1 && nmus==1 && nbm<=2 && njets>=5",lowmet_procs, plot_types);
  pm.Push<Hist1D>(Axis(10,0.,280., "mt", "m_{T} [GeV]", {140.}), "mj14>400 && met>150&&met<200 && nels==1 && nmus==1 && nbm<=2 && njets>=5",lowmet_procs, plot_types);
  // pm.Push<Hist1D>(Axis(6,0.,280., "mt", "m_{T} [GeV]", {140.}), "met>150&&met<200 && nels==2 && nbm<=2 && njets>=5",lowmet_procs, plot_types);
  // pm.Push<Hist1D>(Axis(6,0.,280., "mt", "m_{T} [GeV]", {140.}), "met>150&&met<200 && nmus==2 && nbm<=2 && njets>=5",lowmet_procs, plot_types);
  // pm.Push<Hist1D>(Axis(6,0.,280., "mt", "m_{T} [GeV]", {140.}), "met>150&&met<200 && nleps==2 && nbm<=2 && njets>=5",lowmet_procs, plot_types);

  pm.Push<Hist1D>(Axis(10,0.,400., "mus_pt[0]", "p_{T} (#mu) [GeV]"),"met>150&&met<200 && nleps==2 && nmus>=1 && nbm<=2 && njets>=5",lowmet_procs, plot_types);
  pm.Push<Hist1D>(Axis(10,0.,400., "els_pt[0]", "p_{T} (e)[GeV]"),"met>150&&met<200 && nleps==2 && nels>=1 && nbm<=2 && njets>=5",lowmet_procs, plot_types);
  pm.Push<Hist1D>(Axis(10,0.,400., "mus_pt[0]", "p_{T} (#mu) [GeV]"),"mj14>400 && met>150&&met<200 && nleps==2 && nmus>=1 && nbm<=2 && njets>=5",lowmet_procs, plot_types);
  pm.Push<Hist1D>(Axis(10,0.,400., "els_pt[0]", "p_{T} (e)[GeV]"),"mj14>400 && met>150&&met<200 && nleps==2 && nels>=1 && nbm<=2 && njets>=5",lowmet_procs, plot_types);

  // ------------------- MET 200-350 ----------------------
  pm.Push<Hist1D>(Axis(13, 25, 1000, "mj14", "M_{J} [GeV]",{250.,400.}),"met>200&&met<350 && nels==1 && nmus==1 && nbm<=2 && njets>=5",himet_procs, plot_types);
  // pm.Push<Hist1D>(Axis(13, 25, 1000, "mj14", "M_{J} [GeV]",{250.,400.}),"met>200&&met<350 && nels==2 && nbm<=2 && njets>=5",himet_procs, plot_types);
  // pm.Push<Hist1D>(Axis(13, 25, 1000, "mj14", "M_{J} [GeV]",{250.,400.}),"met>200&&met<350 && nmus==2 && nbm<=2 && njets>=5",himet_procs, plot_types);
  // pm.Push<Hist1D>(Axis(13, 25, 1000, "mj14", "M_{J} [GeV]",{250.,400.}),"met>200&&met<350 && nleps==2 && nbm<=2 && njets>=5",himet_procs, plot_types);

  pm.Push<Hist1D>(Axis(10,0.,280., "mt", "m_{T} [GeV]",{140.}),"met>200&&met<350 && nels==1 && nmus==1 && nbm<=2 && njets>=5",himet_procs, plot_types);
  pm.Push<Hist1D>(Axis(10,0.,280., "mt", "m_{T} [GeV]",{140.}),"mj14>400 && met>200&&met<350 && nels==1 && nmus==1 && nbm<=2 && njets>=5",himet_procs, plot_types);
  // pm.Push<Hist1D>(Axis(6,0.,280., "mt", "m_{T} [GeV]",{140.}),"met>200&&met<350 && nels==2 && nbm<=2 && njets>=5",himet_procs, plot_types);
  // pm.Push<Hist1D>(Axis(6,0.,280., "mt", "m_{T} [GeV]",{140.}),"met>200&&met<350 && nmus==2 && nbm<=2 && njets>=5",himet_procs, plot_types);
  // pm.Push<Hist1D>(Axis(6,0.,280., "mt", "m_{T} [GeV]",{140.}),"met>200&&met<350 && nleps==2 && nbm<=2 && njets>=5",himet_procs, plot_types);

  pm.Push<Hist1D>(Axis(10,0.,400., "mus_pt[0]", "p_{T} (#mu) [GeV]"),"met>200&&met<350 && nleps==2 && nmus>=1 && nbm<=2 && njets>=5",himet_procs, plot_types);
  pm.Push<Hist1D>(Axis(10,0.,400., "els_pt[0]", "p_{T} (e)[GeV]"),"met>200&&met<350 && nleps==2 && nels>=1 && nbm<=2 && njets>=5",himet_procs, plot_types);
  pm.Push<Hist1D>(Axis(10,0.,400., "mus_pt[0]", "p_{T} (#mu) [GeV]"),"mj14>400 && met>200&&met<350 && nleps==2 && nmus>=1 && nbm<=2 && njets>=5",himet_procs, plot_types);
  pm.Push<Hist1D>(Axis(10,0.,400., "els_pt[0]", "p_{T} (e)[GeV]"),"mj14>400 && met>200&&met<350 && nleps==2 && nels>=1 && nbm<=2 && njets>=5",himet_procs, plot_types);

  // ------------------- MET 150-500 ----------------------
  pm.Push<Hist1D>(Axis(12,200.,500., "met", "MET[GeV]",{350.}),"met>200 && nels==1 && nmus==1 && nbm<=2 && njets>=5",allmet_procs, plot_types);
  pm.Push<Hist1D>(Axis(12,200.,500., "met", "MET[GeV]",{350.}),"mj14>400 && met>200 && nels==1 && nmus==1 && nbm<=2 && njets>=5",allmet_procs, plot_types);
  // pm.Push<Hist1D>(Axis(12,200.,500., "met", "MET[GeV]",{350.}),"met>200 && nleps==2 && nels>=1 && nbm<=2 && njets>=5",allmet_procs, plot_types);
  // pm.Push<Hist1D>(Axis(12,200.,500., "met", "MET[GeV]",{350.}),"met>200 && nleps==2 && nmus>=1 && nbm<=2 && njets>=5",allmet_procs, plot_types);

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
