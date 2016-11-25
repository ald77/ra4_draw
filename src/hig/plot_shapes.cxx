///// plot_closure: Compares 2b, 3b, and 4b Higgs distributions

#include <fstream>
#include <iostream>
#include <vector>
#include <chrono>

#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>

#include "TError.h" // Controls error level reporting

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/hist1d.hpp"
#include "core/palette.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  // pick only specific breakdowns to reduce output; possible options: 
  // {"procs","all_btag", "tt_ntrub","tt_btag", "tt_1lntrub","tt_1lbtag","vjets_2lntrub"};
  vector<string> proclist = {"all_btag"};
}

int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info);
  vector<PlotOpt> all_plot_types = {lin_shapes_info};
  Palette colors("txt/colors.txt", "default");


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////        Skims        //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/");

  set<string> filetags = {"*_TTJets*Lept*.root", "*_TTJets_HT*.root", "*_TTZ*.root", "*_TTW*.root",
                          "*_TTGJets*.root", "*_ttHJetTobb*.root","*_TTTT*.root",
                          "*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root",
                          "*_ST_*.root",
                          // "*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root",
                          "*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                          "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"};
  set<string> allfiles = attach_folder(foldermc, filetags);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////     Defining cuts   //////////////////////////////////////////
  // Cuts in baseline speed up the yield finding
  string baseline = "pass && stitch && weight<1 && njets>=4 && njets<=5";
  NamedFunc basefunc(baseline);

  ////// Nb cuts
  string c_2b="nbt==2&&nbm==2";
  string c_3b="nbt>=2&&nbm==3&&nbl==3";
  string c_4b="nbt>=2&&nbm>=3&&nbl>=4";

  NamedFunc ntrub("ntrub",[](const Baby &b) -> NamedFunc::ScalarType{
    int tmp_ntrub(0);
    for (unsigned i(0); i<b.jets_pt()->size(); i++){
      if (!b.jets_h1()->at(i) && !b.jets_h2()->at(i)) continue;
      if (b.jets_hflavor()->at(i)==5) tmp_ntrub++;
    }
    return tmp_ntrub;
  });

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  map<string, vector<shared_ptr<Process> > > procs;

  // comparison by process
  //---------------------------------------------------
  // procs["procs"] = vector<shared_ptr<Process> >();
  // procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (2b)", Process::Type::background, kBlack,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root", 
  //   foldermc+"*_TTZ*.root", foldermc+"*_TTW*.root", foldermc+"*_TTGJets*.root"
  // }, baseline+"&& nvleps==0 &&"+c_2b+"&& ntruleps>=1"));
  // procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (3b)", Process::Type::background, colors("tt_1l"),
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root", 
  //   foldermc+"*_TTZ*.root", foldermc+"*_TTW*.root", foldermc+"*_TTGJets*.root"
  // }, baseline+"&& nvleps==0 &&"+c_3b+"&& ntruleps>=1"));
  // procs["procs"].push_back(Process::MakeShared<Baby_full>("V+jets (3b)", Process::Type::background, kOrange+1,
  //   {foldermc+"*_ZJet*.root", foldermc+"*_WJetsToLNu*.root"}, 
  //   baseline+"&& nvleps==0 &&"+c_3b));
  // procs["procs"].push_back(Process::MakeShared<Baby_full>("Single t (3b)", Process::Type::background, colors("single_t"),
  //   {foldermc+"*_ST_*.root"}, 
  //   baseline+"&& nvleps==0 &&"+c_3b));
  // procs["procs"].push_back(Process::MakeShared<Baby_full>("Other (3b)", Process::Type::background, kPink-2,
  //   {foldermc+"*DYJetsToLL*.root", foldermc+"*_TTTT*.root", foldermc+"*_ttHJetTobb*.root", 
  //   foldermc+"*_WH_HToBB*.root", foldermc+"*_ZH_HToBB*.root", foldermc+"*_WWTo*.root", 
  //   foldermc+"*_WZ*.root", foldermc+"*_ZZ_*.root",
  //   foldermc+"*QCD_HT*0_Tune*.root", foldermc+"*QCD_HT*Inf_Tune*.root"},
    // baseline+"&& nvleps==0 &&"+c_3b));

  // comparison by b-tag category for ttbar
  //---------------------------------------------------
  procs["all_btag"] = vector<shared_ptr<Process> >();
  procs["all_btag"].push_back(Process::MakeShared<Baby_full>("All bkg (2b)", Process::Type::background, kBlack,
                              allfiles, baseline+"&& nvleps==0 &&"+c_2b));
  procs["all_btag"].push_back(Process::MakeShared<Baby_full>("All bkg (3b)", Process::Type::background, kAzure+1,
                              allfiles, baseline+"&& nvleps==0 &&"+c_3b));
  procs["all_btag"].push_back(Process::MakeShared<Baby_full>("All bkg (4b)", Process::Type::background, kPink+2,
                              allfiles, baseline+"&& nvleps==0 &&"+c_4b));

  // comparison by b-tag category for ttbar
  //---------------------------------------------------
  // procs["tt_btag"] = vector<shared_ptr<Process> >();
  // procs["tt_btag"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (2b)", Process::Type::background, kBlack,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"
  // }, baseline+"&& nvleps==0 && ntruleps>=1 &&"+c_2b));
  // procs["tt_btag"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (3b)", Process::Type::background, kAzure+1,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"
  // }, baseline+"&& nvleps==0 && ntruleps>=1 &&"+c_3b));
  // procs["tt_btag"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (4b)", Process::Type::background, kPink+2,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"
  // }, baseline+"&& nvleps==0 && ntruleps>=1 &&"+c_4b));

  // comparison by # true b-tags category for ttbar
  //---------------------------------------------------
  // procs["tt_ntrub"] = vector<shared_ptr<Process> >();
  // procs["tt_ntrub"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (2 true b)", Process::Type::background, kOrange-4,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"
  // }, basefunc && "ntruleps>=1" && ntrub==2));
  // procs["tt_ntrub"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (3 true b)", Process::Type::background, kTeal-8,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"
  // }, basefunc && "ntruleps>=1" && ntrub==3));
  // procs["tt_ntrub"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (#geq4 true b)", Process::Type::background, kAzure-4,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"
  // }, basefunc && "ntruleps>=1" && ntrub>=4));

  // comparison 0 to 1 lep by b-tag category for ttbar
  //---------------------------------------------------
  // procs["tt_1lbtag"] = vector<shared_ptr<Process> >();
  // procs["tt_1lbtag"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (2b)", Process::Type::background, kBlack,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"
  // }, baseline+"&& nvleps==0 && ntruleps>=1 &&"+c_2b));
  // procs["tt_1lbtag"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (3b)", Process::Type::background, kAzure+1,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"
  // }, baseline+"&& nvleps==0 && ntruleps>=1 &&"+c_3b));
  // procs["tt_1lbtag"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (4b)", Process::Type::background, kPink+2,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"
  // }, baseline+"&& nvleps==0 && ntruleps>=1 &&"+c_4b));
  // procs["tt_1lbtag"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (1l, 2b)", Process::Type::background, kBlack,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"
  // }, baseline+"&& nvleps==1 && ntruleps>=1 &&"+c_2b)); procs["tt_1lbtag"].back()->SetLineStyle(2);
  // procs["tt_1lbtag"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (1l, 3b)", Process::Type::background, kAzure+1,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"
  // }, baseline+"&& nvleps==1 && ntruleps>=1 &&"+c_3b)); procs["tt_1lbtag"].back()->SetLineStyle(2);
  // procs["tt_1lbtag"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (1l, 4b)", Process::Type::background, kPink+2,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"
  // }, baseline+"&& nvleps==1 && ntruleps>=1 &&"+c_4b)); procs["tt_1lbtag"].back()->SetLineStyle(2);

  // comparison by # true b-tags category for ttbar
  //---------------------------------------------------
  // procs["tt_1lntrub"] = vector<shared_ptr<Process> >();
  // procs["tt_1lntrub"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (2 true b)", Process::Type::background, kOrange-4,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"
  // }, basefunc && "nvleps==0 && ntruleps>=1" && ntrub==2));
  // procs["tt_1lntrub"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (3 true b)", Process::Type::background, kTeal-8,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"
  // }, basefunc && "nvleps==0 && ntruleps>=1" && ntrub==3));
  // procs["tt_1lntrub"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (#geq4 true b)", Process::Type::background, kAzure-4,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"
  // }, basefunc && "nvleps==0 && ntruleps>=1" && ntrub>=4));
  // procs["tt_1lntrub"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (1l, 2 true b)", Process::Type::background, kOrange-4,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"
  // }, basefunc && "nvleps==1 && ntruleps>=1" && ntrub==2)); procs["tt_1lntrub"].back()->SetLineStyle(2);
  // procs["tt_1lntrub"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (1l, 3 true b)", Process::Type::background, kTeal-8,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"
  // }, basefunc && "nvleps==1 && ntruleps>=1" && ntrub==3)); procs["tt_1lntrub"].back()->SetLineStyle(2);
  // procs["tt_1lntrub"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X (1l, #geq4 true b)", Process::Type::background, kAzure-4,
  //   {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"
  // }, basefunc && "nvleps==1 && ntruleps>=1" && ntrub>=4)); procs["tt_1lntrub"].back()->SetLineStyle(2);

  // comparison by # true b-tags category for vjets
  //---------------------------------------------------
  // procs["vjets_2lntrub"] = vector<shared_ptr<Process> >();
  // procs["vjets_2lntrub"].push_back(Process::MakeShared<Baby_full>("V+jets (2 true b)", Process::Type::background, kOrange,
  //   {foldermc+"*_ZJet*.root", foldermc+"*_WJetsToLNu*.root", foldermc+"*_DYJetsToLL*.root"
  // }, basefunc && "nvleps==0" && ntrub==2));
  // procs["vjets_2lntrub"].push_back(Process::MakeShared<Baby_full>("V+jets (3 true b)", Process::Type::background, kOrange+7,
  //   {foldermc+"*_ZJet*.root", foldermc+"*_WJetsToLNu*.root", foldermc+"*_DYJetsToLL*.root"
  // }, basefunc && "nvleps==0" && ntrub==3));
  // procs["vjets_2lntrub"].push_back(Process::MakeShared<Baby_full>("V+jets (#geq4 true b)", Process::Type::background, kRed+1,
  //   {foldermc+"*_ZJet*.root", foldermc+"*_WJetsToLNu*.root", foldermc+"*_DYJetsToLL*.root"
  // }, basefunc && "nvleps==0" && ntrub>=4));
  // procs["vjets_2lntrub"].push_back(Process::MakeShared<Baby_full>("V+jets (2l, 2 true b)", Process::Type::background, kOrange,
  //   {foldermc+"*_ZJet*.root", foldermc+"*_WJetsToLNu*.root", foldermc+"*_DYJetsToLL*.root"
  // }, basefunc && "nvleps==2 && (mumu_m*(mumu_m>0)+elel_m*(elel_m>0))>80&&(mumu_m*(mumu_m>0)+elel_m*(elel_m>0))<100" && ntrub==2)); procs["vjets_2lntrub"].back()->SetLineStyle(2);
  // procs["vjets_2lntrub"].push_back(Process::MakeShared<Baby_full>("V+jets (2l, 3 true b)", Process::Type::background, kOrange+7,
  //   {foldermc+"*_ZJet*.root", foldermc+"*_WJetsToLNu*.root", foldermc+"*_DYJetsToLL*.root"
  // }, basefunc && "nvleps==2 && (mumu_m*(mumu_m>0)+elel_m*(elel_m>0))>80&&(mumu_m*(mumu_m>0)+elel_m*(elel_m>0))<100" && ntrub==3)); procs["vjets_2lntrub"].back()->SetLineStyle(2);
  // procs["vjets_2lntrub"].push_back(Process::MakeShared<Baby_full>("V+jets (2l, #geq4 true b)", Process::Type::background, kRed+1,
  //   {foldermc+"*_ZJet*.root", foldermc+"*_WJetsToLNu*.root", foldermc+"*_DYJetsToLL*.root"
  // }, basefunc && "nvleps==2 && (mumu_m*(mumu_m>0)+elel_m*(elel_m>0))>80&&(mumu_m*(mumu_m>0)+elel_m*(elel_m>0))<100" && ntrub>=4)); procs["vjets_2lntrub"].back()->SetLineStyle(2);

  vector<TString> metcuts;
  metcuts.push_back("met>150");
  metcuts.push_back("met>150 && met<=200");
  metcuts.push_back("met>200 && met<=300");
  metcuts.push_back("met>300");

  vector<TString> higcuts;
  higcuts.push_back("1");
  higcuts.push_back("ntks==0 && !low_dphi && hig_drmax<2.2");

  //// Adding plots
  PlotMaker pm;

  // TString base_pretty(baseline); base_pretty.ReplaceAll("pass &&","").ReplaceAll("stitch &&","").ReplaceAll("&& weight<1","");

  for (auto &ipr: proclist){
    if (procs.find(ipr)==procs.end()) {cout<<"Breakdown "<<ipr<<" not coded in."<<endl; exit(0);}
    for (auto &imet: metcuts){
      for (auto &ihig: higcuts){
        pm.Push<Hist1D>(Axis(11,0,110,"hig_dm", "#Delta m_{jj} [GeV]", {40.}),
                        baseline+"&&"+ihig+"&&"+imet, procs[ipr], all_plot_types).Tag(ipr);
        pm.Push<Hist1D>(Axis(12,0,240,"hig_am", "<m_{jj}> [GeV]", {100., 140.}),
                        baseline+"&&"+ihig+"&&"+imet, procs[ipr], all_plot_types).Tag(ipr);
        if (!ihig.Contains("hig_drmax"))
          pm.Push<Hist1D>(Axis(25,0,5,"hig_drmax", "#Delta R_{max}", {2.2}),
                        baseline+"&&"+ihig+"&&"+imet, procs[ipr], all_plot_types).Tag(ipr);
      }
    }
  }

  pm.min_print_ = true;
  pm.MakePlots(40.);

  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Making plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
} // main

