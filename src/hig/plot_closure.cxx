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

using namespace std;
using namespace PlotOptTypes;

int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/");

  // Cuts in baseline speed up the yield finding
  string baseline = "pass && stitch && hig_drmax<2.2&&ntks==0&&njets>=4&&njets<=5&&!low_dphi&&nvleps==0";
  ////// Nb cuts
  string c_2b="nbt==2&&nbm==2";
  string c_3b="nbt>=2&&nbm==3&&nbl==3";
  string c_4b="nbt>=2&&nbm>=3&&nbl>=4";

  auto proc_tt_2b = Process::MakeShared<Baby_full>("t#bar{t}, 2b", Process::Type::background, kBlue,
    {foldermc+"*_TTJets*.root"}, baseline+"&&"+c_2b);
  auto proc_tt_3b = Process::MakeShared<Baby_full>("t#bar{t}, 3b", Process::Type::background, kRed+1,
    {foldermc+"*_TTJets*.root"}, baseline+"&&"+c_3b);
  auto proc_tt_4b = Process::MakeShared<Baby_full>("t#bar{t}, 4b", Process::Type::background, kGreen+1,
    {foldermc+"*_TTJets*.root"}, baseline+"&&"+c_4b);


  //// All processes combined in bkg
  vector<string> vnames_bkg({"_TTJets", "_WJetsToLNu", "_TTW", "_TTZ", "DYJetsToLL", 
	"_ZJet", "_ttHJetTobb", "_TTGJets", "_TTTT", 
	"_WH_HToBB", "_ZH_HToBB", "_WWTo", "_WZ", "_ZZ_"});//, 
	//"QCD_HT*to1000_Tune", "QCD_HT*to1500_Tune", "QCD_HT*to2000_Tune", "QCD_HT*Inf_Tune"});
  set<string> names_bkg;
  for(auto name : vnames_bkg)
    names_bkg.insert(name = foldermc + "*" + name + "*.root");

  auto proc_bkg_2b = Process::MakeShared<Baby_full>("All bkg (no QCD), 2b", Process::Type::background, kBlue,
    names_bkg, baseline+"&&"+c_2b);
  auto proc_bkg_3b = Process::MakeShared<Baby_full>("All bkg (no QCD), 3b", Process::Type::background, kOrange+1,
    names_bkg, baseline+"&&"+c_3b);
  auto proc_bkg_4b = Process::MakeShared<Baby_full>("All bkg (no QCD), 4b", Process::Type::background, kGreen+2,
    names_bkg, baseline+"&&"+c_4b);


  vector<shared_ptr<Process> > procs_tt = {proc_tt_2b, proc_tt_3b, proc_tt_4b};
  vector<shared_ptr<Process> > procs_bkg = {proc_bkg_2b, proc_bkg_3b, proc_bkg_4b};

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

  //// Adding plots
  PlotMaker pm;

  NamedFunc cuts = "nvleps==0&&met>100&&!low_dphi&&hig_drmax<2.2&&ntks==0&&njets>=4&&njets<=5";
  pm.Push<Hist1D>(Axis(22,0,110,"hig_dm", "#Deltam_{jj} [GeV]", {40.}),
                  cuts, procs_tt, all_plot_types).Tag("tt");
  pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m_{jj}> [GeV]", {100., 140.}),
                  cuts, procs_tt, all_plot_types).Tag("tt");
  pm.Push<Hist1D>(Axis(22,0,110,"hig_dm", "#Deltam_{jj} [GeV]", {40.}),
                  cuts, procs_bkg, all_plot_types).Tag("bkg");
  pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m_{jj}> [GeV]", {100., 140.}),
                  cuts, procs_bkg, all_plot_types).Tag("bkg");

  pm.min_print_ = true;
  pm.MakePlots(40.);

  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Making plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
} // main

