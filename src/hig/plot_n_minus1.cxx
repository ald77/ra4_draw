// study differences in mjj as a function of b-cat for single vs dilepton ttbar

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/plot_opt.hpp"
#include "core/functions.hpp"
#include "hig/hig_functions.hpp"

using namespace std;
using namespace PlotOptTypes;

void GetOptions(int argc, char *argv[]);

namespace{
  string sample = "tt";
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  Palette colors("txt/colors.txt", "default");


  /////////////////// PLOT STYLES //////////////////////////////////////
  PlotOpt lin_norm_info("txt/plot_styles.txt", "CMSPaper");
  lin_norm_info.Title(TitleType::info)   
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::linear)
    .Stack(StackType::data_norm);
  PlotOpt log_norm_info = lin_norm_info.YAxis(YAxisType::log);
  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};

  PlotOpt log_norm = lin_norm_info.YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.7).Bottom(BottomType::off);
  PlotOpt lin_norm = lin_norm_info.YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::off);
  vector<PlotOpt> plt_norm = {lin_norm, log_norm};

  PlotOpt lin_shapes = lin_norm.Stack(StackType::shapes).Bottom(BottomType::ratio);
  vector<PlotOpt> plt_shapes = {lin_shapes};



  //////////////////////////////////// PROCESSES /////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////// Higgsino //////////////////////////////////////////////
  //string folderhigmc = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higlep1/";
  string folderhigmc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/";
  string folderhigdata = bfolder+"/cms2r0/babymaker/babies/2017_01_27/data/merged_higdata_higlep1/";
  string foldersig(bfolder+"/cms2r0/babymaker/babies/2016_08_10/TChiHH/merged_higmc_higloose/");

  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({"*_TTJets*Lept*.root", "*_TTJets_HT*.root", "*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*_ttHJetTobb*.root","*_TTTT*.root"});
  mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
  mctags["singlet"] = set<string>({"*_ST_*.root"});
  mctags["qcd"]     = set<string>({"*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  set<string> allmctags;
  for (auto &iset: mctags) {
      allmctags.insert(iset.second.begin(), iset.second.end());
  }

  NamedFunc wgt = "weight"* Higfuncs::eff_higtrig;
  NamedFunc base_func = "pass && met/met_calo<5 && nbt>=2 && nvleps==0&&ntks==0&&!low_dphi&&met>150&&weight<.1";

  vector<shared_ptr<Process> > procs_hig;
  procs_hig.push_back(Process::MakeShared<Baby_full>("t#bar{t}+X", 
    Process::Type::background, colors("tt_1l"),    attach_folder(folderhigmc,mctags["ttx"]),     base_func&&"stitch"));
  procs_hig.push_back(Process::MakeShared<Baby_full>("V+jets",     
    Process::Type::background, kOrange+1,          attach_folder(folderhigmc,mctags["vjets"]),   base_func&&"stitch"));
  procs_hig.push_back(Process::MakeShared<Baby_full>("Single t",   
    Process::Type::background, colors("single_t"), attach_folder(folderhigmc,mctags["singlet"]), base_func&&"stitch"));
  procs_hig.push_back(Process::MakeShared<Baby_full>("QCD",        
    Process::Type::background, colors("other"),    attach_folder(folderhigmc,mctags["qcd"]),     base_func&&"stitch")); 
  procs_hig.push_back(Process::MakeShared<Baby_full>("Other",      
    Process::Type::background, kGreen+1,           attach_folder(folderhigmc,mctags["other"]),   base_func&&"stitch"));      

  vector<string> sigm = {"225","400","700"}; 
  vector<int> sig_colors = {kGreen, kRed, kBlue}; // need sigm.size() >= sig_colors.size()
  for (unsigned isig(0); isig<sigm.size(); isig++)
    procs_hig.push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
			sig_colors[isig], {foldersig+"*TChiHH_mGluino-"+sigm[isig]+"*.root"}, base_func));

  // procs_hig.push_back(Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
  //   {folderhigdata+"*root"},  Higfuncs::trig_hig>0. && base_func)); 

  ///////////////////////////////////////////////////////////////////////////////////////////////////


  string cuts = "nbm>=2";
  PlotMaker pm;
  cuts = "nbm>=3&&nbl>=4&&met>150&&hig_drmax<2.2&&hig_dm<40&& njets>=4 && njets<=5 ";
  pm.Push<Hist1D>(Axis(50, 0., 250., "hig_am", "#LTm#GT [GeV]", {100., 140.}),cuts, procs_hig, plt_norm).Tag("hig");
  cuts = "nbm>=3&&nbl>=4&&met>150&&hig_drmax<2.2&&hig_dm<40";
  pm.Push<Hist1D>(Axis(6, 3.5, 9.5, "njets", "N_{jets}", {5.5}),cuts, procs_hig, plt_norm).Tag("hig");

  cuts = "nbm>=3&&nbl>=4&&met>150&&hig_drmax<2.2&& njets>=4 && njets<=5 ";
  pm.Push<Hist1D>(Axis(30, 0., 150., "hig_dm", "#Deltam [GeV]", {40.}),cuts, procs_hig, plt_norm).Tag("hig");

  cuts = "nbm>=3&&nbl>=4&&met>150&&hig_drmax<2.2";
  pm.Push<Hist1D>(Axis(9, 150., 600., "met", "E^{miss}_{T} [GeV]", {150., 200., 300., 450.}),cuts, procs_hig, plt_norm).Tag("hig");
  // cuts = "mt<100 && met>150 && nleps==11";
  // pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nbm", "CSV N_{b,M}"),cuts, procs_hig, plt_norm_info).Tag("hig");
  // pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nbt", "CSV N_{b,T}"),cuts, procs_hig, plt_norm_info).Tag("hig");
  // pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nbdm", "DeepCSV N_{b,M}"),cuts, procs_hig, plt_norm_info).Tag("hig");
  // pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nbdt", "DeepCSV N_{b,T}"),cuts, procs_hig, plt_norm_info).Tag("hig");

  pm.min_print_ = true;
  pm.MakePlots(36.8);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
} // main

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"sample", required_argument, 0, 's'},    
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      sample = optarg;
      break;
    case 0:
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
