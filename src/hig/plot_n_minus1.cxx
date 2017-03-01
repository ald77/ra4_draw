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
using namespace Higfuncs;

void GetOptions(int argc, char *argv[]);

namespace{
  float lumi = 35.9;
  vector<string> sigm = {"225","400","700"}; 
  vector<int> sig_colors = {kGreen+1, kRed, kBlue}; 
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
  string folderhigmc = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higloose/";
  // to plot njets
  folderhigmc = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higlooser/";
  string folderhigdata = bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_higdata_higlep1/";
  string foldersig(bfolder+"/cms2r0/babymaker/babies/2017_01_27/TChiHH/merged_higmc_unskimmed/");

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

  string c_ps = "stitch && pass && pass_ra2_badmu && met/met_calo<5 && weight<1";

  vector<shared_ptr<Process> > procs_hig;
  procs_hig.push_back(Process::MakeShared<Baby_full>("t#bar{t}+X", 
    Process::Type::background, colors("tt_1l"),    attach_folder(folderhigmc,mctags["ttx"]),    c_ps));
  procs_hig.push_back(Process::MakeShared<Baby_full>("V+jets",     
    Process::Type::background, kOrange+1,          attach_folder(folderhigmc,mctags["vjets"]),  c_ps));
  procs_hig.push_back(Process::MakeShared<Baby_full>("Single t",   
    Process::Type::background, colors("single_t"), attach_folder(folderhigmc,mctags["singlet"]),c_ps));
  procs_hig.push_back(Process::MakeShared<Baby_full>("QCD",        
    Process::Type::background, colors("other"),    attach_folder(folderhigmc,mctags["qcd"]),    c_ps)); 
  procs_hig.push_back(Process::MakeShared<Baby_full>("Other",      
    Process::Type::background, kGreen+1,           attach_folder(folderhigmc,mctags["other"]),  c_ps));      

  for (unsigned isig(0); isig<sigm.size(); isig++)
    procs_hig.push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
			sig_colors[isig], {foldersig+"*TChiHH_mGluino-"+sigm[isig]+"*.root"}, 
      "pass_goodv&&pass_ecaldeadcell&&pass_hbhe&&pass_hbheiso&&pass_fsmet"));

  // procs_hig.push_back(Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
  //   {folderhigdata+"*root"},  trig_hig>0. && base_func)); 

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  string uu = "&&";
  string baseline = "njets>=4 && nbdt>=2 && nvleps==0 && ntks==0 && !low_dphi && met>150";

  string cnj = "njets<=5";

  string c2b = "nbdt==2 && nbdm==2";
  string c3b = "nbdt>=2 && nbdm==3 && nbdl==3";
  string c4b = "nbdt>=2 && nbdm>=3 && nbdl>=4";
  string c3bp = "nbdt>=2 && nbdm>=3";

  string hdrmax = "hig_drmax<2.2";
  string hdm = "hig_dm<40";
  string htrim = "hig_dm<40 && hig_am<=200";
  string hig = "higd_drmax<=2.2 && higd_am<=200 && higd_dm <= 40 && (higd_am>100 && higd_am<=140)";
  string sbd = "higd_drmax<=2.2 && higd_am<=200 && higd_dm <= 40 && !(higd_am>100 && higd_am<=140)";

  //to avoid warnings for unused variables...
  string dummy=c2b;dummy=c3b;dummy=c4b;dummy=hdrmax;dummy=hdm;dummy=htrim;dummy=hig;dummy=sbd;

  NamedFunc wgt = weight_higd * eff_higtrig;

  PlotMaker pm;

  pm.Push<Hist1D>(Axis(40, 0., 200., "hig_am", "#LTm#GT [GeV]", {100., 140.}),
    baseline+uu+c4b+uu+cnj+uu+hdrmax+uu+hdm, procs_hig, plt_norm).Weight(wgt).Tag("hig");
  
  pm.Push<Hist1D>(Axis(6, 3.5, 9.5, "njets", "N_{jets}", {5.5}),
    baseline+uu+c3bp+uu+hdrmax+uu+hdm+uu+"met>300", procs_hig, plt_norm).Weight(wgt).Tag("hig");

  pm.Push<Hist1D>(Axis(30, 0., 150., "hig_dm", "#Deltam [GeV]", {40.}),
    baseline+uu+c4b+uu+cnj+uu+hdrmax, procs_hig, plt_norm).Weight(wgt).Tag("hig");


  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
} // main

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"lumi", required_argument, 0, 'l'},    
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "l:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'l':
      lumi = atof(optarg);
      break;
    case 0:
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
