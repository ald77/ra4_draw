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

int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  PlotOpt lin_shapes_info("txt/plot_styles.txt", "CMSPaper");
  lin_shapes_info.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::linear)
    .Stack(StackType::shapes);
  vector<PlotOpt> all_shapes = {lin_shapes_info};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/";

  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({"*_TTJets*Lept*.root", "*_TTJets_HT*.root", "*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*_ttHJetTobb*.root","*_TTTT*.root", "*_ST_*.root"});
  mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});

  // Baseline definitions
  NamedFunc wgt = "weight" * Higfuncs::eff_higtrig;
  string base_func("njets>=4 && njets<=5 && nvleps==0 && nbt>=2 && weight<1"); //met/met_calo
  
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_full>("t#bar{t}+X", 
    Process::Type::background, colors("tt_1l"),    attach_folder(foldermc,mctags["ttx"]),     base_func+"&& pass && stitch"));
  procs.push_back(Process::MakeShared<Baby_full>("V+jets",     
    Process::Type::background, kOrange+1,          attach_folder(foldermc,mctags["vjets"]),   base_func+"&& pass && stitch"));
 
  PlotMaker pm;

  vector<string> metcuts;
  metcuts.push_back("met>150");
  metcuts.push_back("met>100 && met<=150");
  metcuts.push_back("met>150 && met<=200");
  metcuts.push_back("met>200 && met<=300");
  metcuts.push_back("met>300");
  
  map<string, string> xcuts; // useful additional cut definitions
  // xcuts["nm1"] = base_func;
  // xcuts["base"] = base_func+"&& hig_dm<=40 && hig_am<200";
  xcuts["drmax"] = base_func+"&& hig_dm<=40 && hig_am<200 && hig_drmax<=2.2";
  
  string tmp_seln;
  for (auto &ixcut: xcuts) {
    for(unsigned imet(0); imet<metcuts.size(); imet++) { 
      pm.Push<Hist1D>(Axis(12,0,240,"hig_am", "<m> [GeV]", {100., 140.}),
        ixcut.second+"&&"+metcuts[imet], 
        procs, all_shapes).Weight(wgt);
    }
  }

  pm.min_print_ = true;
  pm.MakePlots(40.);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

