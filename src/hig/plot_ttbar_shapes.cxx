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

  PlotOpt lin_shapes_info("txt/plot_styles.txt", "CMSPaper");
  lin_shapes_info.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::linear)
    .Stack(StackType::shapes);
  // PlotOpt log_shapes_info = lin_shapes_info.YAxis(YAxisType::log);
  vector<PlotOpt> plt_types = {lin_shapes_info};

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_met100nj4/";
  set<string> allfiles = {foldermc+"*_TTJets*SingleLept*.root",foldermc+"*_TTJets*DiLept*.root"};

  // Baseline definitions
  string baseline("njets>=4 && njets<=5 && met/met_calo<5");

  ////// Nb cuts
  vector<string> nbcuts;
  nbcuts.push_back("nbt==2&&nbm==2");
  nbcuts.push_back("nbt>=2&&nbm==3&&nbl==3");
  nbcuts.push_back("nbt>=2&&nbm>=3&&nbl>=4");

  string samplename = "t#bar{t}+X";
  // vector<int> colors = {kGreen+3, kGreen+1, kOrange, kAzure+1, kBlue+1};
  vector<int> colors = {kOrange, kAzure+1, kBlue+1};
  vector<shared_ptr<Process> > procs = vector<shared_ptr<Process> >();
  for (unsigned inb(0); inb<nbcuts.size(); inb++){
    procs.push_back(Process::MakeShared<Baby_full>(samplename+" ("+RoundNumber(inb+2,0).Data()+"b)", 
      Process::Type::background, colors[inb], allfiles, baseline+"&& pass &&"+nbcuts[inb]));
  }
  vector<int> colors_trub = {kAzure-4, kTeal-8, kOrange-4, kPink+2, kMagenta-1};
  vector<shared_ptr<Process> > procs_trub = vector<shared_ptr<Process> >();
  for (unsigned inb(2); inb<5; inb++){
      procs_trub.push_back(Process::MakeShared<Baby_full>(samplename+" ("+RoundNumber(inb,0).Data()+" B-hadrons)", 
        Process::Type::background, colors_trub[inb], allfiles, Higfuncs::ntrub==inb && baseline+"&& pass"));
  }

  vector<string> metcuts;
  metcuts.push_back("met>100");
  // metcuts.push_back("met>100 && met<=150");
  // metcuts.push_back("met>150 && met<=200");
  // metcuts.push_back("met>200 && met<=300");
  // metcuts.push_back("met>300");
  
  vector<string> xcuts;
  xcuts.push_back("1");
  xcuts.push_back("hig_drmax<2.2");
  xcuts.push_back("hig_dm<40 && hig_am<200");
  xcuts.push_back("hig_dm<40 && hig_am<200 && hig_drmax<2.2");

  PlotMaker pm;
  for (unsigned imet(0); imet<metcuts.size(); imet++){
    for (unsigned ic(0); ic<xcuts.size(); ic++){
      pm.Push<Hist1D>(Axis(10,0,200,"hig_am", "<m> [GeV]", {100., 140.}),
        baseline+"&& ntrutaush==1 && nvleps==0 &&"+xcuts[ic]+"&&"+metcuts[imet], procs, plt_types).Tag(sample+"0l_shape_bcats");
      pm.Push<Hist1D>(Axis(10,0,200,"hig_am", "<m> [GeV]", {100., 140.}),
        baseline+"&& ntrutaush==1 && nvleps==0 &&"+xcuts[ic]+"&&"+metcuts[imet], procs_trub, plt_types).Tag(sample+"0l_shape_trub");
      pm.Push<Hist1D>(Axis(10,0,200,"hig_am", "<m> [GeV]", {100., 140.}),
        baseline+"&& ntruleps==1 && nleps==1 &&"+xcuts[ic]+"&&"+metcuts[imet], procs, plt_types).Tag(sample+"1l_shape_bcats");
      pm.Push<Hist1D>(Axis(10,0,200,"hig_am", "<m> [GeV]", {100., 140.}),
        baseline+"&& ntruleps==1 && nleps==1 &&"+xcuts[ic]+"&&"+metcuts[imet], procs_trub, plt_types).Tag(sample+"1l_shape_trub");
      pm.Push<Hist1D>(Axis(10,0,200,"hig_am", "<m> [GeV]", {100., 140.}),
        baseline+"&& ntruleps==2 && nleps==2 &&"+xcuts[ic]+"&&"+metcuts[imet], procs, plt_types).Tag(sample+"2l_shape_bcats");
      pm.Push<Hist1D>(Axis(10,0,200,"hig_am", "<m> [GeV]", {100., 140.}),
        baseline+"&& ntruleps==2 && nleps==2 &&"+xcuts[ic]+"&&"+metcuts[imet], procs_trub, plt_types).Tag(sample+"2l_shape_trub");
    }
  }

  pm.min_print_ = true;
  pm.MakePlots(40.);

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
