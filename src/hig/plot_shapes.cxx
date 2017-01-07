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
	bool do_allbkg = false;
  string sample = "search";
  string json = "full";
  bool unblind = false;
  bool do_note = false;
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

  string foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/";
  if (sample=="ttbar") foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higlep1/";
  if (sample=="zll") foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higlep2/";
  if (sample=="qcd") foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higqcd/";
  string folderdata = bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_higloose/";
  if (sample=="ttbar") folderdata = bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_higlep1/";
  if (sample=="zll") folderdata = bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_higlep2/";
  if (sample=="qcd") folderdata = bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_higqcd/";

  set<string> alltags; 
  if (sample=="ttbar" || sample=="search") alltags = {"*_TTJets*Lept*.root", "*_TTJets_HT*.root",
                                      "*_TTZ*.root", "*_TTW*.root", "*_TTGJets*.root", 
                                      "*_ttHJetTobb*.root","*_TTTT*.root"};
  if (sample=="zll") alltags = {"*DYJetsToLL*.root"};
  if (sample=="qcd") alltags = {//"*QCD_HT100to200_Tune*", "*QCD_HT200to300_Tune*",
                                //"*QCD_HT300to500_Tune*", 
                                 "*QCD_HT500to700_Tune*",
                                 "*QCD_HT700to1000_Tune*", "*QCD_HT1000to1500_Tune*", 
                                 "*QCD_HT1500to2000_Tune*", "*QCD_HT2000toInf_Tune*"};
  if (do_allbkg) { // don't include QCD MC unless in QCD control sample
    if(sample=="qcd") alltags = {"*_TTJets*Lept*.root", "*_TTJets_HT*.root", 
            "*_TTZ*.root", "*_TTW*.root", "*_TTGJets*.root", "*_ttHJetTobb*.root","*_TTTT*.root",
            "*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root", "*_ST_*.root",
            "*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root",
            "*_WH_HToBB*.root", "*_ZH_HToBB*.root", "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"};
    else  alltags = {"*_TTJets*Lept*.root", "*_TTJets_HT*.root", 
            "*_TTZ*.root", "*_TTW*.root", "*_TTGJets*.root", "*_ttHJetTobb*.root","*_TTTT*.root",
            "*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root", "*_ST_*.root",
            "*_WH_HToBB*.root", "*_ZH_HToBB*.root", "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"};
  }
  set<string> allfiles = attach_folder(foldermc,alltags);

  // Baseline definitions
  string baseline("njets>=4 && njets<=5 && met/met_calo<5"); //met/met_calo
  // zll skim: ((elel_m>80&&elel_m<100)||(mumu_m>80&&mumu_m<100)) && 
  // nleps==2 && nleps>=1 && Max$(leps_pt)>30 && njets>=4&&njets<=5
  if (sample=="zll") baseline = baseline+"&& nleps==2 && met<50";
  // qcd skim - met>150 && nvleps==0 && (njets==4||njets==5)
  if (sample=="qcd") baseline = baseline+"&& nvleps==0 && low_dphi";
  // ttbar skim - met>100 && nleps==1 && (njets==4||njets==5) && nbm>=2
  if (sample=="ttbar") baseline = baseline+"&& nleps==1 && mt<100";
  // search skim - met>100 && nvleps==0 && (njets==4||njets==5) && nbm>=2
  if (sample=="search") baseline = baseline+"&& nvleps==0";
    

  ////// Nb cuts
  vector<string> nbcuts;
  unsigned firstnb = 0;
  nbcuts.push_back("nbm==0");
  nbcuts.push_back("nbm==1");
  // nbcuts.push_back("nbt==2&&nbm==2");
  if (sample=="ttbar" || sample=="search") {
    firstnb = 2;
    nbcuts.push_back("nbt>=2&&nbm==3&&nbl==3");
    nbcuts.push_back("nbt>=2&&nbm>=3&&nbl>=4");
  } else if (sample=="qcd") {
    nbcuts.push_back("nbt>=2&&nbm>=3");
  }

  string samplename = "t#bar{t}+X";
  if (sample=="qcd") samplename = "QCD";
  if (sample=="zll") samplename = "Z#rightarrow ll";
  if (do_allbkg) samplename = "All bkg.";

  vector<int> colors = {kGreen+3, kGreen+1, kOrange, kAzure+1, kBlue+1};

  vector<shared_ptr<Process> > procs = vector<shared_ptr<Process> >();
  for (unsigned inb(firstnb); inb<nbcuts.size(); inb++){
    // if (sample=="qcd" && inb==nbcuts.size()-1) continue;
    procs.push_back(Process::MakeShared<Baby_full>(samplename+" ("+RoundNumber(inb,0).Data()+"b)", 
      Process::Type::background, colors[inb], allfiles, baseline+"&& pass && pass_ra2_badmu && stitch &&"+nbcuts[inb]));
  }
  vector<int> colors_trub = {kAzure-4, kTeal-8, kOrange-4, kPink+2, kMagenta-1};
  vector<shared_ptr<Process> > procs_trub = vector<shared_ptr<Process> >();
  for (unsigned inb(firstnb); inb<nbcuts.size(); inb++){
    if ((sample=="zll" || sample=="qcd") && inb==nbcuts.size()-1) { // merge 4b into 3b
      procs_trub.push_back(Process::MakeShared<Baby_full>(samplename+" (#geq"+RoundNumber(inb,0).Data()+" B-hadrons)", 
        Process::Type::background, colors_trub[inb], allfiles, Higfuncs::ntrub>=inb && baseline+"&& pass && pass_ra2_badmu && stitch"));
    } else {
      procs_trub.push_back(Process::MakeShared<Baby_full>(samplename+" ("+RoundNumber(inb,0).Data()+" B-hadrons)", 
        Process::Type::background, colors_trub[inb], allfiles, Higfuncs::ntrub==inb && baseline+"&& pass && pass_ra2_badmu && stitch"));
    }
  }

  string jsonCuts = "json4p0";
  float lumi = 4.3;
  if(json=="12p9"){
    lumi = 12.9;
    jsonCuts = "json12p9";
  } else if (json=="full"){
    lumi = 36.2;
    jsonCuts = "1";
  }
  string lumi_s=RoundNumber(lumi,1).Data();

  vector<vector<string>> combos;
  int color_data = kBlue-7;
  if (sample=="zll") { // do 0b vs 1b
    color_data = kOrange+1;
    combos.push_back({"nbm==1","nbm==0"});
    combos.push_back({"nbt==2&&nbm==2","nbm==1"});
  } else if (sample=="qcd") { // do 0b vs 1b and 2b vs 3+b
    color_data = kOrange;
    combos.push_back({"nbm==1","nbm==0"});
    combos.push_back({"nbt>=2&&nbm>=3","nbt==2&&nbm==2"});
    combos.push_back({"nbt==2&&nbm==2","nbm==1"});
  } else if (sample=="search" || sample=="ttbar") {
    if (json=="4p0") { // at low lumi, do 2b vs 3+b
      combos.push_back({"nbt>=2&&nbm>=3","nbt==2&&nbm==2"});
    } else { // at higher lumi, do 2b vs 3b and 2b vs 4b
      combos.push_back({"nbt>=2&&nbm==3&&nbl==3","nbt==2&&nbm==2"});
      combos.push_back({"nbt>=2&&nbm>=3&&nbl>=4","nbt==2&&nbm==2"});
    }
  }

  vector<vector<shared_ptr<Process> >> procs_data;
  for (auto &icomb: combos){
    procs_data.push_back(vector<shared_ptr<Process> >());
    vector<string> combos_label = {"2b","#geq 3b"};
    if (Contains(icomb[0],"nbm==1")) combos_label = {"1b","0b"};
    if (Contains(icomb[0],"nbl>=4")) combos_label = {"2b","4b"};
    if (Contains(icomb[0],"nbl==3")) combos_label = {"2b","3b"};
    if (Contains(icomb[0],"nbt==2")) combos_label = {"2b","1b"};
    string tmpcuts = "pass && pass_ra2_badmu && met/met_calo<5 &&"+jsonCuts+"&&"+icomb[0]; //if not cast here, it crashes
    procs_data.back().push_back(Process::MakeShared<Baby_full>(combos_label[0]+" Data "+lumi_s+" fb^{-1}", 
      Process::Type::data, kBlack, {folderdata+"*root"}, Higfuncs::trig_hig && tmpcuts));
    tmpcuts = "pass && pass_ra2_badmu && met/met_calo<5 &&"+jsonCuts+"&&"+icomb[1];
    procs_data.back().push_back(Process::MakeShared<Baby_full>(combos_label[1]+" Data "+lumi_s+" fb^{-1}", 
      Process::Type::background, kBlack, {folderdata+"*root"}, Higfuncs::trig_hig && tmpcuts));
    procs_data.back().back()->SetFillColor(color_data);
    procs_data.back().back()->SetLineColor(color_data);
    procs_data.back().back()->SetLineWidth(2);
  }


  string metcut = "met>150";
  if (sample=="zll") metcut = "(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))>0";
  else if (sample=="ttbar") metcut = "met>100";

  vector<string> xcuts;
  if (!do_note) xcuts.push_back(metcut);
  xcuts.push_back(metcut+"&& hig_dm<40 && hig_am<200");
  xcuts.push_back(metcut+"&& hig_dm<40 && hig_am<200 && hig_drmax<2.2");

  vector<string> scuts; //additional sample specific options
  scuts.push_back("1");
  if (!do_note) {
    if (sample=="qcd") scuts.push_back("ntks==0");
    if (sample=="search") scuts.push_back("ntks==0 && !low_dphi");
  }

  PlotMaker pm;

  for (unsigned is(0); is<scuts.size(); is++){
    for (unsigned ic(0); ic<xcuts.size(); ic++){
      if (sample!="search" || unblind) {
        for (unsigned i(0); i<combos.size(); i++)
          pm.Push<Hist1D>(Axis(10,0,200,"hig_am", "<m> [GeV]", {100., 140.}),
            baseline+"&&"+xcuts[ic]+"&&"+scuts[is], procs_data[i], plt_types).Tag(sample+"_datavdata"+to_string(i));
      }

      pm.Push<Hist1D>(Axis(10,0,200,"hig_am", "<m> [GeV]", {100., 140.}),
        baseline+"&&"+xcuts[ic]+"&&"+scuts[is], procs, plt_types).Tag(sample+"_shape_bcats");

      pm.Push<Hist1D>(Axis(10,0,200,"hig_am", "<m> [GeV]", {100., 140.}),
        baseline+"&&"+xcuts[ic]+"&&"+scuts[is], procs_trub, plt_types).Tag(sample+"_shape_trub");

      if (sample=="ttbar") {
        pm.Push<Hist1D>(Axis(10,0,200,"hig_am", "<m> [GeV]", {100., 140.}),
          baseline+"&& ntruleps==1 &&"+xcuts[ic]+"&&"+scuts[is], procs, plt_types).Tag(sample+"1l_shape_bcats");
        pm.Push<Hist1D>(Axis(10,0,200,"hig_am", "<m> [GeV]", {100., 140.}),
          baseline+"&& ntruleps==1 &&"+xcuts[ic]+"&&"+scuts[is], procs_trub, plt_types).Tag(sample+"1l_shape_trub");
        pm.Push<Hist1D>(Axis(10,0,200,"hig_am", "<m> [GeV]", {100., 140.}),
          baseline+"&& ntruleps==2 &&"+xcuts[ic]+"&&"+scuts[is], procs, plt_types).Tag(sample+"2l_shape_bcats");
        pm.Push<Hist1D>(Axis(10,0,200,"hig_am", "<m> [GeV]", {100., 140.}),
          baseline+"&& ntruleps==2 &&"+xcuts[ic]+"&&"+scuts[is], procs_trub, plt_types).Tag(sample+"2l_shape_trub");
      }
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
      {"json", required_argument, 0, 'j'},
      {"sample", required_argument, 0, 's'},    
      {"unblind", no_argument, 0, 'u'},   
      {"note", no_argument, 0, 0},    
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "js:tau", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'u':
      unblind = true;
      break;
    case 'a':
      do_allbkg = true;
      break;
    case 's':
      sample = optarg;
      break;
    case 'j':
      json = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "note"){
        do_note = true;
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
        exit(1);
      }
      break;

    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
