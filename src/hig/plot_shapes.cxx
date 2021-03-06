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
  bool paper = true;
  string sample = "search";
  string json = "full";
  bool unblind = true;
  bool note = true;
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
    .Stack(StackType::shapes)
    .LegendColumns(3);
  if (paper) lin_shapes_info = lin_shapes_info.Title(TitleType::simulation);
  // PlotOpt log_shapes_info = lin_shapes_info.YAxis(YAxisType::log);
  vector<PlotOpt> plt_types = {lin_shapes_info};

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higloose/";
  if (sample=="ttbar") foldermc = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higlep1/";
  if (sample=="zll") foldermc = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higlep2/";
  if (sample=="qcd") foldermc = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higqcd/";
  string folderdata = bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_higdata_higloose/";
  if (sample=="ttbar") folderdata = bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_higdata_higlep1/";
  if (sample=="zll") folderdata = bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_higdata_higlep2/";
  if (sample=="qcd") folderdata = bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_higdata_higqcd/";

  set<string> alltags; 
  if (sample=="ttbar" || sample=="search") alltags = {"*TTJets_*Lept*",
                                      "*_TTZ*.root", "*_TTW*.root", "*_TTGJets*.root", 
                                      "*ttHTobb*.root","*_TTTT*.root"};
  if (sample=="zll") alltags = {"*DYJetsToLL*.root"};
  if (sample=="qcd") alltags = {//"*QCD_HT100to200_Tune*", "*QCD_HT200to300_Tune*",
                                //"*QCD_HT300to500_Tune*", 
                                 "*QCD_HT500to700_Tune*",
                                 "*QCD_HT700to1000_Tune*", "*QCD_HT1000to1500_Tune*", 
                                 "*QCD_HT1500to2000_Tune*", "*QCD_HT2000toInf_Tune*"};
  if (do_allbkg) { // don't include QCD MC unless in QCD control sample
    if(sample=="qcd") alltags = {"*TTJets_*Lept*", 
                                 "*_TTZ*.root", "*_TTW*.root", "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root",
                                 "*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root", "*_ST_*.root",
                                 "*QCD_HT100to200_Tune*", "*QCD_HT200to300_Tune*",
                                 "*QCD_HT300to500_Tune*", 
                                 "*QCD_HT500to700_Tune*",
                                 "*QCD_HT700to1000_Tune*", "*QCD_HT1000to1500_Tune*", 
                                 "*QCD_HT1500to2000_Tune*", "*QCD_HT2000toInf_Tune*"
            "*_WH_HToBB*.root", "*_ZH_HToBB*.root", "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"};
    else  alltags = {"*TTJets_*Lept*",
            "*_TTZ*.root", "*_TTW*.root", "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root",
            "*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root", "*_ST_*.root",
            "*_WH_HToBB*.root", "*_ZH_HToBB*.root", "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"};
  }
  set<string> allfiles = attach_folder(foldermc,alltags);

  // Baseline definitions
  string baseline("njets>=4 && njets<=5"); 
  // zll skim: ((elel_m>80&&elel_m<100)||(mumu_m>80&&mumu_m<100)) && 
  // nleps==2 && nleps>=1 && Max$(leps_pt)>30 && njets>=4&&njets<=5
  if (sample=="zll") baseline = baseline+"&& nleps==2 && met<50";
  // qcd skim - met>150 && nvleps==0 && (njets==4||njets==5)
  if (sample=="qcd") baseline = baseline+"&& nvleps==0 && low_dphi";
  // ttbar skim - nleps==1 && (njets==4||njets==5) && nbm>=2
  if (sample=="ttbar") baseline = baseline+"&& nleps==1 && mt<100";
  // search skim - met>100 && nvleps==0 && (njets==4||njets==5) && nbm>=2
  if (sample=="search") baseline = baseline+"&& nvleps==0";
    

  ////// Nb cuts
  vector<string> nbcuts;
  unsigned firstnb = 0;
  nbcuts.push_back("nbdm==0");
  nbcuts.push_back("nbdm==1");
  nbcuts.push_back("nbdt==2&&nbdm==2");
  if (sample=="ttbar" || sample=="search") {
    firstnb = 2;
    nbcuts.push_back("nbdt>=2&&nbdm==3&&nbdl==3");
    nbcuts.push_back("nbdt>=2&&nbdm>=3&&nbdl>=4");
  } 

  string samplename = "t#bar{t}+X";
  if (sample=="qcd") samplename = "QCD";
  if (sample=="zll") samplename = "Z#rightarrow ll";
  if (do_allbkg) samplename = "Bkg.";

  vector<int> colors = {kGreen+3, kGreen+1, kOrange, kAzure+1, kBlue+1};

  // Cuts applied to all processes but not shown in plot title
  string cutsProcs = "pass && pass_ra2_badmu && met/met_calo<5"; 

  vector<shared_ptr<Process> > procs = vector<shared_ptr<Process> >();
  for (unsigned inb(firstnb); inb<nbcuts.size(); inb++){
    // if (sample=="qcd" && inb==nbcuts.size()-1) continue;
    procs.push_back(Process::MakeShared<Baby_full>(samplename+" "+RoundNumber(inb,0).Data()+"b", 
      Process::Type::background, colors[inb], allfiles, baseline +"&&stitch_met&&" + cutsProcs +"&&"+ nbcuts[inb]));
  }
  vector<int> colors_trub = {kAzure-4, kTeal-8, kOrange-4, kPink+2, kMagenta-1};
  vector<shared_ptr<Process> > procs_trub = vector<shared_ptr<Process> >();
  for (unsigned inb(firstnb); inb<nbcuts.size(); inb++){
    if ((sample=="zll" || sample=="qcd") && inb==nbcuts.size()-1) { // merge 4b into 3b
      procs_trub.push_back(Process::MakeShared<Baby_full>(samplename+" #geq"+RoundNumber(inb,0).Data()+" B-hadrons", 
        Process::Type::background, colors_trub[inb], allfiles, Higfuncs::ntrub>=inb &&baseline+"&&stitch_met&&"+cutsProcs));
    } else {
      procs_trub.push_back(Process::MakeShared<Baby_full>(samplename+" "+RoundNumber(inb,0).Data()+" B-hadrons", 
        Process::Type::background, colors_trub[inb], allfiles, Higfuncs::ntrub==inb && baseline +"&&stitch_met&&"+cutsProcs));
    }
  }

  string jsonCuts = "json4p0";
  float lumi = 4.3;
  if(json=="12p9"){
    lumi = 12.9;
    jsonCuts = "json12p9";
  } else if (json=="full"){
    lumi = 35.9;
    jsonCuts = "1";
  }
  string lumi_s=RoundNumber(lumi,1).Data();

  vector<vector<string>> combos, combos_labels;
  int color_data = kBlue-7;
  if (sample=="zll") { // do 0b vs 1b
    color_data = kOrange+1;
    combos.push_back({"nbdm==0","nbdm==1"});
    combos.push_back({"nbdm==1","nbdt==2&&nbdm==2"});
  } else if (sample=="qcd") { // do 0b vs 1b and 2b vs 3+b
    color_data = kOrange;
    combos.push_back({"nbdm==0","nbdm==1"});
    combos.push_back({"nbdm==1","nbdt==2&&nbdm==2"});
    combos.push_back({"nbdt==2&&nbdm==2","nbdt>=2&&nbdm==3"});
    combos.push_back({"nbdt==2&&nbdm==2", "nbdt>=2&&nbdm>=3&&nbdl>=4"});
  } else if (sample=="search" || sample=="ttbar") {
    if (json=="4p0") { // at low lumi, do 2b vs 3+b
      combos.push_back({"nbdt==2&&nbdm==2","nbdt>=2&&nbdm>=3"});
    } else { // at higher lumi, do 2b vs 3b and 2b vs 4b
      combos.push_back({"nbdt==2&&nbdm==2", "nbdt>=2&&nbdm==3&&nbdl==3"});
      combos.push_back({"nbdt==2&&nbdm==2", "nbdt>=2&&nbdm>=3&&nbdl>=4"});
    }
  }
  for (unsigned ind(0); ind<combos.size(); ind++){
    vector<string> label;
    for(unsigned lab=0; lab<2; lab++){
      if (Contains(combos[ind][lab],"nbdm==0")) label.push_back("0b");
      if (Contains(combos[ind][lab],"nbdm==1")) label.push_back("1b");
      if (Contains(combos[ind][lab],"nbdm==2")) label.push_back("2b");
      if (Contains(combos[ind][lab],"nbdm==3")) label.push_back("3b");
      if (Contains(combos[ind][lab],"nbdm>=3") && !Contains(combos[ind][lab],"nbdl>=4")) label.push_back("#geq 3b");
      if (Contains(combos[ind][lab],"nbdl>=4")) label.push_back("4b");
    }
    combos_labels.push_back(label);
  }

  vector<vector<shared_ptr<Process> >> procs_data;
  for (unsigned ind(0); ind<combos.size(); ind++){
    vector<string> icomb = combos[ind], ilab = combos_labels[ind];
    procs_data.push_back(vector<shared_ptr<Process> >());
    string tmpcuts = cutsProcs + " &&"+jsonCuts+"&&"+icomb[0]; //if not cast here, it crashes
    procs_data.back().push_back(Process::MakeShared<Baby_full>(ilab[0]+" Data "+lumi_s+" fb^{-1}", 
							       Process::Type::background, kBlack, {folderdata+"*root"}, 
							       Higfuncs::trig_hig && tmpcuts));
    procs_data.back().back()->SetFillColor(color_data);
    procs_data.back().back()->SetLineColor(color_data);
    procs_data.back().back()->SetLineWidth(2);

    tmpcuts = cutsProcs + " &&"+jsonCuts+"&&"+icomb[1];
    procs_data.back().push_back(Process::MakeShared<Baby_full>(ilab[1]+" Data "+lumi_s+" fb^{-1}", 
							       Process::Type::data, kBlack, 
      {folderdata+"*root"}, Higfuncs::trig_hig && tmpcuts));
  }


  string metcut = "met>150";
  if (sample=="zll") metcut = "(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))>0";
  else if (sample=="ttbar") metcut = "1";

  vector<string> xcuts;
  if (!note) xcuts.push_back(metcut);
  //xcuts.push_back(metcut+"&& higd_dm<40 && higd_am<200");
  xcuts.push_back(metcut+"&& higd_dm<40 && higd_am<200 && higd_drmax<2.2");

  vector<string> scuts; //additional sample specific options
  scuts.push_back("1");
  // if (sample=="qcd") scuts.push_back("ntks==0");
  if (sample=="search") scuts.push_back("ntks==0 && !low_dphi");
  if (sample=="ttbar") scuts.push_back("!low_dphi");
  

  PlotMaker pm;
  NamedFunc wgt = Higfuncs::weight_higd*Higfuncs::eff_higtrig;

  for (unsigned is(0); is<scuts.size(); is++){
    for (unsigned ic(0); ic<xcuts.size(); ic++){
      if (sample!="search" || unblind) {
        for (unsigned i(0); i<combos.size(); i++)
          pm.Push<Hist1D>(Axis(10,0,200,"higd_am", "#LTm#GT [GeV]", {100., 140.}),
            baseline+"&&"+xcuts[ic]+"&&"+scuts[is], procs_data[i], plt_types)
	    .Tag(sample+"_datavdata"+to_string(i)).RatioTitle(combos_labels[i][1],combos_labels[i][0]);
      }

      pm.Push<Hist1D>(Axis(10,0,200,"higd_am", "#LTm#GT [GeV]", {100., 140.}),
        baseline+"&&"+xcuts[ic]+"&&"+scuts[is], procs, plt_types).Weight(wgt).Tag(sample+"_shape_bcats")
        .RatioTitle("Bkg. nb","Bkg. 2b");

      pm.Push<Hist1D>(Axis(10,0,200,"higd_am", "#LTm#GT [GeV]", {100., 140.}),
        baseline+"&&"+xcuts[ic]+"&&"+scuts[is], procs_trub, plt_types).Weight(wgt).Tag(sample+"_shape_trub");

      if (sample=="ttbar" && !note) {
        pm.Push<Hist1D>(Axis(10,0,200,"higd_am", "#LTm#GT [GeV]", {100., 140.}),
          "njets>=4 && njets<=5 && ntruleps==1 &&"+xcuts[ic]+"&&"+scuts[is], procs, plt_types).Weight(wgt).Tag(sample+"1l_shape_bcats");
        pm.Push<Hist1D>(Axis(10,0,200,"higd_am", "#LTm#GT [GeV]", {100., 140.}),
          "njets>=4 && njets<=5 && ntruleps==1 &&"+xcuts[ic]+"&&"+scuts[is], procs_trub, plt_types).Weight(wgt).Tag(sample+"1l_shape_trub");
        pm.Push<Hist1D>(Axis(10,0,200,"higd_am", "#LTm#GT [GeV]", {100., 140.}),
          "njets>=4 && njets<=5 && ntruleps==2 &&"+xcuts[ic]+"&&"+scuts[is], procs, plt_types).Weight(wgt).Tag(sample+"2l_shape_bcats");
        pm.Push<Hist1D>(Axis(10,0,200,"higd_am", "#LTm#GT [GeV]", {100., 140.}),
          "njets>=4 && njets<=5 && ntruleps==2 &&"+xcuts[ic]+"&&"+scuts[is], procs_trub, plt_types).Weight(wgt).Tag(sample+"2l_shape_trub");
      }
    }
  }

  pm.min_print_ = true;
  pm.MakePlots(35.9);

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
        note = true;
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
