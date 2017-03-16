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
  bool do_ra4 = true;
  bool do_pie = false;
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  Palette colors("txt/colors.txt", "default");
  string lsp = "{#lower[-0.1]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}}";
  string t1t_s = "#scale[0.95]{#tilde{g}#kern[0.2]{#tilde{g}}, #tilde{g}#rightarrowt#kern[0.18]{#bar{t}}#kern[0.18]"
    +lsp;


  /////////////////// PLOT STYLES //////////////////////////////////////
  PlotOpt lin_norm_info("txt/plot_styles.txt", "CMSPaper");
  lin_norm_info.Title(TitleType::info)   
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::linear)
    .Stack(StackType::data_norm);
  PlotOpt log_norm_info = lin_norm_info.YAxis(YAxisType::log);
  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};

  PlotOpt log_norm = lin_norm_info.YAxis(YAxisType::log).Title(TitleType::simulation_preliminary).LogMinimum(.7).Bottom(BottomType::off).ShowBackgroundError(false);
  PlotOpt lin_norm = lin_norm_info.YAxis(YAxisType::linear).Title(TitleType::simulation_preliminary).Bottom(BottomType::off).ShowBackgroundError(false);
  vector<PlotOpt> plt_norm = {lin_norm, log_norm};

  PlotOpt lin_shapes = lin_norm.Stack(StackType::shapes).Bottom(BottomType::ratio);
  vector<PlotOpt> plt_shapes = {lin_shapes};



  //////////////////////////////////// PROCESSES /////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  // Baseline definitions
  string baseline("nleps==1  && nveto==0 && njets>=6 && nbm>=2 && met/met_calo<5. && pass");

  string foldermc = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_mcbase_stdnj5/";
  string folderdata = bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_database_stdnj5/";
  string tag ="";
  auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets*Lept*"+tag+"*"}, baseline+"&&ntruleps<=1&&stitch_met");
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {foldermc+"*_TTJets*Lept*"+tag+"*"}, baseline+"&&ntruleps>=2&&stitch_met");
  auto wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {foldermc+"*_WJetsToLNu*"+tag+"*"},baseline+"&&stitch");
  auto single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {foldermc+"*_ST_*"+tag+"*"}, baseline);
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {foldermc+"*_TTWJets*"+tag+"*", foldermc+"*_TTGJets*"+tag+"*", foldermc+"*_TTZTo*"+tag+"*"}, baseline);
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {foldermc+"*DYJetsToLL*"+tag+"*", foldermc+"*_QCD_HT*"+tag+"*",
        foldermc+"*_ZJet*"+tag+"*", foldermc+"*_WWTo*"+tag+"*",
        foldermc+"*ggZH_HToBB*"+tag+"*", foldermc+"*ttHJetTobb*"+tag+"*",
        foldermc+"*_TTTT_*"+tag+"*",
        foldermc+"*_WH_HToBB*"+tag+"*", foldermc+"*_WZTo*"+tag+"*",
        foldermc+"*_ZH_HToBB*"+tag+"*", foldermc+"_ZZ_*"+tag+"*"}, baseline);

  auto t1tttt_nc = Process::MakeShared<Baby_full>(t1t_s+"}(1800,100)", Process::Type::signal, colors("t1tttt"),
    {"/net/cms29/cms29r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/*SMS-T1tttt_mGluino-1800_mLSP-100_*.root"}, baseline);
  auto t1tttt_c = Process::MakeShared<Baby_full>(t1t_s+"}(1400,1000)", Process::Type::signal, colors("t1tttt"),
    {"/net/cms29/cms29r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/*SMS-T1tttt_mGluino-1400_mLSP-1000*.root"}, baseline);
  t1tttt_c->SetLineStyle(2);

  auto data = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {folderdata+"*.root"},baseline+"&&trig_ra4");

  vector<shared_ptr<Process> > procs = {t1tttt_nc, t1tttt_c, tt1l, tt2l, wjets, single_t, ttv, other};
  vector<shared_ptr<Process> > procs_tt = {tt2l, tt1l};



  string cuts = "nbm>=2";
  PlotMaker pm;
  if(do_ra4){
    cuts = "nbm>=2&&mt<=140&&met>350";
    pm.Push<Hist1D>(Axis(13, 25., 1000., "mj14", "M_{J} [GeV]", {250,400.}),cuts, procs, plt_norm).Tag("ra4");
    cuts = "nbm>=2&&mt>140&&met>350";
    pm.Push<Hist1D>(Axis(13, 25., 1000., "mj14", "M_{J} [GeV]", {250,400.}),cuts, procs, plt_norm).Tag("ra4");
    // cuts = "nbm>=2&&mt<=140";
    // pm.Push<Hist1D>(Axis(13, 25., 1000., "mj14", "M_{J} [GeV]", {250,400.}),cuts, procs, plt_norm).Tag("ra4");
    // cuts = "nbm>=2&&mt>140";
    // pm.Push<Hist1D>(Axis(13, 25., 1000., "mj14", "M_{J} [GeV]", {250,400.}),cuts, procs, plt_norm).Tag("ra4");
    cuts = "nbm>=2&&mj14>250&&met>350";
    pm.Push<Hist1D>(Axis(14, 0., 280., "mt", "m_{T} [GeV]", {140.}),cuts, procs, plt_norm).Tag("ra4");
    // cuts = "mj14>250";
    // pm.Push<Hist1D>(Axis(7, 5.5, 12.5, "njets", "N_{jets}", {8.5}),cuts, procs_tt, plt_shapes).Tag("ra4");
  } // do_ra4


  ///////////////////// PIE CHARTS /////////////////// (Not updated since Manuel's 02_22 seminar macro)
  vector<TString> metcuts;
  metcuts.push_back("met>100 && met<=150");
  metcuts.push_back("met>150 && met<=200");
  metcuts.push_back("met>200 && met<=350");
  metcuts.push_back("met>350 && met<=500");
  metcuts.push_back("met>500");
  //  metcuts.push_back("met>300");

  vector<TString> nbcuts;
  nbcuts.push_back("nbm>=2");

  vector<TString> njcuts;
  njcuts.push_back("njets>=6");

  vector<TString> mtcuts({"mt<=140", "mt>140"});
  //vector<TString> mjcuts({"mj14>250&&mj14<=400", "mj14>400"});
  vector<TString> mjcuts({"mj14>250"});

  vector<TString> piecuts;
  vector<TableRow> table_cuts;

  //// nleps = 1
  for(auto &imet: metcuts) 
    for(auto &inb: nbcuts) 
      for(auto &inj: njcuts) 
	for(auto &imj: mjcuts) 
	  for(auto &imt: mtcuts) {
	    piecuts.push_back("nleps==1 && nveto==0 && "+imet+"&&"+inb+"&&"+inj+"&&"+imt+"&&"+imj);
	  }

  for(size_t icut=0; icut<piecuts.size(); icut++)
    table_cuts.push_back(TableRow("$"+CodeToLatex(piecuts[icut].Data())+"$", piecuts[icut].Data()));  
  if(do_pie) pm.Push<Table>("chart_full",  table_cuts, procs, true, true, true, false);

  pm.min_print_ = true;
  pm.MakePlots(35.9);

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
