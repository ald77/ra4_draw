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
  bool do_ra4 = false;
  bool do_hig = true;
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

  PlotOpt log_norm = lin_norm_info.YAxis(YAxisType::log).Title(TitleType::simulation).LogMinimum(.7).Bottom(BottomType::off);
  PlotOpt lin_norm = lin_norm_info.YAxis(YAxisType::linear).Title(TitleType::simulation).Bottom(BottomType::off);
  vector<PlotOpt> plt_norm = {lin_norm, log_norm};

  PlotOpt lin_shapes = lin_norm.Stack(StackType::shapes).Bottom(BottomType::ratio);
  vector<PlotOpt> plt_shapes = {lin_shapes};



  //////////////////////////////////// PROCESSES /////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  // Baseline definitions
  string baseline("nleps==1  && nveto==0 && njets>=6 && nbm>=2 && met/met_calo<5 && pass");

  string tag = "metG200"; tag="";
  string foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_met100_stdnj5/";
  string folderdata = bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_database_standard/";

  auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets*Lept*"+tag+"*", foldermc+"*_TTJets_HT*"+tag+"*"}, baseline+"&&ntruleps<=1&&stitch");
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {foldermc+"*_TTJets*Lept*"+tag+"*", foldermc+"*_TTJets_HT*"+tag+"*"}, baseline+"&&ntruleps>=2&&stitch");
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

  auto t1tttt_nc = Process::MakeShared<Baby_full>(t1t_s+"}(1500,100)", Process::Type::signal, colors("t1tttt"),
    {foldermc+"*SMS-T1tttt_mGluino-1500_mLSP-100*"+tag+"*"}, baseline);
  auto t1tttt_c = Process::MakeShared<Baby_full>(t1t_s+"}(1200,800)", Process::Type::signal, colors("t1tttt"),
    {foldermc+"*SMS-T1tttt_mGluino-1200_mLSP-800*"+tag+"*"}, baseline);
  t1tttt_c->SetLineStyle(2);

  auto data = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {folderdata+"*.root"},baseline+"&&trig_ra4");

  vector<shared_ptr<Process> > procs = {t1tttt_nc, t1tttt_c, tt1l, tt2l, wjets, single_t, ttv, other};
  vector<shared_ptr<Process> > procs_tt = {tt2l, tt1l};

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
  NamedFunc base_func = "pass && met/met_calo<5 && njets>=4 && njets<=5 && nbt>=2 && nvleps==0&&ntks==0&&!low_dphi&&met>150&&hig_dm<40&&weight<.1";

  vector<shared_ptr<Process> > procs_hig;
  vector<string> sigm = {"400"}; 
  vector<int> sig_colors = {kRed}; // need sigm.size() >= sig_colors.size()
  for (unsigned isig(0); isig<sigm.size(); isig++)
    procs_hig.push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)", Process::Type::background, 
        sig_colors[isig], {foldersig+"*TChiHH_mGluino-"+sigm[isig]+"*.root"}, base_func));
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

  // vector<string> sigm = {"225","400","700"}; 
  // vector<int> sig_colors = {kGreen, kRed, kBlue}; // need sigm.size() >= sig_colors.size()
  // for (unsigned isig(0); isig<sigm.size(); isig++)
  //   procs_hig.push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
  //       sig_colors[isig], {foldersig+"*TChiHH_mGluino-"+sigm[isig]+"*.root"}, base_func));

  // procs_hig.push_back(Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
  //   {folderhigdata+"*root"},  Higfuncs::trig_hig>0. && base_func)); 

  ///////////////////////////////////////////////////////////////////////////////////////////////////


  string cuts = "nbm>=2";
  PlotMaker pm;
  if(do_ra4){
    cuts = "nbm>=2&&mt<=140&&mj14<1000&&met>300";
    pm.Push<Hist1D>(Axis(13, 25., 1000., "mj14", "M_{J} [GeV]", {250,400.}),cuts, procs, plt_norm).Tag("ra4");
    cuts = "nbm>=2&&mt>140&&mj14<1000&&met>300";
    pm.Push<Hist1D>(Axis(13, 25., 1000., "mj14", "M_{J} [GeV]", {250,400.}),cuts, procs, plt_norm).Tag("ra4");
    // cuts = "nbm>=2&&mt<=140";
    // pm.Push<Hist1D>(Axis(13, 25., 1000., "mj14", "M_{J} [GeV]", {250,400.}),cuts, procs, plt_norm).Tag("ra4");
    // cuts = "nbm>=2&&mt>140";
    // pm.Push<Hist1D>(Axis(13, 25., 1000., "mj14", "M_{J} [GeV]", {250,400.}),cuts, procs, plt_norm).Tag("ra4");
    cuts = "nbm>=2&&mj14>250&&met>300&&mt<280";
    pm.Push<Hist1D>(Axis(14, 0., 280., "mt", "m_{T} [GeV]", {140.}),cuts, procs, plt_norm).Tag("ra4");
    // cuts = "mj14>250";
    // pm.Push<Hist1D>(Axis(7, 5.5, 12.5, "njets", "N_{jets}", {8.5}),cuts, procs_tt, plt_shapes).Tag("ra4");
  } // do_ra4

  if(do_hig){
    cuts = "nbm>=3&&nbl>=4&&met>150&&hig_drmax<2.2";
    pm.Push<Hist1D>(Axis(50, 0., 250., "hig_am", "#LTm#GT [GeV]", {100., 140.}),cuts, procs_hig, plt_norm).Tag("hig");
    cuts = "nbm>=3&&nbl>=4&&met>200&&hig_drmax<2.2";
    pm.Push<Hist1D>(Axis(50, 0., 250., "hig_am", "#LTm#GT [GeV]", {100., 140.}),cuts, procs_hig, plt_norm).Tag("hig");
    cuts = "nbm>=3&&nbl>=4&&met>300&&hig_drmax<2.2";
    pm.Push<Hist1D>(Axis(50, 0., 250., "hig_am", "#LTm#GT [GeV]", {100., 140.}),cuts, procs_hig, plt_norm).Tag("hig");

    cuts = "nbm>=3&&nbl>=4&&met>150&&hig_drmax<2.2";
    pm.Push<Hist1D>(Axis(9, 150., 600., "met", "E^{miss}_{T} [GeV]", {150., 200., 300., 450.}),cuts, procs_hig, plt_norm).Tag("hig");
    // cuts = "mt<100 && met>150 && nleps==11";
    // pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nbm", "CSV N_{b,M}"),cuts, procs_hig, plt_norm_info).Tag("hig");
    // pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nbt", "CSV N_{b,T}"),cuts, procs_hig, plt_norm_info).Tag("hig");
    // pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nbdm", "DeepCSV N_{b,M}"),cuts, procs_hig, plt_norm_info).Tag("hig");
    // pm.Push<Hist1D>(Axis(6, -0.5, 5.5, "nbdt", "DeepCSV N_{b,T}"),cuts, procs_hig, plt_norm_info).Tag("hig");
  } // do_hig


  ///////////////////// PIE CHARTS ///////////////////
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
