///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

#include "TError.h" // Controls error level reporting

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/plot_opt.hpp"

namespace{
  //bool do_met150 = true;
}

using namespace std;

int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_stdnj5/");
  //if(do_met150) foldermc = (bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_met150/");
  Palette colors("txt/colors.txt", "default");

  // Cuts in baseline speed up the yield finding
  string baseline = "mj14>250 && nleps>=1 && ht>500 && met>150 && pass && njets>=5 && weight<1"; // Excluding one QCD event
  //if(do_met150)  baseline = "mj14>250 && nleps>=1 && ht>500 && met>150 && pass && njets>=5";

  auto proc_tt1l = Process::MakeShared<Baby_full>("t#bar{t} (l)", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets_HT*.root"},
    baseline+" && stitch && ntruleps==1");
  auto proc_tt2l = Process::MakeShared<Baby_full>("t#bar{t} (ll)", Process::Type::background, colors("tt_2l"),
    {foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"},
    baseline+" && stitch && ntruleps==2 && ntrutaush==0");
  auto proc_ttltau = Process::MakeShared<Baby_full>("t#bar{t} (#tau_{h}l)", Process::Type::background, colors("tt_ltau"),
    {foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"},
    baseline+" && stitch && ntruleps==2 && ntrutaush>=1");
  auto proc_wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {foldermc+"*_WJetsToLNu*.root"}, baseline+" && stitch");
  auto proc_single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {foldermc+"*_ST_*.root"}, baseline);
  auto proc_ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {foldermc+"*_TTWJets*.root", foldermc+"*_TTZ*.root"}, baseline);
  auto proc_other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {foldermc+"*DYJetsToLL*.root",foldermc+"*QCD_HT*0_Tune*.root",foldermc+"*QCD_HT*Inf_Tune*.root",
        foldermc+"*_ZJet*.root",foldermc+"*_ttHJetTobb*.root",
        foldermc+"*_TTGJets*.root",foldermc+"*_TTTT*.root",
        foldermc+"*_WH_HToBB*.root",foldermc+"*_ZH_HToBB*.root",
        foldermc+"*_WWTo*.root",foldermc+"*_WZ*.root",foldermc+"*_ZZ_*.root"},
    baseline+" && stitch");

  // auto proc_tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l", Process::Type::background, colors("tt_1l"),
  //   {foldermc+"*_TTJets_Tune*.root"},
  //   baseline+" && ntruleps==1");
  // auto proc_tt2l = Process::MakeShared<Baby_full>("t#bar{t} 2l", Process::Type::background, colors("tt_2l"),
  //   {foldermc+"*_TTJets_Tune**.root"},
  //   baseline+" && ntruleps==2");



  vector<shared_ptr<Process> > all_procs = {proc_tt1l, proc_tt2l, proc_ttltau, proc_wjets, proc_single_t, proc_ttv, proc_other};
  //vector<shared_ptr<Process> > all_procs = {proc_tt1l, proc_tt2l};

  ////// MET cuts
  TString c_vlowmet = "met>150 && met<=200";
  TString c_lowmet  = "met>200 && met<=350";
  TString c_midmet  = "met>350 && met<=500";
  TString c_higmet  = "met>500";
  vector<TString> metcuts({c_vlowmet, c_lowmet, c_midmet, c_higmet});
  //  if(do_met150) metcuts = vector<TString>({c_vlowmet});

  ////// Nb cuts
  TString c_lownb = "nbm==1";
  TString c_midnb = "nbm==2";
  TString c_hignb = "nbm>=3";
  vector<TString> nbcuts({c_lownb, c_midnb, c_hignb});
  //vector<TString> nbcuts({"nbm>=2"});

  ////// Njets cuts
  TString c_lownj = "njets>=6 && njets<=8";
  TString c_hignj = "njets>=9";
  TString c_nj5   = "njets==5";
  vector<TString> njcuts({c_nj5, c_lownj, c_hignj});
  //vector<TString> njcuts({"njets>=6"});

  ////// mT cuts
  TString c_lowmt = "mt<=140";
  TString c_higmt = "mt>140";
  vector<TString> mtcuts({c_lowmt, c_higmt});

  // Adding nleps==1 cuts
  vector<TString> cuts;
  for(size_t imet=0; imet<metcuts.size(); imet++)
    for(size_t inb=0; inb<nbcuts.size(); inb++)
      for(size_t inj=0; inj<njcuts.size(); inj++)
	for(size_t imt=0; imt<mtcuts.size(); imt++){
	  cuts.push_back("nleps==1 && nveto==0 && "+metcuts[imet]+"&&"+nbcuts[inb]+"&&"+njcuts[inj]+"&&"+mtcuts[imt]);
	  cuts.push_back("nleps==1 && nveto==1 && "+metcuts[imet]+"&&"+nbcuts[inb]+"&&"+njcuts[inj]+"&&"+mtcuts[imt]);
	}

  // Adding nleps==2 cuts
  njcuts = vector<TString>({"njets>=5 && njets<=7", "njets>=8"});
  for(size_t imet=0; imet<metcuts.size(); imet++)
    for(size_t inj=0; inj<njcuts.size(); inj++)
      cuts.push_back("nleps==2 && "+metcuts[imet]+"&&"+njcuts[inj]);
 
  vector<TableRow> table_cuts;
  for(size_t icut=0; icut<cuts.size(); icut++){
    table_cuts.push_back(TableRow("$"+CodeToLatex(cuts[icut].Data())+"$", cuts[icut].Data()));
  }

  PlotMaker pm;
  pm.Push<Table>("chart",  table_cuts, all_procs, true, true, true);
  pm.min_print_ = true;
  pm.MakePlots(40.);

  time(&endtime);
  cout<<endl<<"Making "<<table_cuts.size()<<" piecharts took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
