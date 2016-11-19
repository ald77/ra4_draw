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
  float lumi = 35.;
}

using namespace std;

NamedFunc nb_tru("nb_tru",[](const Baby &b) -> NamedFunc::ScalarType{
  int nbtru(0);
  for (unsigned i(0); i<b.jets_pt()->size(); i++){
    if (!b.jets_h1()->at(i) && !b.jets_h2()->at(i)) continue;
    if (b.jets_hflavor()->at(i)==5) nbtru++;
  }
  return nbtru;
});
  
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

  string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/");

  // Cuts in baseline speed up the yield finding
  string preseln = "pass && stitch && nvleps==0 && ntks==0 && !low_dphi && njets>=4 && njets<=5"; //keep order to allow replacement with "preseln" in filename & tex
  string baseline = preseln + "&& hig_drmax<2.2";

  Palette colors("txt/colors.txt", "default");

  string ntupletag = "";
  set<string> allfiles = {foldermc+"*_TTJets*Lept*"+ntupletag+"*.root", foldermc+"*_TTJets_HT*"+ntupletag+"*.root",
         foldermc+"*_WJetsToLNu*"+ntupletag+"*.root",foldermc+"*_ST_*"+ntupletag+"*.root",
         foldermc+"*_TTW*"+ntupletag+"*.root",foldermc+"*_TTZ*"+ntupletag+"*.root",
         foldermc+"*_TTGJets*"+ntupletag+"*.root",foldermc+"*_TTTT*"+ntupletag+"*.root",
         foldermc+"*QCD_HT*Inf_Tune*"+ntupletag+"*.root", foldermc+"*QCD_HT*0_Tune*"+ntupletag+"*.root",
         foldermc+"*DYJetsToLL*"+ntupletag+"*.root",
         foldermc+"*_ZJet*"+ntupletag+"*.root",foldermc+"*_ttHJetTobb*"+ntupletag+"*.root",
         foldermc+"*_WH_HToBB*"+ntupletag+"*.root",foldermc+"*_ZH_HToBB*"+ntupletag+"*.root",
         foldermc+"*_WWTo*"+ntupletag+"*.root",foldermc+"*_WZ*"+ntupletag+"*.root",foldermc+"*_ZZ_*"+ntupletag+"*.root"
       };

  // allfiles = set<string>({foldermc+"*_TTJets_Tune*"});

  map<string, vector<shared_ptr<Process> > > procs;

  procs["procs"] = vector<shared_ptr<Process> >();
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t} (l)", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets_HT*.root"},
    baseline+" && ntruleps==1"));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t} (ll)", Process::Type::background, colors("tt_2l"),
    {foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"},
    baseline+" && ntruleps==2"));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("Z+jets", Process::Type::background, kOrange+1,
    {foldermc+"*_ZJet*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {foldermc+"*_WJetsToLNu*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {foldermc+"*_ST_*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t}+W/Z/#gamma", Process::Type::background, kBlue-2,
    {foldermc+"*_TTZ*.root", foldermc+"*_TTW*.root", foldermc+"*_TTGJets*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("QCD", Process::Type::background, colors("other"),
  {foldermc+"*QCD_HT*0_Tune*.root",
    foldermc+"*QCD_HT*Inf_Tune*.root"},
    baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("Other", Process::Type::background, kPink-2,
    {foldermc+"*DYJetsToLL*.root",
    foldermc+"*_TTTT*.root",
    foldermc+"*_ttHJetTobb*.root",
    foldermc+"*_WH_HToBB*.root",
    foldermc+"*_ZH_HToBB*.root",
    foldermc+"*_WWTo*.root",
    foldermc+"*_WZ*.root",
    foldermc+"*_ZZ_*.root"},
    baseline));

  NamedFunc base_func(baseline);
  procs["cats"] = vector<shared_ptr<Process> >();
  procs["cats"].push_back(Process::MakeShared<Baby_full>
          ("#leq 1 B-hadron", Process::Type::background, kPink+2,
           allfiles, base_func && nb_tru<=1));
  procs["cats"].push_back(Process::MakeShared<Baby_full>
  			  ("2 B-hadrons", Process::Type::background, kOrange-4,
  			   allfiles, base_func && nb_tru==2));
  procs["cats"].push_back(Process::MakeShared<Baby_full>
  			  ("3 B-hadrons", Process::Type::background, kTeal-8, 
           allfiles, base_func &&  nb_tru==3));
  procs["cats"].push_back(Process::MakeShared<Baby_full>
  			  ("#geq 4 B-hadrons", Process::Type::background, kAzure-4, 
  			   allfiles, base_func && nb_tru>=4));

  PlotMaker pm;

  vector<TString> metcuts;
  // metcuts.push_back("met>100 && met<=150");
  metcuts.push_back("met>150 && met<=200");
  metcuts.push_back("met>200 && met<=300");
  metcuts.push_back("met>300 && met<=500");
  metcuts.push_back("met>500");

  vector<TString> nbcuts;
  nbcuts.push_back("nbt<=1&&nbm==2");
  nbcuts.push_back("nbt==2&&nbm==2");
  nbcuts.push_back("nbt>=2&&nbm==3&&nbl==3");
  nbcuts.push_back("nbt>=2&&nbm>=3&&nbl>=4");

  vector<TString> regs;
  regs.push_back("hig_am>100 && hig_am<=140 && hig_dm <= 40");
  regs.push_back("(hig_am<=100 || hig_am>140 || hig_dm > 40)");

  vector<TString> cuts;
  vector<TableRow> table_cuts;

  for(auto &imet: metcuts) {
    for(auto &inb: nbcuts) {
      for (auto &ireg: regs) {
        cuts.push_back(baseline+"&&"+imet+"&&"+inb+"&&"+ireg);
      }
    }
  }
  
  for(size_t icut=0; icut<cuts.size(); icut++)
    table_cuts.push_back(TableRow("$"+CodeToLatex(cuts[icut].Data())+"$", cuts[icut].Data()));  
  for(auto &ipr: procs) 
    pm.Push<Table>("chart_"+ipr.first,  table_cuts, ipr.second, true, true, true, true);

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making "<<table_cuts.size()<<" piecharts took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
