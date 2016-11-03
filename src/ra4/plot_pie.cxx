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
  bool do_allplots = false;
}

using namespace std;

NamedFunc offshellw("offshellw",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (abs(b.mc_id()->at(i))!=24) continue;
      if (b.mc_mass()->at(i) > 140.) {
        return 1;
      }
    }
    return 0;
  });

//// Lowers the njets requirements for dilepton bins
TString lowerNjets(TString &cut){
  TString lowcut = cut;
  for(int nj=6; nj<=9; nj++){
    TString nj_s = ""; nj_s += nj;
    TString njlo_s = ""; njlo_s += nj-1;
    lowcut.ReplaceAll(nj_s, njlo_s);
  }
  return lowcut;
}

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

  string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_met100_stdnj5/");
  //if(do_met150) foldermc = (bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_met150/");
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

  // Cuts in baseline speed up the yield finding
  string baseline = "pass && stitch && mj14>250 && nleps>=1 && st>500 && met>100 && njets>=5 && weight<1"; // Excluding one QCD event

  map<string, vector<shared_ptr<Process> > > procs;


  procs["procs"] = vector<shared_ptr<Process> >();
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t} (l)", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets_HT*.root"},
    baseline+" && ntruleps==1"));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t} (ll)", Process::Type::background, colors("tt_2l"),
    {foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"},
    baseline+" && ntruleps==2 && ntrutaush==0"));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t} (#tau_{h}l)", Process::Type::background, colors("tt_ltau"),
    {foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"},
    baseline+" && ntruleps==2 && ntrutaush>=1"));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {foldermc+"*_WJetsToLNu*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {foldermc+"*_ST_*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t}Z", Process::Type::background, kOrange+7,
    {foldermc+"*_TTZ*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t}W", Process::Type::background, colors("ttv"),
    {foldermc+"*_TTW*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("QCD", Process::Type::background, colors("other"),
  {foldermc+"*QCD_HT*0_Tune*.root",
    foldermc+"*QCD_HT*Inf_Tune*.root"},
    baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("tt#gamma", Process::Type::background, kMagenta+3,
    {foldermc+"*_TTGJets*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("tttt", Process::Type::background, kBlue+2,
    {foldermc+"*_TTTT*.root"}, baseline));
  procs["procs"].push_back(Process::MakeShared<Baby_full>("Other", Process::Type::background, kPink-2,
    {foldermc+"*DYJetsToLL*.root",
    foldermc+"*_ZJet*.root",
    foldermc+"*_ttHJetTobb*.root",
    foldermc+"*_WH_HToBB*.root",
    foldermc+"*_ZH_HToBB*.root",
    foldermc+"*_WWTo*.root",
    foldermc+"*_WZ*.root",
    foldermc+"*_ZZ_*.root"},
    baseline));


  NamedFunc multNeu = "(type==5000 || type==13000 || type==15000 || type==16000)";
  NamedFunc multNeu2l = "(ntruleps>=2 || (ntruleps<=1&&(type==5000 || type==13000 || type==15000 || type==16000)))";
  procs["cats"] = vector<shared_ptr<Process> >();
  procs["cats"].push_back(Process::MakeShared<Baby_full>
  			  ("#geq2#nu^{prompt}", Process::Type::background, kCyan-3,
  			   allfiles, baseline && multNeu2l));
  procs["cats"].push_back(Process::MakeShared<Baby_full>
  			  ("#leq1#kern[.1]{#nu^{pr.}}, mismeas.", Process::Type::background, kRed-4, allfiles, 
  			   baseline && "ntruleps<=1 && mt_tru<=140" && !multNeu));
  procs["cats"].push_back(Process::MakeShared<Baby_full>
  			  ("#leq1#kern[.1]{#nu^{pr.}}, #geq1#kern[.1]{#nu}^{non-prompt}", Process::Type::background, kGreen-3, 
  			   allfiles, baseline && "ntruleps<=1 && mt_tru>140" && !multNeu && offshellw==0.));
  procs["cats"].push_back(Process::MakeShared<Baby_full>
  			  ("#leq1#kern[.1]{#nu^{pr.}}, off-shell W", Process::Type::background, kOrange,
  			   allfiles, baseline && "ntruleps<=1 && mt_tru>140" && !multNeu && offshellw>0.));

  PlotMaker pm;

  vector<TString> metcuts;
  //metcuts.push_back("met>100");
  metcuts.push_back("met>100 && met<=150");
  metcuts.push_back("met>150 && met<=200");
  metcuts.push_back("met>200 && met<=350");
  metcuts.push_back("met>350 && met<=500");
  metcuts.push_back("met>500");
  metcuts.push_back("met>200");

  vector<TString> nbcuts;
  //nbcuts.push_back("nbm==0");
  nbcuts.push_back("nbm>=1");
  if(do_allplots){
    nbcuts.push_back("nbm==1");
    nbcuts.push_back("nbm==2");
    nbcuts.push_back("nbm>=3");
  }

  vector<TString> njcuts;
  njcuts.push_back("njets==5");
  //njcuts.push_back("njets>=5");
  njcuts.push_back("njets>=6");
  njcuts.push_back("njets>=6 && njets<=8");
  njcuts.push_back("njets>=9");
  
  vector<TString> mtcuts({"mt<=140", "mt>140"});
  
  // Adding nleps==1 cuts
  vector<TString> cuts;
  vector<TableRow> table_cuts;

  //// nleps = 1
  for(auto &imet: metcuts) 
    for(auto &inb: nbcuts) 
      for(auto &inj: njcuts) 
	for(auto &imt: mtcuts) {
	  cuts.push_back("nleps==1 && nveto==0 && "+imet+"&&"+inb+"&&"+inj+"&&"+imt);
	  cuts.push_back("nleps==1 && nveto==1 && "+imet+"&&"+inb+"&&"+inj+"&&"+imt);
	}
  


  //// nleps = 2 and ll+lv
  for(auto &imet: metcuts) 
    for(size_t ind=0; ind<njcuts.size(); ind++) {
      if(njcuts[ind]=="njets==5") continue;
      cuts.push_back("nleps==2 && nbm<=2 && "+imet+"&&"+lowerNjets(njcuts[ind]));
      if(do_allplots) 
	cuts.push_back("((nleps==2 && nbm<=2 && "+lowerNjets(njcuts[ind])+") || (nleps==1 && nveto==1 && nbm>=1 && mt>140 && "
		       +njcuts[ind]+")) && "+imet);
    }

  for(size_t icut=0; icut<cuts.size(); icut++)
    table_cuts.push_back(TableRow("$"+CodeToLatex(cuts[icut].Data())+"$", cuts[icut].Data()));  
  for(auto &ipr: procs) 
    pm.Push<Table>("chart_"+ipr.first,  table_cuts, ipr.second, true, true, true, false);

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making "<<table_cuts.size()<<" piecharts took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
