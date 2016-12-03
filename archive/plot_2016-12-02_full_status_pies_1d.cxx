#include <cmath>
#include <stdio.h>
#include <chrono>

#include "TError.h"
#include "TVector2.h"

#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist1d.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "core/table.hpp"
#include "core/slide_maker.hpp"

using namespace std;
using namespace PlotOptTypes;


NamedFunc max_b_pt("max_b_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float maxPt=-999.;
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (abs(b.mc_id()->at(i))!=5) continue;
      if(b.mc_pt()->at(i) > maxPt) maxPt = b.mc_pt()->at(i);
    }
    return maxPt;
  });

NamedFunc max_t_pt("max_t_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float maxPt=-999.;
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (abs(b.mc_id()->at(i))!=6) continue;
      if(b.mc_pt()->at(i) > maxPt) maxPt = b.mc_pt()->at(i);
    }
    return maxPt;
  });


int main(){
  gErrorIgnoreLevel = 6000;

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  double lumi = 32.2;
  string bfolder("");
  string hostname(execute("echo $HOSTNAME"));
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string folder(bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_met100_stdnj5/");

  Palette colors("txt/colors.txt", "default");

  string ntupletag = "metG200";
  ntupletag = "";

  set<string>t1ncfiles({"*T1tttt*1500*"+ntupletag+"*.root"});
  set<string>t1cfiles({"*T1tttt*1200*"+ntupletag+"*.root"});


  set<string>ttfiles({"*_TTJets*Lept*"+ntupletag+"*.root", "*_TTJets_HT*"+ntupletag+"*.root"});
  set<string>wfiles({"*_WJetsToLNu*"+ntupletag+"*.root"});
  set<string>stfiles({"*_ST_*"+ntupletag+"*.root"});
  set<string>ofiles({"*_WH_HToBB*"+ntupletag+"*.root", "*_ZH_HToBB*"+ntupletag+"*.root","*_WWTo*"+ntupletag+"*.root", "*_WZ*"+ntupletag+"*.root", "*_ZZ_*"+ntupletag+"*.root", "*_TTZ*"+ntupletag+"*.root", 
	"*_TTW*"+ntupletag+"*.root","*_TTGJets*"+ntupletag+"*.root", "*_ttHJetTobb*"+ntupletag+"*.root","*_TTTT*"+ntupletag+"*.root", "*QCD_HT*0_Tune*"+ntupletag+"*.root", "*QCD_HT*Inf_Tune*"+ntupletag+"*.root", 
	"*_ZJet*"+ntupletag+"*.root", "*DYJetsToLL*"+ntupletag+"*.root"});


  string pscuts = "pass&&stitch";

  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_full>("t#bar{t}(1l), m_{T}<140", Process::Type::background, colors("tt_1l"),
						 attach_folder(folder,ttfiles), pscuts+"&& ntruleps==1&&mt<140&&mj14<850"));
  procs.push_back(Process::MakeShared<Baby_full>("t#bar{t}(2l), m_{T}>140", Process::Type::background, colors("tt_2l"),
						 attach_folder(folder,ttfiles), pscuts+"&& ntruleps==2&&mt>140&&mj14<850"));
  procs.push_back(Process::MakeShared<Baby_full>("T1tttt(1500,100)", Process::Type::background, 2,
						 attach_folder(folder,t1ncfiles), pscuts+"&& mj14<850"));
  // procs.push_back(Process::MakeShared<Baby_full>("T1tttt(1200,800)", Process::Type::background, kGreen+1,
  // 						attach_folder(folder,t1cfiles), pscuts+"&& mj14<850"));


  vector<shared_ptr<Process> > proc_pies;
  proc_pies.push_back(Process::MakeShared<Baby_full>("t#bar{t}(1l)", Process::Type::background, colors("tt_1l"),
						     attach_folder(folder,ttfiles), pscuts+"&& ntruleps==1"));
  proc_pies.push_back(Process::MakeShared<Baby_full>("t#bar{t}(2l)", Process::Type::background, colors("tt_2l"),
						     attach_folder(folder,ttfiles), pscuts+"&& ntruleps==2"));
  proc_pies.push_back(Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
						     attach_folder(folder,wfiles), pscuts));
  proc_pies.push_back(Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
						     attach_folder(folder,stfiles), pscuts));
  proc_pies.push_back(Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
						     attach_folder(folder,ofiles), pscuts));

  
  PlotOpt linear_shapes("txt/plot_styles.txt", "CMSPaper");
  linear_shapes.Title(TitleType::info)
    .YAxis(YAxisType::linear)
    .Bottom(BottomType::ratio)
    .Stack(StackType::shapes)
    .RatioMaximum(2.8);
  vector<PlotOpt> plot_types = {linear_shapes};


  PlotMaker pm;


  ///// Plotting 1D
  pm.Push<Hist1D>(Axis(15,100.,850.,"mj14","M_{J} [GeV]",{250, 400}),"nleps==1&&st>500&&met>200&&njets>=6&&nbm>=1",
		  procs, plot_types);


  ///// Plotting pies
  vector<TString> metcuts;
  metcuts.push_back("met>100 && met<=150");
  metcuts.push_back("met>150 && met<=200");
  metcuts.push_back("met>200 && met<=350");
  metcuts.push_back("met>350 && met<=500");
  metcuts.push_back("met>500");
  metcuts.push_back("met>200");

  vector<TString> nbcuts;
  nbcuts.push_back("nbm>=1");
  //nbcuts.push_back("nbm>=2");

  vector<TString> njcuts;
  njcuts.push_back("njets>=6");

  vector<TString> mtcuts({"mt<=140", "mt>140"});
  vector<TString> mjcuts({"mj14>250", "mj14>250&&mj14<=400", "mj14>400"});

  vector<TString> cuts;
  vector<TableRow> table_cuts;

  //// nleps = 1
  for(auto &imet: metcuts) 
    for(auto &inb: nbcuts) 
      for(auto &inj: njcuts) 
	for(auto &imj: mjcuts) 
	  for(auto &imt: mtcuts) {
	    cuts.push_back("nleps==1 && nveto==0 && "+imet+"&&"+inb+"&&"+inj+"&&"+imt+"&&"+imj);
	  }
  for(size_t icut=0; icut<cuts.size(); icut++)
    table_cuts.push_back(TableRow("$"+CodeToLatex(cuts[icut].Data())+"$", cuts[icut].Data()));  
  pm.Push<Table>("chart_full",  table_cuts, proc_pies, true, true, true, false);

  pm.min_print_ = true;
  pm.MakePlots(lumi);


  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Making "<<pm.Figures().size()<<" plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
}

