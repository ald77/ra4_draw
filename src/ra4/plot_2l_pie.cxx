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

NamedFunc N_fail_acceptance("N_fail_acceptance",[](const Baby &b) -> NamedFunc::ScalarType{
    int nfail=0;
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if ((abs(b.mc_id()->at(i))==11 && (b.mc_pt()->at(i)<20 || b.mc_eta()->at(i)>2.5)) ||(abs(b.mc_id()->at(i))==13 && (b.mc_pt()->at(i)<20 || b.mc_eta()->at(i)>2.4)) ) nfail++;
      
    }
    return nfail;
  });

NamedFunc N_fail_iso("N_fail_iso",[](const Baby &b) -> NamedFunc::ScalarType{
    //(Sum$(els_sigid&&els_tm&&els_miniso>0.1)+Sum$(mus_sigid&&mus_tm&&mus_miniso>0.2))>=1
    int nfail=0;
    for (unsigned i(0); i<b.els_pt()->size(); i++){
      if (b.els_sigid()->at(i)&&b.els_tm()->at(i)&&b.els_miniso()->at(i)>0.1) nfail++;
    }
    for (unsigned i(0); i<b.mus_pt()->size(); i++){
      if (b.mus_sigid()->at(i)&&b.mus_tm()->at(i)&&b.mus_miniso()->at(i)>0.2) nfail++;
    }
    return nfail;
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

  string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/unskimmed/");
  
  Palette colors("txt/colors.txt", "default");

  string ntupletag = "";
  set<string> allfiles = {foldermc+"*_TTJets_DiLept*"+ntupletag+"*.root", foldermc+"*_TTJets_HT*"+ntupletag+"*.root",
       };

  // allfiles = set<string>({foldermc+"*_TTJets_Tune*"});

  // Cuts in baseline speed up the yield finding
  string baseline = "pass && stitch && mj14>250 && st>500 && met>200 && njets>=6 && nbm >= 1 && ntruleps>=2";

  map<string, vector<shared_ptr<Process> > > procs;
  map<string, vector<shared_ptr<Process> > > procs_no_sel;
  map<string, vector<shared_ptr<Process> > > procs_lost;
 

  // procs["procs"] = vector<shared_ptr<Process> >();
  //procs["procs"].push_back(Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
  //							  allfiles, baseline+" && ntruleps>=2 "));


  procs["nreco_cut"] = vector<shared_ptr<Process> >();
  procs["nreco_cut"].push_back(Process::MakeShared<Baby_full>
  			  ("2 reco leptons", Process::Type::background, kCyan-3,
  			   allfiles, baseline+"&&nleps>=2"));
  procs["nreco_cut"].push_back(Process::MakeShared<Baby_full>
  			  ("1 reco lepton", Process::Type::background, kRed-4, allfiles, 
  			   baseline + "&&nleps==1"));
  procs["nreco_cut"].push_back(Process::MakeShared<Baby_full>
  			  ("0 reco leptons", Process::Type::background, kGreen-3, 
  			   allfiles, baseline+"&&nleps==0"));

  procs_no_sel["nreco"] = vector<shared_ptr<Process> >();
  procs_no_sel["nreco"].push_back(Process::MakeShared<Baby_full>
  			  ("2 reco leptons", Process::Type::background, kCyan-3,
  			   allfiles, "nleps>=2"));
  procs_no_sel["nreco"].push_back(Process::MakeShared<Baby_full>
  			  ("1 reco lepton", Process::Type::background, kRed-4, 
			   allfiles, "nleps==1"));
  procs_no_sel["nreco"].push_back(Process::MakeShared<Baby_full>
  			  ("0 reco leptons", Process::Type::background, kGreen-3, 
  			   allfiles, "nleps==0"));
 

  string baseline_2l = "pass && stitch && mj14>250 && nleps==1 && ntruleps>=2 && st>500 && met>200 && njets>=6 && nbm >= 1";
  procs_lost["cats"] = vector<shared_ptr<Process> >();
  procs_lost["cats"].push_back(Process::MakeShared<Baby_full>
  			  ("Is a hadronic tau", Process::Type::background, kCyan-3,
  			   allfiles, baseline_2l+"&& ntrutaush>0"));
  procs_lost["cats"].push_back(Process::MakeShared<Baby_full>
  			  ("Fails acceptance", Process::Type::background, kRed-4, allfiles, 
  			   baseline_2l+"&&ntrutaush==0"&& N_fail_acceptance > 0.));
  procs_lost["cats"].push_back(Process::MakeShared<Baby_full>
  			  ("Passes acceptance and ID but fails isolation", Process::Type::background, kGreen-3, 
  			   allfiles, baseline+"&&ntrutaush==0" && N_fail_acceptance==0. && N_fail_iso > 0.));
  procs_lost["cats"].push_back(Process::MakeShared<Baby_full>
  			  ("Passes acceptance but fails ID", Process::Type::background, kOrange, 
  			   allfiles, baseline+"&&ntrutaush==0" &&  N_fail_acceptance==0. && N_fail_iso == 0.)); 



  PlotMaker pm;

  vector<TString> metcuts;
  /*  //metcuts.push_back("met>100");
  metcuts.push_back("met>100 && met<=150");
  metcuts.push_back("met>150 && met<=200");*/
  metcuts.push_back("met>200 && met<=350");
  metcuts.push_back("met>350 && met<=500");
  metcuts.push_back("met>500");
  metcuts.push_back("met>200");

  vector<TString> nbcuts;
  //nbcuts.push_back("nbm==0");
  nbcuts.push_back("nbm>=1");
  if(false){
    nbcuts.push_back("nbm==1");
    nbcuts.push_back("nbm==2");
    nbcuts.push_back("nbm>=3");
  }

  vector<TString> njcuts;
  njcuts.push_back("njets>=6");
  njcuts.push_back("njets>=6 && njets<=8");
  njcuts.push_back("njets>=9");
  
  vector<TString> mtcuts({"1","mt<=140", "mt>140"});
  
  // Adding nleps==1 cuts
  vector<TString> cuts;
  vector<TableRow> table_cuts;
  vector<TableRow> table_cuts_no_sel;
  //// nleps = 1
  for(auto &imet: metcuts) 
    for(auto &inb: nbcuts) 
      for(auto &inj: njcuts) 
	for(auto &imt: mtcuts) {
	  cuts.push_back(imet+"&&"+inb+"&&"+inj+"&&"+imt);
	  //cuts.push_back("nleps==1 && "+imet+"&&"+inb+"&&"+inj+"&&"+imt);
	  cuts.push_back("nveto==0 && "+imet+"&&"+inb+"&&"+inj+"&&"+imt);
	  cuts.push_back("nveto==1 && "+imet+"&&"+inb+"&&"+inj+"&&"+imt);
	}
  
  

  table_cuts_no_sel.push_back(TableRow("$"+CodeToLatex("1")+"$", "1"));  

  for(size_t icut=0; icut<cuts.size(); icut++)
    table_cuts.push_back(TableRow("$"+CodeToLatex(cuts[icut].Data())+"$", cuts[icut].Data()));  

  for(auto &ipr: procs) 
    pm.Push<Table>("chart_"+ipr.first,  table_cuts_no_sel, ipr.second, true, true, true, false);
  for(auto &ipr: procs_no_sel) 
    pm.Push<Table>("chart_"+ipr.first,  table_cuts_no_sel, ipr.second, true, true, true, false);
  for(auto &ipr: procs_lost) 
    pm.Push<Table>("chart_"+ipr.first,  table_cuts, ipr.second, true, true, true, false);

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making "<<table_cuts.size()<<" piecharts took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
