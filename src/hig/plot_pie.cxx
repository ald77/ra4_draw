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

using namespace std;

namespace{
  float lumi = 36.;
  bool do_procs = true;
  bool do_cats_ntrub = false;
  vector<string> selns= {"nom","qcd","onelep","dilep"};
  // vector<string> selns= {"qcd","dilep"};
  enum proc_types{ttx, vjets, singlet, qcd, other, nprocs};
}
  
int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  NamedFunc nb_tru("nb_tru",[](const Baby &b) -> NamedFunc::ScalarType{
    int nbtru(0);
    for (unsigned i(0); i<b.jets_pt()->size(); i++){
      if (!b.jets_h1()->at(i) && !b.jets_h2()->at(i)) continue;
      if (b.jets_hflavor()->at(i)==5) nbtru++;
    }
    return nbtru;
  });

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  map<string, string> folders; 
  folders["nom"] = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/";
  folders["qcd"] = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/";
  folders["onelep"] = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_hig_nlep1/";
  folders["dilep"] = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_hig_nlep2/";

  // Cuts in baseline speed up the yield finding
  map<string, string> baseline;
  baseline["nom"]   = "pass && stitch && njets>=4 && njets<=5 && nvleps==0 && ntks==0 && !low_dphi";
  baseline["qcd"]   = "pass && stitch && njets>=4 && njets<=5 && nvleps==0 && ntks==0 && low_dphi";
  baseline["onelep"] = "pass && stitch && njets>=4 && njets<=5 && nleps==1 && mt<=140";
  //z-mass window already in the skim
  baseline["dilep"] = "pass && stitch && njets>=4 && njets<=5 && nleps==2 && (mumu_m*(mumu_m>0)+elel_m*(elel_m>0))>80&&(mumu_m*(mumu_m>0)+elel_m*(elel_m>0))<100";

  Palette colors("txt/colors.txt", "default");

  vector<set<string>> files; files.resize(nprocs);
  files[ttx]     = set<string>({"*_TTJets*Lept*.root", "*_TTJets_HT*.root", "*_TTZ*.root", "*_TTW*.root",
                                "*_TTGJets*.root", "*_ttHJetTobb*.root","*_TTTT*.root"});
  files[vjets]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
  files[singlet] = set<string>({"*_ST_*.root"});
  files[qcd]     = set<string>({"*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root"});
  files[other]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});

  map<string, vector<shared_ptr<Process> > > procs;
  if (do_procs) {  
    for (auto &iseln: selns){
      procs[iseln] = vector<shared_ptr<Process> >();
      procs[iseln].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X", Process::Type::background, colors("tt_1l"),
        attach_folder(folders[iseln], files[ttx]), baseline[iseln]));
      procs[iseln].push_back(Process::MakeShared<Baby_full>("V+jets", Process::Type::background, kOrange+1,
        attach_folder(folders[iseln], files[vjets]), baseline[iseln]));
      procs[iseln].push_back(Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
        attach_folder(folders[iseln], files[singlet]), baseline[iseln]));
      procs[iseln].push_back(Process::MakeShared<Baby_full>("QCD", Process::Type::background, colors("other"),
        attach_folder(folders[iseln], files[qcd]), baseline[iseln]));
      procs[iseln].push_back(Process::MakeShared<Baby_full>("Other", Process::Type::background, kGreen+1,
        attach_folder(folders[iseln], files[other]), baseline[iseln]));
    }
  } 
  if (do_cats_ntrub) {
    set<string> allfiles;
    for (auto &iset: files) 
      allfiles.insert(iset.begin(), iset.end());
    for (auto &iseln: selns){
      set<string> full_fset = attach_folder(folders[iseln], allfiles);
      NamedFunc base_func(baseline[iseln]);
      procs["cats"+iseln] = vector<shared_ptr<Process> >();
      procs["cats"+iseln].push_back(Process::MakeShared<Baby_full>
              ("#leq 1 B-hadron", Process::Type::background, kPink+2,
               full_fset, base_func && nb_tru<=1));
      procs["cats"+iseln].push_back(Process::MakeShared<Baby_full>
      			  ("2 B-hadrons", Process::Type::background, kOrange-4,
      			   full_fset, base_func && nb_tru==2));
      procs["cats"+iseln].push_back(Process::MakeShared<Baby_full>
      			  ("3 B-hadrons", Process::Type::background, kTeal-8, 
               full_fset, base_func &&  nb_tru==3));
      procs["cats"+iseln].push_back(Process::MakeShared<Baby_full>
      			  ("#geq 4 B-hadrons", Process::Type::background, kAzure-4, 
      			   full_fset, base_func && nb_tru>=4));
    }
  }

  PlotMaker pm;

  map<string, vector<string> > metcuts;
  metcuts["nom"] = {"met>150&&met<=200", "met>200&&met<=300","met>300"};
  metcuts["qcd"] = metcuts["nom"];
  metcuts["onelep"] = metcuts["nom"];
  metcuts["dilep"] = {"(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))>50"};
  

  vector<TString> nbcuts;
  nbcuts.push_back("nbt==2&&nbm==2");
  nbcuts.push_back("nbt>=2&&nbm==3&&nbl==3");
  nbcuts.push_back("nbt>=2&&nbm>=3&&nbl>=4");

  vector<TString> regs;
  regs.push_back("1");
  regs.push_back("hig_am>100 && hig_am<=140 && hig_dm <= 40");
  regs.push_back("(hig_am<=100 || hig_am>140 || hig_dm > 40)");

  map<string, string> xcuts; // additional cut options, done for each selection
  xcuts["1"] = "1";
  xcuts["hdrmax"] = "hig_drmax<2.2";

  for (auto &iseln: selns) {
    for (auto &ixcut: xcuts) {
      vector<TString> cuts;
      vector<TableRow> table_cuts;
      for(auto &imet: metcuts[iseln]) { 
        for(auto &inb: nbcuts) {
          for (auto &ireg: regs) {
            cuts.push_back(baseline[iseln]+"&&"+imet+"&&"+inb+"&&"+ireg+"&&"+ixcut.second);
          }

        }
      }
      for(size_t icut=0; icut<cuts.size(); icut++)
        table_cuts.push_back(TableRow("$"+CodeToLatex(cuts[icut].Data())+"$", cuts[icut].Data()));  
      if (do_procs) pm.Push<Table>("chart_"+iseln+"_"+ixcut.first,  table_cuts, procs[iseln], true, true, true, false);
      if (do_cats_ntrub) 
        pm.Push<Table>("chartcats_"+iseln+"_"+ixcut.first,  table_cuts, procs["cats"+iseln], true, true, true, false);
    }
  }

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making piecharts took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
