///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

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

namespace{
  float lumi = 36.2;
  // choose processes to include, options are: "ttx", "vjets", "singlet", "qcd", "other", "ttonly"
  set<string> proc_types = {"ttx", "vjets", "singlet", "qcd", "other"}; // for default data/MC
  // set<string> proc_types = {}; // to make signal plots only
  // set<string> proc_types = {"ttonly"};
  // signal points to include and their colors
  //for signal plots only
  vector<string> sigm = {"400"};//,"400","750","1000"}; 
  // vector<int> sig_colors = {kMagenta+2 , kGreen, kRed, kBlue, kAzure+10}; // need sigm.size() >= sig_colors.size()
}
  
int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  set<string> skims;
  skims = {"lep0"};

  map<string, string> foldermc; //ordered in number of leptons
  foldermc["lep0"] = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/";
  foldermc["lep1"] = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higlep1/";
  foldermc["lep2"] = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higlep2/";

  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({"*_TTJets*Lept*.root", "*_TTJets_HT*.root", "*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*_ttHJetTobb*.root","*_TTTT*.root"});
  mctags["ttonly"]  = set<string>({"*_TTJets*Lept*.root", "*_TTJets_HT*.root"});
  mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
  mctags["singlet"] = set<string>({"*_ST_*.root"});
  mctags["qcd"]     = set<string>({"*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  set<string> allmctags;
  for (auto &iset: mctags) {
    if (proc_types.find(iset.first)!=proc_types.end()) 
      allmctags.insert(iset.second.begin(), iset.second.end());
  }
  string foldersig = bfolder+"/cms2r0/babymaker/babies/2016_08_10/TChiHH/merged_higmc_unskimmed/";
  foldersig += "*TChiHH_mGluino-";

  string c_ps = "pass && stitch";

  map<string, vector<shared_ptr<Process> > > procs;
  for (auto &iskim: skims){
    procs[iskim] = vector<shared_ptr<Process> >();
    if (proc_types.find("other")!=proc_types.end())       
      procs[iskim].push_back(Process::MakeShared<Baby_full>("Other", Process::Type::background, kGreen+1,
							    attach_folder(foldermc[iskim], mctags["other"]),c_ps));
    if (proc_types.find("singlet")!=proc_types.end())       
      procs[iskim].push_back(Process::MakeShared<Baby_full>("Single t", Process::Type::background, 1,
							    attach_folder(foldermc[iskim],mctags["singlet"]),c_ps));
    if (proc_types.find("qcd")!=proc_types.end())       
      procs[iskim].push_back(Process::MakeShared<Baby_full>("QCD", Process::Type::background, colors("other"),
							    attach_folder(foldermc[iskim], mctags["qcd"]),c_ps)); 
    if (proc_types.find("vjets")!=proc_types.end())       
      procs[iskim].push_back(Process::MakeShared<Baby_full>("V+jets", Process::Type::background, kOrange+1,
							    attach_folder(foldermc[iskim],mctags["vjets"]),c_ps));
    if (proc_types.find("ttx")!=proc_types.end()) 
      procs[iskim].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
							    attach_folder(foldermc[iskim], mctags["ttx"]),c_ps));
    if (proc_types.find("ttonly")!=proc_types.end()) 
      procs[iskim].push_back(Process::MakeShared<Baby_full>("t#bar{t}", Process::Type::background, colors("tt_1l"),
							    attach_folder(foldermc[iskim],mctags["ttonly"]),c_ps));
    if (iskim == "lep0") {
      procs["sig"+iskim] = vector<shared_ptr<Process> >(procs[iskim]);
      if (proc_types.size()==0) { // have to pretend signal is background, otherwise crashes
       for (unsigned isig(0); isig<sigm.size(); isig++)
         procs["sig"+iskim].push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)",
          Process::Type::background, 2,
          {foldersig+sigm[isig]+"*.root"}, "1"));
      } else {
       for (unsigned isig(0); isig<sigm.size(); isig++)
         procs["sig"+iskim].push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)", 
								      Process::Type::signal, 2,
								      {foldersig+sigm[isig]+"*.root"}, "1"));
      }
    }
  }

 

  PlotMaker pm;

  //just to be pretty... already in skim...
  string njcut = "njets>=4 && njets<=5";
  string c_2b = "nbt==2&&nbm==2";
  string c_3b = "nbt>=2&&nbm==3&&nbl==3";
  string c_4b = "nbt>=2&&nbm>=3&&nbl>=4";

  map<string, string> xcuts; // useful additional cut definitions
  xcuts["fullhig"] = "hig_drmax<=2.2 && hig_am<=200 && hig_dm <= 40 && (hig_am>100 && hig_am<=140)";
  xcuts["fullsbd"] = "hig_drmax<=2.2 && hig_am<=200 && hig_dm <= 40 && !(hig_am>100 && hig_am<=140)";

  // string baseline = "pass && pass_ra2_badmu && stitch && nvleps==0 && met>150 && met/met_calo<5";
  string baseline = "pass_ra2_badmu && met/met_calo<5 && nvleps==0 && ntks==0 && met>150 && !(low_dphi)";
  string himet = "400";

  //        Cutflow table
  //-------------------------------- 
  NamedFunc wgt = "weight" * Higfuncs::eff_higtrig;
  pm.Push<Table>("cutflow", vector<TableRow>{
  TableRow("HIG, 2b, $150<E_{T}^{miss}\\leq$200",    baseline + " && met<=200 &&"                 +c_2b+"&&"+njcut+"&&"+xcuts["fullhig"],0,0, wgt),
  TableRow("HIG, 2b, $200<E_{T}^{miss}\\leq$300",    baseline + " && met>200 && met<=300 &&"      +c_2b+"&&"+njcut+"&&"+xcuts["fullhig"],0,0, wgt),
  TableRow("HIG, 2b, $300<E_{T}^{miss}\\leq$"+himet, baseline + " && met>300 && met<="+himet+" &&"+c_2b+"&&"+njcut+"&&"+xcuts["fullhig"],0,0, wgt),
  TableRow("HIG, 2b, $E_{T}^{miss}>$"+himet,         baseline + " && met>"+himet+" &&"            +c_2b+"&&"+njcut+"&&"+xcuts["fullhig"],0,1, wgt),

  TableRow("HIG, 3b, $150<E_{T}^{miss}\\leq$200",    baseline + " && met<=200 &&"                 +c_3b+"&&"+njcut+"&&"+xcuts["fullhig"],0,0, wgt),
  TableRow("HIG, 3b, $200<E_{T}^{miss}\\leq$300",    baseline + " && met>200 && met<=300 &&"      +c_3b+"&&"+njcut+"&&"+xcuts["fullhig"],0,0, wgt),
  TableRow("HIG, 3b, $300<E_{T}^{miss}\\leq$"+himet, baseline + " && met>300 && met<="+himet+" &&"+c_3b+"&&"+njcut+"&&"+xcuts["fullhig"],0,0, wgt),
  TableRow("HIG, 3b, $E_{T}^{miss}>$"+himet,         baseline + " && met>"+himet+" &&"            +c_3b+"&&"+njcut+"&&"+xcuts["fullhig"],0,1, wgt),

  TableRow("HIG, 4b, $150<E_{T}^{miss}\\leq$200",    baseline + " && met<=200 &&"                 +c_4b+"&&"+njcut+"&&"+xcuts["fullhig"],0,0, wgt),
  TableRow("HIG, 4b, $200<E_{T}^{miss}\\leq$300",    baseline + " && met>200 && met<=300 &&"      +c_4b+"&&"+njcut+"&&"+xcuts["fullhig"],0,0, wgt),
  TableRow("HIG, 4b, $300<E_{T}^{miss}\\leq$"+himet, baseline + " && met>300 && met<="+himet+" &&"+c_4b+"&&"+njcut+"&&"+xcuts["fullhig"],0,0, wgt),
  TableRow("HIG, 4b, $E_{T}^{miss}>$"+himet,         baseline + " && met>"+himet+" &&"            +c_4b+"&&"+njcut+"&&"+xcuts["fullhig"],0,2, wgt),

  TableRow("SBD, 2b, $150<E_{T}^{miss}\\leq$200",    baseline + " && met<=200 &&"                 +c_2b+"&&"+njcut+"&&"+xcuts["fullsbd"],0,0, wgt),
  TableRow("SBD, 2b, $200<E_{T}^{miss}\\leq$300",    baseline + " && met>200 && met<=300 &&"      +c_2b+"&&"+njcut+"&&"+xcuts["fullsbd"],0,0, wgt),
  TableRow("SBD, 2b, $300<E_{T}^{miss}\\leq$"+himet, baseline + " && met>300 && met<="+himet+" &&"+c_2b+"&&"+njcut+"&&"+xcuts["fullsbd"],0,0, wgt),
  TableRow("SBD, 2b, $E_{T}^{miss}>$"+himet,         baseline + " && met>"+himet+" &&"            +c_2b+"&&"+njcut+"&&"+xcuts["fullsbd"],0,1, wgt),

  TableRow("SBD, 3b, $150<E_{T}^{miss}\\leq$200",    baseline + " && met<=200 &&"                 +c_3b+"&&"+njcut+"&&"+xcuts["fullsbd"],0,0, wgt),
  TableRow("SBD, 3b, $200<E_{T}^{miss}\\leq$300",    baseline + " && met>200 && met<=300 &&"      +c_3b+"&&"+njcut+"&&"+xcuts["fullsbd"],0,0, wgt),
  TableRow("SBD, 3b, $300<E_{T}^{miss}\\leq$"+himet, baseline + " && met>300 && met<="+himet+" &&"+c_3b+"&&"+njcut+"&&"+xcuts["fullsbd"],0,0, wgt),
  TableRow("SBD, 3b, $E_{T}^{miss}>$"+himet,         baseline + " && met>"+himet+" &&"            +c_3b+"&&"+njcut+"&&"+xcuts["fullsbd"],0,1, wgt),

  TableRow("SBD, 4b, $150<E_{T}^{miss}\\leq$200",    baseline + " && met<=200 &&"                 +c_4b+"&&"+njcut+"&&"+xcuts["fullsbd"],0,0, wgt),
  TableRow("SBD, 4b, $200<E_{T}^{miss}\\leq$300",    baseline + " && met>200 && met<=300 &&"      +c_4b+"&&"+njcut+"&&"+xcuts["fullsbd"],0,0, wgt),
  TableRow("SBD, 4b, $300<E_{T}^{miss}\\leq$"+himet, baseline + " && met>300 && met<="+himet+" &&"+c_4b+"&&"+njcut+"&&"+xcuts["fullsbd"],0,0, wgt),
  TableRow("SBD, 4b, $E_{T}^{miss}>$"+himet,         baseline + " && met>"+himet+" &&"            +c_4b+"&&"+njcut+"&&"+xcuts["fullsbd"],0,2, wgt),

	},procs["siglep0"],0);


  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
