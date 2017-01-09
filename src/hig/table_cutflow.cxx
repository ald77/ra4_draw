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
  vector<string> sigm = {"225","400","700"}; 
}
  
int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/";
  string foldersig = bfolder+"/cms2r0/babymaker/babies/2016_08_10/TCHiHH/merged_higmc_unskimmed/";
  foldersig += "*TChiHH_mGluino-";

  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({"*_TTJets*Lept*.root", "*_TTJets_HT*.root", "*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*_ttHJetTobb*.root","*_TTTT*.root"});
  mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
  mctags["singlet"] = set<string>({"*_ST_*.root"});
  mctags["qcd"]     = set<string>({"*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  
  string c_ps = "pass && stitch";
  vector<shared_ptr<Process> > procs = vector<shared_ptr<Process> >();
  procs.push_back(Process::MakeShared<Baby_full>("Other", Process::Type::background, 1,
					    attach_folder(foldermc, mctags["other"]),c_ps));
  procs.push_back(Process::MakeShared<Baby_full>("Single t", Process::Type::background, 1,
					    attach_folder(foldermc,mctags["singlet"]),c_ps));
  procs.push_back(Process::MakeShared<Baby_full>("QCD", Process::Type::background, 1,
					    attach_folder(foldermc, mctags["qcd"]),c_ps)); 
  procs.push_back(Process::MakeShared<Baby_full>("V+jets", Process::Type::background, 1,
					    attach_folder(foldermc,mctags["vjets"]),c_ps));
  procs.push_back(Process::MakeShared<Baby_full>("t#bar{t}+X", Process::Type::background,1,
					    attach_folder(foldermc, mctags["ttx"]),c_ps));
  for (unsigned isig(0); isig<sigm.size(); isig++)
    procs.push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)", 
      Process::Type::signal, 1, {foldersig+sigm[isig]+"*.root"}, "1"));

  string c_2b = "nbt==2&&nbm==2";
  string c_3b = "nbt>=2&&nbm==3&&nbl==3";
  string c_4b = "nbt>=2&&nbm>=3&&nbl>=4";
  string c_ge3b = "nbt>=2&&nbm>=3";

  string hig = "hig_drmax<=2.2 && hig_am<=200 && hig_dm <= 40 && (hig_am>100 && hig_am<=140)";
  string sbd = "hig_drmax<=2.2 && hig_am<=200 && hig_dm <= 40 && !(hig_am>100 && hig_am<=140)";

  string baseline = "pass_ra2_badmu && met/met_calo<5 && nvleps==0 && njets>=4 && njets<=5";
  string sigonly = "type>100e3";
  NamedFunc wgt = "weight" * Higfuncs::eff_higtrig;

  string ncols = to_string(procs.size()+2);  
  PlotMaker pm;
  pm.Push<Table>("cutflow", vector<TableRow>{
  TableRow("No selection                ", 
    sigonly,0,0, wgt),
  // TableRow("Filters                     ", 
  //   sigonly+"&& pass_ra2_badmu && met/met_calo<5",0,0, wgt),
  TableRow("$0\\ell$, $\\text{4-5 jets}$  ", 
    baseline+"&&"+sigonly,0,0, wgt),
  TableRow("$N_{\\text{b,T}}\\geq 2$      ", 
    baseline+"&&"+sigonly + " && nbt>=2",0,0, wgt),
	TableRow("$E_{T}^{miss} > 150$", 
		baseline + " && met>150 && nbt>=2",0,0, wgt),
	TableRow("Track veto", 
		baseline + " && ntks==0 && met>150 && nbt>=2",0,0, wgt),
	TableRow("$\\Delta\\phi_{1,2}>0.5,\\Delta\\phi_{3,4}>0.3$",        
		baseline + " && ntks==0 && met>150 && nbt>=2   && !low_dphi",0,0, wgt),
	TableRow("$|\\Delta m| < 40$",     
		baseline + " && ntks==0 && met>150 && nbt>=2   && !low_dphi && hig_dm<=40",0,0, wgt),
	TableRow("$\\Delta R_{\\text{max}} < 2.2$",                    
		baseline +"  && ntks==0 && met>150 && nbt>=2   && !low_dphi && hig_dm<=40 && hig_drmax<=2.2",0,1,wgt),

  TableRow("\\multicolumn{"+ncols+"}{c}{HIG: $100<\\left< m \\right>\\leq140$}\\\\%", 
    "met>1e6",0,1, wgt),

  TableRow("HIG", 
    baseline + " && ntks==0 && met>150 && nbt>=2   && !low_dphi &&"+hig,0,0, wgt),
	TableRow("3b + 4b", 
    baseline +"  && ntks==0 && met>150 &&"+c_ge3b+"&& !low_dphi &&"+hig,0,0,wgt),
  TableRow("4b", 
    baseline + " && ntks==0 && met>150 &&"+c_4b+"  && !low_dphi &&"+hig,0,0, wgt),
  TableRow("$E_{T}^{miss}>200$", 
    baseline + " && ntks==0 && met>200 &&"+c_4b+"  && !low_dphi &&"+hig,0,0,wgt),
  TableRow("$E_{T}^{miss}>300$", 
    baseline + " && ntks==0 && met>300 &&"+c_4b+"  && !low_dphi &&"+hig,0,0,wgt),
  TableRow("$E_{T}^{miss}>450$", 
    baseline + " && ntks==0 && met>450 &&"+c_4b+"  && !low_dphi &&"+hig,0,1,wgt),

  TableRow("\\multicolumn{"+ncols+"}{c}{SBD: $\\left< m\\right> <100$ or $140<\\left< m\\right>\\leq200$}\\\\%", 
    "met>1e6",0,1, wgt),

  TableRow("SBD", 
    baseline + " && ntks==0 && met>150 && nbt>=2   && !low_dphi &&"+sbd,0,0, wgt),
  TableRow("3b + 4b", 
    baseline +"  && ntks==0 && met>150 &&"+c_ge3b+"&& !low_dphi &&"+sbd,0,0,wgt),
  TableRow("4b", 
    baseline + " && ntks==0 && met>150 &&"+c_4b+"  && !low_dphi &&"+sbd,0,0, wgt),
  TableRow("$E_{T}^{miss}>200$", 
    baseline + " && ntks==0 && met>200 &&"+c_4b+"  && !low_dphi &&"+sbd,0,0,wgt),
  TableRow("$E_{T}^{miss}>300$", 
    baseline + " && ntks==0 && met>300 &&"+c_4b+"  && !low_dphi &&"+sbd,0,0,wgt),
  TableRow("$E_{T}^{miss}>450$", 
    baseline + " && ntks==0 && met>450 &&"+c_4b+"  && !low_dphi &&"+sbd,0,0,wgt),

	},procs,0);


  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
