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
  float lumi = 35.9;
  bool deep = true;
  // vector<string> sigm = {}; 
  vector<string> sigm = {"225","700"}; 
}
  
int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc(bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higloose/");
  string foldersig(bfolder+"/cms2r0/babymaker/babies/2017_02_26/TChiHH/merged_higmc_unskimmed/");

  map<string, set<string>> mctags; 
  mctags["all"]     = set<string>({"*TTJets_*Lep*", "*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root",
                                     "*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root",
                                     "*_ST_*.root",
                                     "*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root",
                                     "*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  
  string c_ps = "pass && stitch_met";
  vector<shared_ptr<Process> > procs = vector<shared_ptr<Process> >();
  procs.push_back(Process::MakeShared<Baby_full>("All bkg.", Process::Type::background, 1,
					    attach_folder(foldermc, mctags["all"]),c_ps));
  for (unsigned isig(0); isig<sigm.size(); isig++)
    procs.push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)", 
      Process::Type::signal, 1, {foldersig+"*TChiHH_mGluino-"+sigm[isig]+"*.root"}, "1"));


  string c_2bt = "nbt>=2";
  string c_2b = "nbt==2&&nbm==2";
  string c_3b = "nbt>=2&&nbm==3&&nbl==3";
  string c_4b = "nbt>=2&&nbm>=3&&nbl>=4";
  string c_ge3b = "nbt>=2&&nbm>=3";

  string c_hig_dm = "hig_dm<=40";
  string c_drmax = "hig_drmax<=2.2";
  string hig = "hig_drmax<=2.2 && hig_am<=200 && hig_dm <= 40 && (hig_am>100 && hig_am<=140)";
  string sbd = "hig_drmax<=2.2 && hig_am<=200 && hig_dm <= 40 && !(hig_am>100 && hig_am<=140)";
  // string trim = "hig_drmax<=2.2 && hig_am<=200 && hig_dm <= 40";
  string trim = "hig_drmax<=2.2 && hig_dm <= 40";

  NamedFunc wgt = Higfuncs::weight_hig * Higfuncs::eff_higtrig;

  if (deep) {
    c_2bt  = "nbdt>=2";
    c_2b   = "nbdt==2&&nbdm==2";
    c_3b   = "nbdt>=2&&nbdm==3&&nbdl==3";
    c_4b   = "nbdt>=2&&nbdm>=3&&nbdl>=4";
    c_ge3b = "nbdt>=2&&nbdm>=3";

    c_hig_dm = "higd_dm<=40";
    c_drmax  = "higd_drmax<=2.2";
    hig = "higd_drmax<=2.2 && higd_am<=200 && higd_dm <= 40 && (higd_am>100 && higd_am<=140)";
    sbd = "higd_drmax<=2.2 && higd_am<=200 && higd_dm <= 40 && !(higd_am>100 && higd_am<=140)";
    // trim = "higd_drmax<=2.2 && higd_am<=200 && higd_dm <= 40";
    trim = "higd_drmax<=2.2 && higd_dm <= 40";

    wgt = Higfuncs::weight_higd * Higfuncs::eff_higtrig;
  }

  string baseline = "pass_ra2_badmu && met/met_calo<5 && nvleps==0 && njets>=4 && njets<=5";
  string sigonly = "type>100e3";

  string ncols = to_string(procs.size()+2);  
  PlotMaker pm;
  string name = "csv";
  if (deep) name = "deepcsv";
  pm.Push<Table>(name, vector<TableRow>{
  TableRow("No pre-selection, $\\geq$ 2b", 
    sigonly+"&&"+c_2bt,0,0, Higfuncs::weight_higd * "1/eff_jetid"),
  TableRow("No pre-selection, $\\geq$ 3b",  
    sigonly+"&&"+c_ge3b,0,0, Higfuncs::weight_higd * "1/eff_jetid"),
  TableRow("No pre-selection, 4b", 
    sigonly+"&&"+c_4b,0,0, Higfuncs::weight_higd * "1/eff_jetid"),
  TableRow("Baseline, $p_{\\rm T}^{\\rm miss}>150$, $\\geq$ 2b", 
    baseline + " && ntks==0 && met>150 &&"+c_2bt+" && !low_dphi && "+trim,0,0, wgt),
	TableRow("Baseline, $p_{\\rm T}^{\\rm miss}>150$, $\\geq$ 3b", 
    baseline +"  && ntks==0 && met>150 &&"+c_ge3b+"&& !low_dphi && "+trim,0,0,wgt),
  TableRow("Baseline, $p_{\\rm T}^{\\rm miss}>150$, 4b", 
    baseline + " && ntks==0 && met>150 &&"+c_4b+"  && !low_dphi && "+trim,0,0, wgt),
  TableRow("Baseline, $p_{\\rm T}^{\\rm miss}>300$, $\\geq$ 2b", 
    baseline + " && ntks==0 && met>300 &&"+c_2bt+" && !low_dphi && "+trim,0,0, wgt),
  TableRow("Baseline, $p_{\\rm T}^{\\rm miss}>300$, $\\geq$ 3b", 
    baseline +"  && ntks==0 && met>300 &&"+c_ge3b+"&& !low_dphi && "+trim,0,0,wgt),
  TableRow("Baseline, $p_{\\rm T}^{\\rm miss}>300$, 4b", 
    baseline + " && ntks==0 && met>300 &&"+c_4b+"  && !low_dphi && "+trim,0,0, wgt),
	},procs,0);


  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
