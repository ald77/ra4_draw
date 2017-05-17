///// table_preds: Makes piecharts

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
#include "core/event_scan.hpp"
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
  float lumi = 35.9;
  bool doSignal = true;
  bool csv = false;
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc(bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higtight/");
  string folderdata(bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_higdata_higloose/");
  string foldersig(bfolder+"/cms2r0/babymaker/babies/2017_02_26/TChiHH/merged_higmc_higtight/");

  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({"*TTJets_*Lept*", "*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
  mctags["qcd"]     = set<string>({"*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root", "*_ST_*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});

  string c_ps = "pass && stitch_met";

  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_full>("Other", Process::Type::background, kGreen+1,
							    attach_folder(foldermc, mctags["other"]),c_ps));
  procs.push_back(Process::MakeShared<Baby_full>("QCD", Process::Type::background, colors("other"),
							    attach_folder(foldermc, mctags["qcd"]),c_ps)); 
  procs.push_back(Process::MakeShared<Baby_full>("V+jets", Process::Type::background, kOrange+1,
							    attach_folder(foldermc,mctags["vjets"]),c_ps));
  procs.push_back(Process::MakeShared<Baby_full>("t#bar{t}+X", Process::Type::background,colors("tt_1l"),
                  attach_folder(foldermc, mctags["ttx"]),c_ps));
  procs.push_back(Process::MakeShared<Baby_full>("Data", Process::Type::data, 1,
                  {folderdata+"*root"}, Higfuncs::trig_hig>0. && "pass"));

  auto data = Process::MakeShared<Baby_full>("Data", Process::Type::data, 1,
                  {folderdata+"*root"}, Higfuncs::trig_hig>0. && "pass");
  procs.push_back(data);

  if (doSignal) {
    vector<string> sigm({"225","400", "700"});
    for (unsigned isig(0); isig<sigm.size(); isig++)
      procs.push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)", 
        Process::Type::signal, 1, {foldersig+"*TChiHH_mGluino-"+sigm[isig]+"*.root"}, "pass_goodv&&pass_ecaldeadcell&&pass_hbhe&&pass_hbheiso&&pass_fsmet"));
  }
 

  string filters = "pass_ra2_badmu && met/met_calo<5";
  string baseline = filters+"&& njets>=4 && njets<=5 && nvleps==0 && ntks==0 && !low_dphi";

  string c_2b = "nbdt==2&&nbdm==2";
  string c_3b = "nbdt>=2&&nbdm==3&&nbdl==3";
  string c_4b = "nbdt>=2&&nbdm>=3&&nbdl>=4";
  string hig = "higd_drmax<=2.2 && higd_am<=200 && higd_dm <= 40 && (higd_am>100 && higd_am<=140)";
  string sbd = "higd_drmax<=2.2 && higd_am<=200 && higd_dm <= 40 && !(higd_am>100 && higd_am<=140)";
  NamedFunc wgt = Higfuncs::weight_higd * Higfuncs::eff_higtrig;//*(1+Higfuncs::wgt_syst_ttx);// * Higfuncs::wgt_comp;
  if (csv) {
    c_2b = "nbt==2&&nbm==2";
    c_3b = "nbt>=2&&nbm==3&&nbl==3";
    c_4b = "nbt>=2&&nbm>=3&&nbl>=4";
    hig = "hig_drmax<=2.2 && hig_am<=200 && hig_dm <= 40 && (hig_am>100 && hig_am<=140)";
    sbd = "hig_drmax<=2.2 && hig_am<=200 && hig_dm <= 40 && !(hig_am>100 && hig_am<=140)";
    wgt = Higfuncs::weight_hig * Higfuncs::eff_higtrig;  
  }

  //        Cutflow table
  //-------------------------------- 
  PlotMaker pm;
  string tabname = "regions";
  if (csv) tabname +="_csv";
  pm.Push<Table>(tabname, vector<TableRow>{
  // TableRow("SBD, 2b", baseline + " && met>150 && met<=200 &&" +c_2b+"&&"+sbd,0,0, wgt),
  // TableRow("SBD, 3b", baseline + " && met>150 && met<=200 &&" +c_3b+"&&"+sbd,0,0, wgt),
  // TableRow("SBD, 4b", baseline + " && met>150 && met<=200 &&" +c_4b+"&&"+sbd,0,0, wgt),
  // TableRow("HIG, 2b", baseline + " && met>150 && met<=200 &&" +c_2b+"&&"+hig,0,1, wgt),
  // TableRow("HIG, 3b", baseline + " && met>150 && met<=200 &&" +c_3b+"&&"+hig,0,1, wgt),
  // TableRow("HIG, 4b", baseline + " && met>150 && met<=200 &&" +c_4b+"&&"+hig,0,1, wgt),
  
  // TableRow("$200\\leq E_{T}^{miss}<300$"),
  // TableRow("SBD, 2b", baseline + " && met>200 && met<=300 &&" +c_2b+"&&"+sbd,0,0, wgt),
  // TableRow("SBD, 3b", baseline + " && met>200 && met<=300 &&" +c_3b+"&&"+sbd,0,0, wgt),
  // TableRow("SBD, 4b", baseline + " && met>200 && met<=300 &&" +c_4b+"&&"+sbd,0,0, wgt),
  // TableRow("HIG, 2b", baseline + " && met>200 && met<=300 &&" +c_2b+"&&"+hig,0,1, wgt),
  // TableRow("HIG, 3b", baseline + " && met>200 && met<=300 &&" +c_3b+"&&"+hig,0,1, wgt),
  // TableRow("HIG, 4b", baseline + " && met>200 && met<=300 &&" +c_4b+"&&"+hig,0,1, wgt),
  
  // TableRow("$300\\leq E_{T}^{miss}<450$"),
  // TableRow("SBD, 2b", baseline + " && met>300 && met<=450 &&" +c_2b+"&&"+sbd,0,0, wgt),
  // TableRow("SBD, 3b", baseline + " && met>300 && met<=450 &&" +c_3b+"&&"+sbd,0,0, wgt),
  // TableRow("SBD, 4b", baseline + " && met>300 && met<=450 &&" +c_4b+"&&"+sbd,0,0, wgt),
  // TableRow("HIG, 2b", baseline + " && met>300 && met<=450 &&" +c_2b+"&&"+hig,0,1, wgt),
  // TableRow("HIG, 3b", baseline + " && met>300 && met<=450 &&" +c_3b+"&&"+hig,0,1, wgt),
  // TableRow("HIG, 4b", baseline + " && met>300 && met<=450 &&" +c_4b+"&&"+hig,0,1, wgt),
  
  // TableRow("$E_{T}^{miss}>450$"),
  // TableRow("SBD, 2b", baseline + " && met>450 &&"             +c_2b+"&&"+sbd,0,0, wgt),
  // TableRow("SBD, 3b", baseline + " && met>450 &&"             +c_3b+"&&"+sbd,0,0, wgt),
  // TableRow("SBD, 4b", baseline + " && met>450 &&"             +c_4b+"&&"+sbd,0,0, wgt),
  // TableRow("HIG, 2b", baseline + " && met>450 &&"             +c_2b+"&&"+hig,0,1, wgt),
  // TableRow("HIG, 3b", baseline + " && met>450 &&"             +c_3b+"&&"+hig,0,1, wgt),
  // TableRow("HIG, 4b", baseline + " && met>450 &&"             +c_4b+"&&"+hig,0,1, wgt),

  TableRow("$150\\leq E_{T}^{miss}<200$"),
  TableRow("SBD, 2b", baseline + " && met>150 && met<=200 &&" +c_2b+"&&"+sbd,0,0, wgt),
  TableRow("HIG, 2b", baseline + " && met>150 && met<=200 &&" +c_2b+"&&"+hig,0,1, wgt),
  TableRow("SBD, 3b", baseline + " && met>150 && met<=200 &&" +c_3b+"&&"+sbd,0,0, wgt),
  TableRow("HIG, 3b", baseline + " && met>150 && met<=200 &&" +c_3b+"&&"+hig,0,1, wgt),
  TableRow("SBD, 4b", baseline + " && met>150 && met<=200 &&" +c_4b+"&&"+sbd,0,0, wgt),
  TableRow("HIG, 4b", baseline + " && met>150 && met<=200 &&" +c_4b+"&&"+hig,0,1, wgt),
  
  TableRow("$200\\leq E_{T}^{miss}<300$"),
  TableRow("SBD, 2b", baseline + " && met>200 && met<=300 &&" +c_2b+"&&"+sbd,0,0, wgt),
  TableRow("HIG, 2b", baseline + " && met>200 && met<=300 &&" +c_2b+"&&"+hig,0,1, wgt),
  TableRow("SBD, 3b", baseline + " && met>200 && met<=300 &&" +c_3b+"&&"+sbd,0,0, wgt),
  TableRow("HIG, 3b", baseline + " && met>200 && met<=300 &&" +c_3b+"&&"+hig,0,1, wgt),
  TableRow("SBD, 4b", baseline + " && met>200 && met<=300 &&" +c_4b+"&&"+sbd,0,0, wgt),
  TableRow("HIG, 4b", baseline + " && met>200 && met<=300 &&" +c_4b+"&&"+hig,0,1, wgt),
  
  TableRow("$300\\leq E_{T}^{miss}<450$"),
  TableRow("SBD, 2b", baseline + " && met>300 && met<=450 &&" +c_2b+"&&"+sbd,0,0, wgt),
  TableRow("HIG, 2b", baseline + " && met>300 && met<=450 &&" +c_2b+"&&"+hig,0,1, wgt),
  TableRow("SBD, 3b", baseline + " && met>300 && met<=450 &&" +c_3b+"&&"+sbd,0,0, wgt),
  TableRow("HIG, 3b", baseline + " && met>300 && met<=450 &&" +c_3b+"&&"+hig,0,1, wgt),
  TableRow("SBD, 4b", baseline + " && met>300 && met<=450 &&" +c_4b+"&&"+sbd,0,0, wgt),
  TableRow("HIG, 4b", baseline + " && met>300 && met<=450 &&" +c_4b+"&&"+hig,0,1, wgt),
  
  TableRow("$E_{T}^{miss}>450$"),
  TableRow("SBD, 2b", baseline + " && met>450 &&"             +c_2b+"&&"+sbd,0,0, wgt),
  TableRow("HIG, 2b", baseline + " && met>450 &&"             +c_2b+"&&"+hig,0,1, wgt),
  TableRow("SBD, 3b", baseline + " && met>450 &&"             +c_3b+"&&"+sbd,0,0, wgt),
  TableRow("HIG, 3b", baseline + " && met>450 &&"             +c_3b+"&&"+hig,0,1, wgt),
  TableRow("SBD, 4b", baseline + " && met>450 &&"             +c_4b+"&&"+sbd,0,0, wgt),
  TableRow("HIG, 4b", baseline + " && met>450 &&"             +c_4b+"&&"+hig,0,1, wgt),
	},procs,0);

  string listname = "eventlist";
  if (csv) listname +="_csv";

  pm.Push<EventScan>(listname, baseline && "met>150 && nbdt>=2 && nbdm>=3" && hig, 
                     vector<NamedFunc>{"run", "lumiblock", "event"},
                     vector<shared_ptr<Process> >{data});

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making cutflow took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}


void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"no_signal", no_argument, 0, 'n'},    
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s:n", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'n':
      doSignal = false;
      break;
    case 0:
      // optname = long_options[option_index].name;
      // if(optname == "note"){
      //   note = true;
      // }else{
      //   printf("Bad option! Found option name %s\n", optname.c_str());
      //   exit(1);
      // }
      // break;

    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
