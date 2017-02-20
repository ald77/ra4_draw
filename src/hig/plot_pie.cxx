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
#include "core/plot_maker.hpp"
#include "core/slide_maker.hpp"
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
  //fixme:simplify options
  float lumi = 35.9;
  // options "zll", "qcd", "ttbar", "search"
  string sample = "search";
  bool do_trim = true;
  bool split_higsbd = false;
  bool note = true;
}
  
int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  NamedFunc ntrub("ntrub",[](const Baby &b) -> NamedFunc::ScalarType{
    int tmp_ntrub(0);
    for (unsigned i(0); i<b.jets_pt()->size(); i++){
      if (!b.jets_h1d()->at(i) && !b.jets_h2d()->at(i)) continue;
      if (b.jets_hflavor()->at(i)==5) tmp_ntrub++;
    }
    return tmp_ntrub;
  });

  PlotOpt lin_shapes_info("txt/plot_styles.txt", "CMSPaper");
  lin_shapes_info.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::linear)
    .Stack(StackType::shapes);
  PlotOpt log_shapes_info = lin_shapes_info;
  log_shapes_info.YAxis(YAxisType::log);
  vector<PlotOpt> plt_types = {lin_shapes_info, log_shapes_info};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc(bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higloose/");
  if (sample=="ttbar") foldermc = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higlep1/";
  if (sample=="zll") foldermc = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higlep2/";
  if (sample=="qcd") foldermc = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higqcd/";

  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({"*TTJets_*Lept*", "*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*ttHTobb*.root","*_TTTT*.root"});
  mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
  mctags["singlet"] = set<string>({"*_ST_*.root"});
  mctags["qcd"]     = set<string>({"*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  set<string> allmctags;
  for (auto &iset: mctags) {
      allmctags.insert(iset.second.begin(), iset.second.end());
  }

  // Baseline definitions
  NamedFunc wgt = Higfuncs::weight_higd * Higfuncs::eff_higtrig;
  NamedFunc base_func("pass && pass_ra2_badmu && stitch_met && met/met_calo<5 && njets>=4 && njets<=5");
  if (do_trim) base_func = base_func && "higd_dm<=40 && higd_am<=200";
  if (sample=="zll")    base_func = base_func && "met<50";
  if (sample=="qcd")    base_func = base_func && "ntks==0 && low_dphi && met>=150";
  if (sample=="ttbar")  base_func = base_func && "mt<100";
  if (sample=="search") base_func = base_func && "nveto==0 && ntks==0 && !low_dphi && met>=150";

  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_full>("t#bar{t}+X", 
    Process::Type::background, colors("tt_1l"),    attach_folder(foldermc,mctags["ttx"]),     base_func));
  procs.push_back(Process::MakeShared<Baby_full>("V+jets",     
    Process::Type::background, kOrange+1,          attach_folder(foldermc,mctags["vjets"]),   base_func));
  procs.push_back(Process::MakeShared<Baby_full>("Single t",   
    Process::Type::background, colors("single_t"), attach_folder(foldermc,mctags["singlet"]), base_func));
  procs.push_back(Process::MakeShared<Baby_full>("QCD",        
    Process::Type::background, colors("other"),    attach_folder(foldermc,mctags["qcd"]),     base_func)); 
  procs.push_back(Process::MakeShared<Baby_full>("Other",      
    Process::Type::background, kGreen+1,           attach_folder(foldermc,mctags["other"]),   base_func));      

  set<string> allfiles = attach_folder(foldermc, allmctags);
  vector<shared_ptr<Process> > procs_ntrub;
  procs_ntrub.push_back(Process::MakeShared<Baby_full>
    ("0 B-hadron",       Process::Type::background, kAzure-4,   allfiles, base_func && ntrub<1));
  procs_ntrub.push_back(Process::MakeShared<Baby_full>
    ("1 B-hadron",       Process::Type::background, kTeal-8,    allfiles, base_func && ntrub==1));
  procs_ntrub.push_back(Process::MakeShared<Baby_full>
    ("2 B-hadrons",      Process::Type::background, kOrange-4,  allfiles, base_func && ntrub==2));
  procs_ntrub.push_back(Process::MakeShared<Baby_full>
    ("3 B-hadrons",      Process::Type::background, kPink+2,    allfiles, base_func && ntrub==3));
  procs_ntrub.push_back(Process::MakeShared<Baby_full>
    ("4 B-hadrons", Process::Type::background, kMagenta-1, allfiles, base_func && ntrub==4));

  vector<string> xcuts; // useful additional cut definitions
  if (!note) xcuts.push_back("1");
  xcuts.push_back("higd_drmax<=2.2");
  
  PlotMaker pm;
  SlideMaker sm("slide_pies_"+sample+".tex","1610");
  vector<string> pnames;
  vector<TableRow> table_cuts;

  // pie charts for the "nb - additional cuts" plane
  //-------------------------------------------------
  vector<string> nbcuts;
  if (sample=="zll" || sample=="qcd") {
    nbcuts.push_back("nbdm==0");
    nbcuts.push_back("nbdm==1");
  }
  nbcuts.push_back("nbdt==2&&nbdm==2");
  if (sample=="ttbar" || sample=="search" || !note) {
    nbcuts.push_back("nbdt>=2&&nbdm==3&&nbdl==3");
    nbcuts.push_back("nbdt>=2&&nbdm>=3&&nbdl>=4");
  }

  string tag = "ntrub";
  for(auto &ixcut: xcuts) {
    for(auto &inb: nbcuts) {
      string icut = inb+"&&"+ixcut;
      table_cuts.push_back(TableRow("", icut, 0, 0, wgt));  
      pnames.push_back("pie_"+sample+"_"+tag+"_"+CodeToPlainText(icut)+"_perc_lumi"+RoundNumber(lumi,0).Data()+".pdf");
    }
  }
  pm.Push<Table>(sample+"_"+tag,  table_cuts, procs_ntrub, true, true, true);
  sm.AddSlide(pnames, nbcuts.size(), "X-axis: Number of b-tags, Y-axis: additional cuts");  
  //push the table with the same cuts but different procs, so also have to change the pie chart names
  if (sample=="search" || !note) {
    pm.Push<Table>(sample+"_procs",  table_cuts, procs, true, true, true);
    sm.AddSlideWithReplace(tag,"procs", pnames, nbcuts.size(), "X-axis: Number of b-tags, Y-axis: additional cuts");
  }
  
  // pie charts for the "met - nb cuts" plane, one slide per set of additional cuts
  //--------------------------------------------------------------------------------
  vector<string> metcuts;
  string metdef = "met";
  if (sample=="zll") metdef = "(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))";
  if (sample=="ttbar" || sample=="zll"){
    metcuts.push_back(metdef+">0&&"+metdef+"<=75");
    metcuts.push_back(metdef+">75&&"+metdef+"<=150");
  }
  metcuts.push_back(metdef+">150&&"+metdef+"<=200");
  metcuts.push_back(metdef+">200&&"+metdef+"<=300");
  // if (sample=="search") {
    metcuts.push_back(metdef+">300&&"+metdef+"<=450");
    metcuts.push_back(metdef+">450");
  // } else {
  //   metcuts.push_back(metdef+">300");
  // }

  table_cuts.clear();
  for(auto &ixcut: xcuts) {
    pnames.clear();
    string slide_ttl = "Additional cuts: $"+CodeToLatex(ixcut)+"$";
    if (ixcut=="1") slide_ttl = "No additional cuts";
    for(auto &inb: nbcuts) {
      for(auto &imet: metcuts) {
        string icut = inb+"&&"+imet+"&&"+ixcut;
        if (split_higsbd && !note) {
          table_cuts.push_back(TableRow("", icut+"&&!(higd_am>100&&higd_am<=140)", 0, 0, wgt));  
          pnames.push_back("pie_"+sample+"_procs_"+CodeToPlainText(icut+"&&!(higd_am>100&&higd_am<=140)")+
            "_perc_lumi"+RoundNumber(lumi,0).Data()+".pdf");
          table_cuts.push_back(TableRow("", icut+"&&(higd_am>100&&higd_am<=140)", 0, 0, wgt));  
          pnames.push_back("pie_"+sample+"_procs_"+CodeToPlainText(icut+"&&(higd_am>100&&higd_am<=140)")+
            "_perc_lumi"+RoundNumber(lumi,0).Data()+".pdf");
        } else {
          // procs
          table_cuts.push_back(TableRow("", icut, 0, 0, wgt));  
          pnames.push_back("pie_"+sample+"_procs_"+CodeToPlainText(icut)+"_perc_lumi"+RoundNumber(lumi,0).Data()+".pdf");
        }
      }
    }
    if (split_higsbd) sm.AddSlide(pnames, metcuts.size()*2, slide_ttl);
    else sm.AddSlide(pnames, metcuts.size(), slide_ttl);

    if (sample=="ttbar" || !note){
      if (!split_higsbd) sm.AddSlideWithReplace("procs", tag, pnames, metcuts.size(), slide_ttl);
    }
  }
  pm.Push<Table>(sample+"_procs",  table_cuts, procs, true, true, true);
  if (!split_higsbd && (sample=="ttbar" || !note))
      pm.Push<Table>(sample+"_"+tag,  table_cuts, procs_ntrub, true, true, true);

  pm.min_print_ = true;
  pm.multithreaded_ = true;

  pm.MakePlots(lumi);
  sm.Close();

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"bcats", required_argument, 0, 'b'},
      {"lumi", required_argument, 0, 'l'},    // Luminosity to normalize MC with (no data)
      {"sample", required_argument, 0, 's'},    // Which sample to use: standard, met150, 2015 data
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s:l:gn:d:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'l':
      lumi = atof(optarg);
      break;
    case 's':
      sample = optarg;
      break;
    case 0:
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
