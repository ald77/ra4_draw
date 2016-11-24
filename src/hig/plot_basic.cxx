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

using namespace std;
using namespace PlotOptTypes;

namespace{
  float lumi = 36.;
  bool do_cats_ntrub = false;
  // turn on/off control region plots
  bool do_cr_pies = false;
  bool do_cr_plots = false;
  // simplify selection
  bool do_onemet = true; // only met>150
  bool do_ge3b = false;  // do btag cats: 2b and 3b+
  bool do_23vs4b = true; // combine 2b with 3b, keep 4b separate
  // choose processes to include, options are: "ttx", "vjets", "singlet", "qcd", "other", "ttonly"
  // default bkg includes: {"ttx", "vjets", "singlet", "qcd", "other"}
  set<string> proc_types = {"ttx", "vjets", "singlet", "qcd", "other"};
  // set<string> proc_types = {"ttonly"};
  // signal points to include and their colors
  vector<string> sigm = {"225","300","400","700"}; 
  vector<int> sig_colors = {kRed, kGreen, kOrange+3, kBlue}; // need sigm.size() >= sig_colors.size()
}
  
int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  NamedFunc ntrub("ntrub",[](const Baby &b) -> NamedFunc::ScalarType{
    int tmp_ntrub(0);
    for (unsigned i(0); i<b.jets_pt()->size(); i++){
      if (!b.jets_h1()->at(i) && !b.jets_h2()->at(i)) continue;
      if (b.jets_hflavor()->at(i)==5) tmp_ntrub++;
    }
    return tmp_ntrub;
  });

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::off)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info);
  vector<PlotOpt> all_plot_types = {lin_lumi_info, log_lumi_info};

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  set<string> skims;
  if (do_cr_pies || do_cr_plots) skims = {"lep0","lep1","lep2"};
  else skims = {"lep0"};

  map<string, string> foldermc; //ordered in number of leptons
  foldermc["lep0"] = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/";
  foldermc["lep1"] = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_hig_nlep1/";
  foldermc["lep2"] = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_hig_nlep2/";

  string foldersig = bfolder+"/cms2r0/babymaker/babies/2016_08_10/TCHiHH/merged_higmc_higloose/";

  // Cuts in baseline speed up the yield finding
  map<string, string> baseline;
  baseline["lep0"] = "pass && stitch && weight<1 && nvleps==0";
  baseline["lep1"] = "pass && stitch && weight<1 && nleps==1 && mt<=140";
  baseline["lep2"] = "pass && stitch && weight<1 && nleps==2 && (mumu_m*(mumu_m>0)+elel_m*(elel_m>0))>80&&(mumu_m*(mumu_m>0)+elel_m*(elel_m>0))<100";

  Palette colors("txt/colors.txt", "default");

  map<string, set<string>> filetags; 
  filetags["ttx"]     = set<string>({"*_TTJets*Lept*.root", "*_TTJets_HT*.root", "*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*_ttHJetTobb*.root","*_TTTT*.root"});
  filetags["ttonly"]  = set<string>({"*_TTJets*Lept*.root", "*_TTJets_HT*.root"});
  filetags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
  filetags["singlet"] = set<string>({"*_ST_*.root"});
  filetags["qcd"]     = set<string>({"*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root"});
  filetags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  set<string> allfiletags;
  for (auto &iset: filetags) {
    if (proc_types.find(iset.first)!=proc_types.end()) 
      allfiletags.insert(iset.second.begin(), iset.second.end());
  }

  map<string, vector<shared_ptr<Process> > > procs;
  for (auto &iskim: skims){
    procs[iskim] = vector<shared_ptr<Process> >();
    if (proc_types.find("ttx")!=proc_types.end()) 
      procs[iskim].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X", Process::Type::background, colors("tt_1l"),
        attach_folder(foldermc[iskim], filetags["ttx"]), baseline[iskim]));
    if (proc_types.find("ttonly")!=proc_types.end()) 
      procs[iskim].push_back(Process::MakeShared<Baby_full>("t#bar{t}", Process::Type::background, colors("tt_1l"),
        attach_folder(foldermc[iskim], filetags["ttonly"]), baseline[iskim]));
    if (proc_types.find("vjets")!=proc_types.end())       
      procs[iskim].push_back(Process::MakeShared<Baby_full>("V+jets", Process::Type::background, kOrange+1,
        attach_folder(foldermc[iskim], filetags["vjets"]), baseline[iskim]));
    if (proc_types.find("singlet")!=proc_types.end())       
      procs[iskim].push_back(Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
        attach_folder(foldermc[iskim], filetags["singlet"]), baseline[iskim]));
    if (proc_types.find("qcd")!=proc_types.end())       
      procs[iskim].push_back(Process::MakeShared<Baby_full>("QCD", Process::Type::background, colors("other"),
        attach_folder(foldermc[iskim], filetags["qcd"]), baseline[iskim]));
    if (proc_types.find("other")!=proc_types.end())       
      procs[iskim].push_back(Process::MakeShared<Baby_full>("Other", Process::Type::background, kGreen+1,
        attach_folder(foldermc[iskim], filetags["other"]), baseline[iskim]));
    if (iskim == "lep0") {
      procs["sig"+iskim] = vector<shared_ptr<Process> >(procs[iskim]);
      for (unsigned isig(0); isig<sigm.size(); isig++)
        procs["sig"+iskim].push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
          sig_colors[isig], {foldersig+"*TChiHH_mGluino-"+sigm[isig]+"*.root"}, baseline[iskim]));
    }
  }
  if (do_cats_ntrub) {
    for (auto &iskim: skims){
      set<string> allfiles = attach_folder(foldermc[iskim], allfiletags);
      NamedFunc base_func(baseline[iskim]);
      procs["bcat_"+iskim] = vector<shared_ptr<Process> >();
      procs["bcat_"+iskim].push_back(Process::MakeShared<Baby_full>
              ("#leq 1 B-hadron", Process::Type::background, kPink+2,
               allfiles, base_func && ntrub<=1));
      procs["bcat_"+iskim].push_back(Process::MakeShared<Baby_full>
      			  ("2 B-hadrons", Process::Type::background, kOrange-4,
      			   allfiles, base_func && ntrub==2));
      procs["bcat_"+iskim].push_back(Process::MakeShared<Baby_full>
      			  ("3 B-hadrons", Process::Type::background, kTeal-8, 
               allfiles, base_func &&  ntrub==3));
      procs["bcat_"+iskim].push_back(Process::MakeShared<Baby_full>
      			  ("#geq 4 B-hadrons", Process::Type::background, kAzure-4, 
      			   allfiles, base_func && ntrub>=4));
    }
  }

  PlotMaker pm;

  //just to be pretty... already in skim...
  string njcut = "njets>=4 && njets<=5";

  map<string, string> selns;
  selns["base"]   = "pass && stitch && nvleps==0 && ntks==0 && !low_dphi";
  selns["notrkphi"] = "pass && stitch && nvleps==0";
  if (do_cr_plots || do_cr_pies) {
    selns["qcd"]   = "pass && stitch && nvleps==0 && ntks==0 && low_dphi";
    selns["lep1"] = baseline["lep1"];
    selns["lep2"] = baseline["lep2"];
  }

  map<string, vector<string> > metcuts;
  if (do_onemet) metcuts["base"] = {"met>150"};
  else metcuts["base"] = {"met>150&&met<=200", "met>200&&met<=300","met>300"};
  metcuts["qcd"] = metcuts["base"];
  metcuts["notrkphi"] = metcuts["base"];
  metcuts["lep1"] = metcuts["base"];
  metcuts["lep2"] = {"(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))>50"};
  
  vector<string> nbcuts;
  if (do_23vs4b) {
    nbcuts.push_back("((nbt==2&&nbm==2)||(nbt>=2&&nbm==3&&nbl==3))");
    nbcuts.push_back("nbt>=2&&nbm>=3&&nbl>=4");
  } else {
    nbcuts.push_back("nbt==2&&nbm==2");
    if (do_ge3b) {
      nbcuts.push_back("nbt>=2&&nbm>=3");
    } else {
      nbcuts.push_back("nbt>=2&&nbm==3&&nbl==3");
      nbcuts.push_back("nbt>=2&&nbm>=3&&nbl>=4");
    }
  }

  map<string, string> xcuts; // useful additional cut definitions
  xcuts["drmax"] = "hig_drmax<2.2";
  xcuts["hig"] = "hig_am>100 && hig_am<=140 && hig_dm <= 40";
  xcuts["fullhig"] = "hig_am>100 && hig_am<=140 && hig_dm <= 40 && hig_drmax<2.2";
  xcuts["sbd"] = "(hig_am<=100 || hig_am>140 || hig_dm > 40)";

  //     N-1 and other 1D distributions
  //----------------------------------------
  for (auto &iseln: selns) {
    //decide which is the relevant set of procs
    string iskimcat = "lep0";
    if (iseln.first=="lep1" || iseln.first=="lep2") iskimcat = iseln.first;
    string iproc = iskimcat=="lep0" ? ("sig"+iskimcat) : iskimcat; //want to include signal for selections with 0 leptons
    //for the moment skip control regions
    if (!do_cr_plots && (iseln.first=="lep1" || iseln.first=="lep2" || iseln.first=="qcd")) continue;
    for(auto &imet: metcuts[iseln.first]) { 
      for(auto &inb: nbcuts) {  
        if (iseln.first!= "notrkphi") {
          pm.Push<Hist1D>(Axis(15,0,150,"hig_dm", "#Deltam [GeV]", {40.}),
            iseln.second+"&&"+njcut+"&&"+imet+"&&"+inb+"&&hig_am>100 && hig_am<=140 && hig_drmax<2.2", 
            procs[iproc], all_plot_types).Tag(iseln.first);
          pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m> [GeV]", {100., 140.}),
            iseln.second+"&&"+njcut+"&&"+imet+"&&"+inb+"&&hig_dm<=40 && hig_drmax<2.2", 
            procs[iproc], all_plot_types).Tag(iseln.first);
          pm.Push<Hist1D>(Axis(20,0,4,"hig_drmax", "#DeltaR_{max}", {2.2}),
            iseln.second+"&&"+njcut+"&&"+imet+"&&"+inb+"&&hig_am>100 && hig_am<=140 && hig_dm <= 40", 
            procs[iproc], all_plot_types).Tag(iseln.first);
          string tmp_seln = iseln.second+"&&"+njcut+"&&"+imet+"&&"+inb+"&&" + xcuts["fullhig"];
          pm.Push<Hist1D>(Axis(15,0,600,"jets_pt[0]", "Jet 1 p_{T} [GeV]"), tmp_seln, procs[iproc], all_plot_types).Tag(iseln.first);
          pm.Push<Hist1D>(Axis(17,0,340,"jets_pt[1]", "Jet 2 p_{T} [GeV]"), tmp_seln, procs[iproc], all_plot_types).Tag(iseln.first);
          pm.Push<Hist1D>(Axis(12,0,240,"jets_pt[2]", "Jet 3 p_{T} [GeV]"), tmp_seln, procs[iproc], all_plot_types).Tag(iseln.first);
          pm.Push<Hist1D>(Axis(12,0,240,"jets_pt[3]", "Jet 4 p_{T} [GeV]"), tmp_seln, procs[iproc], all_plot_types).Tag(iseln.first);
          tmp_seln = iseln.second+"&&"+njcut+"&& met>100 &&"+inb+"&&" + xcuts["fullhig"];
          pm.Push<Hist1D>(Axis(10,100,600,"met", "E_{T}^{miss} [GeV]"), tmp_seln, procs[iproc], all_plot_types).Tag(iseln.first);
          tmp_seln = iseln.second+"&&"+njcut+"&&"+imet+"&& nbm>=2 &&" + xcuts["fullhig"];
          pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbl", "N_{b}^{L}"), tmp_seln, procs[iproc], all_plot_types).Tag(iseln.first);
          pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbm", "N_{b}^{M}"), tmp_seln, procs[iproc], all_plot_types).Tag(iseln.first);
          pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbt", "N_{b}^{T}"), tmp_seln, procs[iproc], all_plot_types).Tag(iseln.first);
        } else if (iseln.first == "notrkphi") {
          string tmp_seln = iseln.second+"&& !low_dphi &&"+njcut+"&&"+imet+"&&"+inb;
          pm.Push<Hist1D>(Axis(5,-0.5,4.5,"ntks", "N_{tks}"),
            tmp_seln, procs[iproc], all_plot_types).Tag(iseln.first);
          tmp_seln = iseln.second+"&& ntks==0 &&"+njcut+"&&"+imet+"&&"+inb;
          pm.Push<Hist1D>(Axis(16,0.,3.2,"dphi2", "#Delta#phi_{2}"),
            tmp_seln+"&& dphi1>0.5", procs[iproc], all_plot_types).Tag(iseln.first);
          pm.Push<Hist1D>(Axis(16,0.,3.2,"dphi3", "#Delta#phi_{3}"),
            tmp_seln+"&& dphi1>0.5 && dphi2>0.5", procs[iproc], all_plot_types).Tag(iseln.first);
          pm.Push<Hist1D>(Axis(16,0.,3.2,"dphi4", "#Delta#phi_{4}"),
            tmp_seln+"&& dphi1>0.5 && dphi2>0.5 && dphi3>0.3", procs[iproc], all_plot_types).Tag(iseln.first);
        }
      }
    }
  }

  //        Pie charts
  //--------------------------------
  for (auto &iseln: selns) {
    //decide which is the relevant set of procs
    string iproc = "lep0";
    if (iseln.first=="lep1" || iseln.first=="lep2") iproc = iseln.first;
    if (!do_cr_pies && (iseln.first=="lep1" || iseln.first=="lep2" || iseln.first=="qcd")) continue;
    for (auto &ixcut: {"drmax"}) {
      vector<TString> cuts;
      vector<TableRow> table_cuts;
      for(auto &imet: metcuts[iseln.first]) {
        for(auto &inb: nbcuts) {
          for (auto &ireg: {"hig","sbd"}) {
            cuts.push_back(iseln.second+"&&"+njcut+"&&"+imet+"&&"+inb+"&&"+xcuts[ireg]+"&&"+xcuts[ixcut]);
          }
        }
      }
      for(size_t icut=0; icut<cuts.size(); icut++)
        table_cuts.push_back(TableRow("$"+CodeToLatex(cuts[icut].Data())+"$", cuts[icut].Data()));  
      pm.Push<Table>("chart_"+iseln.first+"_"+ixcut,  table_cuts, procs[iproc], true, true, true, false);
      if (do_cats_ntrub) 
        pm.Push<Table>("chartcats_"+iseln.first+"_"+ixcut,  table_cuts, procs["cats"+iproc], true, true, true, false);
    }
  }

  //        Cutflow table
  //--------------------------------
  pm.Push<Table>("cutflow", vector<TableRow>{
    TableRow("$E_{T}^{miss} > 150$, $\\text{2M b-tags}$, $\\text{4 or 5 jets}$, $0\\ell$", 
      "met>150 && nbm>=2 && njets>=4 && njets<=5"),
    TableRow("$E_{T}^{miss} > 150$, $\\text{2T b-tags}$, $\\text{4 or 5 jets}$, $0\\ell$", 
      "met>150 && nbt>=2 && njets>=4 && njets<=5"),
    TableRow("$\\Delta\\phi_{\\text{min}}$", 
      "met>150 && nbt>=2 && njets>=4 && njets<=5 && !low_dphi"),
    TableRow("$\\Delta R_{\\text{max}} < 2.2$", 
      "met>150 && nbt>=2 && njets>=4 && njets<=5 && !low_dphi && hig_drmax<2.2" ,0,1),
    TableRow("$\\Delta m < 40$", 
      "met>150 && nbt>=2 && njets>=4 && njets<=5 && !low_dphi && hig_drmax<2.2 && hig_dm<=40"),
    TableRow("$\\left< m \\right> \\in (100,140)$", 
      "met>150 && nbt>=2 && njets>=4 && njets<=5 && !low_dphi && hig_drmax<2.2 && hig_dm<=40 && hig_am>100 && hig_am<=140",0,1),
    TableRow("3b", 
      "met>150 && nbt>=2&&nbm==3&&nbl==3 && njets>=4 && njets<=5 && !low_dphi && hig_drmax<2.2 && hig_dm<=40 && hig_am>100 && hig_am<=140"),
    TableRow("4b", 
      "met>150 && nbt>=2&&nbm>=3&&nbl>=4 && njets>=4 && njets<=5 && !low_dphi && hig_drmax<2.2 && hig_dm<=40 && hig_am>100 && hig_am<=140"),
    TableRow("3b, $E_{T}^{miss}>$200", 
      "met>200 && nbt>=2&&nbm==3&&nbl==3 && njets>=4 && njets<=5 && !low_dphi && hig_drmax<2.2 && hig_dm<=40 && hig_am>100 && hig_am<=140"),
    TableRow("4b, $E_{T}^{miss}>$200", 
      "met>200 && nbt>=2&&nbm>=3&&nbl>=4 && njets>=4 && njets<=5 && !low_dphi && hig_drmax<2.2 && hig_dm<=40 && hig_am>100 && hig_am<=140"),
    TableRow("3b, $E_{T}^{miss}>$300", 
      "met>300 && nbt>=2&&nbm==3&&nbl==3 && njets>=4 && njets<=5 && !low_dphi && hig_drmax<2.2 && hig_dm<=40 && hig_am>100 && hig_am<=140"),
    TableRow("4b, $E_{T}^{miss}>$300", 
      "met>300 && nbt>=2&&nbm>=3&&nbl>=4 && njets>=4 && njets<=5 && !low_dphi && hig_drmax<2.2 && hig_dm<=40 && hig_am>100 && hig_am<=140"),
  },procs["siglep0"],0);

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
