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

using namespace std;
using namespace PlotOptTypes;

namespace{
  float lumi = 36.2;
  bool do_data = false;
  bool do_cats_ntrub = false;
  bool do_loose_seln = false; // removes track veto and delta phi requirement
  //output
  bool do_stackplots = false;
  bool do_cflow = false;
  bool do_pies = true;
  bool do_shapes = false;
  bool do_cr_pies = false;
  bool do_cr_stackplots = false;
  // simplify selection
  bool do_metg150_metg200 = false; // only integrated met>150 and met > 200
  bool do_ge3b = false;  // do btag cats: 2b and 3b+
  bool do_23vs4b = false; // combine 2b with 3b, keep 4b separate
  // choose processes to include, options are: "ttx", "vjets", "singlet", "qcd", "other", "ttonly"
  set<string> proc_types = {"ttx", "vjets", "singlet", "qcd", "other"}; // for default data/MC
  // set<string> proc_types = {}; // to make signal plots only
  // set<string> proc_types = {"ttonly"};
  // signal points to include and their colors
  vector<string> sigm = {"225","350","700"}; 
  vector<int> sig_colors = {kGreen, kRed, kBlue}; // need sigm.size() >= sig_colors.size()
  //for signal plots only
  // vector<string> sigm = {"175","225","350","700","1000"}; 
  // vector<int> sig_colors = {kMagenta+2 , kGreen, kRed, kBlue, kAzure+10}; // need sigm.size() >= sig_colors.size()
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
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::ratio);
  vector<PlotOpt> all_plot_types = {lin_lumi_info, log_lumi_info};
  vector<PlotOpt> all_shapes = {lin_shapes_info, log_shapes_info};

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  set<string> skims;
  if (do_cr_pies || do_cr_stackplots) skims = {"lep0","lep1","lep2"};
  else skims = {"lep0"};

  map<string, string> foldermc; //ordered in number of leptons
  foldermc["lep0"] = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/";
  foldermc["lep1"] = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_hig_nlep1/";
  foldermc["lep2"] = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_hig_nlep2/";

  string foldersig = bfolder+"/cms2r0/babymaker/babies/2016_08_10/TCHiHH/merged_higmc_unskimmed/";

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

  set<string> allfiletags_noqcd;
  for (auto &iset: filetags) {
    if (iset.first == "qcd") continue;
    if (proc_types.find(iset.first)!=proc_types.end()) 
      allfiletags_noqcd.insert(iset.second.begin(), iset.second.end());
  }

  map<string, vector<shared_ptr<Process> > > procs;
  if (do_stackplots || do_pies) {
    for (auto &iskim: skims){
      procs[iskim] = vector<shared_ptr<Process> >();
      if (proc_types.find("ttx")!=proc_types.end()) 
        procs[iskim].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X", Process::Type::background, colors("tt_1l"),
          attach_folder(foldermc[iskim], filetags["ttx"]), "pass && stitch"));
      if (proc_types.find("ttonly")!=proc_types.end()) 
        procs[iskim].push_back(Process::MakeShared<Baby_full>("t#bar{t}", Process::Type::background, colors("tt_1l"),
          attach_folder(foldermc[iskim], filetags["ttonly"]), "pass && stitch"));
      if (proc_types.find("vjets")!=proc_types.end())       
        procs[iskim].push_back(Process::MakeShared<Baby_full>("V+jets", Process::Type::background, kOrange+1,
          attach_folder(foldermc[iskim], filetags["vjets"]), "pass && stitch"));
      if (proc_types.find("singlet")!=proc_types.end())       
        procs[iskim].push_back(Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
          attach_folder(foldermc[iskim], filetags["singlet"]), "pass && stitch"));
      if (proc_types.find("qcd")!=proc_types.end())       
        procs[iskim].push_back(Process::MakeShared<Baby_full>("QCD", Process::Type::background, colors("other"),
          attach_folder(foldermc[iskim], filetags["qcd"]), "pass && stitch")); // && weight<1
      if (proc_types.find("other")!=proc_types.end())       
        procs[iskim].push_back(Process::MakeShared<Baby_full>("Other", Process::Type::background, kGreen+1,
          attach_folder(foldermc[iskim], filetags["other"]), "pass && stitch"));
      if (iskim == "lep0") {
        procs["sig"+iskim] = vector<shared_ptr<Process> >(procs[iskim]);
        if (proc_types.size()==0) { // have to pretend signal is background, otherwise crashes
          for (unsigned isig(0); isig<sigm.size(); isig++)
            procs["sig"+iskim].push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)", Process::Type::background, 
              sig_colors[isig], {foldersig+"*TChiHH_mGluino-"+sigm[isig]+"*.root"}, "1"));
        } else {
          for (unsigned isig(0); isig<sigm.size(); isig++)
            procs["sig"+iskim].push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
              sig_colors[isig], {foldersig+"*TChiHH_mGluino-"+sigm[isig]+"*.root"}, "1"));
        }
      }
    }
  }

  if (do_cats_ntrub && proc_types.size()>0) {
    for (auto &iskim: skims){
      set<string> allfiles = attach_folder(foldermc[iskim], allfiletags);
      NamedFunc base_func("pass && stitch");
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

  if (do_shapes){
    set<string> allfiles = attach_folder(foldermc["lep0"], allfiletags); 
    procs["shapes"] = vector<shared_ptr<Process> >();
    procs["shapes"].push_back(Process::MakeShared<Baby_full>("2b all bkg", Process::Type::background, kBlack,
      allfiles, "pass && stitch && weight<0.2" && Functions::hig_nb==2));
    procs["shapes"].push_back(Process::MakeShared<Baby_full>("3b all bkg", Process::Type::background, kAzure+1,
      allfiles, "pass && stitch && weight<0.2" && Functions::hig_nb==3));
    procs["shapes"].push_back(Process::MakeShared<Baby_full>("4b all bkg", Process::Type::background, kPink+2,
      allfiles, "pass && stitch && weight<0.2" && Functions::hig_nb==4));
  }

  PlotMaker pm;

  //just to be pretty... already in skim...
  string njcut = "njets>=4 && njets<=5";

  map<string, string> selns;
  selns["base"]   = "nvleps==0 && ntks==0 && !low_dphi";
  if (do_loose_seln) selns["notrkphi"] = "nvleps==0";
  if (do_cr_stackplots || do_cr_pies) {
    selns["qcd"]   = "nvleps==0 && ntks==0 && low_dphi";
    selns["lep1"] = "nleps==1 && mt<=100";
    selns["lep2"] = "nleps==2 && (mumu_m*(mumu_m>0)+elel_m*(elel_m>0))>80&&(mumu_m*(mumu_m>0)+elel_m*(elel_m>0))<100";
  }

  map<string, vector<string> > metcuts;
  if (do_metg150_metg200) metcuts["base"] = {"met>150","met>200"};
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
  xcuts["drmax"] = "hig_drmax<=2.2";
  xcuts["hig"] = "hig_am>100 && hig_am<=140 && hig_dm <= 40";
  xcuts["fullhig"] = "hig_am>100 && hig_am<=140 && hig_dm <= 40 && hig_drmax<=2.2";
  xcuts["sbd"] = "(hig_am<=100 || (hig_am>140 && hig_am<=200)) && hig_dm <= 40";
  xcuts["fullsbd"] = "(hig_am<=100 || (hig_am>140 && hig_am<=200)) && hig_dm <= 40 && hig_drmax<=2.2";

  //    Signal plots with no (min) selection
  //----------------------------------------------
  if (sigm.size()>4) {
    pm.Push<Hist1D>(Axis(30,0,150,"hig_dm", "#Deltam [GeV]", {40.}), "njets>=4", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m> [GeV]", {100., 140.}), "njets>=4", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m> [GeV]", {100., 140.}), "njets>=4  && hig_dm<=40", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(30,0,150,"hig_dm", "#Deltam [GeV]", {40.}), "njets>=4 && nbt>=2&&nbm>=3&&nbl>=4", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m> [GeV]", {100., 140.}), "njets>=4 && nbt>=2&&nbm>=3&&nbl>=4", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(40,0,4,"hig_drmax", "#DeltaR_{max}", {2.2}), "njets>=4", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(20,0,4,"hig_drmax", "#DeltaR_{max}", {2.2}), "njets>=4 && met>200", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(30,0,600,"jets_pt[0]", "Jet 1 p_{T} [GeV]"), "njets>=1", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(20,0,400,"jets_pt[1]", "Jet 2 p_{T} [GeV]"), "njets>=2", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(24,0,240,"jets_pt[2]", "Jet 3 p_{T} [GeV]"), "njets>=3", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(24,0,240,"jets_pt[3]", "Jet 4 p_{T} [GeV]"), "njets>=4", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(30,0,600,"jets_pt[0]", "Jet 1 p_{T} [GeV]"), "njets>=4", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(20,0,400,"jets_pt[1]", "Jet 2 p_{T} [GeV]"), "njets>=4", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(24,0,240,"jets_pt[2]", "Jet 3 p_{T} [GeV]"), "njets>=4", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(24,0,240,"jets_pt[3]", "Jet 4 p_{T} [GeV]"), "njets>=4", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbl", "N_{b}^{L}"), "1", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbm", "N_{b}^{M}"), "1", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbt", "N_{b}^{T}"), "1", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbl", "N_{b}^{L}"), "njets>=4", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbm", "N_{b}^{M}"), "njets>=4", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbt", "N_{b}^{T}"), "njets>=4", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(8,-0.5,7.5,"njets", "N_{jets}"), "1", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(26,0,1300,"met", "E_{T}^{miss} [GeV]"), "1", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(26,0,1300,"met", "E_{T}^{miss} [GeV]"), "njets>=4", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(26,0,1300,"ht", "H_{T} [GeV]"), "1", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
    pm.Push<Hist1D>(Axis(26,0,1300,"ht", "H_{T} [GeV]"), "njets>=4", procs["siglep0"], 
      vector<PlotOpt>({lin_shapes_info, log_shapes_info})).Tag("sig");
  }

  //     N-1 and other 1D distributions
  //----------------------------------------
  if (do_stackplots) {
    for (auto &iseln: selns) {
      NamedFunc wgt = "weight" * Functions::eff_mettrig;
      if (!do_data) wgt = "weight";
      if (proc_types.size()==0) continue;
      //decide which is the relevant set of procs
      string iskimcat = "lep0";
      if (iseln.first=="lep1" || iseln.first=="lep2") iskimcat = iseln.first;
      string iproc = iskimcat=="lep0" ? ("sig"+iskimcat) : iskimcat; //want to include signal for selections with 0 leptons
      //for the moment skip control regions
      if (!do_cr_stackplots && (iseln.first=="lep1" || iseln.first=="lep2" || iseln.first=="qcd")) continue;
      for(unsigned imet(0); imet<metcuts[iseln.first].size(); imet++) { 
        for(auto &inb: nbcuts) {  
          if (iseln.first!= "notrkphi") {
            pm.Push<Hist1D>(Axis(15,0,150,"hig_dm", "#Deltam [GeV]", {40.}),
              iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&&"+inb+"&&hig_am>100 && hig_am<=140 && hig_drmax<=2.2", 
              procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m> [GeV]", {100., 140.}),
              iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&&"+inb+"&&hig_dm<=40 && hig_drmax<=2.2", 
              procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(20,0,4,"hig_drmax", "#DeltaR_{max}", {2.2}),
              iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&&"+inb+"&&hig_am>100 && hig_am<=140 && hig_dm <= 40", 
              procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            string tmp_seln = iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&&"+inb+"&&" + xcuts["fullhig"];
            pm.Push<Hist1D>(Axis(15,0,600,"jets_pt[0]", "Jet 1 p_{T} [GeV]"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(17,0,340,"jets_pt[1]", "Jet 2 p_{T} [GeV]"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(12,0,240,"jets_pt[2]", "Jet 3 p_{T} [GeV]"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(12,0,240,"jets_pt[3]", "Jet 4 p_{T} [GeV]"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            tmp_seln = iseln.second+"&&"+njcut+"&& met>100 &&"+inb+"&&" + xcuts["fullhig"];
            pm.Push<Hist1D>(Axis(10,100,600,"met", "E_{T}^{miss} [GeV]",{150,200,300}), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            tmp_seln = iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&& nbm>=2 &&" + xcuts["fullhig"];
            pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbl", "N_{b}^{L}"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbm", "N_{b}^{M}"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbt", "N_{b}^{T}"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            tmp_seln = iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&& nbm>=2 && hig_dm<=40";
            pm.Push<Hist1D>(Axis(6,-0.5,5.5,Functions::hig_nb, "b-tag category"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            tmp_seln = iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&& nbm>=2 && hig_dm<=40 && hig_drmax<=2.2";
            pm.Push<Hist1D>(Axis(6,-0.5,5.5,Functions::hig_nb, "b-tag category"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
          } else if (iseln.first == "notrkphi") {
            if (imet>0) continue; 
            string tmp_seln = iseln.second+"&& ntks==0 && !low_dphi &&"+njcut+"&& met>100 &&"+inb;
            pm.Push<Hist1D>(Axis(10,100,600,"met", "E_{T}^{miss} [GeV]",{150,200,300}), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            tmp_seln = iseln.second+"&& !low_dphi &&"+njcut+"&& met>150 &&"+inb;
            pm.Push<Hist1D>(Axis(5,-0.5,4.5,"ntks", "N_{tks}"),
              tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            tmp_seln = iseln.second+"&& ntks==0 &&"+njcut+"&& met>150 &&"+inb;
            pm.Push<Hist1D>(Axis(32,0.,3.2,"dphi2", "#Delta#phi_{2}",{0.5}),
              tmp_seln+"&& dphi1>0.5", procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(32,0.,3.2,"dphi3", "#Delta#phi_{3}",{0.3}),
              tmp_seln+"&& dphi1>0.5 && dphi2>0.5", procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(32,0.,3.2,"dphi4", "#Delta#phi_{4}",{0.3}),
              tmp_seln+"&& dphi1>0.5 && dphi2>0.5 && dphi3>0.3", procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
          }
        }
      }
    }
  }

  //        Pie charts
  //--------------------------------
  if (do_pies) { 
    NamedFunc wgt = "weight" * Functions::eff_mettrig;
    for (auto &iseln: selns) {
      if (proc_types.size()==0) continue;
      //decide which is the relevant set of procs
      string iproc = "lep0";
      if (iseln.first=="lep1" || iseln.first=="lep2") iproc = iseln.first;
      if (!do_cr_pies && (iseln.first=="lep1" || iseln.first=="lep2" || iseln.first=="qcd")) continue;
      vector<TString> cuts;
      vector<TableRow> table_cuts;
      for(auto &imet: metcuts[iseln.first]) {
        for(auto &inb: nbcuts) {
          for (auto &ireg: {"fullhig","fullsbd"}) {
            cuts.push_back(iseln.second+"&&"+njcut+"&&"+imet+"&&"+inb+"&&"+xcuts[ireg]);
          }
        }
      }
      for(size_t icut=0; icut<cuts.size(); icut++)
        table_cuts.push_back(TableRow("$"+CodeToLatex(cuts[icut].Data())+"$", cuts[icut].Data(), 0, 0, wgt));  
      pm.Push<Table>("chart_"+iseln.first,  table_cuts, procs[iproc], true, true, true, false);
      if (do_cats_ntrub) 
        pm.Push<Table>("chartcats_"+iseln.first,  table_cuts, procs["bcat_"+iproc], true, true, true, false);
    }
  }

  //        Cutflow table
  //-------------------------------- 
  if (do_cflow && proc_types.size()>0) {
    NamedFunc wgt = "weight" * Functions::eff_mettrig;
    pm.Push<Table>("cutflow", vector<TableRow>{
      TableRow("$E_{T}^{miss} > 150$, 2 CSVM, $\\text{4-5 jets}$, $0\\ell/$tk", 
        "nvleps==0 && ntks==0 && met>150 && nbm>=2 &&"+njcut,0,0, wgt),
      TableRow("2 CSVT", 
        "nvleps==0 && ntks==0 && met>150 && nbt>=2 &&"+njcut,0,0, wgt),
      TableRow("$\\Delta\\phi_{1,2}>0.5,\\Delta\\phi_{3,4}>0.3$",        
        "nvleps==0 && ntks==0 && met>150 && nbt>=2 &&"+njcut+"&& !low_dphi",0,0, wgt),
      TableRow("$\\Delta m < 40$",     
        "nvleps==0 && ntks==0 && met>150 && nbt>=2 &&"+njcut+"&& !low_dphi && hig_dm<=40",0,0, wgt),
      TableRow("$\\Delta R_{\\text{max}} < 2.2$",                    
        "nvleps==0 && ntks==0 && met>150 && nbt>=2 &&"+njcut+"&& !low_dphi && hig_drmax<=2.2 && hig_dm<=40",0,0, wgt),
      TableRow("$\\left< m \\right> \\in (100,140)$", 
        "nvleps==0 && ntks==0 && met>150 && nbt>=2 &&"+njcut+"&& !low_dphi &&"+xcuts["fullhig"],0,0, wgt),
      TableRow("HIG, 3b, $150<E_{T}^{miss}\\leq$200", 
        "nvleps==0 && ntks==0 && met>150 && met<=200 && nbt>=2&&nbm==3&&nbl==3 &&"+njcut+"&& !low_dphi &&"+xcuts["fullhig"],1,0, wgt),
      TableRow("SBD, 3b, $150<E_{T}^{miss}\\leq$200", 
        "nvleps==0 && ntks==0 && met>150 && met<=200 && nbt>=2&&nbm==3&&nbl==3 &&"+njcut+"&& !low_dphi &&"+xcuts["fullsbd"],0,0, wgt),
      TableRow("HIG, 4b, $150<E_{T}^{miss}\\leq$200", 
        "nvleps==0 && ntks==0 && met>150 && met<=200 && nbt>=2&&nbm>=3&&nbl>=4 &&"+njcut+"&& !low_dphi &&"+xcuts["fullhig"],0,0, wgt),
      TableRow("SBD, 4b, $150<E_{T}^{miss}\\leq$200", 
        "nvleps==0 && ntks==0 && met>150 && met<=200 && nbt>=2&&nbm>=3&&nbl>=4 &&"+njcut+"&& !low_dphi &&"+xcuts["fullsbd"],0,0, wgt),
      TableRow("HIG, 3b, $200<E_{T}^{miss}\\leq$300", 
        "nvleps==0 && ntks==0 && met>200 && met<=300 && nbt>=2&&nbm==3&&nbl==3 &&"+njcut+"&& !low_dphi &&"+xcuts["fullhig"],0,0, wgt),
      TableRow("SBD, 3b, $200<E_{T}^{miss}\\leq$300", 
        "nvleps==0 && ntks==0 && met>200 && met<=300 && nbt>=2&&nbm==3&&nbl==3 &&"+njcut+"&& !low_dphi &&"+xcuts["fullsbd"],0,0, wgt),
      TableRow("HIG, 4b, $200<E_{T}^{miss}\\leq$300", 
        "nvleps==0 && ntks==0 && met>200 && met<=300 && nbt>=2&&nbm>=3&&nbl>=4 &&"+njcut+"&& !low_dphi &&"+xcuts["fullhig"],0,0, wgt),
      TableRow("SBD, 4b, $200<E_{T}^{miss}\\leq$300", 
        "nvleps==0 && ntks==0 && met>200 && met<=300 && nbt>=2&&nbm>=3&&nbl>=4 &&"+njcut+"&& !low_dphi &&"+xcuts["fullsbd"],0,0, wgt),
      TableRow("HIG, 3b, $E_{T}^{miss}>$300",         
        "nvleps==0 && ntks==0 && met>300 && nbt>=2&&nbm==3&&nbl==3 &&"+njcut+"&& !low_dphi &&"+xcuts["fullhig"],0,0, wgt),
      TableRow("SBD, 3b, $E_{T}^{miss}>$300",         
        "nvleps==0 && ntks==0 && met>300 && nbt>=2&&nbm==3&&nbl==3 &&"+njcut+"&& !low_dphi &&"+xcuts["fullsbd"],0,0, wgt),
      TableRow("HIG, 4b, $E_{T}^{miss}>$300",         
        "nvleps==0 && ntks==0 && met>300 && nbt>=2&&nbm>=3&&nbl>=4 &&"+njcut+"&& !low_dphi &&"+xcuts["fullhig"],0,0, wgt),
      TableRow("SBD, 4b, $E_{T}^{miss}>$300",         
        "nvleps==0 && ntks==0 && met>300 && nbt>=2&&nbm>=3&&nbl>=4 &&"+njcut+"&& !low_dphi &&"+xcuts["fullsbd"],0,0, wgt),
    },procs["siglep0"],0);
  }

  //    Background closure - shapes
  //-----------------------------------
  if (do_shapes) {
    string tmp_seln = "nvleps==0 && ntks==0 && !low_dphi";
    if (do_loose_seln) tmp_seln = "nvleps==0";
    for(auto &imet: metcuts["base"]) { 
      pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m> [GeV]", {100., 140.}),
              tmp_seln+"&&"+njcut+"&&"+imet, 
              procs["shapes"], all_shapes).Tag("base");
      pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m> [GeV]", {100., 140.}),
              tmp_seln+"&&"+njcut+"&&"+imet+"&&hig_dm<=40", 
              procs["shapes"], all_shapes).Tag("base");
      pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m> [GeV]", {100., 140.}),
              tmp_seln+"&&"+njcut+"&&"+imet+"&&hig_dm<=40 && hig_drmax<=2.2", 
              procs["shapes"], all_shapes).Tag("base");      
    }
  }

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
