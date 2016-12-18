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
  //fixme:simplify options
  float lumi = 36.2;
  // float lumi = 4.3;
  string json = "json4p0";
  bool do_data = false;
  bool do_data_shapes = false;
  bool unblind = false;
  bool do_cats_ntrub = false;
  bool do_loose_seln = false; // removes track veto and delta phi requirement
  //output
  bool do_stackplots = true;
  bool do_cflow = false;
  bool do_shapes = false;
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
  if (do_data) {
    log_lumi_info = log_lumi().Title(TitleType::info).Bottom(BottomType::ratio);
    lin_lumi_info = lin_lumi().Title(TitleType::info).Bottom(BottomType::ratio);
  }
  vector<PlotOpt> all_plot_types = {lin_lumi_info, log_lumi_info};
  vector<PlotOpt> all_shapes = {lin_shapes_info};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  set<string> skims;
  if (do_cr_stackplots || do_data_shapes) skims = {"lep0","lep1","lep2"}; 
  else skims = {"lep0"};

  map<string, string> foldermc; //ordered in number of leptons
  foldermc["lep0"] = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/";
  foldermc["lep1"] = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higlep1/";
  foldermc["lep2"] = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_hig_nlep2/";

  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({"*_TTJets*Lept*.root", "*_TTJets_HT*.root", "*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*_ttHJetTobb*.root","*_TTTT*.root"});
  mctags["ttonly"]  = set<string>({"*_TTJets_*Lept*.root", "*_TTJets_HT*.root"});
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

  map<string, string> folderdata; //ordered in number of leptons
  folderdata["lep0"] = bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_higloose/";
  folderdata["lep1"] = bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_higlep1/";
  folderdata["lep2"] = bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_higlep2/";

  string foldersig = bfolder+"/cms2r0/babymaker/babies/2016_08_10/TCHiHH/merged_higmc_unskimmed/";

  map<string, vector<shared_ptr<Process> > > procs;
  if (do_stackplots || do_cr_stackplots) {
    for (auto &iskim: skims){
      procs[iskim] = vector<shared_ptr<Process> >();
      if (proc_types.find("ttx")!=proc_types.end()) 
        procs[iskim].push_back(Process::MakeShared<Baby_full>("t#bar{t}+X", Process::Type::background, colors("tt_1l"),
          attach_folder(foldermc[iskim], mctags["ttx"]), "pass && stitch"));
      if (proc_types.find("ttonly")!=proc_types.end()) 
        procs[iskim].push_back(Process::MakeShared<Baby_full>("t#bar{t}", Process::Type::background, colors("tt_1l"),
          attach_folder(foldermc[iskim], mctags["ttonly"]), "pass && stitch"));
      if (proc_types.find("vjets")!=proc_types.end())       
        procs[iskim].push_back(Process::MakeShared<Baby_full>("V+jets", Process::Type::background, kOrange+1,
          attach_folder(foldermc[iskim], mctags["vjets"]), "pass && stitch"));
      if (proc_types.find("singlet")!=proc_types.end())       
        procs[iskim].push_back(Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
          attach_folder(foldermc[iskim], mctags["singlet"]), "pass && stitch"));
      if (proc_types.find("qcd")!=proc_types.end())       
        procs[iskim].push_back(Process::MakeShared<Baby_full>("QCD", Process::Type::background, colors("other"),
          attach_folder(foldermc[iskim], mctags["qcd"]), "pass && stitch")); 
      if (proc_types.find("other")!=proc_types.end())       
        procs[iskim].push_back(Process::MakeShared<Baby_full>("Other", Process::Type::background, kGreen+1,
          attach_folder(foldermc[iskim], mctags["other"]), "pass && stitch"));
      if (do_data) {
        procs[iskim].push_back(Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
          {folderdata[iskim]+"*RunB*root"}, "trig_ra4 && pass &&"+json)); 
          // {folderdata[iskim]+"*RunB*root"}, "(trig[13]||trig[33]||trig[14]||trig[15]||trig[30]||trig[31]) && pass &&"+json)); 
      }
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
      set<string> allfiles = attach_folder(foldermc[iskim], allmctags);
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
    for (auto &iskim: skims){
      set<string> allfiles = attach_folder(foldermc[iskim], allmctags); 
      procs["shapes"+iskim] = vector<shared_ptr<Process> >();
      procs["shapes"+iskim].push_back(Process::MakeShared<Baby_full>("2b all bkg", Process::Type::background, kBlack,
        allfiles, "pass && stitch && weight<0.2" && Functions::hig_nb==2));
      procs["shapes"+iskim].push_back(Process::MakeShared<Baby_full>("3b all bkg", Process::Type::background, kAzure+1,
        allfiles, "pass && stitch && weight<0.2" && Functions::hig_nb==3));
      procs["shapes"+iskim].push_back(Process::MakeShared<Baby_full>("4b all bkg", Process::Type::background, kPink+2,
        allfiles, "pass && stitch && weight<0.2" && Functions::hig_nb==4));
    }
  }

  if (do_data_shapes){
    for (auto &iskim: skims){
      procs["data_shapes"+iskim] = vector<shared_ptr<Process> >();
      procs["data_shapes"+iskim].push_back(Process::MakeShared<Baby_full>("#geq 3b data", Process::Type::data, kBlack,
        {folderdata[iskim]+"*RunB*root"}, "trig_ra4 && pass" && json && Functions::hig_nb>=3));
      procs["data_shapes"+iskim].push_back(Process::MakeShared<Baby_full>("2b data", Process::Type::background, kBlack,
        {folderdata[iskim]+"*RunB*root"}, "trig_ra4 && pass" && json && Functions::hig_nb==2));
      procs["data_shapes"+iskim].back()->SetFillColor(kBlue-7);
      procs["data_shapes"+iskim].back()->SetLineColor(kBlue-7);
      procs["data_shapes"+iskim].back()->SetLineWidth(2);
    }
  }

  PlotMaker pm;

  //just to be pretty... already in skim...
  string njcut = "njets>=4 && njets<=5";

  map<string, string> selns;
  selns["base"]   = "nvleps==0 && ntks==0 && !low_dphi";
  if (do_loose_seln) selns["notrkphi"] = "nvleps==0";
  if (do_cr_stackplots || do_data_shapes || do_shapes) { 
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
  if (do_data) {
    nbcuts.push_back("nbt==2&&nbm==2");
    nbcuts.push_back("nbt>=2&&nbm>=3");
  } else {
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

  // Temporary funcs for other b-tag category options
  //-----------------------------------------------------
  const NamedFunc hig_nb_mmmm("hig_nb_mmmm",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.nbm()>=2) return min(4,b.nbm());
    else return 0;
  });
  const NamedFunc hig_nb_ttll("hig_nb_ttll",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.nbt()==2 && b.nbl()==2) return 2;
    else if (b.nbt()>=2 && b.nbl()==3) return 3;
    else if (b.nbt()>=2 && b.nbl()>=4) return 4;
    else return 0;
  });
  const NamedFunc hig_nb_tmml("hig_nb_tmml",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.nbt()>=1 && b.nbm()==2) return 2;
    else if (b.nbt()>=1 && b.nbm()==3 && b.nbl()==3) return 3;
    else if (b.nbt()>=1 && b.nbm()>=3 && b.nbl()>=4) return 4;
    else return 0;
  });
  vector<NamedFunc> btag_xopts = {hig_nb_mmmm, hig_nb_tmml, hig_nb_ttll};

  //     N-1 and other 1D distributions
  //----------------------------------------
  if ((do_stackplots || do_cr_stackplots) && !do_data) {
    for (auto &iseln: selns) {
      NamedFunc wgt = "weight" * Functions::eff_higtrig;
      if (!do_data) wgt = "weight";
      //decide which is the relevant set of procs
      string iskimcat = "lep0";
      if (iseln.first=="lep1" || iseln.first=="lep2") iskimcat = iseln.first;
      string iproc = iskimcat=="lep0" ? ("sig"+iskimcat) : iskimcat; //want to include signal for selections with 0 leptons
      //for the moment skip control regions
      if (!do_cr_stackplots && (iseln.first=="lep1" || iseln.first=="lep2" || iseln.first=="qcd")) continue;
      for(unsigned imet(0); imet<metcuts[iseln.first].size(); imet++) { 
        if (iseln.first!= "notrkphi") {
          string tmp_seln = iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&& nbm>=2 &&" + xcuts["fullhig"];
          pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbl", "N_{b}^{L}"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
          pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbm", "N_{b}^{M}"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
          pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbt", "N_{b}^{T}"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
          if (iseln.first=="base"){
            tmp_seln = iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&& nbm>=2 && hig_dm<=40 && hig_drmax<=2.2 && hig_am<=200";
            pm.Push<Hist1D>(Axis(5,0.5,5.5,Functions::hig_nb, "b-tag category (TTML)"), tmp_seln && Functions::hig_nb>0., procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(5,0.5,5.5,hig_nb_ttll, "b-tag category (TTLL)"), tmp_seln && hig_nb_ttll>0., procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(5,0.5,5.5,hig_nb_tmml, "b-tag category (TMML)"), tmp_seln && hig_nb_tmml>0., procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(5,0.5,5.5,hig_nb_mmmm, "b-tag category (MMMM)"), tmp_seln && hig_nb_mmmm>0., procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
          }
        }
        for(unsigned inb(0); inb<nbcuts.size(); inb++) {
          if (iseln.first!= "notrkphi") {
            pm.Push<Hist1D>(Axis(15,0,150,"hig_dm", "#Deltam [GeV]", {40.}),
              iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&&"+nbcuts[inb]+"&&hig_am>100 && hig_am<=140 && hig_drmax<=2.2", 
              procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m> [GeV]", {100., 140.}),
              iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&&"+nbcuts[inb]+"&&hig_dm<=40 && hig_drmax<=2.2", 
              procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(20,0,4,"hig_drmax", "#DeltaR_{max}", {2.2}),
              iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&&"+nbcuts[inb]+"&&hig_am>100 && hig_am<=140 && hig_dm <= 40", 
              procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            string tmp_seln = iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&&"+nbcuts[inb]+"&&" + xcuts["fullhig"];
            pm.Push<Hist1D>(Axis(15,0,600,"jets_pt[0]", "Jet 1 p_{T} [GeV]"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(17,0,340,"jets_pt[1]", "Jet 2 p_{T} [GeV]"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(12,0,240,"jets_pt[2]", "Jet 3 p_{T} [GeV]"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(12,0,240,"jets_pt[3]", "Jet 4 p_{T} [GeV]"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            tmp_seln = iseln.second+"&&"+njcut+"&& met>100 &&"+nbcuts[inb]+"&&" + xcuts["fullhig"];
            pm.Push<Hist1D>(Axis(10,100,600,"met", "E_{T}^{miss} [GeV]",{150,200,300}), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
          } else if (iseln.first == "notrkphi") {
            if (imet>0) continue; 
            string tmp_seln = iseln.second+"&& ntks==0 && !low_dphi &&"+njcut+"&& met>100 &&"+nbcuts[inb];
            pm.Push<Hist1D>(Axis(10,100,600,"met", "E_{T}^{miss} [GeV]",{150,200,300}), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            tmp_seln = iseln.second+"&& !low_dphi &&"+njcut+"&& met>150 &&"+nbcuts[inb];
            pm.Push<Hist1D>(Axis(5,-0.5,4.5,"ntks", "N_{tks}"),
              tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            tmp_seln = iseln.second+"&& ntks==0 &&"+njcut+"&& met>150 &&"+nbcuts[inb];
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

  //    1D distributions with data
  //----------------------------------------
  if ((do_stackplots || do_cr_stackplots) && do_data) {
    for (auto &iseln: selns) {
      NamedFunc wgt = "weight" * Functions::eff_higtrig;
      //decide which is the relevant set of procs
      string iskimcat = "lep0";
      if (iseln.first=="lep1" || iseln.first=="lep2") iskimcat = iseln.first;
      string iproc = iskimcat=="lep0" ? ("sig"+iskimcat) : iskimcat; //want to include signal for selections with 0 leptons
      //for the moment skip control regions
      if (!do_cr_stackplots && (iseln.first=="lep1" || iseln.first=="lep2" || iseln.first=="qcd")) continue;
      for(unsigned imet(0); imet<metcuts[iseln.first].size(); imet++) { 
        for(unsigned inb(0); inb<nbcuts.size(); inb++) {
          if (iskimcat=="lep0" && !unblind && inb>0) continue;         
          if (iseln.first!= "notrkphi") {
            pm.Push<Hist1D>(Axis(15,0,150,"hig_dm", "#Deltam [GeV]", {40.}),
              iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&&"+nbcuts[inb], 
              procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m> [GeV]", {100., 140.}),
              iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&&"+nbcuts[inb]+"&&hig_dm<=40 && hig_drmax<=2.2", 
              procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m> [GeV]", {100., 140.}),
              iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&&"+nbcuts[inb]+"&&hig_dm<=40", 
              procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(20,0,4,"hig_drmax", "#DeltaR_{max}", {2.2}),
              iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&&"+nbcuts[inb]+"&&hig_dm<=40", 
              procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            string tmp_seln = iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&&"+nbcuts[inb];
            pm.Push<Hist1D>(Axis(15,0,600,"jets_pt[0]", "Jet 1 p_{T} [GeV]"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(17,0,340,"jets_pt[1]", "Jet 2 p_{T} [GeV]"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(12,0,240,"jets_pt[2]", "Jet 3 p_{T} [GeV]"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            pm.Push<Hist1D>(Axis(12,0,240,"jets_pt[3]", "Jet 4 p_{T} [GeV]"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            tmp_seln = iseln.second+"&&"+njcut+"&& met>150 &&"+nbcuts[inb];
            pm.Push<Hist1D>(Axis(18,150,600,"met", "E_{T}^{miss} [GeV]",{150,200,300}), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            if (iseln.first=="lep1"){
              pm.Push<Hist1D>(Axis(18,150,600,"met", "E_{T}^{miss} [GeV]",{150,200,300}), tmp_seln+"&&nels==1", procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
              pm.Push<Hist1D>(Axis(18,150,600,"met", "E_{T}^{miss} [GeV]",{150,200,300}), tmp_seln+"&&nmus==1", procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            }
            tmp_seln = iseln.second+"&&"+njcut+"&&"+metcuts[iseln.first][imet]+"&&"+nbcuts[inb];
            pm.Push<Hist1D>(Axis(20,0,2000,"ht", "H_{T} [GeV]"), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
          } else if (iseln.first == "notrkphi") {
            if (imet>0) continue; 
            string tmp_seln = iseln.second+"&& ntks==0 && !low_dphi &&"+njcut+"&& met>100 &&"+nbcuts[inb];
            pm.Push<Hist1D>(Axis(10,100,600,"met", "E_{T}^{miss} [GeV]",{150,200,300}), tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            tmp_seln = iseln.second+"&& !low_dphi &&"+njcut+"&& met>150 &&"+nbcuts[inb];
            pm.Push<Hist1D>(Axis(5,-0.5,4.5,"ntks", "N_{tks}"),
              tmp_seln, procs[iproc], all_plot_types).Weight(wgt).Tag(iseln.first);
            tmp_seln = iseln.second+"&& ntks==0 &&"+njcut+"&& met>150 &&"+nbcuts[inb];
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

  //        Cutflow table
  //-------------------------------- 
  if (do_cflow && proc_types.size()>0) {
    NamedFunc wgt = "weight" * Functions::eff_higtrig;
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
  if (do_shapes || do_data_shapes) {
    for(auto &iseln: selns) {
      string iskimcat = "lep0";
      if (iseln.first=="lep1" || iseln.first=="lep2") iskimcat = iseln.first;
      string iproc = do_data_shapes ? "data_shapes"+iskimcat : "shapes"+iskimcat; //fixme: want both opts
      if (iseln.first=="base" && !unblind && do_data_shapes) continue;
      for(auto &imet: {"met>=150"}){//metcuts[iseln.first]) { 
        pm.Push<Hist1D>(Axis(12,0,240,"hig_am", "<m> [GeV]", {100., 140.}),
                iseln.second+"&&"+njcut+"&&"+imet, 
                procs[iproc], all_shapes).Tag("datavsdata_"+iseln.first); //fixme: tag for mc
        pm.Push<Hist1D>(Axis(12,0,240,"hig_am", "<m> [GeV]", {100., 140.}),
                iseln.second+"&&"+njcut+"&&"+imet+"&&hig_dm<=40", 
                procs[iproc], all_shapes).Tag("datavsdata_"+iseln.first);
        pm.Push<Hist1D>(Axis(10,0,200,"hig_am", "<m> [GeV]", {100., 140.}),
                iseln.second+"&&"+njcut+"&&"+imet+"&&hig_dm<=40 && hig_am<=200", 
                procs[iproc], all_shapes).Tag("datavsdata_"+iseln.first);
        pm.Push<Hist1D>(Axis(12,0,240,"hig_am", "<m> [GeV]", {100., 140.}),
                iseln.second+"&&"+njcut+"&&"+imet+"&&hig_dm<=40 && hig_drmax<=2.2", 
                procs[iproc], all_shapes).Tag("datavsdata_"+iseln.first);      
      }
    }
  }

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
