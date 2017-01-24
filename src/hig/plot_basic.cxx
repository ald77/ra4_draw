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
  float lumi = 4.3;
  string sample = "search";
  string json = "1";
  bool do_data = true;
  bool unblind = false;
  bool do_loose = false; // removes track veto and delta phi requirement for the search region to make dphi "N-1" plots
  // simplify selection
  bool do_metint = true; // only integrated met>150 and met > 200
  bool do_ge3b = false;  // do btag cats: 2b and 3b+
  bool subtr_ttx = false;
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
  
int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);
  if (json=="json4p0") lumi = 4.3;
  else if (json=="1") lumi = 36.2; 
  else {
    cout<<"Json "<<json<<" has not been implemented!!"<<endl;
    exit(0);
  }
  // protect from silly options
  if (sample=="ttbar" || sample=="search" || !do_data) subtr_ttx = false;

  time_t begtime, endtime;
  time(&begtime);

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
  vector<PlotOpt> linplot = {lin_lumi_info};
  PlotOpt lin_lumi_info_print = lin_lumi().Title(TitleType::info).Bottom(BottomType::ratio).PrintVals(true);
  vector<PlotOpt> linplotprint = {lin_lumi_info_print};
  PlotOpt log_lumi_info_print = log_lumi().Title(TitleType::info).Bottom(BottomType::ratio).PrintVals(true);
  vector<PlotOpt> logplotprint = {log_lumi_info_print};
  vector<PlotOpt> logplot = {log_lumi_info};
  vector<PlotOpt> all_shapes = {lin_shapes_info};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/";
  if (sample=="ttbar") foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higlep1met0/";
  if (sample=="zll") foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higlep2/";
  if (sample=="qcd") foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higqcd/";
  string folderdata(bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_higloose/");
  if (sample=="ttbar") folderdata = bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_higlep1met0/";
  if (sample=="zll") folderdata = bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_higlep2/";
  if (sample=="qcd") folderdata = bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_higqcd/";
  string foldersig(bfolder+"/cms2r0/babymaker/babies/2016_08_10/TChiHH/merged_higmc_higloose/");

  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({"*_TTJets*Lept*.root", "*_TTJets_HT*.root", "*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*_ttHJetTobb*.root","*_TTTT*.root"});
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
  NamedFunc wgt_subtr("wgt_subtr",[&](const Baby &b){
    return Higfuncs::wgt_subtr_ttx(b, json);
  });
  NamedFunc wgt = "weight" * Higfuncs::eff_higtrig;
  if (subtr_ttx) wgt *= wgt_subtr;

  string base_func("njets>=4 && njets<=5 && met/met_calo<5"); //met/met_calo
  // zll skim: ((elel_m>80&&elel_m<100)||(mumu_m>80&&mumu_m<100)) && 
  // nleps==2 && Max$(leps_pt)>40 && (njets==4||njets==5)
  if (sample=="zll") base_func = base_func+"&& nleps==2 && met<50";
  // qcd skim - met>150 && nvleps==0 && (njets==4||njets==5)
  if (sample=="qcd") base_func = base_func+"&& ntks==0 && nvleps==0 && low_dphi";
  // ttbar skim - met>100 && nleps==1 && (njets==4||njets==5) && nbt>=2
  if (sample=="ttbar") base_func = base_func+"&& nleps==1 && mt<100";
  // search skim - met>100 && nvleps==0 && (njets==4||njets==5) && nbt>=2
  if (sample=="search") {
    if (do_loose) base_func = base_func+"&& nvleps==0";
    else base_func = base_func+"&& nvleps==0 && ntks==0 && !low_dphi";
  } 

  vector<shared_ptr<Process> > procs;
  if (!subtr_ttx) 
    procs.push_back(Process::MakeShared<Baby_full>("t#bar{t}+X", 
      Process::Type::background, colors("tt_1l"),    attach_folder(foldermc,mctags["ttx"]),     
      base_func+"&& pass && pass_ra2_badmu && stitch"));
  procs.push_back(Process::MakeShared<Baby_full>("V+jets",     
    Process::Type::background, kOrange+1,          attach_folder(foldermc,mctags["vjets"]),   
    base_func+"&& pass && pass_ra2_badmu && stitch"));
  procs.push_back(Process::MakeShared<Baby_full>("Single t",   
    Process::Type::background, colors("single_t"), attach_folder(foldermc,mctags["singlet"]), 
    base_func+"&& pass && pass_ra2_badmu && stitch"));
  procs.push_back(Process::MakeShared<Baby_full>("QCD",        
    Process::Type::background, colors("other"),    attach_folder(foldermc,mctags["qcd"]),     
    base_func+"&& pass && pass_ra2_badmu && stitch")); 
  procs.push_back(Process::MakeShared<Baby_full>("Other",      
    Process::Type::background, kGreen+1,           attach_folder(foldermc,mctags["other"]),   
    base_func+"&& pass && pass_ra2_badmu && stitch"));      

  if (do_data) {
    if (subtr_ttx) {
      set<string> tt_data_files = attach_folder(foldermc, mctags["ttx"]);
      tt_data_files.insert(folderdata+"*root");
      // data or MC specific cuts will be applied via the weight
      procs.push_back(Process::MakeShared<Baby_full>("Data - t#bar{t}X", Process::Type::data, kBlack,
        tt_data_files, "pass && pass_ra2_badmu &&"+base_func)); 
    } else {
      procs.push_back(Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
        {folderdata+"*root"},  Higfuncs::trig_hig>0. && " pass && pass_ra2_badmu &&"+json+"&&"+base_func)); 
    }
  }
  if (sample == "search") {
    for (unsigned isig(0); isig<sigm.size(); isig++)
      procs.push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
        sig_colors[isig], {foldersig+"*TChiHH_mGluino-"+sigm[isig]+"*.root"}, base_func +" && pass_ra2_badmu"));
  }

  PlotMaker pm;

  vector<string> metcuts;
  string metdef = "met";
  if (sample=="zll") metdef = "(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))";
  if (do_metint) {
    if (sample=="zll") metcuts.push_back(metdef+">0");
    else if (sample=="ttbar") metcuts.push_back(metdef+">0");
    metcuts.push_back(metdef+">150");
  } else {
    metcuts.push_back(metdef+">150&&"+metdef+"<=200");
    metcuts.push_back(metdef+">200&&"+metdef+"<=300");
    metcuts.push_back(metdef+">300&&"+metdef+"<=450");
    metcuts.push_back(metdef+">450");
  }
  
  vector<string> nbcuts;
  if (sample=="zll" || sample=="qcd") {
    nbcuts.push_back("nbm==0");
    nbcuts.push_back("nbm==1");
  } 
  nbcuts.push_back("nbt==2&&nbm==2");
  if (sample!="search" || unblind || !do_data) {
    if (do_ge3b) {
      nbcuts.push_back("nbt>=2&&nbm>=3");
    } else {
      nbcuts.push_back("nbt>=2&&nbm==3&&nbl==3");
      nbcuts.push_back("nbt>=2&&nbm>=3&&nbl>=4");
    }
  }

  map<string, string> xcuts; // useful additional cut definitions
  xcuts["nm1"] = base_func;
  xcuts["base"] = base_func+"&& hig_dm<=40 && hig_am<200";
  // xcuts["drmax"] = base_func+"&& hig_dm<=40 && hig_am<200 && hig_drmax<=2.2";
  // xcuts["hig"] = "hig_am>100 && hig_am<=140 && hig_dm <= 40 && hig_drmax<=2.2";
  // xcuts["sbd"] = "(hig_am<=100 || (hig_am>140 && hig_am<=200)) && hig_dm <= 40 && hig_drmax<=2.2";
 
  string tmp_seln = base_func;
  // nb plots integrated in MET
  if (metcuts.size()>0) tmp_seln += "&&"+metcuts[0];
  // pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbl", "N_{b}^{L}"), tmp_seln, procs, linplot).Weight(wgt).Tag(sample);
  // pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbm", "N_{b}^{M}"), tmp_seln, procs, linplot).Weight(wgt).Tag(sample);
  // pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbt", "N_{b}^{T}"), tmp_seln, procs, linplot).Weight(wgt).Tag(sample);
  if (sample=="search") {
    if (!do_data || unblind) {
      pm.Push<Hist1D>(Axis(5,0.5,5.5,Higfuncs::hig_nb, "b-tag category (TTML)"), 
        tmp_seln && Higfuncs::hig_nb>0., procs, linplot).Weight(wgt).Tag(sample);
      pm.Push<Hist1D>(Axis(5,0.5,5.5,Higfuncs::hig_nb_ttll, "b-tag category (TTLL)"), 
        tmp_seln && Higfuncs::hig_nb_ttll>0., procs, linplot).Weight(wgt).Tag(sample);
      pm.Push<Hist1D>(Axis(5,0.5,5.5,Higfuncs::hig_nb_tmml, "b-tag category (TMML)"), 
        tmp_seln && Higfuncs::hig_nb_tmml>0., procs, linplot).Weight(wgt).Tag(sample);
      pm.Push<Hist1D>(Axis(5,0.5,5.5,Higfuncs::hig_nb_mmmm, "b-tag category (MMMM)"), 
        tmp_seln && Higfuncs::hig_nb_mmmm>0., procs, linplot).Weight(wgt).Tag(sample);
    }
  } else { 
    pm.Push<Hist1D>(Axis(6,-0.5,5.5,Higfuncs::hig_nb_extended, "Extended b-tag categories (TTML)"), 
      tmp_seln && Higfuncs::hig_nb_extended<6, procs, linplot).Weight(wgt).Tag(sample);
    pm.Push<Hist1D>(Axis(5,0.5,5.5,Higfuncs::hig_nb, "b-tag categories (TTML)"), 
      tmp_seln && Higfuncs::hig_nb>0., procs, linplotprint).Weight(wgt).Tag(sample); 
    pm.Push<Hist1D>(Axis(5,0.5,5.5,Higfuncs::hig_nb, "b-tag categories (TTML)"), 
      tmp_seln && Higfuncs::hig_nb>0. && "hig_dm<=40 && hig_am<=200", procs, linplot).Weight(wgt).Tag(sample);            
    pm.Push<Hist1D>(Axis(5,0.5,5.5,Higfuncs::hig_nb, "b-tag categories (TTML)"), 
      tmp_seln && Higfuncs::hig_nb>0. && "hig_dm<=40 && hig_am<=200 && hig_drmax<=2.2", procs, linplot).Weight(wgt).Tag(sample);    
    cout<<tmp_seln<<endl;
  }
  if (sample!="search") {
    vector<double> metbins = {150, 200, 300, 600};
    vector<double> metbins_ext = {0,50, 100, 150, 200, 300, 600};
    if (sample=="zll") pm.Push<Hist1D>(Axis(metbins,metdef, "p_{T}^{Z} [GeV]",{150,200,300}), 
      tmp_seln && "nbt>=2" && metdef+">150", procs, logplotprint).Weight(wgt).Tag(sample+"_norm");
    else {
      if (sample=="ttbar") 
        pm.Push<Hist1D>(Axis(metbins_ext,"met", "E_{T}^{miss} [GeV]",{150,200,300}), 
          tmp_seln && "nbt>=2", procs, logplotprint).Weight(wgt).Tag(sample+"_normext");
        pm.Push<Hist1D>(Axis(metbins,"met", "E_{T}^{miss} [GeV]",{150,200,300}), 
          tmp_seln && "nbt>=2" && metdef+">"+to_string(metbins[0]), procs, logplotprint).Weight(wgt).Tag(sample+"_norm");
    }
  }


  if (!subtr_ttx) {
    for (auto &ixcut: xcuts) {
      tmp_seln = base_func;
      for(unsigned imet(0); imet<metcuts.size(); imet++) { 
        tmp_seln = metcuts[imet];
        for(unsigned inb(0); inb<nbcuts.size(); inb++) {
          if (ixcut.first=="nm1") { 
            pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m> [GeV]", {100., 140.}),
              ixcut.second+"&&"+metcuts[imet]+"&&"+nbcuts[inb]+"&&hig_dm<40", 
              procs, linplot).Weight(wgt).Tag(sample);
            pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m> [GeV]", {100., 140.}),
              ixcut.second+"&&"+metcuts[imet]+"&&"+nbcuts[inb]+"&&hig_dm<40 && hig_drmax<=2.2", 
              procs, linplot).Weight(wgt).Tag(sample);
            tmp_seln = ixcut.second+"&&"+metcuts[imet]+"&&"+nbcuts[inb];
            pm.Push<Hist1D>(Axis(15,0,150,"hig_dm", "#Deltam [GeV]", {40.}), 
              tmp_seln, procs, linplot).Weight(wgt).Tag(sample);
          }
          tmp_seln = ixcut.second+"&&"+metcuts[imet]+"&&"+nbcuts[inb];
          pm.Push<Hist1D>(Axis(20,0,2000,"ht", "H_{T} [GeV]"), 
            tmp_seln, procs, logplot).Weight(wgt).Tag(sample);
          if (imet==0) {
            if (sample=="zll") pm.Push<Hist1D>(Axis(24,0,600,metdef, "p_{T}^{Z} [GeV]",{150,200,300}), 
              tmp_seln, procs, logplot).Weight(wgt).Tag(sample);
            else pm.Push<Hist1D>(Axis(24,0,600,"met", "E_{T}^{miss} [GeV]",{150,200,300}), 
              tmp_seln, procs, logplot).Weight(wgt).Tag(sample);
          }  
          if (ixcut.first=="base") // do only with the trimmed selection
            pm.Push<Hist1D>(Axis(20,0,4,"hig_drmax", "#DeltaR_{max}", {2.2}),
            ixcut.second+"&&"+metcuts[imet]+"&&"+nbcuts[inb], 
            procs, linplot).Weight(wgt).Tag(sample);
          // pm.Push<Hist1D>(Axis(15,0,600,"jets_pt[0]", "Jet 1 p_{T} [GeV]"), 
          //   tmp_seln, procs, linplot).Weight(wgt).Tag(sample);
          // pm.Push<Hist1D>(Axis(17,0,340,"jets_pt[1]", "Jet 2 p_{T} [GeV]"), 
          //   tmp_seln, procs, linplot).Weight(wgt).Tag(sample);
          // pm.Push<Hist1D>(Axis(12,0,240,"jets_pt[2]", "Jet 3 p_{T} [GeV]"), 
          //   tmp_seln, procs, linplot).Weight(wgt).Tag(sample);
          // pm.Push<Hist1D>(Axis(12,0,240,"jets_pt[3]", "Jet 4 p_{T} [GeV]"), 
          //   tmp_seln, procs, linplot).Weight(wgt).Tag(sample);
          // if (sample=="ttbar" && !do_note){
          //   tmp_seln = ixcut.second+"&&"+metcuts[imet]+"&&"+nbcuts[inb];
          //   pm.Push<Hist1D>(Axis(18,150,600,"met", "E_{T}^{miss} [GeV]",{150,200,300}), 
          //     tmp_seln+"&&nels==1", procs, linplot).Weight(wgt).Tag(sample);
          //   pm.Push<Hist1D>(Axis(18,150,600,"met", "E_{T}^{miss} [GeV]",{150,200,300}), 
          //     tmp_seln+"&&nmus==1", procs, linplot).Weight(wgt).Tag(sample);
          // }
          
            // else if (sample=="search" && do_loose) {
          //   if (imet>0) continue; 
          //   tmp_seln = ixcut.second+"&& ntks==0 && !low_dphi && met>100 &&"+nbcuts[inb];
          //   pm.Push<Hist1D>(Axis(10,100,600,"met", "E_{T}^{miss} [GeV]",{150,200,300}), tmp_seln, procs, linplot).Weight(wgt).Tag(sample);
          //   tmp_seln = ixcut.second+"&& !low_dphi && met>150 &&"+nbcuts[inb];
          //   pm.Push<Hist1D>(Axis(5,-0.5,4.5,"ntks", "N_{tks}"),
          //     tmp_seln, procs, linplot).Weight(wgt).Tag(sample);
          //   tmp_seln = ixcut.second+"&& ntks==0 && met>150 &&"+nbcuts[inb];
          //   pm.Push<Hist1D>(Axis(32,0.,3.2,"dphi2", "#Delta#phi_{2}",{0.5}),
          //     tmp_seln+"&& dphi1>0.5", procs, linplot).Weight(wgt).Tag(sample);
          //   pm.Push<Hist1D>(Axis(32,0.,3.2,"dphi3", "#Delta#phi_{3}",{0.3}),
          //     tmp_seln+"&& dphi1>0.5 && dphi2>0.5", procs, linplot).Weight(wgt).Tag(sample);
          //   pm.Push<Hist1D>(Axis(32,0.,3.2,"dphi4", "#Delta#phi_{4}",{0.3}),
          //     tmp_seln+"&& dphi1>0.5 && dphi2>0.5 && dphi3>0.3", procs, linplot).Weight(wgt).Tag(sample);
          // }
        }
      }
    }
  }

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"lumi", required_argument, 0, 'l'},    // Luminosity to normalize MC with (no data)
      {"sample", required_argument, 0, 's'},    // Which sample to use: standard, met150, 2015 data
      {"subttx", no_argument, 0, 0},             
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s:l:", long_options, &option_index);
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
      optname = long_options[option_index].name;
      if(optname == "subttx"){
        subtr_ttx = true;
      }else {
        printf("Bad option! Found option name %s\n", optname.c_str());
        exit(1);
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
