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
  // float lumi = 36.2;
  string sample = "search";
  float lumi = 4.3;
  string json = "json4p0";
  bool do_note = true;
  bool do_data = true;
  bool unblind = false;
  bool do_loose = false; // removes track veto and delta phi requirement for the search region to make dphi "N-1" plots
  // simplify selection
  bool do_metg150_metg200 = true; // only integrated met>150 and met > 200
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
  
int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

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
  vector<PlotOpt> all_plot_types = {lin_lumi_info};//, log_lumi_info};
  vector<PlotOpt> all_shapes = {lin_shapes_info};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/";
  if (sample=="ttbar") foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higlep1/";
  if (sample=="zll") foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_nj4zcandl40/";
  if (sample=="qcd") foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higqcd/";
  string folderdata(bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_higloose/");
  if (sample=="ttbar") folderdata = bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_higlep1/";
  if (sample=="zll") folderdata = bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_nj4zcandl40/";
  if (sample=="qcd") folderdata = bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_nl0nj4met150/";
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
  NamedFunc wgt = "weight" * Higfuncs::eff_higtrig;
  string base_func("njets>=4 && njets<=5 && met/met_calo<5"); //met/met_calo
  // zll skim: ((elel_m>80&&elel_m<100)||(mumu_m>80&&mumu_m<100)) && 
  // nleps==2 && Max$(leps_pt)>40
  if (sample=="zll") base_func = base_func+"&& nleps==2 && met<50";
  // qcd skim - met>150 && nvleps==0 && (njets==4||njets==5)
  if (sample=="qcd") base_func = base_func+"&& nvleps==0 && ntks==0 && low_dphi";
  // ttbar skim - met>100 && nleps==1 && (njets==4||njets==5) && nbm>=2
  if (sample=="ttbar") base_func = base_func+"&& nleps==1 && mt<100";
  // search skim - met>100 && nvleps==0 && (njets==4||njets==5) && nbm>=2
  if (sample=="search") {
    if (do_loose) base_func = base_func+"&& nvleps==0";
    else base_func = base_func+"&& nvleps==0 && ntks==0 && !low_dphi";
  } 

  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_full>("t#bar{t}+X", 
    Process::Type::background, colors("tt_1l"),    attach_folder(foldermc,mctags["ttx"]),     base_func+"&& pass && stitch"));
  procs.push_back(Process::MakeShared<Baby_full>("V+jets",     
    Process::Type::background, kOrange+1,          attach_folder(foldermc,mctags["vjets"]),   base_func+"&& pass && stitch"));
  procs.push_back(Process::MakeShared<Baby_full>("Single t",   
    Process::Type::background, colors("single_t"), attach_folder(foldermc,mctags["singlet"]), base_func+"&& pass && stitch"));
  procs.push_back(Process::MakeShared<Baby_full>("QCD",        
    Process::Type::background, colors("other"),    attach_folder(foldermc,mctags["qcd"]),     base_func+"&& pass && stitch")); 
  procs.push_back(Process::MakeShared<Baby_full>("Other",      
    Process::Type::background, kGreen+1,           attach_folder(foldermc,mctags["other"]),   base_func+"&& pass && stitch"));      

  if (do_data) {
    procs.push_back(Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
      {folderdata+"*RunB*root"},  Higfuncs::trig_hig>0. && " pass &&"+json+"&&"+base_func)); 
  }
  // need to modify base_func to use this
  // if (sample == "search") {
  //   for (unsigned isig(0); isig<sigm.size(); isig++)
  //     procs.push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
  //       sig_colors[isig], {foldersig+"*TChiHH_mGluino-"+sigm[isig]+"*.root"}, ));
  // }

  PlotMaker pm;

  vector<string> metcuts;
  string metdef = "met";
  if (sample=="zll") metdef = "(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))";
  if (do_metg150_metg200) {
    metcuts.push_back(metdef+">150");
    // metcuts.push_back(metdef+">200");
  } else {
    metcuts.push_back(metdef+">150&&"+metdef+"<=200");
    metcuts.push_back(metdef+">200&&"+metdef+"<=300");
    metcuts.push_back(metdef+">300");
  }
  
  vector<string> nbcuts;
  if (sample=="zll" || sample=="qcd") {
    nbcuts.push_back("nbm==0");
    nbcuts.push_back("nbm==1");
  } else {
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
  }

  map<string, string> xcuts; // useful additional cut definitions
  xcuts["nm1"] = base_func;
  xcuts["base"] = base_func+"&& hig_dm<=40 && hig_am<200";
  // xcuts["drmax"] = base_func+"&& hig_dm<=40 && hig_am<200 && hig_drmax<=2.2";
  // xcuts["hig"] = "hig_am>100 && hig_am<=140 && hig_dm <= 40 && hig_drmax<=2.2";
  // xcuts["sbd"] = "(hig_am<=100 || (hig_am>140 && hig_am<=200)) && hig_dm <= 40 && hig_drmax<=2.2";


  // Temporary funcs for other b-tag category options
  //-----------------------------------------------------
  // const NamedFunc hig_nb_mmmm("hig_nb_mmmm",[](const Baby &b) -> NamedFunc::ScalarType{
  //   if (b.nbm()>=2) return min(4,b.nbm());
  //   else return 0;
  // });
  // const NamedFunc hig_nb_ttll("hig_nb_ttll",[](const Baby &b) -> NamedFunc::ScalarType{
  //   if (b.nbt()==2 && b.nbl()==2) return 2;
  //   else if (b.nbt()>=2 && b.nbl()==3) return 3;
  //   else if (b.nbt()>=2 && b.nbl()>=4) return 4;
  //   else return 0;
  // });
  // const NamedFunc hig_nb_tmml("hig_nb_tmml",[](const Baby &b) -> NamedFunc::ScalarType{
  //   if (b.nbt()>=1 && b.nbm()==2) return 2;
  //   else if (b.nbt()>=1 && b.nbm()==3 && b.nbl()==3) return 3;
  //   else if (b.nbt()>=1 && b.nbm()>=3 && b.nbl()>=4) return 4;
  //   else return 0;
  // });
  // vector<NamedFunc> btag_xopts = {hig_nb_mmmm, hig_nb_tmml, hig_nb_ttll};

  //    1D distributions with data
  //----------------------------------------
  if (do_data) {
    for (auto &ixcut: xcuts) {
      for(unsigned imet(0); imet<metcuts.size(); imet++) { 
        for(unsigned inb(0); inb<nbcuts.size(); inb++) {
          if (sample=="search" && !unblind && inb>0) continue;         
          if (!do_loose) {
            if (ixcut.first=="nm1") { // do only in the loosest selection
              pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m> [GeV]", {100., 140.}),
              ixcut.second+"&&"+metcuts[imet]+"&&"+nbcuts[inb]+"&&hig_dm<40", 
              procs, all_plot_types).Weight(wgt).Tag(sample);
              pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m> [GeV]", {100., 140.}),
              ixcut.second+"&&"+metcuts[imet]+"&&"+nbcuts[inb]+"&&hig_dm<40 && hig_drmax<=2.2", 
              procs, all_plot_types).Weight(wgt).Tag(sample);
              string tmp_seln = ixcut.second+"&&"+metcuts[imet]+"&&"+nbcuts[inb];
              if (!do_note || sample=="search") {
                pm.Push<Hist1D>(Axis(20,0,2000,"ht", "H_{T} [GeV]"), tmp_seln, procs, all_plot_types).Weight(wgt).Tag(sample);
                pm.Push<Hist1D>(Axis(18,150,600,"met", "E_{T}^{miss} [GeV]",{150,200,300}), tmp_seln, procs, all_plot_types).Weight(wgt).Tag(sample);
                pm.Push<Hist1D>(Axis(15,0,150,"hig_dm", "#Deltam [GeV]", {40.}), tmp_seln, procs, all_plot_types).Weight(wgt).Tag(sample);
              }
            }
            if (ixcut.first=="base") // do only with the trimmed selection
              pm.Push<Hist1D>(Axis(20,0,4,"hig_drmax", "#DeltaR_{max}", {2.2}),
              ixcut.second+"&&"+metcuts[imet]+"&&"+nbcuts[inb], 
              procs, all_plot_types).Weight(wgt).Tag(sample);
            // pm.Push<Hist1D>(Axis(15,0,600,"jets_pt[0]", "Jet 1 p_{T} [GeV]"), tmp_seln, procs, all_plot_types).Weight(wgt).Tag(sample);
            // pm.Push<Hist1D>(Axis(17,0,340,"jets_pt[1]", "Jet 2 p_{T} [GeV]"), tmp_seln, procs, all_plot_types).Weight(wgt).Tag(sample);
            // pm.Push<Hist1D>(Axis(12,0,240,"jets_pt[2]", "Jet 3 p_{T} [GeV]"), tmp_seln, procs, all_plot_types).Weight(wgt).Tag(sample);
            // pm.Push<Hist1D>(Axis(12,0,240,"jets_pt[3]", "Jet 4 p_{T} [GeV]"), tmp_seln, procs, all_plot_types).Weight(wgt).Tag(sample);
            // if (sample=="ttbar"){
            //   string tmp_seln = ixcut.second+"&&"+metcuts[imet]+"&&"+nbcuts[inb];
            //   pm.Push<Hist1D>(Axis(18,150,600,"met", "E_{T}^{miss} [GeV]",{150,200,300}), tmp_seln+"&&nels==1", procs, all_plot_types).Weight(wgt).Tag(sample);
            //   pm.Push<Hist1D>(Axis(18,150,600,"met", "E_{T}^{miss} [GeV]",{150,200,300}), tmp_seln+"&&nmus==1", procs, all_plot_types).Weight(wgt).Tag(sample);
            // }
          } else if (sample=="search" && do_loose) {
            // if (imet>0) continue; 
            // string tmp_seln = ixcut.second+"&& ntks==0 && !low_dphi && met>100 &&"+nbcuts[inb];
            // pm.Push<Hist1D>(Axis(10,100,600,"met", "E_{T}^{miss} [GeV]",{150,200,300}), tmp_seln, procs, all_plot_types).Weight(wgt).Tag(sample);
            // tmp_seln = ixcut.second+"&& !low_dphi && met>150 &&"+nbcuts[inb];
            // pm.Push<Hist1D>(Axis(5,-0.5,4.5,"ntks", "N_{tks}"),
            //   tmp_seln, procs, all_plot_types).Weight(wgt).Tag(sample);
            // tmp_seln = ixcut.second+"&& ntks==0 && met>150 &&"+nbcuts[inb];
            // pm.Push<Hist1D>(Axis(32,0.,3.2,"dphi2", "#Delta#phi_{2}",{0.5}),
            //   tmp_seln+"&& dphi1>0.5", procs, all_plot_types).Weight(wgt).Tag(sample);
            // pm.Push<Hist1D>(Axis(32,0.,3.2,"dphi3", "#Delta#phi_{3}",{0.3}),
            //   tmp_seln+"&& dphi1>0.5 && dphi2>0.5", procs, all_plot_types).Weight(wgt).Tag(sample);
            // pm.Push<Hist1D>(Axis(32,0.,3.2,"dphi4", "#Delta#phi_{4}",{0.3}),
            //   tmp_seln+"&& dphi1>0.5 && dphi2>0.5 && dphi3>0.3", procs, all_plot_types).Weight(wgt).Tag(sample);
          }
        }
      }
    }
  }

  //     N-1 and other 1D MC-only distributions
  //----------------------------------------
  else {
    wgt = "weight";
    for(unsigned imet(0); imet<metcuts.size(); imet++) { 
      // if (!do_loose) {
      //   string tmp_seln = metcuts[imet]+"&& nbm>=2 &&" + xcuts["hig"];
      //   pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbl", "N_{b}^{L}"), tmp_seln, procs, all_plot_types).Weight(wgt).Tag(sample);
      //   pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbm", "N_{b}^{M}"), tmp_seln, procs, all_plot_types).Weight(wgt).Tag(sample);
      //   pm.Push<Hist1D>(Axis(6,0.5,6.5,"nbt", "N_{b}^{T}"), tmp_seln, procs, all_plot_types).Weight(wgt).Tag(sample);
      //   tmp_seln = metcuts[imet]+"&& nbm>=2 && hig_dm<=40 && hig_drmax<=2.2 && hig_am<=200";
      //   pm.Push<Hist1D>(Axis(5,0.5,5.5,Functions::hig_nb, "b-tag category (TTML)"), tmp_seln && Functions::hig_nb>0., procs, all_plot_types).Weight(wgt).Tag(sample);
      //   pm.Push<Hist1D>(Axis(5,0.5,5.5,hig_nb_ttll, "b-tag category (TTLL)"), tmp_seln && hig_nb_ttll>0., procs, all_plot_types).Weight(wgt).Tag(sample);
      //   pm.Push<Hist1D>(Axis(5,0.5,5.5,hig_nb_tmml, "b-tag category (TMML)"), tmp_seln && hig_nb_tmml>0., procs, all_plot_types).Weight(wgt).Tag(sample);
      //   pm.Push<Hist1D>(Axis(5,0.5,5.5,hig_nb_mmmm, "b-tag category (MMMM)"), tmp_seln && hig_nb_mmmm>0., procs, all_plot_types).Weight(wgt).Tag(sample);
      // }
      for(unsigned inb(0); inb<nbcuts.size(); inb++) {
        if (!do_loose) {
          // pm.Push<Hist1D>(Axis(15,0,150,"hig_dm", "#Deltam [GeV]", {40.}),
          //   metcuts[imet]+"&&"+nbcuts[inb]+"&&hig_am>100 && hig_am<=140 && hig_drmax<=2.2", 
          //   procs, all_plot_types).Weight(wgt).Tag(sample);
          pm.Push<Hist1D>(Axis(25,0,250,"hig_am", "<m> [GeV]", {100., 140.}),
            metcuts[imet]+"&&"+nbcuts[inb]+"&&hig_dm<=40 && hig_drmax<=2.2", 
            procs, all_plot_types).Weight(wgt).Tag(sample);
          pm.Push<Hist1D>(Axis(20,0,4,"hig_drmax", "#DeltaR_{max}", {2.2}),
            metcuts[imet]+"&&"+nbcuts[inb]+"&&hig_am>100 && hig_am<=140 && hig_dm <= 40", 
            procs, all_plot_types).Weight(wgt).Tag(sample);
          // string tmp_seln = metcuts[imet]+"&&"+nbcuts[inb]+"&&" + xcuts["hig"];
          // pm.Push<Hist1D>(Axis(15,0,600,"jets_pt[0]", "Jet 1 p_{T} [GeV]"), tmp_seln, procs, all_plot_types).Weight(wgt).Tag(sample);
          // pm.Push<Hist1D>(Axis(17,0,340,"jets_pt[1]", "Jet 2 p_{T} [GeV]"), tmp_seln, procs, all_plot_types).Weight(wgt).Tag(sample);
          // pm.Push<Hist1D>(Axis(12,0,240,"jets_pt[2]", "Jet 3 p_{T} [GeV]"), tmp_seln, procs, all_plot_types).Weight(wgt).Tag(sample);
          // pm.Push<Hist1D>(Axis(12,0,240,"jets_pt[3]", "Jet 4 p_{T} [GeV]"), tmp_seln, procs, all_plot_types).Weight(wgt).Tag(sample);
          pm.Push<Hist1D>(Axis(10,100,600,"met", "E_{T}^{miss} [GeV]",{150,200,300}), metcuts[imet]+"&&"+nbcuts[inb], procs, all_plot_types).Weight(wgt).Tag(sample);
        } else if (sample=="search" && do_loose) {
          // if (imet>0) continue; 
          // string tmp_seln = "ntks==0 && !low_dphi && met>100 &&"+nbcuts[inb];
          // pm.Push<Hist1D>(Axis(10,100,600,"met", "E_{T}^{miss} [GeV]",{150,200,300}), tmp_seln, procs, all_plot_types).Weight(wgt).Tag(sample);
          // tmp_seln = "!low_dphi && met>150 &&"+nbcuts[inb];
          // pm.Push<Hist1D>(Axis(5,-0.5,4.5,"ntks", "N_{tks}"),
          //   tmp_seln, procs, all_plot_types).Weight(wgt).Tag(sample);
          // tmp_seln = "ntks==0 && met>150 &&"+nbcuts[inb];
          // pm.Push<Hist1D>(Axis(32,0.,3.2,"dphi2", "#Delta#phi_{2}",{0.5}),
          //   tmp_seln+"&& dphi1>0.5", procs, all_plot_types).Weight(wgt).Tag(sample);
          // pm.Push<Hist1D>(Axis(32,0.,3.2,"dphi3", "#Delta#phi_{3}",{0.3}),
          //   tmp_seln+"&& dphi1>0.5 && dphi2>0.5", procs, all_plot_types).Weight(wgt).Tag(sample);
          // pm.Push<Hist1D>(Axis(32,0.,3.2,"dphi4", "#Delta#phi_{4}",{0.3}),
          //   tmp_seln+"&& dphi1>0.5 && dphi2>0.5 && dphi3>0.3", procs, all_plot_types).Weight(wgt).Tag(sample);
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
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
