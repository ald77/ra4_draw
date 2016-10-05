#include <cmath>
#include <stdio.h>
#include <chrono>

#include "TError.h"
#include "TVector2.h"

#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist1d.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "core/table.hpp"
#include "core/slide_maker.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  bool do_metbins = true;       
  bool do_met150 = true;        
  bool do_met100 = true;        
}

int main(){
  gErrorIgnoreLevel = 6000;

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  double lumi = 2.6;
  string bfolder("");
  string hostname(execute("echo $HOSTNAME"));
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string ntupletag="metG200"; 
  if(do_met150 && do_met100) ntupletag="";
  string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_met100_stdnj5/");

  Palette colors("txt/colors.txt", "default");


  NamedFunc baseline = "stitch && nleps==1 && nveto==0 && st>500 && met>100 && njets>=5 && nbm>=1 && weight<1";

  set<string> allfiles = {foldermc+"*_TTJets*Lept*"+ntupletag+"*.root", foldermc+"*_TTJets_HT*"+ntupletag+"*.root",
			  foldermc+"*_WJetsToLNu*"+ntupletag+"*.root",foldermc+"*_ST_*"+ntupletag+"*.root",
			  foldermc+"*_TTW*"+ntupletag+"*.root",foldermc+"*_TTZ*"+ntupletag+"*.root",
			  foldermc+"*_TTGJets*"+ntupletag+"*.root",foldermc+"*_TTTT*"+ntupletag+"*.root",
			  foldermc+"*QCD_HT*Inf_Tune*"+ntupletag+"*.root", 
			  foldermc+"*QCD_HT*0_Tune*"+ntupletag+"*.root",foldermc+"*DYJetsToLL*"+ntupletag+"*.root",
			  foldermc+"*_ZJet*"+ntupletag+"*.root",foldermc+"*_ttHJetTobb*"+ntupletag+"*.root",
			  foldermc+"*_WH_HToBB*"+ntupletag+"*.root",foldermc+"*_ZH_HToBB*"+ntupletag+"*.root",
			  foldermc+"*_WWTo*"+ntupletag+"*.root",foldermc+"*_WZ*"+ntupletag+"*.root",
			  foldermc+"*_ZZ_*"+ntupletag+"*.root"
  };

  set<string> ttfiles = {foldermc+"*_TTJets*Lept*"+ntupletag+"*.root", foldermc+"*_TTJets_HT*"+ntupletag+"*.root"};
  set<string> wfiles = {foldermc+"*_WJetsToLNu*"+ntupletag+"*.root"};
  // set<string> stfiles = {foldermc+"*_ST_*"+ntupletag+"*.root"};
  set<string> otherfiles = {foldermc+"*_ST_*"+ntupletag+"*.root", foldermc+"*_TTW*"+ntupletag+"*.root",foldermc+"*_TTZ*"+ntupletag+"*.root",
			  foldermc+"*_TTGJets*"+ntupletag+"*.root",foldermc+"*_TTTT*"+ntupletag+"*.root",
			  foldermc+"*QCD_HT*Inf_Tune*"+ntupletag+"*.root", 
			  foldermc+"*QCD_HT*0_Tune*"+ntupletag+"*.root",foldermc+"*DYJetsToLL*"+ntupletag+"*.root",
			  foldermc+"*_ZJet*"+ntupletag+"*.root",foldermc+"*_ttHJetTobb*"+ntupletag+"*.root",
			  foldermc+"*_WH_HToBB*"+ntupletag+"*.root",foldermc+"*_ZH_HToBB*"+ntupletag+"*.root",
			  foldermc+"*_WWTo*"+ntupletag+"*.root",foldermc+"*_WZ*"+ntupletag+"*.root",
			  foldermc+"*_ZZ_*"+ntupletag+"*.root"};

  auto proc_2l = Process::MakeShared<Baby_full>("tru 2l", Process::Type::background, kPink-2,
						allfiles, baseline && "stitch && ntruleps>=2");
  auto proc_1lgood = Process::MakeShared<Baby_full>("tru 1l m_{T}-m_{T}^{tru}#leq60", Process::Type::background, kPink-6,
						    allfiles, baseline && "stitch && ntruleps==1 && mt-mt_tru<=60");
  auto proc_1lbad = Process::MakeShared<Baby_full>("tru 1l m_{T}-m_{T}^{tru}>60", Process::Type::background, kBlue+1,
						   allfiles, baseline && "stitch && ntruleps==1 && mt-mt_tru>60");

  auto proc_2l_himt = Process::MakeShared<Baby_full>("True 2l, m_{T}>140", Process::Type::background,kCyan-3, //kGreen+1,
						     allfiles, baseline && "stitch && ntruleps>=2 && mt>140");
  auto proc_1lgood_himt = Process::MakeShared<Baby_full>("tru 1l well meas. m_{T}>140", Process::Type::background, kGreen-3,
						    allfiles, baseline && "stitch && ntruleps==1 && mt_tru>=140 && mt>140");
  auto proc_1lbad_himt = Process::MakeShared<Baby_full>("tru 1l mismeas. m_{T}>140", Process::Type::background, kRed,
						   allfiles, baseline && "stitch && ntruleps==1 && mt_tru<=140 && mt>140");
  auto proc_0l_himt = Process::MakeShared<Baby_full>("True 0l, m_{T}>140", Process::Type::background,colors("other"),
						     allfiles, baseline && "stitch && ntruleps==0 && mt>140");

  //// Low mT
  auto proc_lomt = Process::MakeShared<Baby_full>("All bkg. m_{T}#leq140", Process::Type::background, 1,
						  allfiles, baseline && "stitch && mt<=140");
  auto proc_tt_lomt = Process::MakeShared<Baby_full>("t#bar{t} m_{T}#leq140", Process::Type::background, 1,
						ttfiles, baseline && "stitch && mt<=140 && ntruleps==1");

  //// High mT 0l, 1l, 2l
  auto proc_1l_himt = Process::MakeShared<Baby_full>("True 1l, m_{T}>140", Process::Type::background,colors("tt_1l"),
						     allfiles, baseline && "stitch && ntruleps==1 && mt>140");
  auto proc_01l_himt = Process::MakeShared<Baby_full>("True 0-1l, m_{T}>140", Process::Type::background,kRed+1,
						    allfiles, baseline && "stitch && ntruleps<=1 && mt>140");

  //// High mT contributions
  auto proc_tt2l_himt = Process::MakeShared<Baby_full>("t#bar{t} 2l, m_{T}>140", Process::Type::background, 
						       colors("tt_2l"), ttfiles, 
						       baseline && "stitch && ntruleps>=2 && mt>140");
  auto proc_tt1l_himt = Process::MakeShared<Baby_full>("t#bar{t} 1l, m_{T}>140", Process::Type::background, 
						       colors("tt_1l"), ttfiles, 
						       baseline && "stitch && ntruleps==1 && mt>140");
  auto proc_w_himt = Process::MakeShared<Baby_full>("W+jets, m_{T}>140", Process::Type::background, 
						    colors("wjets"), wfiles, 
						    baseline && "stitch && mt>140");
  // auto proc_st_himt = Process::MakeShared<Baby_full>("Single t, m_{T}>140", Process::Type::background, 
  // 						    colors("single_t"), stfiles, 
  // 						    baseline && "stitch && mt>140");
  auto proc_other_himt = Process::MakeShared<Baby_full>("Other, m_{T}>140", Process::Type::background, 
							colors("other"), otherfiles, 
							baseline && "stitch && mt>140");
  

  //vector<shared_ptr<Process> > all_procs = {proc_1lgood, proc_2l, proc_1lbad, proc_0l};
  vector<shared_ptr<Process> > goodbad_procs = {proc_lomt, proc_1lbad_himt, proc_1lgood_himt, proc_2l_himt};
  vector<shared_ptr<Process> > goodbad0_procs = {proc_lomt, proc_0l_himt, proc_1lbad_himt, proc_1lgood_himt, proc_2l_himt};
  vector<shared_ptr<Process> > proc_procs = {proc_lomt, proc_other_himt, 
					    proc_w_himt, proc_tt1l_himt, proc_tt2l_himt};
  vector<shared_ptr<Process> > lep_procs = {proc_lomt, proc_01l_himt, proc_2l_himt};
  
  PlotOpt log_shapes("txt/plot_styles.txt", "CMSPaper");
  log_shapes.Title(TitleType::info)
  .Bottom(BottomType::ratio)
  .YAxis(YAxisType::log)
  .Stack(StackType::shapes)
  .RatioMaximum(2.8);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  vector<PlotOpt> plot_types = {lin_shapes};


  vector<NamedFunc> metbins;
  if (do_met100) metbins.push_back(NamedFunc("met>100&&met<=150"));
  if (do_met150) metbins.push_back(NamedFunc("met>150&&met<=200"));
  if (do_metbins) {
    metbins.push_back(NamedFunc("met>200&&met<=350"));
    metbins.push_back(NamedFunc("met>350&&met<=500"));
    metbins.push_back(NamedFunc("met>500"));
  } else {
    metbins.push_back(NamedFunc("met>200"));
  } 
  vector<string> nobjbins;
  nobjbins = vector<string>();
  // nobjbins.push_back("njets==5");
  // nobjbins.push_back("njets>=6 && njets<=8");
  // nobjbins.push_back("njets>=9");
  nobjbins.push_back("njets>=6");
  nobjbins.push_back("njets>=9");
  
  int Nplots=0;
  PlotMaker pm;
  for (auto &imet: metbins) {
    for (auto &inobj: nobjbins) {
      NamedFunc icut = imet && inobj;
      pm.Push<Hist1D>(Axis(10, 100., 850., "mj14", "M_{J} [GeV]", {250., 400.}),
		      icut, goodbad_procs, plot_types).Tag("thin_goodbad");
      pm.Push<Hist1D>(Axis(10, 100., 850., "mj14", "M_{J} [GeV]", {250., 400.}),
		      icut, goodbad0_procs, plot_types).Tag("thin_goodbad0");
      pm.Push<Hist1D>(Axis(10, 100., 850., "mj14", "M_{J} [GeV]", {250., 400.}),
		      icut, proc_procs, plot_types).Tag("thin_proc");
      pm.Push<Hist1D>(Axis(10, 100., 850., "mj14", "M_{J} [GeV]", {250., 400.}),
		      icut, lep_procs, plot_types).Tag("thin_lep");

      pm.Push<Hist1D>(Axis(5, 100., 850., "mj14", "M_{J} [GeV]", {250., 400.}),
		      icut, goodbad_procs, plot_types).Tag("fat_goodbad");
      pm.Push<Hist1D>(Axis(5, 100., 850., "mj14", "M_{J} [GeV]", {250., 400.}),
		      icut, goodbad0_procs, plot_types).Tag("fat_goodbad0");
      pm.Push<Hist1D>(Axis(5, 100., 850., "mj14", "M_{J} [GeV]", {250., 400.}),
		      icut, proc_procs, plot_types).Tag("fat_proc");
      pm.Push<Hist1D>(Axis(5, 100., 850., "mj14", "M_{J} [GeV]", {250., 400.}),
		      icut, lep_procs, plot_types).Tag("fat_lep");

      // pm.Push<Hist1D>(Axis(10, 100., 850., "mj14", "M_{J} [GeV]", {250., 400.}),
      // 		      icut, lep_procs, plot_types).Tag("lep");
      // pm.Push<Hist1D>(Axis(12, 0., 420., "mt", "m_{T} [GeV]", {140.}),
      // 		      icut, all_procs, plot_types).Tag("mt");
      Nplots += 1;
    } // Loop over nobj bins
  } // Loop over met bins
  pm.min_print_ = true;
  pm.MakePlots(lumi);



  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Making "<<Nplots<<" plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
}

