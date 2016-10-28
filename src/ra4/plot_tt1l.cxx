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


NamedFunc max_b_pt("max_b_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float maxPt=-999.;
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (abs(b.mc_id()->at(i))!=5) continue;
      if(b.mc_pt()->at(i) > maxPt) maxPt = b.mc_pt()->at(i);
    }
    return maxPt;
  });

NamedFunc max_t_pt("max_t_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float maxPt=-999.;
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (abs(b.mc_id()->at(i))!=6) continue;
      if(b.mc_pt()->at(i) > maxPt) maxPt = b.mc_pt()->at(i);
    }
    return maxPt;
  });


int main(){
  gErrorIgnoreLevel = 6000;

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  double lumi = 2.6;
  string bfolder("");
  string hostname(execute("echo $HOSTNAME"));
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string ntupletag=""; 
  string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_met100_stdnj5/");

  Palette colors("txt/colors.txt", "default");


  NamedFunc baseline = "stitch && nleps==1 && nveto==0 && st>500 && met>100 && njets>=5 && nbm>=1 && weight<1";

  set<string> ttfiles = {foldermc+"*_TTJets*Lept*"+ntupletag+"*.root", foldermc+"*_TTJets_HT*"+ntupletag+"*.root"};

  /// Study of ttbar 1l
  auto proc_tt1l_lomt = Process::MakeShared<Baby_full>("t#bar{t} 1l, m_{T}#leq140", Process::Type::background, 
						       colors("tt_1l"), ttfiles, 
						       baseline && "ntruleps<=1 && mt<=140");
  auto proc_ttltau = Process::MakeShared<Baby_full>("t#bar{t} l#tau_{h}, m_{T}>140", Process::Type::background, 
						    kBlue-6, ttfiles, 
						    baseline && "ntruels+ntrumus+ntrutausl==1 && ntrutaush==1 && mt>140");
  auto proc_tt2l = Process::MakeShared<Baby_full>("t#bar{t} 2l, m_{T}>140", Process::Type::background, 
						  colors("tt_2l"), ttfiles, 
						  baseline && "ntruels+ntrumus+ntrutausl==2 && mt>140");
  auto proc_tt1l_ghimt = Process::MakeShared<Baby_full>("t#bar{t} 1l, m_{T}>140, m_{T}^{tru}>140", Process::Type::background, 
						       kGreen-3, ttfiles, 
						       baseline && "ntruleps<=1 && mt>140&&mt_tru>140");
  auto proc_tt1l_bhimt = Process::MakeShared<Baby_full>("t#bar{t} 1l, m_{T}>140, m_{T}^{tru}#leq140", Process::Type::background, 
						       kRed-4, ttfiles, 
						       baseline && "ntruleps<=1 && mt>140&&mt_tru<140");

  

  vector<shared_ptr<Process> > tt_procs = {proc_tt1l_lomt, proc_tt2l, proc_ttltau, proc_tt1l_ghimt, proc_tt1l_bhimt};
  vector<shared_ptr<Process> > tt1l_procs = {proc_tt1l_lomt, proc_tt1l_ghimt, proc_tt1l_bhimt};
  
  PlotOpt log_shapes("txt/plot_styles.txt", "CMSPaper");
  log_shapes.Title(TitleType::info)
  .Bottom(BottomType::ratio)
  .YAxis(YAxisType::log)
  .Stack(StackType::shapes)
  .RatioMaximum(2.8);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  vector<PlotOpt> plot_types = {lin_shapes};


  PlotMaker pm;


  pm.Push<Hist1D>(Axis(20,100.,850.,"mj14","M_{J} [GeV]",{250., 400.}),"met>100", tt_procs, plot_types).Tag("tt1l");

  pm.Push<Hist1D>(Axis(35,0.,700.,max_b_pt,"Max b-quark p_{T} [GeV]"),"met>100", tt1l_procs, plot_types).Tag("tt1l");
  pm.Push<Hist1D>(Axis(30,-0.5,1.,"(met-met_tru)/met","(MET-MET^{tru})/MET"),"met>100", 
		  tt1l_procs, plot_types).Tag("tt1l");
  pm.Push<Hist1D>(Axis(60,-100.,200.,"met-met_tru","MET-MET^{tru} [GeV]"),"met>100", tt1l_procs, plot_types).Tag("tt1l");


  pm.Push<Hist1D>(Axis(30,300.,1500.,"ht","H_{T} [GeV]"),"met>100", tt1l_procs, plot_types).Tag("tt1l");

  pm.Push<Hist1D>(Axis(30,500.,1700.,"st","S_{T} [GeV]"),"met>100", tt1l_procs, plot_types).Tag("tt1l");
  pm.Push<Hist1D>(Axis(30,300.,1500.,"ht_tru","H_{T}^{tru} [GeV]"),"met>100", tt1l_procs, plot_types).Tag("tt1l");
  pm.Push<Hist1D>(Axis(48,-800.,400.,"ht-ht_tru","H_{T}-H_{T}^{tru} [GeV]"),"met>100", tt1l_procs, plot_types).Tag("tt1l");
  pm.Push<Hist1D>(Axis(40,-100.,300.,"mt-mt_tru","m_{T}-m_{T}^{tru} [GeV]"),"met>100", tt1l_procs, plot_types).Tag("tt1l");

  pm.Push<Hist1D>(Axis(30,-1.,5.,"(met-met_tru)/met_tru","(MET-MET^{tru})/MET^{tru}"),"met>100", 
		  tt1l_procs, plot_types).Tag("tt1l");
  pm.Push<Hist1D>(Axis(30,200.,800.,"met","MET [GeV]"),"met>200", tt1l_procs, plot_types).Tag("tt1l");
  pm.Push<Hist1D>(Axis(30,300.,1500.,"m_tt","m_{t#bar{t}} [GeV]"),"met>100", tt1l_procs, plot_types).Tag("tt1l");
  pm.Push<Hist1D>(Axis(6,-0.5,5.5,"nisr","N_{jets}^{ISR}"),"met>100", tt1l_procs, plot_types).Tag("tt1l");
  pm.Push<Hist1D>(Axis(40,0.,200.,"leps_pt[0]","Lepton p_{T} [GeV]"),"met>100", tt1l_procs, plot_types).Tag("tt1l");
  pm.Push<Hist1D>(Axis(50,0.,500.,"isr_tru_pt","t#bar{t} true p_{T} [GeV]"),"met>100", tt1l_procs, plot_types).Tag("tt1l");
  pm.Push<Hist1D>(Axis(35,0.,700.,max_t_pt,"Max t-quark p_{T} [GeV]"),"met>100", tt1l_procs, plot_types).Tag("tt1l");


  pm.min_print_ = true;
  pm.MakePlots(lumi);


  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Making "<<pm.Figures().size()<<" plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
}

