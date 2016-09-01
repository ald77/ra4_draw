#include <cmath>
#include <stdio.h>
#include <chrono>

#include "TError.h"
#include "TVector2.h"

#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/histo_stack.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

int main(){
  gErrorIgnoreLevel = 6000;

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  double lumi = 2.6;
  string bfolder("");
  string hostname(execute("echo $HOSTNAME"));
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  Palette colors("txt/colors.txt", "default");

  //string ntupletag="_metG200";
  string ntupletag="";
  //string folder_mc(bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_stdnj5/");
  string folder_mc(bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_nl1st500met150nj5/");

  //speed up if met 150-200 not needed
  //if (!do_met150) ntupletag = "_metG150";

  set<string> files_tt({folder_mc+"*_TTJets*Lept*"+ntupletag+"*.root", 
	folder_mc+"*_TTJets_HT*"+ntupletag+"*.root"});
  set<string> files_wjets({folder_mc+"*_WJetsToLNu*"+ntupletag+"*.root"});
  set<string> files_st({folder_mc+"*_ST_*"+ntupletag+"*.root"});
  set<string> files_other({
    folder_mc+"*DYJetsToLL*"+ntupletag+"*.root", folder_mc+"*_QCD_HT*"+ntupletag+"*.root", 
    folder_mc+"*_ZJet*"+ntupletag+"*.root", folder_mc+"*_WWTo*"+ntupletag+"*.root", 
    folder_mc+"*ggZH_HToBB*"+ntupletag+"*.root", folder_mc+"*ttHJetTobb*"+ntupletag+"*.root",
    folder_mc+"*_TTGJets*"+ntupletag+"*.root", folder_mc+"*_TTTT_*"+ntupletag+"*.root", 
    folder_mc+"*_TTWJets*"+ntupletag+"*.root", folder_mc+"*_TTZTo*"+ntupletag+"*.root",
    folder_mc+"*_WH_HToBB*"+ntupletag+"*.root", folder_mc+"*_WZTo*"+ntupletag+"*.root",
    folder_mc+"*_ZH_HToBB*"+ntupletag+"*.root", folder_mc+"_ZZ_*"+ntupletag+"*.root"});
  
  set<string> files_nontt(files_other);
  files_nontt.insert(files_wjets.begin(), files_wjets.end());
  files_nontt.insert(files_st.begin(), files_st.end());

  set<string> files_all(files_nontt);
  files_all.insert(files_tt.begin(), files_tt.end());

  PlotOpt log_shapes("txt/plot_styles.txt", "CMSPaper");
  log_shapes.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::shapes)
    .RatioMaximum(2.4);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  vector<PlotOpt> plot_types = {lin_shapes};

  NamedFunc baseline_1l = "stitch && nleps==1 && nveto==0 && nbm>=1 && weight<1";
  NamedFunc baseline_2l = "stitch && nleps==2";
  NamedFunc baseline_lveto = "stitch && nleps==1 && nveto==1 && nbm>=1 && mt>140";

  auto tt1l_lowmt = Process::MakeShared<Baby_full>("t#bar{t}(1l), m_{T} #leq 140", 
    Process::Type::background, colors("tt_1l"), files_tt, baseline_1l && "ntruleps>=1 && mt<=140");
  auto tt1l_highmt = Process::MakeShared<Baby_full>("t#bar{t}(1l), m_{T} > 140", 
    Process::Type::background, colors("tt_1l"), files_tt, baseline_1l && "ntruleps>=1 && mt>140");
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t}(2l)", 
    Process::Type::background, kGreen+2, files_tt, baseline_2l && "ntruleps>=1");
  auto ttlveto = Process::MakeShared<Baby_full>("t#bar{t}(lv), m_{T} > 140", 
    Process::Type::background, kOrange+1, files_tt, baseline_lveto && "ntruleps>=1");

  auto tt1l_highmt_2trul = Process::MakeShared<Baby_full>("t#bar{t}(1l), m_{T} > 140, tru: 2l", 
    Process::Type::background, colors("tt_2l"), files_tt,
    baseline_1l && "ntruleps>=1 && ntruels+ntrumus+ntrutausl==2 && mt>140");
  auto tt1l_highmt_1trul_1trutauh = Process::MakeShared<Baby_full>("t#bar{t}(1l), m_{T} > 140, tru: l#tau_{h}", 
    Process::Type::background, kBlue-6, files_tt,
    baseline_1l && "ntruleps>=1 && ntruels+ntrumus+ntrutausl==1 && ntrutaush==1 && mt>140");
  auto tt1l_highmt_1trul = Process::MakeShared<Baby_full>("t#bar{t}(1l), m_{T} > 140, tru: 1l", 
    Process::Type::background, kTeal-6, files_tt,
    baseline_1l && "ntruleps>=1 && ntruels+ntrumus+ntrutausl==1 && ntrutaush==0 && mt>140");

  auto bkg1l_lowmt = Process::MakeShared<Baby_full>("Tot. bkgd (1l), m_{T}#leq140", 
    Process::Type::background, kGray+3, files_all, baseline_1l && "mt<=140");
  bkg1l_lowmt->SetLineStyle(2);
  auto bkg1l_highmt = Process::MakeShared<Baby_full>("Tot. bkgd (1l), m_{T}>140", 
    Process::Type::background, kGray+3, files_all, baseline_1l && "mt>140");
  auto bkg2l = Process::MakeShared<Baby_full>("Tot. bkgd (2l)", 
    Process::Type::background, kGreen-5, files_all, baseline_2l);
  auto bkglveto = Process::MakeShared<Baby_full>("Tot. bkgd (lv), m_{T} > 140", 
    Process::Type::background, kOrange+1, files_all, baseline_lveto);

  auto wjets1l_lowmt = Process::MakeShared<Baby_full>("W+jets (1l), m_{T}#leq140", 
    Process::Type::background, colors("wjets"), files_wjets, baseline_1l && "mt<=140");
  //wjets1l_lowmt->SetLineStyle(2);
  auto st1l_lowmt = Process::MakeShared<Baby_full>("Single t (1l), m_{T}#leq140", 
    Process::Type::background, colors("single_t"), files_st, baseline_1l && "mt<=140");
  //st1l_lowmt->SetLineStyle(2);
  auto other1l_lowmt = Process::MakeShared<Baby_full>("Other (1l), m_{T}#leq140", 
    Process::Type::background, colors("other"), files_other, baseline_1l && "mt<=140");
  //other1l_lowmt->SetLineStyle(2);
  auto nontt2l = Process::MakeShared<Baby_full>("Non-t#bar{t} (2l)", 
    Process::Type::background, colors("other"), files_nontt, baseline_2l);
  auto nontt1l_highmt = Process::MakeShared<Baby_full>("Non-t#bar{t} (1l), m_{T}>140", 
    Process::Type::background, colors("other"), files_nontt, baseline_1l && "mt>140");

  map<string, vector<shared_ptr<Process> > > procs;
  procs["totbkg"] = vector<shared_ptr<Process> >({bkg1l_lowmt, bkglveto, bkg2l, bkg1l_highmt});
  procs["dilep"] = vector<shared_ptr<Process> >({tt1l_lowmt, nontt2l, ttlveto, tt2l});
  procs["lowmt"] = vector<shared_ptr<Process> >({tt1l_lowmt, other1l_lowmt, st1l_lowmt, wjets1l_lowmt});
  procs["highmt"] = vector<shared_ptr<Process> >({tt1l_lowmt, tt1l_highmt_2trul, 
                                                  tt1l_highmt_1trul_1trutauh, tt1l_highmt_1trul, nontt1l_highmt});

  vector<NamedFunc> htopt;
  htopt.push_back(NamedFunc("ht"));
  //htopt.push_back(NamedFunc("st"));
  htopt.push_back(NamedFunc("st", [](const Baby &b) -> NamedFunc::ScalarType{
    float st = b.ht();
    for (const auto &pt: *(b.leps_pt())) st += pt; 
    return st;
  }));
  // htopt.push_back(NamedFunc("ht1l_stmin2l", [](const Baby &b) -> NamedFunc::ScalarType{
  //   float ht_proxy = b.ht();
  //   if (b.nleps()==2) ht_proxy =  b.ht()+b.leps_pt()->at(1);
  //   return ht_proxy;
  // }));
  // htopt.push_back(NamedFunc("ht1l_stmax2l", [](const Baby &b) -> NamedFunc::ScalarType{
  //   float ht_proxy = b.ht();
  //   if (b.nleps()==2) ht_proxy =  b.ht()+b.leps_pt()->at(0);
  //   return ht_proxy;
  // }));
  // htopt.push_back(NamedFunc("ht1l_stave2l", [](const Baby &b) -> NamedFunc::ScalarType{
  //   float ht_proxy = b.ht();
  //   if (b.nleps()==2) ht_proxy =  b.ht()+(b.leps_pt()->at(0)+b.leps_pt()->at(1))/2.;
  //   return ht_proxy;
  // }));

  int Nplots=0;
  PlotMaker pm;
  for (auto mj_lep: {true}){
    vector<string> nobj_cuts;
    if (mj_lep){
      nobj_cuts.push_back("met>200&&njets+nleps==6");
      nobj_cuts.push_back("met>200&&njets+nleps>=7 && njets+nleps<=9");
      nobj_cuts.push_back("met>200&&njets+nleps>=10");
      nobj_cuts.push_back("met>150&&met<=200&&njets+nleps==6");
      nobj_cuts.push_back("met>150&&met<=200&&njets+nleps>=7 && njets+nleps<=9");
      nobj_cuts.push_back("met>150&&met<=200&&njets+nleps>=10");
    } else {
      nobj_cuts.push_back("njets>=6 && njets<=8");
      nobj_cuts.push_back("njets>=9");
    }

    string var("mj14"), xtitle("M_{J} [GeV]");
    //if (!mj_lep) {var = "mj14"; xtitle = "M_{J} (no lep) [GeV]";}
    for (auto &iht: htopt){
      for (auto &iproc:procs){
        for (auto inobj: nobj_cuts) {
	  if(iproc.first == "dilep" && Contains(inobj,"==6")) continue; 
          pm.Push<HistoStack>(HistoDef(iproc.first+"_"+iht.PlainName(), 
            10, 100., 850., var, xtitle, inobj && iht>500, "weight", {250.,400.}), iproc.second, plot_types);
	  Nplots++;
        }
      }
    }
  }
  pm.min_print_ = true;
  pm.MakePlots(lumi);

  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Making "<<Nplots<<" plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
}
