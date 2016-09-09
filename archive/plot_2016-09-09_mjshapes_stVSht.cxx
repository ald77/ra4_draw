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
#include "core/table.hpp"
#include "core/slide_maker.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  bool do_metbins = false;       
  bool do_met150 = true;  
  bool do_njbins = false;    
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

  Palette colors("txt/colors.txt", "default");

  //if folder name contains "reclustered", plots will also be made for MJ_nolep
  string folder_mc(bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_stdnj5/");
  string ntupletag = "*.root";
  if (!do_met150) ntupletag = "*metG200*.root";

  set<string> files_tt({folder_mc+"*_TTJets*Lept"+ntupletag, folder_mc+"*_TTJets_HT"+ntupletag});
  set<string> files_wjets({folder_mc+"*_WJetsToLNu"+ntupletag});
  set<string> files_st({folder_mc+"*_ST_"+ntupletag});
  set<string> files_other({
    folder_mc+"*DYJetsToLL"+ntupletag, 
    folder_mc+"*_QCD_HT*00_Tune"+ntupletag, folder_mc+"*_QCD_HT*Inf_Tune"+ntupletag,
    folder_mc+"*_ZJet"+ntupletag, folder_mc+"*_WWTo"+ntupletag, 
    folder_mc+"*ggZH_HToBB"+ntupletag, folder_mc+"*ttHJetTobb"+ntupletag,
    folder_mc+"*_TTGJets"+ntupletag, folder_mc+"*_TTTT_"+ntupletag, 
    folder_mc+"*_TTWJets"+ntupletag, folder_mc+"*_TTZTo"+ntupletag,
    folder_mc+"*_WH_HToBB"+ntupletag, folder_mc+"*_WZTo"+ntupletag,
    folder_mc+"*_ZH_HToBB"+ntupletag, folder_mc+"_ZZ_"+ntupletag});
  
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
  tt1l_lowmt->SetLineStyle(2);
  auto tt1l_highmt = Process::MakeShared<Baby_full>("t#bar{t}(1l), m_{T} > 140", 
    Process::Type::background, colors("tt_1l"), files_tt, baseline_1l && "ntruleps>=1 && mt>140");
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t}(2l)", 
    Process::Type::background, kGreen+2, files_tt, baseline_2l && "ntruleps>=1");
  auto ttlveto = Process::MakeShared<Baby_full>("t#bar{t}(lv), m_{T} > 140", 
    Process::Type::background, kOrange+1, files_tt, baseline_lveto && "ntruleps>=1");
  auto tt2lveto = Process::MakeShared<Baby_full>("t#bar{t}(2l or 1l+trk)", 
    Process::Type::background, kPink+2, files_tt, (baseline_lveto || baseline_2l) && "ntruleps>=1");

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
  // wjets1l_lowmt->SetLineStyle(2);
  auto st1l_lowmt = Process::MakeShared<Baby_full>("Single t (1l), m_{T}#leq140", 
    Process::Type::background, colors("single_t"), files_st, baseline_1l && "mt<=140");
  // st1l_lowmt->SetLineStyle(2);
  auto other1l_lowmt = Process::MakeShared<Baby_full>("Other (1l), m_{T}#leq140", 
    Process::Type::background, colors("other"), files_other, baseline_1l && "mt<=140");
  // other1l_lowmt->SetLineStyle(2);
  auto nontt2l = Process::MakeShared<Baby_full>("Non-t#bar{t} (2l)", 
    Process::Type::background, colors("other"), files_nontt, baseline_2l);
  auto nontt1l_highmt = Process::MakeShared<Baby_full>("Non-t#bar{t} (1l), m_{T}>140", 
    Process::Type::background, colors("other"), files_nontt, baseline_1l && "mt>140");

  map<string, vector<shared_ptr<Process> > > procs, procs_pie, procs_isr;
  // compare the shapes of the low-mT components: tt, st, w, other
  procs["lowmt"] = vector<shared_ptr<Process> >({tt1l_lowmt, other1l_lowmt, st1l_lowmt, wjets1l_lowmt});
  procs_pie["lowmt"] = vector<shared_ptr<Process> >({tt1l_lowmt, other1l_lowmt, st1l_lowmt, wjets1l_lowmt});
  // compare the tt@low-mT shape to the high-mT components: tru 2l, tru 1l+tauh, non-tt 
  procs["highmt"] = vector<shared_ptr<Process> >({tt1l_lowmt, tt1l_highmt_2trul, 
                                                  tt1l_highmt_1trul_1trutauh, tt1l_highmt_1trul, nontt1l_highmt});
  procs_pie["highmt"] = vector<shared_ptr<Process> >({tt1l_highmt_2trul, 
                                                  tt1l_highmt_1trul_1trutauh, tt1l_highmt_1trul, nontt1l_highmt});
  // compare the tt@low-mT shape to the dilepton test pieces: tt 2l, non-tt 2l and lveto
  procs["dilep"] = vector<shared_ptr<Process> >({tt1l_lowmt, nontt2l, ttlveto, tt2l});
  procs_pie["dilep"] = vector<shared_ptr<Process> >({nontt2l, ttlveto, tt2l});
  // compare full background: low-mT, high-mT, 2l and lveto 
  procs["totbkg"] = vector<shared_ptr<Process> >({bkg1l_lowmt, bkg2l, bkglveto, bkg1l_highmt});
  procs_pie["totbkg"] = vector<shared_ptr<Process> >({bkg2l, bkglveto, bkg1l_highmt});

  procs["2lonly"] = vector<shared_ptr<Process> >({tt1l_lowmt, tt1l_highmt, tt2l, tt2lveto});

  vector<NamedFunc> metbins;
  if (do_met150) metbins.push_back(NamedFunc("met>150&&met<=200"));
  if (do_metbins) {
    metbins.push_back(NamedFunc("met>200&&met<=350"));
    metbins.push_back(NamedFunc("met>350&&met<=500"));
    metbins.push_back(NamedFunc("met>500"));
  } else {
    metbins.push_back(NamedFunc("met>150"));
    metbins.push_back(NamedFunc("met>200"));
  } 

  vector<string> nobjbins;
  if (do_njbins) {
    nobjbins.push_back("njets+nleps==6");
    nobjbins.push_back("njets+nleps>=7 && njets+nleps<=9");
    nobjbins.push_back("njets+nleps>=10");
  } else {
    nobjbins.push_back("njets+nleps>=7");
  }

  PlotMaker pm;
  for (auto &iht: {"ht>500","st>500"}){
    for (auto &iproc:procs){
      vector<TableRow> table_cuts;
      for (auto &imet: metbins) {
        for (auto &inobj: nobjbins) {
          // no dilepton test for njets==5
          if(iproc.first == "dilep" && Contains(inobj,"==6")) continue; 
          vector<shared_ptr<Process> > procs_tmp = iproc.second;
          //similarly, remove 2l for njets==5 category
          if (iproc.first == "totbkg" && Contains(inobj,"==6")) 
            procs_tmp = vector<shared_ptr<Process> >({bkg1l_lowmt, bkglveto, bkg1l_highmt});
          // histograms
          NamedFunc icut = imet && inobj && iht;
          pm.Push<HistoStack>(HistoDef(iproc.first, 
            15, 100., 700., "mj14", "M_{J} [GeV]", icut, "weight", {250.,400.}), procs_tmp, plot_types);
        } // loop over nobj bins
      } // loop over met bins
    } // loop over proc sets
  }
  pm.min_print_ = true;
  pm.MakePlots(lumi, "mjshapes");

  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Making plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
}

