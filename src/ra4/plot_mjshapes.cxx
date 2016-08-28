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
  //beware setting all to true, results in many...many...plots... 
  bool do_metbins = false;       
  bool do_met150 = true;        
  bool do_htopts = false;        // if false, does Jae's proposal only; if true, produces a pdf with comparisons
  bool compare_mj_nolep = true;  
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
  string folder_mc(bfolder+"/cms2r0/babymaker/babies/reclustered/2016_06_14/mc/merged_nl1st500met150nj5/");
  //speed up if met 150-200 not needed
  if (!do_met150) folder_mc = bfolder+"/cms2r0/babymaker/babies/reclustered/2016_06_14/mc/merged_nleps1met200nj5/";

  set<string> files_tt({folder_mc+"*_TTJets*Lept*.root", folder_mc+"*_TTJets_HT*.root"});
  set<string> files_wjets({folder_mc+"*_WJetsToLNu*.root"});
  set<string> files_st({folder_mc+"*_ST_*.root"});
  set<string> files_other({
    folder_mc+"*DYJetsToLL*.root", folder_mc+"*_QCD_HT*.root", 
    folder_mc+"*_ZJet*.root", folder_mc+"*_WWTo*.root", 
    folder_mc+"*ggZH_HToBB*.root", folder_mc+"*ttHJetTobb*.root",
    folder_mc+"*_TTGJets*.root", folder_mc+"*_TTTT_*.root", 
    folder_mc+"*_TTWJets*.root", folder_mc+"*_TTZTo*.root",
    folder_mc+"*_WH_HToBB*.root", folder_mc+"*_WZTo*.root",
    folder_mc+"*_ZH_HToBB*.root", folder_mc+"_ZZ_*.root"});
  
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
  procs_isr["highmt"] = vector<shared_ptr<Process> >({tt1l_lowmt, tt1l_highmt_2trul, 
                                                  tt1l_highmt_1trul_1trutauh, tt1l_highmt_1trul});
  procs_pie["highmt"] = vector<shared_ptr<Process> >({tt1l_highmt_2trul, 
                                                  tt1l_highmt_1trul_1trutauh, tt1l_highmt_1trul, nontt1l_highmt});
  // compare the tt@low-mT shape to the dilepton test pieces: tt 2l, non-tt 2l and lveto
  procs["dilep"] = vector<shared_ptr<Process> >({tt1l_lowmt, nontt2l, ttlveto, tt2l});
  procs_isr["dilep"] = vector<shared_ptr<Process> >({tt1l_lowmt, ttlveto, tt2l});
  procs_pie["dilep"] = vector<shared_ptr<Process> >({nontt2l, ttlveto, tt2l});
  // compare full background: low-mT, high-mT, 2l and lveto 
  procs["totbkg"] = vector<shared_ptr<Process> >({bkg1l_lowmt, bkg2l, bkglveto, bkg1l_highmt});
  procs_pie["totbkg"] = vector<shared_ptr<Process> >({bkg2l, bkglveto, bkg1l_highmt});

  vector<NamedFunc> htopt;
  if (do_htopts) {
    htopt.push_back(NamedFunc("ht"));
    htopt.push_back(NamedFunc("st", [](const Baby &b) -> NamedFunc::ScalarType{
      float st = b.ht();
      for (const auto &pt: *(b.leps_pt())) st += pt; 
      return st;
    }));
    htopt.push_back(NamedFunc("ht1l_stmin2l", [](const Baby &b) -> NamedFunc::ScalarType{
      float ht_proxy = b.ht();
      if (b.nleps()==2) ht_proxy =  b.ht()+b.leps_pt()->at(1);
      return ht_proxy;
    }));
    htopt.push_back(NamedFunc("ht1l_stave2l", [](const Baby &b) -> NamedFunc::ScalarType{
      float ht_proxy = b.ht();
      if (b.nleps()==2) ht_proxy =  b.ht()+(b.leps_pt()->at(0)+b.leps_pt()->at(1))/2.;
      return ht_proxy;
    }));
  }
  htopt.push_back(NamedFunc("ht1l_stmax2l", [](const Baby &b) -> NamedFunc::ScalarType{
    float ht_proxy = b.ht();
    if (b.nleps()==2) ht_proxy =  b.ht()+b.leps_pt()->at(0);
    return ht_proxy;
  }));

  vector<NamedFunc> metbins;
  if (do_met150) metbins.push_back(NamedFunc("met>150&&met<=200"));
  if (do_metbins) {
    metbins.push_back(NamedFunc("met>200&&met<=350"));
    metbins.push_back(NamedFunc("met>350&&met<=500"));
    metbins.push_back(NamedFunc("met>500"));
  } else {
    metbins.push_back(NamedFunc("met>200"));
  } 

  NamedFunc ave_toppt("ave_toppt",[](const Baby &b) -> NamedFunc::ScalarType{
    float toppt = 0;
    for (size_t imc(0); imc<b.mc_pt()->size(); imc++){
      if (fabs(b.mc_id()->at(imc))==6 && b.mc_status()->at(imc)==62) toppt += b.mc_pt()->at(imc);
    }
    return toppt/2.;
  });

  map<bool, vector<string> > nobjbins;
  // with leptons clustered
  nobjbins[true] = vector<string>();
  nobjbins[true].push_back("njets+nleps==6");
  nobjbins[true].push_back("njets+nleps>=7 && njets+nleps<=9");
  nobjbins[true].push_back("njets+nleps>=10");
  // without leptons clustered
  nobjbins[false] = vector<string>();
  nobjbins[false].push_back("njets==5");
  nobjbins[false].push_back("njets>=6 && njets<=8");
  nobjbins[false].push_back("njets>=9");
  
  set<string> pie_names, plot_names;
  int Nplots(0);
  PlotMaker pm;
  vector<bool> rcopts = {true, false}; // plot both 
  if (!Contains(folder_mc,"reclustered") || !compare_mj_nolep) rcopts = {true};
  for (auto mj_lep: {true, false}){
    string xtitle = mj_lep ? "M_{J} [GeV]" : "M_{J} (no lep) [GeV]";
    string var = (Contains(folder_mc,"reclustered") && mj_lep) ? "mj14_original" : "mj14";
    for (auto &iproc:procs){
      vector<TableRow> table_cuts;
      for (auto &imet: metbins) {
        for (auto &inobj: nobjbins[mj_lep]) {
          for (auto &iht: htopt){
            // no dilepton test for njets==5
            if(iproc.first == "dilep" && Contains(inobj,"==6")) continue; 
            vector<shared_ptr<Process> > procs_tmp = iproc.second;
            //similarly, remove 2l for njets==5 category
            if (iproc.first == "totbkg" && Contains(inobj,"==6")) 
              procs_tmp = vector<shared_ptr<Process> >({bkg1l_lowmt, bkglveto, bkg1l_highmt});
            // histograms
            NamedFunc icut = imet && inobj && iht>500;
            pm.Push<HistoStack>(HistoDef(iproc.first+"_"+iht.PlainName(), 
              10, 100., 850., var, xtitle, icut, "weight", {250.,400.}), procs_tmp, plot_types);
            plot_names.insert(pm.GetLast<HistoStack>()->definition_.Name()+"_OPT_"+plot_types[0].TypeString()+".pdf");

            if (procs_isr.find(iproc.first)!=procs_isr.end()) {
              pm.Push<HistoStack>(HistoDef(iproc.first+"_"+iht.PlainName(), 
                10, 0., 800., ave_toppt, "Ave. top p_{T} [GeV]", icut && var+">250", "weight"), procs_isr[iproc.first], plot_types);
              plot_names.insert(pm.GetLast<HistoStack>()->definition_.Name()+"_OPT_"+plot_types[0].TypeString()+".pdf");
              pm.Push<HistoStack>(HistoDef(iproc.first+"_"+iht.PlainName(), 
                10, 0., 800., "isr_tru_pt", "True ISR p_{T} [GeV]", icut && var+">250", "weight"), procs_isr[iproc.first], plot_types);
              plot_names.insert(pm.GetLast<HistoStack>()->definition_.Name()+"_OPT_"+plot_types[0].TypeString()+".pdf");
            }
            // pies
            table_cuts.push_back(TableRow("",icut));
            pie_names.insert("pie_"+iproc.first+"_"+icut.PlainName()+"_perc_lumi"+RoundNumber(lumi,0).Data()+".pdf");
            Nplots++;
          } // loop over ht options
        } // loop over nobj bins
      } // loop over met bins
      pm.Push<Table>(iproc.first,  table_cuts, procs_pie[iproc.first], false, true, true);
    } // loop over proc sets
  } // loop over cluster vs not-cluster leptons
  pm.min_print_ = true;
  pm.MakePlots(lumi);

  // slides comparing closure for MJ_lep and MJ_nolep
  if (Contains(folder_mc,"reclustered") && compare_mj_nolep) {
    SlideMaker sm_mjnolep("mjlep_vs_nolep.tex","1610");
    for (auto &iht: htopt){
      for (auto &imet: metbins) {
        for (size_t inobj(0); inobj<nobjbins[true].size(); inobj++) {
          vector<string> pnames;
          for (auto mj_lep: {true, false}) {
            for (auto &iproc: procs){
              string proc = iproc.first;
              if(proc == "dilep" && (Contains(nobjbins[mj_lep][inobj],"==6") || Contains(nobjbins[mj_lep][inobj],"==5"))) continue; 
              NamedFunc icut = imet && nobjbins[mj_lep][inobj] && iht>500;
              string var = mj_lep ? "mj14_original" : "mj14";
              string iname = proc+"_"+iht.PlainName()+"_VAR_"+var+"_CUT_"+icut.PlainName() 
                             +"_WGT_weight_OPT_"+plot_types[0].TypeString()+".pdf";
              if (plot_names.find(iname)!=plot_names.end()) pnames.push_back(iname);
            }
          }
          sm_mjnolep.AddSlide(pnames,4);
        }
      }
    }
    sm_mjnolep.Close();
  }
  
  // Make slides to show comparison between the different modified HT options
  if (do_htopts){
    for (auto mj_lep: {true, false}){
      string sname = mj_lep ? "plots_mj.tex" : "plots_mj_nolep.tex";
      SlideMaker sm_htopt(sname);
      string var = (Contains(folder_mc,"reclustered") && mj_lep) ? "mj14_original" : "mj14";
      for (auto &iproc: procs){
        string proc = iproc.first;
        for (auto &imet: metbins) {
          for (auto &inobj: nobjbins[mj_lep]) {
            vector<string> pnames(15,"");
            for (size_t iht(0); iht<htopt.size(); iht++){
              NamedFunc icut = imet && inobj && htopt[iht]>500;
              string iname = proc+"_"+htopt[iht].PlainName()+"_VAR_"+var + "_CUT_" + icut.PlainName() 
                             + "_WGT_weight_OPT_"+plot_types[0].TypeString()+".pdf";
              if (plot_names.find(iname)!=plot_names.end()){
                pnames[iht] = "pie_"+proc+"_"+icut.PlainName()+"_perc_lumi"+RoundNumber(lumi,0).Data()+".pdf";
                pnames[htopt.size()+iht] = iname;
              }
              icut = icut && var+">250";
              iname = proc+"_"+htopt[iht].PlainName()+"_VAR_isr_tru_pt_CUT_" + icut.PlainName() 
                             + "_WGT_weight_OPT_"+plot_types[0].TypeString()+".pdf";
              if (plot_names.find(iname)!=plot_names.end()) pnames[2*htopt.size()+iht] = iname;
            }
            sm_htopt.AddSlide(pnames,5);
          }
        }
      }
      sm_htopt.Close();
    }
  }

  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Making "<<Nplots<<" plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
}

