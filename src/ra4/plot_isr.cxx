#include <cmath>
#include <algorithm>

#include "TError.h"
#include "TVector2.h"
#include "TString.h"

#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/histo_stack.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace {
  const string isrtype = "ttisr";
  bool do_tt1l = false;
  double lumi = 4.34;
  bool single_thread = false;

  double CSVMedium = 0.800;
}

void addSlices(PlotMaker &pm, const vector<double> slices, NamedFunc svar,
               const vector<double> xbins, NamedFunc xvar, string xlabel,
               const NamedFunc &baseline, const NamedFunc &weight,
               const vector<shared_ptr<Process> > &proc,
               const vector<PlotOpt> &plot_types, int tag_digits=0);

bool isRelIsoEl(const Baby &b, size_t iel);
NamedFunc::ScalarType nRelIsoEls(const Baby &b);
NamedFunc::ScalarType maxRelIsoElsPt(const Baby &b);

bool isRelIsoMu(const Baby &b, size_t imu);
NamedFunc::ScalarType nRelIsoMus(const Baby &b);
NamedFunc::ScalarType maxRelIsoMusPt(const Baby &b);

bool isGoodJet(const Baby &b, size_t ijet);
NamedFunc::VectorType isrJetsPt(const Baby &b, float ptThresh=30.);
NamedFunc::ScalarType isrSystemPt(const Baby &b);

NamedFunc::ScalarType nJetsWeights_ttisr(const Baby &b, bool use_baby_nisr);
NamedFunc::ScalarType nJetsWeights_visr(const Baby &b);

int main(){
  gErrorIgnoreLevel = 6000;

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  Palette colors("txt/colors.txt", "default");

  //// Processes for ISR skims
  string dir_mc_isr = bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_"+isrtype+"/";
  auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {dir_mc_isr+"*_TTJets*SingleLept*.root"}, "ntruleps<=1");
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {dir_mc_isr+"*_TTJets*DiLept*.root"}, "ntruleps>=2");
  auto single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {dir_mc_isr+"*_ST_*.root"});
  auto dyjets = Process::MakeShared<Baby_full>("DY+jets", Process::Type::background, kOrange-3,
    //{dir_mc_isr+"*DYJetsToLL_M-50_Tu*.root"});
    {dir_mc_isr+"*DYJetsToLL_M-50_*.root"},"stitch"); // Inclusive + HT-binned DY
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {dir_mc_isr+"*_TTWJets*.root", dir_mc_isr+"*_TTZTo*.root", dir_mc_isr+"*_TTGJets*.root"});
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {dir_mc_isr+"*_WJetsToLNu*.root",
        // dir_mc_isr+"*QCD_HT*.root", // Has some ugly high-weight events
        dir_mc_isr+"*_ZJet*.root", dir_mc_isr+"*_WWTo*.root",
        dir_mc_isr+"*ggZH_HToBB*.root", dir_mc_isr+"*ttHJetTobb*.root",
        dir_mc_isr+"*_TTTT_*.root",
        dir_mc_isr+"*_WH_HToBB*.root", dir_mc_isr+"*_WZTo*.root",
        dir_mc_isr+"*_ZH_HToBB*.root", dir_mc_isr+"_ZZ_*.root"});

  //// For "honest" closure
  string dir_mc_isr_alg = bfolder+"/cms2r0/babymaker/babies/2016_07_18/mc/merged_"+isrtype+"/";
  auto tt1l_alg = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {dir_mc_isr_alg+"*_TTJets*SingleLept*.root"}, "ntruleps<=1");
  auto tt2l_alg = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {dir_mc_isr_alg+"*_TTJets*DiLept*.root"}, "ntruleps>=2");

  //// Only used for W+jets
  auto wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {dir_mc_isr+"*_WJetsToLNu*.root"},
    "stitch");
  auto other_w = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {dir_mc_isr+"*DYJetsToLL_M-50_Tu*.root",dir_mc_isr+"*QCD_HT*.root",
        dir_mc_isr+"*_ZJet*.root", dir_mc_isr+"*_WWTo*.root",
        dir_mc_isr+"*ggZH_HToBB*.root", dir_mc_isr+"*ttHJetTobb*.root",
        dir_mc_isr+"*_TTTT_*.root",
        dir_mc_isr+"*_WH_HToBB*.root", dir_mc_isr+"*_WZTo*.root",
        dir_mc_isr+"*_ZH_HToBB*.root", dir_mc_isr+"_ZZ_*.root"});

  string dir_data_isr = bfolder+"/cms2r0/babymaker/babies/2016_06_26/data/skim_"+isrtype+"/";
  string lumi_label = RoundNumber(lumi,1).Data();
  auto data = Process::MakeShared<Baby_full>("Data "+lumi_label+" fb^{-1}", Process::Type::data, kBlack,
    {dir_data_isr+"*.root"},
    "pass && (trig[19]||trig[23])");

  vector<shared_ptr<Process> > procs;
  if (isrtype=="zisr") procs = {data, dyjets, tt2l, tt1l, single_t, ttv, other};
  else if (isrtype=="ttisr") procs = {data, tt2l, tt1l, dyjets, single_t, ttv, other};
  else if (isrtype=="wisr") procs = {data, wjets, tt1l, tt2l, single_t, ttv, other_w};
  else {cout<<isrtype<<" not supported, exiting"<<endl<<endl; return 0;}

  vector<shared_ptr<Process> > procs_alg = {data, tt2l_alg, tt1l_alg, dyjets, single_t, ttv, other};

  //// Processes for 1l ttbar closure
  string dir_mc_std(bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_standard/");
  string dir_data_std(bfolder+"/cms2r0/babymaker/babies/2016_06_26/data/merged_standard/");

  auto std_data = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {dir_data_std+"/*.root"},
    "(trig[4]||trig[8]||trig[13]||trig[33])");

  auto std_tt1l = Process::MakeShared<Baby_full>("tt 1lep", Process::Type::background, colors("tt_1l"),
    {dir_mc_std+"*_TTJets*SingleLept*.root"},
    "ntruleps==1");
  auto std_tt2l = Process::MakeShared<Baby_full>("tt 2lep", Process::Type::background, colors("tt_2l"),
    {dir_mc_std+"*_TTJets*DiLept*.root"},
    "ntruleps==2");
  auto std_wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {dir_mc_std+"*_WJetsToLNu*.root"},
    "stitch");
  auto std_singlet = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {dir_mc_std+"*_ST_*.root"});
  auto std_ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {dir_mc_std+"*_TTWJets*.root", dir_mc_std+"*_TTZTo*.root"});
  auto std_other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {dir_mc_std+"*DYJetsToLL*.root",dir_mc_std+"*QCD_HT*.root",
        dir_mc_std+"*_ZJet*.root",dir_mc_std+"*_ttHJetTobb*.root",
        dir_mc_std+"*_TTGJets*.root",dir_mc_std+"*_TTTT*.root",
        dir_mc_std+"*_WH_HToBB*.root",dir_mc_std+"*_ZH_HToBB*.root",
        dir_mc_std+"*_WWTo*.root",dir_mc_std+"*_WZ*.root",dir_mc_std+"*_ZZ_*.root"},
    "stitch");

  string dir_mc_std_alg = bfolder+"/cms2r0/babymaker/babies/2016_07_18/mc/merged_standard/";
  auto std_tt1l_alg = Process::MakeShared<Baby_full>("tt 1lep", Process::Type::background, colors("tt_1l"),
    {dir_mc_std_alg+"*_TTJets*SingleLept*.root"},
    "ntruleps==1");
  auto std_tt2l_alg = Process::MakeShared<Baby_full>("tt 2lep", Process::Type::background, colors("tt_2l"),
    {dir_mc_std_alg+"*_TTJets*DiLept*.root"},
    "ntruleps==2");

  vector<shared_ptr<Process> > procs_1l = {std_data, std_tt1l, std_tt2l, std_wjets, std_singlet, std_ttv, std_other};
  vector<shared_ptr<Process> > procs_1l_alg = {std_data, std_tt1l_alg, std_tt2l_alg, std_wjets, std_singlet, std_ttv, std_other};

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .Bottom(BottomType::off)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  vector<PlotOpt> plot_types = {log_lumi, lin_lumi};
  vector<PlotOpt> plot_vals = {lin_lumi().PrintVals(true)};
  PlotMaker pm;

  // tt_isr skim def:
  // "nvleps==2 && nleps>=1 && nbm==2 &&
  // max(Max$(mus_pt*(mus_tight&&mus_reliso<.1)),Max$(els_pt*(els_tight&&els_reliso<.1)))>30"
  NamedFunc nreliso_els("nreliso_els",nRelIsoEls);
  NamedFunc nreliso_mus("nreliso_mus",nRelIsoMus);
  NamedFunc baseline = "nleps==2" && nreliso_els+nreliso_mus>=1;
  NamedFunc baseline_w = nreliso_els+nreliso_mus==1 && "ht>200&&met>100&&nbl==0";
  NamedFunc baseline_1l = "nleps==1 && ht>500 && met>200 && nbm>=1 && pass";
  if(isrtype=="wisr") baseline = baseline_w;

  NamedFunc max_reliso_elspt("max_reliso_elspt",maxRelIsoElsPt);
  NamedFunc max_reliso_muspt("max_reliso_muspt",maxRelIsoMusPt);
  NamedFunc isr_jetspt("isr_jetspt",[&](const Baby &b){
      return isrJetsPt(b, 30.);
    });
  NamedFunc nisrjets("nisrjets", [&](const Baby &b){
      return isrJetsPt(b, 30.).size();
    });
  NamedFunc isr_ht("isr_ht", [&](const Baby &b){
      vector<double> jets_pt = isrJetsPt(b, 30.);
      double ht = 0;
      for (auto &jpt: jets_pt) ht += jpt;
      return ht;
    });
  NamedFunc isr_jetspt20("isr_jetspt20",[&](const Baby &b){
      return isrJetsPt(b, 20);
    });
  NamedFunc nisrjets20("nisrjets20", [&](const Baby &b){
      return isrJetsPt(b, 20).size();
    });
  NamedFunc nisrjets50("nisrjets50", [&](const Baby &b){
      return isrJetsPt(b, 50).size();
    });
  NamedFunc nisrjets75("nisrjets75", [&](const Baby &b){
      return isrJetsPt(b, 75).size();
    });
  NamedFunc isr_syspt("isr_syspt", isrSystemPt);

  // definitions for njets in slices of ISR pT
  const vector<double> isr_syspt_slices = {0, 100, 200, 300};
  const vector<double> nisrjet_bins = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5};
  vector<double> nisrjet_bins_vals = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5};
  if(isrtype!="ttisr") nisrjet_bins_vals = nisrjet_bins;

  // definitions for ISR pT in slices of njets
  vector<double> nisrjet_slices = {0,1,2,3,4,5};
  const vector<double> isr_syspt_bins = {0, 50, 100, 150, 200, 300, 400, 600, 800};

  const vector<double> isr_ht_slices = {0,100, 200, 300, 400, 500, 2000};
  const vector<double> isr_ht_bins = {0,100, 200, 300};

  vector<double> ptbins = {30,40,50,75,100,150,200,300,400,600};
  vector<double> ptbins_zoom = {20,25,30,35,40,50,75,100,150,200};

  vector<NamedFunc> weight_opts;
  weight_opts.push_back(NamedFunc("default", [](const Baby &b) -> NamedFunc::ScalarType{return b.weight()/b.eff_trig()/b.w_toppt();}));
  if (isrtype=="ttisr"){
    //reweight TTJets only according to nisr = b.njets()-2 --> perfect closure
    weight_opts.push_back(NamedFunc("w_ttisr_njets", [](const Baby &b) -> NamedFunc::ScalarType{return nJetsWeights_ttisr(b, false);}));
    //reweight TTJets only according to nisr = b.nisr() --> honest closure
    weight_opts.push_back(NamedFunc("w_ttisr_alg", [](const Baby &b) -> NamedFunc::ScalarType{return nJetsWeights_ttisr(b, true);}));
  } else if (isrtype=="zisr"){
    weight_opts.push_back(NamedFunc("w_visr", nJetsWeights_visr));
  }
  for (const auto &iweight: weight_opts){
    vector<shared_ptr<Process> > *iprocs = &procs;
    if(CodeToPlainText(iweight.Name())=="w_ttisr_alg") iprocs = &procs_alg;
    if(do_tt1l) {//// 1l ttbar closure
      if(CodeToPlainText(iweight.Name())=="w_ttisr_alg")
        pm.Push<HistoStack>(HistoDef("tt1l", 13, -0.5, 12.5, "njets", "Number of jets", baseline_1l, iweight), procs_1l_alg, plot_types);
      else
        pm.Push<HistoStack>(HistoDef("tt1l", 13, -0.5, 12.5, "njets", "Number of jets", baseline_1l, iweight), procs_1l, plot_types);
    } else {
      addSlices(pm, isr_syspt_slices, isr_syspt, nisrjet_bins, nisrjets, "ISR jet multiplicity", baseline, iweight, *iprocs, plot_types);
      addSlices(pm, nisrjet_slices, nisrjets, isr_syspt_bins, isr_syspt, "ISR p_{T} [GeV]", baseline, iweight, *iprocs, {log_lumi});
      addSlices(pm, isr_ht_slices, isr_ht, isr_syspt_bins, isr_syspt, "ISR p_{T} [GeV]", baseline, iweight, *iprocs, {log_lumi});
      addSlices(pm, isr_ht_slices, isr_ht, nisrjet_bins, nisrjets, "ISR jet multiplicity", baseline, iweight, *iprocs, {log_lumi});
      addSlices(pm, nisrjet_slices, nisrjets, isr_syspt_bins, isr_ht, "ISR H_{T} [GeV]", baseline, iweight, *iprocs, {log_lumi});

      pm.Push<HistoStack>(HistoDef(isrtype, isr_syspt_bins, isr_syspt, "ISR p_{T} [GeV]", baseline && "ht>300", iweight), *iprocs, vector<PlotOpt>({log_lumi}));
      pm.Push<HistoStack>(HistoDef(isrtype, ptbins, isr_jetspt[0.], "Leading ISR jet p_{T} [GeV]", baseline && nisrjets>0., iweight), *iprocs, plot_types);
      pm.Push<HistoStack>(HistoDef(isrtype, ptbins, isr_jetspt[1], "2^{nd} ISR jet p_{T} [GeV]", baseline && nisrjets>1, iweight), *iprocs, plot_types);
      pm.Push<HistoStack>(HistoDef(isrtype, ptbins, isr_jetspt[2], "3^{rd} ISR jet p_{T} [GeV]", baseline && nisrjets>2, iweight), *iprocs, plot_types);
      pm.Push<HistoStack>(HistoDef(isrtype, ptbins, isr_jetspt[3], "4^{th} ISR jet p_{T} [GeV]", baseline && nisrjets>3, iweight), *iprocs, plot_types);

      pm.Push<HistoStack>(HistoDef(isrtype, ptbins_zoom, isr_jetspt20[0.], "Leading ISR jet p_{T} [GeV]", baseline && nisrjets20>0., iweight), *iprocs, plot_types);
      pm.Push<HistoStack>(HistoDef(isrtype, ptbins_zoom, isr_jetspt20[1], "2^{nd} ISR jet p_{T} [GeV]", baseline && nisrjets20>1, iweight), *iprocs, plot_types);
      pm.Push<HistoStack>(HistoDef(isrtype, ptbins_zoom, isr_jetspt20[2], "3^{rd} ISR jet p_{T} [GeV]", baseline && nisrjets20>2, iweight), *iprocs, plot_types);
      pm.Push<HistoStack>(HistoDef(isrtype, ptbins_zoom, isr_jetspt20[3], "4^{th} ISR jet p_{T} [GeV]", baseline && nisrjets20>3, iweight), *iprocs, plot_types);

      pm.Push<HistoStack>(HistoDef(isrtype, ptbins, max_reliso_elspt, "Leading electron p_{T} [GeV]", baseline && nreliso_els>0., iweight), *iprocs, plot_types);
      pm.Push<HistoStack>(HistoDef(isrtype, ptbins, max_reliso_muspt, "Leading muon p_{T} [GeV]", baseline && nreliso_mus>0., iweight), *iprocs, plot_types);

      pm.Push<HistoStack>(HistoDef(isrtype,20,0.,500., "met", "MET [GeV]", baseline, iweight), *iprocs, plot_types);
      pm.Push<HistoStack>(HistoDef(isrtype,15,0.,1500., "ht", "H_{T} [GeV]", baseline, iweight), *iprocs, plot_types);
      pm.Push<HistoStack>(HistoDef(isrtype,15,0.,1500., "mj14", "M_{J} [GeV]", baseline, iweight), *iprocs, plot_types);
      pm.Push<HistoStack>(HistoDef(isrtype+"_vals",nisrjet_bins_vals, nisrjets, "ISR jet multiplicity", baseline, iweight), *iprocs, plot_vals);
      pm.Push<HistoStack>(HistoDef(isrtype,nisrjet_bins, nisrjets50, "Number of 50 GeV ISR jets", baseline, iweight), *iprocs, plot_types);
      pm.Push<HistoStack>(HistoDef(isrtype,nisrjet_bins, nisrjets75, "Number of 75 GeV ISR jets", baseline, iweight), *iprocs, plot_types);
      if(isrtype=="zisr"){
        pm.Push<HistoStack>(HistoDef(isrtype,nisrjet_bins, nisrjets, "ISR jet multiplicity", baseline && "ht>200", iweight), *iprocs, plot_types);
        pm.Push<HistoStack>(HistoDef(isrtype,isr_syspt_bins, isr_syspt, "ISR p_{T} [GeV]", baseline && "ht>200", iweight), *iprocs, plot_types);
      }
    } // if not wjets_tt1l
  } // Loop over weights

  if (single_thread) pm.multithreaded_ = false;
  pm.MakePlots(lumi);
}

void addSlices(PlotMaker &pm, const vector<double> slices, NamedFunc svar,
               const vector<double> xbins, NamedFunc xvar, string xlabel,
               const NamedFunc &baseline, const NamedFunc &weight,
               const vector<shared_ptr<Process> > &proc,
               const vector<PlotOpt> &plot_types, int tag_digits){

  //add the inclusive version first
  pm.Push<HistoStack>(HistoDef(isrtype+"_incl", xbins, xvar, xlabel, baseline, weight), proc, plot_types);
  for(unsigned i(0); i<slices.size(); i++){
    NamedFunc cut = baseline && svar>=slices[i];
    if (i<(slices.size()-1)) cut = cut && svar<slices[i+1];

    string tag = CodeToPlainText(isrtype+"_"+svar.Name()+RoundNumber(slices[i],tag_digits).Data());
    pm.Push<HistoStack>(HistoDef(tag, xbins, xvar, xlabel, cut, weight), proc, plot_types);
  }
}

bool isRelIsoEl(const Baby &b, size_t iel){
  return iel<b.els_pt()->size()
      && b.els_pt()->at(iel)>30.
      && fabs(b.els_sceta()->at(iel))<2.
      && b.els_tight()->at(iel)
      && b.els_reliso()->at(iel) < 0.1;
}

NamedFunc::ScalarType nRelIsoEls(const Baby &b){
  int nels = 0;
  for (size_t iel(0); iel<b.els_pt()->size(); iel++){
    if (isRelIsoEl(b,iel)) nels++;
  }
  return nels;
}

NamedFunc::ScalarType maxRelIsoElsPt(const Baby &b){
  double max_pt = 0;
  for (size_t iel(0); iel<b.els_pt()->size(); iel++){
    if (isRelIsoEl(b,iel) && b.els_pt()->at(iel)>max_pt) max_pt = b.els_pt()->at(iel);
  }
  return max_pt;
}

bool isRelIsoMu(const Baby &b, size_t imu){
  return imu<b.mus_pt()->size()
      && b.mus_pt()->at(imu)>30.
      && fabs(b.mus_eta()->at(imu))<2.
      && b.mus_tight()->at(imu)
      && b.mus_reliso()->at(imu) < 0.1;
}

NamedFunc::ScalarType nRelIsoMus(const Baby &b){
  int nmus = 0;
  for (size_t imu(0); imu<b.mus_pt()->size(); imu++){
    if (isRelIsoMu(b,imu)) nmus++;
  }
  return nmus;
}

NamedFunc::ScalarType maxRelIsoMusPt(const Baby &b){
  double max_pt = 0;
  for (size_t imu(0); imu<b.mus_pt()->size(); imu++){
    if (isRelIsoMu(b,imu) && b.mus_pt()->at(imu)>max_pt) max_pt = b.mus_pt()->at(imu);
  }
  return max_pt;
}

bool isGoodJet(const Baby &b, size_t ijet){
  return ijet<b.jets_pt()->size()
      && fabs(b.jets_eta()->at(ijet))<2.4
      && !b.jets_islep()->at(ijet);
}

NamedFunc::VectorType isrJetsPt(const Baby &b, float ptThresh){
  vector<double> isr_jetspt;
  for (size_t ijet(0); ijet<b.jets_pt()->size(); ijet++){
    if (!isGoodJet(b, ijet) || b.jets_pt()->at(ijet)<ptThresh) continue;
    if (isrtype=="ttisr" && b.jets_csv()->at(ijet)>CSVMedium) continue;
    isr_jetspt.push_back(b.jets_pt()->at(ijet));
  }
  std::sort(isr_jetspt.begin(), isr_jetspt.end(), std::greater<double>());
  return isr_jetspt;
}

NamedFunc::ScalarType isrSystemPt(const Baby &b){
    if (isrtype=="ttisr") return b.jetsys_nob_pt();
    else return b.jetsys_pt();
}

NamedFunc::ScalarType nJetsWeights_ttisr(const Baby &b, bool use_baby_nisr){
  if (b.ntrupv()<0) return 1.; // Do not reweight Data

  float wgt = b.weight()/b.eff_trig()/b.w_toppt();

  int nisrjets = b.njets();
  if (b.SampleType()==20) {
    if (use_baby_nisr) nisrjets = b.nisr();
    else nisrjets = b.njets() - 2;
  }
  // weights derived in TTJets and applied using the nisr calculation algorithm
  if      (nisrjets==0) return 1.099*wgt; //  +- 0.012
  else if (nisrjets==1) return 0.969*wgt; //  +- 0.014
  else if (nisrjets==2) return 0.870*wgt; //  +- 0.020
  else if (nisrjets==3) return 0.772*wgt; //  +- 0.031
  else if (nisrjets==4) return 0.712*wgt; //  +- 0.051
  else if (nisrjets==5) return 0.661*wgt; //  +- 0.088
  else if (nisrjets>=6) return 0.566*wgt; //  +- 0.133
  else return wgt;
}

NamedFunc::ScalarType nJetsWeights_visr(const Baby &b){
  if (b.ntrupv()<0) return 1.; // Do not reweight Data

  float wgt = b.weight()/b.eff_trig()/b.w_toppt();
  if(b.SampleType()<30 && b.SampleType()>=60) return wgt;

  int nisrjets(b.njets());
  // weights derived in DY+jets
  if      (nisrjets==0) return 0.981*wgt; //  +- 0.001
  else if (nisrjets==1) return 1.071*wgt; //  +- 0.001
  else if (nisrjets==2) return 1.169*wgt; //  +- 0.003
  else if (nisrjets==3) return 1.157*wgt; //  +- 0.007
  else if (nisrjets==4) return 1.014*wgt; //  +- 0.013
  else if (nisrjets==5) return 0.920*wgt; //  +- 0.025
  else if (nisrjets==6) return 0.867*wgt; //  +- 0.048
  else if (nisrjets>=7) return 0.935*wgt; //  +- 0.088
  else return wgt;
}
