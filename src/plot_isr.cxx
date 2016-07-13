#include "plot_dilep_angles.hpp"

#include <cmath>
#include <algorithm>

#include "TError.h"
#include "TVector2.h"
#include "TString.h"

#include "plot_maker.hpp"
#include "plot_opt.hpp"
#include "palette.hpp"
#include "histo_stack.hpp"
#include "event_scan.hpp"
#include "utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace {
  const string isrtype = "ttisr";
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
NamedFunc::VectorType isrJetsPt(const Baby &b);
NamedFunc::VectorType isrJetsPt20(const Baby &b);
NamedFunc::ScalarType isrSystemPt(const Baby &b);

NamedFunc::ScalarType nisrMatch(const Baby &b);
NamedFunc::ScalarType nJetsReweightingTTJets(const Baby &b);


int main(){
  gErrorIgnoreLevel = 6000;
  bool sb = false;
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    sb = true;

  Palette colors("txt/colors.txt", "default");

  string folder_mc = "/cms2r0/babymaker/babies/2016_06_14/mc/merged_"+isrtype+"/";
  if (sb) folder_mc = "/net/cms2"+folder_mc;
  auto tt1l = Proc<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {folder_mc+"*_TTJets*SingleLept*.root"}, "ntruleps<=1");
  auto tt2l = Proc<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {folder_mc+"*_TTJets*DiLept*.root"}, "ntruleps>=2");
  auto single_t = Proc<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {folder_mc+"*_ST_*.root"});
  auto dyjets = Proc<Baby_full>("DY+jets", Process::Type::background, colors("wjets"),
    {folder_mc+"*DYJetsToLL_M-50_Tu*.root"});

  auto ttv = Proc<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {folder_mc+"*_TTWJets*.root", folder_mc+"*_TTZTo*.root", folder_mc+"*_TTGJets*.root"});
  auto other = Proc<Baby_full>("Other", Process::Type::background, colors("other"),
    {{folder_mc+"*_WJetsToLNu*.root"},
        folder_mc+"*_ZJet*.root", folder_mc+"*_WWTo*.root",
        folder_mc+"*ggZH_HToBB*.root", folder_mc+"*ttHJetTobb*.root",
        folder_mc+"*_TTTT_*.root",
        folder_mc+"*_WH_HToBB*.root", folder_mc+"*_WZTo*.root",
        folder_mc+"*_ZH_HToBB*.root", folder_mc+"_ZZ_*.root"});

  string folder_data = "/cms26r0/babymaker/babies/2016_06_26/data/skim_"+isrtype+"/";
  if (sb) folder_data = "/net/cms26" + folder_data;
  string lumi_label = RoundNumber(lumi,1).Data();
  auto data = Proc<Baby_full>("Data "+lumi_label+" fb^{-1}", Process::Type::data, kBlack,
    {folder_data+"*.root"},
    "pass && (trig[19]||trig[23])");

  auto t1tttt_nc = Proc<Baby_full>("T1tttt(1800,200)", Process::Type::signal, colors("t1tttt"),
    {folder_mc+"*SMS-T1tttt_mGluino-1800_mLSP-200*.root"});
  auto t1tttt_c = Proc<Baby_full>("T1tttt(1200,800)", Process::Type::signal, colors("t1tttt"),
    {folder_mc+"*SMS-T1tttt_mGluino-1200_mLSP-800*.root"});
  t1tttt_c->SetLineStyle(2);

  vector<shared_ptr<Process> > procs;
  if (isrtype=="zisr") procs = {data, dyjets, tt2l, tt1l, single_t, ttv, other};
  else procs = {data, tt2l, tt1l, dyjets, single_t, ttv, other};

  //// Processes for 1l ttbar closure
  string bfolder("");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder
  string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_standard/");
  string folderdata(bfolder+"/cms2r0/babymaker/babies/2016_06_26/data/merged_standard/");
  string baseline_1l("nleps==1 && ht>500 && met>200 && nbm>=2 && pass");

  auto proc_data = Proc<Baby_full>("Data", Process::Type::data, kBlack, 
    {folderdata+"/*.root"}, 
    baseline_1l+" && (trig[4]||trig[8]||trig[13]||trig[33])");

  auto proc_tt1l = Proc<Baby_full>("tt 1lep", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets*SingleLept*.root"},
    baseline_1l+" && ntruleps==1");
  auto proc_tt2l = Proc<Baby_full>("tt 2lep", Process::Type::background, colors("tt_2l"),
    {foldermc+"*_TTJets*DiLept*.root"},
    baseline_1l+" && ntruleps==2");
  auto proc_wjets = Proc<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {foldermc+"*_WJetsToLNu*.root"},
    baseline_1l+" && stitch");
  auto proc_singlet = Proc<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {foldermc+"*_ST_*.root"});
  auto proc_ttv = Proc<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {foldermc+"*_TTWJets*.root", foldermc+"*_TTZTo*.root"});
  auto proc_other = Proc<Baby_full>("Other", Process::Type::background, colors("other"),
    {foldermc+"*DYJetsToLL*.root",foldermc+"*QCD_HT*.root",
	foldermc+"*_ZJet*.root",foldermc+"*_ttHJetTobb*.root",
	foldermc+"*_TTGJets*.root",foldermc+"*_TTTT*.root",
	foldermc+"*_WH_HToBB*.root",foldermc+"*_ZH_HToBB*.root",
	foldermc+"*_WWTo*.root",foldermc+"*_WZ*.root",foldermc+"*_ZZ_*.root"},
    baseline_1l+" && stitch");

  vector<shared_ptr<Process> > procs_1l = {proc_data, proc_tt1l, proc_tt2l, proc_wjets, proc_singlet, proc_ttv, proc_other};


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

  PlotMaker pm;

  // tt_isr skim def:
  // "nvleps==2 && nleps>=1 && nbm==2 &&
  // max(Max$(mus_pt*(mus_tight&&mus_reliso<.1)),Max$(els_pt*(els_tight&&els_reliso<.1)))>30"
  NamedFunc nreliso_els("nreliso_els",nRelIsoEls);
  NamedFunc nreliso_mus("nreliso_mus",nRelIsoMus);
  NamedFunc baseline = "nleps==2" && nreliso_els+nreliso_mus>=1;

  NamedFunc max_reliso_elspt("max_reliso_elspt",maxRelIsoElsPt);
  NamedFunc max_reliso_muspt("max_reliso_muspt",maxRelIsoMusPt);
  NamedFunc isr_jetspt("isr_jetspt",isrJetsPt);
  NamedFunc nisrjets("nisrjets", [&](const Baby &b){
    return isrJetsPt(b).size();
  });
  NamedFunc isr_jetspt20("isr_jetspt20",isrJetsPt20);
  NamedFunc nisrjets20("nisrjets20", [&](const Baby &b){
    return isrJetsPt20(b).size();
  });
  NamedFunc isr_syspt("isr_syspt", isrSystemPt);

  // definitions for njets in slices of ISR pT
  const vector<double> isr_syspt_slices = {0, 50, 100, 150, 200, 300, 400};
  const vector<double> nisrjet_bins = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5};
  // definitions for ISR pT in slices of njets
  vector<double> nisrjet_slices = {0,1,2,3,4,5};
  const vector<double> isr_syspt_bins = {0, 50, 100, 150, 200, 300, 400, 600, 800};

  vector<double> ptbins = {30,40,50,75,100,150,200,300,400,600};
  vector<double> ptbins_zoom = {20,25,30,35,40,50,75,100,150,200};

  vector<NamedFunc> weight_opts;
  weight_opts.push_back(NamedFunc("default", [](const Baby &b) -> NamedFunc::ScalarType{return b.weight()/b.eff_trig()/b.w_toppt();}));
  weight_opts.push_back(NamedFunc("w_njets", nJetsReweightingTTJets));

  for (const auto &iweight: weight_opts){
    addSlices(pm, isr_syspt_slices, isr_syspt, nisrjet_bins, nisrjets, "ISR jet multiplicity", baseline, iweight, procs, {lin_lumi});
    addSlices(pm, nisrjet_slices, nisrjets, isr_syspt_bins, isr_syspt, "ISR p_{T} [GeV]", baseline, iweight, procs, {log_lumi});

    pm.Push<HistoStack>(HistoDef(isrtype, ptbins, isr_jetspt[0.], "Leading ISR jet p_{T} [GeV]", baseline && nisrjets>0., iweight), procs, plot_types);
    pm.Push<HistoStack>(HistoDef(isrtype, ptbins, isr_jetspt[1], "2^{nd} ISR jet p_{T} [GeV]", baseline && nisrjets>1, iweight), procs, plot_types);
    pm.Push<HistoStack>(HistoDef(isrtype, ptbins, isr_jetspt[2], "3^{rd} ISR jet p_{T} [GeV]", baseline && nisrjets>2, iweight), procs, plot_types);
    pm.Push<HistoStack>(HistoDef(isrtype, ptbins, isr_jetspt[3], "4^{th} ISR jet p_{T} [GeV]", baseline && nisrjets>3, iweight), procs, plot_types);

    pm.Push<HistoStack>(HistoDef(isrtype, ptbins_zoom, isr_jetspt20[0.], "Leading ISR jet p_{T} [GeV]", baseline && nisrjets20>0., iweight), procs, plot_types);
    pm.Push<HistoStack>(HistoDef(isrtype, ptbins_zoom, isr_jetspt20[1], "2^{nd} ISR jet p_{T} [GeV]", baseline && nisrjets20>1, iweight), procs, plot_types);
    pm.Push<HistoStack>(HistoDef(isrtype, ptbins_zoom, isr_jetspt20[2], "3^{rd} ISR jet p_{T} [GeV]", baseline && nisrjets20>2, iweight), procs, plot_types);
    pm.Push<HistoStack>(HistoDef(isrtype, ptbins_zoom, isr_jetspt20[3], "4^{th} ISR jet p_{T} [GeV]", baseline && nisrjets20>3, iweight), procs, plot_types);
    
    pm.Push<HistoStack>(HistoDef(isrtype, ptbins, max_reliso_elspt, "Leading electron p_{T} [GeV]", baseline && nreliso_els>0., iweight), procs, plot_types);
    pm.Push<HistoStack>(HistoDef(isrtype, ptbins, max_reliso_muspt, "Leading muon p_{T} [GeV]", baseline && nreliso_mus>0., iweight), procs, plot_types);

    pm.Push<HistoStack>(HistoDef(isrtype,20,0.,500., "met", "MET [GeV]", baseline, iweight), procs, plot_types);
    pm.Push<HistoStack>(HistoDef(isrtype,15,0.,1500., "ht", "H_{T} [GeV]", baseline, iweight), procs, plot_types);
    pm.Push<HistoStack>(HistoDef(isrtype,15,0.,1500., "mj14", "M_{J} [GeV]", baseline, iweight), procs, plot_types);

    //// 1l ttbar closure
    pm.Push<HistoStack>(HistoDef(13, -0.5, 12.5, "njets", "Number of jets", baseline_1l, iweight),  procs_1l, plot_types);
  }

  // MC study
  // vector<shared_ptr<Process> > proc_tt = {tt2l};
  if (isrtype=="ttisr") pm.Push<HistoStack>(HistoDef(isrtype,20,-2.,2., "(jetsys_nob_pt-isr_tru_pt)/isr_tru_pt", "(reco - true)/true ISR p_{T}", baseline), vector<shared_ptr<Process> >({tt2l}), vector<PlotOpt>({lin_shapes}));
  else pm.Push<HistoStack>(HistoDef(isrtype,20,-2.,2., "(jetsys_pt-isr_tru_pt)/isr_tru_pt", "(reco - true)/true ISR p_{T}", baseline), vector<shared_ptr<Process> >({dyjets}), vector<PlotOpt>({lin_shapes}));

  pm.Push<HistoStack>(HistoDef(isr_syspt_bins, "isr_tru_pt", "true ISR p_{T}", baseline), vector<shared_ptr<Process> >({tt2l,t1tttt_nc,t1tttt_c}), vector<PlotOpt>({lin_shapes}));

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

    string tag = isrtype+"_"+svar.PrettyName()+RoundNumber(slices[i],tag_digits).Data();
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

NamedFunc::VectorType isrJetsPt(const Baby &b){
  vector<double> isr_jetspt;
  for (size_t ijet(0); ijet<b.jets_pt()->size(); ijet++){
    if (!isGoodJet(b, ijet) || b.jets_pt()->at(ijet)<30.) continue;
    if (isrtype=="ttisr" && b.jets_csv()->at(ijet)>CSVMedium) continue;
    isr_jetspt.push_back(b.jets_pt()->at(ijet));
  }
  std::sort(isr_jetspt.begin(), isr_jetspt.end(), std::greater<double>());
  return isr_jetspt;
}

NamedFunc::VectorType isrJetsPt20(const Baby &b){
  vector<double> isr_jetspt;
  for (size_t ijet(0); ijet<b.jets_pt()->size(); ijet++){
    if (!isGoodJet(b, ijet) || b.jets_pt()->at(ijet)<20.) continue;
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

NamedFunc::ScalarType nisrMatch(const Baby &b){
  int Nisr=0;
  for (size_t ijet(0); ijet<b.jets_pt()->size(); ijet++){
    if(!isGoodJet(b, ijet) || b.jets_pt()->at(ijet)<30) continue;
    bool matched=false;
    for (size_t imc(0); imc<b.mc_pt()->size(); imc++){
      if(b.mc_status()->at(imc)!=23 || abs(b.mc_id()->at(imc))>5) continue;
      if(!(abs(b.mc_mom()->at(imc))==6 || abs(b.mc_mom()->at(imc))==23 || 
     abs(b.mc_mom()->at(imc))==24 || abs(b.mc_mom()->at(imc))==15)) continue; // In our ntuples where all taus come from W
        float dR = deltaR(b.jets_eta()->at(ijet), b.jets_phi()->at(ijet), b.mc_eta()->at(imc), b.mc_phi()->at(imc));
      if(dR<0.4){
        matched = true;
        break;
      }
    } // Loop over MC particles
    if(!matched) Nisr++;
  } // Loop over jets

  return Nisr;
}

NamedFunc::ScalarType nJetsReweightingTTJets(const Baby &b){
  if (b.ntrupv()<0) return 1.;

  int nisrjets(floor(b.njets()+0.5));
  double wgt = b.weight()/b.eff_trig()/b.w_toppt();
  if (nisrjets==0) return 1.10*wgt; //  +- 0.01
  else if (nisrjets==1) return 0.969*wgt; //  +- 0.02
  else if (nisrjets==2) return 0.870*wgt; //  +- 0.05
  else if (nisrjets==3) return 0.772*wgt; //  +- 0.09
  else if (nisrjets==4) return 0.712*wgt; //  +- 0.21
  else if (nisrjets>=5) return 0.635*wgt; //  +- 0.21
  else return wgt;
}


