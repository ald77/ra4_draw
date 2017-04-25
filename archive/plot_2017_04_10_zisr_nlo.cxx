#include <cmath>
#include <algorithm>
#include <getopt.h>

#include "TError.h"
#include "TVector2.h"
#include "TString.h"

#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist1d.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"
#include "hig/hig_functions.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace {
  string isrtype = "zcand";

  bool single_thread = false;
  bool quick = false; 
  bool rewgt = false; 

  float lumi = 35.9;
}

void addSlices(PlotMaker &pm, const vector<double> slices, string svar,
               const vector<double> xbins, string xvar, string xlabel,
               const NamedFunc &baseline, const NamedFunc &weight,
               const vector<shared_ptr<Process> > &proc,
               const vector<PlotOpt> &plot_types, int tag_digits=0);

bool isGoodJet(const Baby &b, size_t ijet);
NamedFunc::VectorType isrJetsPt(const Baby &b, float ptThresh=30.);

NamedFunc::ScalarType ZpTWeights(const Baby &b);

void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  Palette colors("txt/colors.txt", "default");

  //// Processes for ISR skims
  string dir_mc_isr = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_minisrmc_"+isrtype+"/";
  string dir_mc_isr_nlo = bfolder+"/cms2r0/babymaker/babies/2017_04_01/mc/merged_minisrmc_"+isrtype+"/";
  auto ttbar = Process::MakeShared<Baby_full>("t#bar{t}", Process::Type::background, colors("tt_1l"),
    {dir_mc_isr+"*_TTJets_SingleLeptFromT_Tune*.root", dir_mc_isr+"*_TTJets_SingleLeptFromTbar_Tune*.root",
     dir_mc_isr+"*_TTJets_DiLept_Tune*.root"}, "1");
  auto single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {dir_mc_isr+"*_ST_*.root"});
  auto dyjets = Process::MakeShared<Baby_full>("DY+jets (LO)", Process::Type::background, kOrange+1,
     {dir_mc_isr+"*DYJetsToLL_M-50_*.root"},"stitch"); // Inclusive + HT-binned DY
  // auto dyjets = Process::MakeShared<Baby_full>("DY+jets (LO)", Process::Type::background, kOrange+1,
  //    {dir_mc_isr+"*DYJetsToLL_M-50_Tu*.root"},"1"); // Inclusive + HT-binned DY
  auto dyjets_nlo = Process::MakeShared<Baby_full>("DY+jets (NLO)", Process::Type::background, kRed-3,
    {dir_mc_isr_nlo+"*DYJetsToLL*.root"},"(ptll_me<100||type!=6200)"); 
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, kTeal-8,
    {dir_mc_isr+"*_TTWJets*.root", dir_mc_isr+"*_TTZTo*.root", dir_mc_isr+"*_TTGJets*.root"});
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {dir_mc_isr+"*_ZJet*.root", dir_mc_isr+"*_WJetsToLNu*.root",
    dir_mc_isr+"*QCD_HT*0_Tune*.root", dir_mc_isr+"*QCD_HT*Inf_Tune*.root",
    dir_mc_isr+"*ggZH_HToBB*.root", dir_mc_isr+"*_ttHJetTobb*.root",
    dir_mc_isr+"*_TTTT*.root", dir_mc_isr+"*_WWTo*.root",
    dir_mc_isr+"*_WH_HToBB*.root",dir_mc_isr+"*_ZH_HToBB*.root",
    dir_mc_isr+"*_WZ*.root",dir_mc_isr+"*_ZZ_*.root"},"stitch");

  string dir_data_isr = bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_minisrdata_"+isrtype+"/*root";
  string lumi_label = RoundNumber(lumi,1).Data();
  auto data = Process::MakeShared<Baby_full>("Data "+lumi_label+" fb^{-1}", Process::Type::data, kBlack,
    {dir_data_isr}, "pass" && Higfuncs::trig_hig>0.);

  vector<shared_ptr<Process> > procs = {data, dyjets, ttbar, single_t, ttv, other};
  vector<shared_ptr<Process> > procs_nlo = {data, dyjets_nlo, ttbar, single_t, ttv, other};
  vector<vector<shared_ptr<Process> > > vprocs;
  vprocs.push_back(procs);
  vprocs.push_back(procs_nlo);

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  vector<PlotOpt> plot_log = {log_lumi};
  vector<PlotOpt> plot_lin = {lin_lumi};
  vector<PlotOpt> plot_vals = {lin_lumi().PrintVals(true)};
  PlotMaker pm;

  string baseline = "nleps==2&&leps_pt[0]>40&&((elel_m>80&&elel_m<100)||(mumu_m>80&&mumu_m<100))";


  // definitions for njets in slices of ISR pT
  const vector<double> isr_syspt_slices = {0,500};
  const vector<double> isr_syspt_bins = {0, 50, 100, 150, 200, 300, 400, 600, 800, 1000};

  // definitions for ISR pT in slices of njets
  vector<double> nisrjet_slices = {0,6};
  const vector<double> nisrjet_bins = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5,7.5, 8.5, 9.5, 10.5};

  const vector<double> isr_ht_slices = {0,500};

  vector<NamedFunc> weight_opts;
  weight_opts.push_back("weight");
  weight_opts.push_back(NamedFunc("nominal", 
    [](const Baby &b) -> NamedFunc::ScalarType{
      double wgt_ = b.weight(); 
      if (b.type()==6202) wgt_*=0.984;
      else if (b.type()==6203) wgt_*=1.169;
      else if (b.type()==6204) wgt_*=1.267;
      else if (b.type()==6205) wgt_*=1.216;
      else if (b.type()==6100) wgt_*=0.84;
      else if (b.type()==6101) wgt_*=1.061;
      else if (b.type()==6102) wgt_*=0.954;
      else if (b.type()==6103) wgt_*=0.983;
      else if (b.type()==6104) wgt_*=0.974;
      else if (b.type()==6105) wgt_*=0.855;
      else if (b.type()==6106 || b.type()==6107) wgt_*=1.04;
      return wgt_;
    }));

  for (unsigned iproc(0); iproc<vprocs.size(); iproc++ ){
    for (const auto &iweight: weight_opts){
      string tag = isrtype;
      tag+=to_string(iproc);
      //print ratio

      pm.Push<Hist1D>(Axis(11,-0.5,10.5, "njets", "ISR jet multiplicity"), 
        baseline, vprocs[iproc], plot_log).Weight(iweight).Tag(tag);
      pm.Push<Hist1D>(Axis(25,0,1250, "jetsys_pt", "ISR p_{T} [GeV]"), 
        baseline, vprocs[iproc], plot_log).Weight(iweight).Tag(tag);
      pm.Push<Hist1D>(Axis(20,0.,2000., "ht", "H_{T} [GeV]"), 
        baseline, vprocs[iproc], plot_log).Weight(iweight).Tag(tag);

      pm.Push<Hist1D>(Axis(11,-0.5,10.5, "njets", "ISR jet multiplicity"), 
        baseline+"&&ht>500", vprocs[iproc], plot_log).Weight(iweight).Tag(tag);
      pm.Push<Hist1D>(Axis(25,0,1250, "jetsys_pt", "ISR p_{T} [GeV]"), 
        baseline+"&&njets>=6", vprocs[iproc], plot_log).Weight(iweight).Tag(tag);
      pm.Push<Hist1D>(Axis(20,0.,2000., "ht", "H_{T} [GeV]"), 
        baseline+"&&njets>=6", vprocs[iproc], plot_log).Weight(iweight).Tag(tag);

    } // Loop over weights
  }

  if (single_thread) pm.multithreaded_ = false;
  pm.min_print_ = true;
  pm.MakePlots(lumi);
    double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Making plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
} // main

void addSlices(PlotMaker &pm, const vector<double> slices, string svar,
               const vector<double> xbins, string xvar, string xlabel,
               const NamedFunc &baseline, const NamedFunc &weight,
               const vector<shared_ptr<Process> > &proc,
               const vector<PlotOpt> &plot_types, int tag_digits){

  //add the inclusive version first
  pm.Push<Hist1D>(Axis(xbins, xvar, xlabel), baseline, proc, plot_types).Weight(weight).Tag(isrtype+"_incl");
  for(unsigned i(0); i<slices.size(); i++){
    NamedFunc cut = baseline && (svar+">="+RoundNumber(slices[i],tag_digits).Data());
    if (i<(slices.size()-1)) cut = cut && (svar+"<"+RoundNumber(slices[i+1],tag_digits).Data());

    string tag = CodeToPlainText(isrtype+"_"+svar+RoundNumber(slices[i],tag_digits).Data());
    pm.Push<Hist1D>(Axis(xbins, xvar, xlabel), cut, proc, plot_types).Weight(weight).Tag(tag);
  }
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
    isr_jetspt.push_back(b.jets_pt()->at(ijet));
  }
  return isr_jetspt;
}


NamedFunc::ScalarType ZpTWeights(const Baby &b){
  if (b.type()<1000) return 1.; // Do not reweight Data

  float wgt = b.weight();
  if((b.type()>=2000 && b.type()<3000) ||   //wjets
    (b.type()>=6000 && b.type()<7000)) {    //dyjets

  wgt = b.weight()/b.w_isr();
  // int nisrjets(b.njets());
  // weights derived in DY+jets
  // if      (nisrjets==0) return 0.990*wgt; //  +- 0.001
  // else if (nisrjets==1) return 1.038*wgt; //  +- 0.001
  // else if (nisrjets==2) return 1.039*wgt; //  +- 0.003
  // else if (nisrjets==3) return 1.021*wgt; //  +- 0.007
  // else if (nisrjets==4) return 0.945*wgt; //  +- 0.013
  // else if (nisrjets==5) return 0.898*wgt; //  +- 0.025
  // else if (nisrjets==6) return 0.856*wgt; //  +- 0.048
  // else                  return 0.834*wgt; //  +- 0.088

    float sys_pt(b.jetsys_pt());

    // weights derived in DY+jets
    if      (sys_pt<=50)  return 0.993*wgt;
    else if (sys_pt<=100) return 1.045*wgt;
    else if (sys_pt<=150) return 1.171*wgt;
    else if (sys_pt<=200) return 1.142*wgt;
    else if (sys_pt<=300) return 1.050*wgt;
    else if (sys_pt<=400) return 0.993*wgt;
    else if (sys_pt<=600) return 0.906*wgt;
    else                  return 0.778*wgt; 
  } else {
    return wgt;
  }
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"type", required_argument, 0, 't'},  // Method to run on (if you just want one)
      {"quick", no_argument, 0, 0},           // Used inclusive ttbar for quick testing
      {"rewgt", no_argument, 0, 0},           // Used inclusive ttbar for quick testing
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "t:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 't':
      isrtype = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "quick"){
        quick = true;
      }else if(optname == "rewgt"){
        rewgt = true;
      }else{
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