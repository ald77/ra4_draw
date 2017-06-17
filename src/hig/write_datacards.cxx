// Adapted from ra4/src/syscalc_scan.cxx

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>
#include <ctime>
#include <algorithm>
#include <unistd.h> // getopt in Macs
#include <stdlib.h> // atof
#include <getopt.h>

#include "TSystem.h"
#include "TString.h"
#include "TMath.h"
#include "TError.h" // Controls error level reporting

#include "core/named_func.hpp"
#include "core/baby_full.hpp"
#include "core/utilities.hpp"

using namespace std;
namespace {
  bool fake_PU = false;
  float bf = 1.;
  bool incl_nonbb = false;
  bool incl_nonhh = false;
  bool old_cards = false;
  TString luminosity = "35.9";
  TString nom_wgt = "weight*eff_trig"; // nominal weight to use
  TString infolder = "/cms2r0/babymaker/babies/2017_03_17/TChiHH/merged_higsys_higsys/";
  TString infile = "*SMS-TChiHH_mGluino-200_mLSP-1_*.root";
  TString outfolder = ".";

  vector<string> metbins = {"met0","met1","met2"};
  // vector<string> metbins = {"met0","met1","met2","met3"};
  bool do_3bonly = false;
  enum SysType {kConst, kWeight, kSmear, kCorr, kMetSwap, kPU};
  bool nosys = false;

  const vector<double> v_data_npv{6.540e-06, 2.294e-05, 6.322e-05, 8.558e-05, 1.226e-04, 1.642e-04, 1.917e-04, 3.531e-04, 9.657e-04, 2.155e-03, 4.846e-03, 9.862e-03, 1.651e-02, 2.401e-02, 3.217e-02, 4.078e-02, 4.818e-02, 5.324e-02, 5.612e-02, 5.756e-02, 5.841e-02, 5.886e-02, 5.831e-02, 5.649e-02, 5.376e-02, 5.044e-02, 4.667e-02, 4.257e-02, 3.833e-02, 3.406e-02, 2.982e-02, 2.567e-02, 2.169e-02, 1.799e-02, 1.464e-02, 1.170e-02, 9.178e-03, 7.058e-03, 5.306e-03, 3.884e-03, 2.757e-03, 1.890e-03, 1.247e-03, 7.901e-04, 4.795e-04, 2.783e-04, 1.544e-04, 8.181e-05, 4.141e-05, 2.004e-05, 9.307e-06, 4.178e-06, 1.846e-06, 8.350e-07, 4.150e-07, 2.458e-07, 1.779e-07, 1.488e-07, 1.339e-07, 1.238e-07, 1.153e-07, 1.071e-07, 9.899e-08, 9.095e-08, 8.301e-08, 7.527e-08, 6.778e-08, 6.063e-08, 5.387e-08, 4.753e-08, 4.166e-08, 3.627e-08, 3.136e-08, 2.693e-08, 2.297e-08};
  TH1D h_data_npv("h_data_npv", "Data;N_{PV};P(N_{PV})", v_data_npv.size(), -0.5, v_data_npv.size()-0.5);
  TH1D h_mc_npv("h_mc_npv", "MC;N_{PV};P(N_{PV})", v_data_npv.size(), -0.5, v_data_npv.size()-0.5);
  double pu_low = 0.;
  double pu_high = 0.;

  vector<double> global_fit, observed;
  
}

class bindef {
public:
  bindef(TString itag, TString icut): tag(itag), cut(icut){};
  TString tag, cut;
};

class sysdef {
public:
  sysdef(TString ilabel, TString itag, SysType isystype): label(ilabel), tag(itag), sys_type(isystype) {
    v_wgts = vector<TString>();
  }
  // as will appear in the AN latex table
  TString label;
  // as will appear in the file handed to ra4 stats
  TString tag;
  // Is it a const, a weight or does it actually change the analysis variables, e.g. like lumi, b_tag or JEC
  SysType sys_type;
  // if sys_type = kSmear, what is the index in e.g. sys_met, check in babymaker:
  // https://github.com/manuelfs/babymaker/blob/2c0d9b2bde517b0bb129b8b3afffa77a581123e1/bmaker/interface/utilities.hh#L17 
  // if sys_type = kCorr, what is the index in e.g. sys_met, where the shifted *Up* value is stored, assuming Down is Up+1
  size_t shift_index; 
  // if sys_type = kWeight, add all weights to be used
  vector<TString> v_wgts;
  // here, we will store where this systematic begins in the big yields & entires vectors that we get from getYields()
  // from there, indices are order like: nVariations*iBin + iVariation 
  size_t ind;
};

TString nom2sys_bin(TString ibin, size_t shift_index);
TString nom2genmet(TString ibin);
void fillHiggsinoSys(ofstream &fcard);
vector<double> getYields(Baby_full &baby, const NamedFunc &baseline, const vector<NamedFunc> &bincuts,
                         vector<double> &yield, vector<double> &w2, double lumi,
                         bool do_trig = false, const TString &flag = "");
void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing bfanches
  time_t begtime, endtime;
  time(&begtime);
  
  GetOptions(argc, argv);
  gSystem->mkdir(outfolder, kTRUE);
  // if (nosys) infolder.ReplaceAll("higsys_higsys","higmc_higtight");
  // TString infile = "/cms2r0/babymaker/babies/2015_11_27/sms/split_sms/renorm/baby_SMS-T1tttt_mGluino-1500_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15FSPremix-MCRUN2_74_V9_renorm.root";
  string prs = infile.Data();
  int mglu = stoi(prs.substr(prs.find("ino-")+4,prs.find("_mLSP")-prs.find("ino-")-4));
  int mlsp = stoi(prs.substr(prs.find("LSP-")+4,prs.find("_Tune")-prs.find("LSP-")-4));
  cout<<"Working on: mGluino = "<<mglu<<" mLSP = "<<mlsp<<endl;
  string glu_lsp("mGluino-"+to_string(mglu)+"_mLSP-"+to_string(mlsp));
  string model = "TChiHZ";
  if (infolder.Contains("03_17")) model = "TChiHH";

  //------------- SYSTEMATICS DEFINITIONS -----------------------------
  //// tables has a vector of the tables you want to print
  TString baseline("higd_drmax<2.2&&ntks==0&&njets>=4&&njets<=5&&!low_dphi&&nvleps==0&&pass_ra2_badmu&&met/met_calo<5");
  vector<bindef> v_bins;
  TString cut2b="nbdt==2&&nbdm==2", cut3b="nbdt>=2&&nbdm==3&&nbdl==3", cut4b="nbdt>=2&&nbdm>=3&&nbdl>=4";
  TString cuthig="higd_am>100&&higd_am<140&&higd_dm<40";
  TString cutsbd="!(higd_am>100&&higd_am<140)&&higd_dm<40&&higd_am<200";
  for (unsigned imet(0); imet<metbins.size(); imet++) {
    if (metbins[imet]=="met0") {
      v_bins.push_back(bindef("sbd_2b_met0", cut2b+"&&"+cutsbd+"&& met>150&&met<=200"));    global_fit.push_back(1560.1);    observed.push_back(1559);  
      v_bins.push_back(bindef("hig_2b_met0", cut2b+"&&"+cuthig+"&& met>150&&met<=200"));    global_fit.push_back(656.2);     observed.push_back(658);   
      v_bins.push_back(bindef("sbd_3b_met0", cut3b+"&&"+cutsbd+"&& met>150&&met<=200"));    global_fit.push_back(140.3);     observed.push_back(145);   
      v_bins.push_back(bindef("hig_3b_met0", cut3b+"&&"+cuthig+"&& met>150&&met<=200"));    global_fit.push_back(57.7);      observed.push_back(53);    
      if (!do_3bonly) {
        v_bins.push_back(bindef("sbd_4b_met0", cut4b+"&&"+cutsbd+"&& met>150&&met<=200"));    global_fit.push_back(48.1);      observed.push_back(45);    
        v_bins.push_back(bindef("hig_4b_met0", cut4b+"&&"+cuthig+"&& met>150&&met<=200"));    global_fit.push_back(21.9);      observed.push_back(25);    
      } 
    }
    if (metbins[imet]=="met1") {
      v_bins.push_back(bindef("sbd_2b_met1", cut2b+"&&"+cutsbd+"&& met>200&&met<=300"));    global_fit.push_back(588.0);     observed.push_back(585);   
      v_bins.push_back(bindef("hig_2b_met1", cut2b+"&&"+cuthig+"&& met>200&&met<=300"));    global_fit.push_back(333.1);     observed.push_back(336);   
      v_bins.push_back(bindef("sbd_3b_met1", cut3b+"&&"+cutsbd+"&& met>200&&met<=300"));    global_fit.push_back(55.3);      observed.push_back(61);    
      v_bins.push_back(bindef("hig_3b_met1", cut3b+"&&"+cuthig+"&& met>200&&met<=300"));    global_fit.push_back(30.6);      observed.push_back(25);    
      if (!do_3bonly) {
        v_bins.push_back(bindef("sbd_4b_met1", cut4b+"&&"+cutsbd+"&& met>200&&met<=300"));    global_fit.push_back(15.6);      observed.push_back(13);    
        v_bins.push_back(bindef("hig_4b_met1", cut4b+"&&"+cuthig+"&& met>200&&met<=300"));    global_fit.push_back(11.4);      observed.push_back(14);    
      }
    }
    if (metbins[imet]=="met2") {
      v_bins.push_back(bindef("sbd_2b_met2", cut2b+"&&"+cutsbd+"&& met>300&&met<=450"));    global_fit.push_back(72.4);      observed.push_back(74);    
      v_bins.push_back(bindef("hig_2b_met2", cut2b+"&&"+cuthig+"&& met>300&&met<=450"));    global_fit.push_back(40.6);      observed.push_back(39);    
      v_bins.push_back(bindef("sbd_3b_met2", cut3b+"&&"+cutsbd+"&& met>300&&met<=450"));    global_fit.push_back(5.7);       observed.push_back(4);     
      v_bins.push_back(bindef("hig_3b_met2", cut3b+"&&"+cuthig+"&& met>300&&met<=450"));    global_fit.push_back(3.3);       observed.push_back(5);     
      if (!do_3bonly) {
        v_bins.push_back(bindef("sbd_4b_met2", cut4b+"&&"+cutsbd+"&& met>300&&met<=450"));    global_fit.push_back(1.9);       observed.push_back(2);     
        v_bins.push_back(bindef("hig_4b_met2", cut4b+"&&"+cuthig+"&& met>300&&met<=450"));    global_fit.push_back(1.1);       observed.push_back(1);     
      }
    }

    if (metbins[imet]=="met3") {
      v_bins.push_back(bindef("sbd_2b_met3", cut2b+"&&"+cutsbd+"&& met>450"));              global_fit.push_back(5.4);       observed.push_back(5);     
      v_bins.push_back(bindef("hig_2b_met3", cut2b+"&&"+cuthig+"&& met>450"));              global_fit.push_back(4.6);       observed.push_back(5);     
      v_bins.push_back(bindef("sbd_3b_met3", cut3b+"&&"+cutsbd+"&& met>450"));              global_fit.push_back(0.6);       observed.push_back(1);     
      v_bins.push_back(bindef("hig_3b_met3", cut3b+"&&"+cuthig+"&& met>450"));              global_fit.push_back(0.4);       observed.push_back(0);     
      if (!do_3bonly) {
        v_bins.push_back(bindef("sbd_4b_met3", cut4b+"&&"+cutsbd+"&& met>450"));              global_fit.push_back(0.0001);       observed.push_back(0);     
        v_bins.push_back(bindef("hig_4b_met3", cut4b+"&&"+cuthig+"&& met>450"));              global_fit.push_back(0.0001);       observed.push_back(0);     
      }
    }
  }

  // ------------ fits for extrapolating yields including contamination ---------------
  map<TString, vector<double> > hz_fit;
  hz_fit["sbd_2b"] = { 1.79760e+00, -1.42525e-01 };
  hz_fit["hig_2b"] = { 7.47671e-01,  1.83420e-01 };
  hz_fit["sbd_3b"] = { 6.07960e-01, -1.03622e+00 };
  hz_fit["hig_3b"] = { 1.25460e-01, -2.26228e+00 };

  map<TString, vector<double> > zz_fit;
  zz_fit["sbd_2b"] = { 1.81166e+00, -1.42525e-01 };
  zz_fit["hig_2b"] = { 7.85145e-02, -1.19438e+00 };
  zz_fit["sbd_3b"] = { 3.62738e-01, -1.18628e+00 };
  zz_fit["hig_3b"] = { 8.91456e-03, -3.97757e+00 };

  //------------- SYSTEMATICS DEFINITIONS -----------------------------
  vector<sysdef> v_sys;
  // order as they will appear in latex table
  // *Nominal must stay in the first spot!!* (will be skipped in table)
  v_sys.push_back(sysdef("Nominal", "nominal", kWeight)); 
  v_sys.back().v_wgts.push_back("1.");
  nom_wgt = "weight/w_btag*w_bhig_deep*eff_trig"; // nominal weight to use
  if (!nosys) {
    v_sys.push_back(sysdef("Trigger efficiency", "trig_HH", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_trig["+to_string(i)+"]"); // the TChiHH babies have relative mutliplicative sys_trig
    v_sys.push_back(sysdef("B-tag efficiency", "btag_HF_deep_HH", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_bchig_deep["+to_string(i)+"]/w_bhig_deep");
    v_sys.push_back(sysdef("B-tag efficiency FS", "btag_FS_HF_deep_HH", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_fs_bchig_deep["+to_string(i)+"]/w_bhig_deep");
    v_sys.push_back(sysdef("Mistag efficiency", "btag_LF_deep_HH", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_udsghig_deep["+to_string(i)+"]/w_bhig_deep");
    v_sys.push_back(sysdef("Mistag efficiency FS", "btag_FS_LF_deep_HH",kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_fs_udsghig_deep["+to_string(i)+"]/w_bhig_deep");

    v_sys.push_back(sysdef("Gen vs reco MET FS", "MET",kMetSwap));
    
    v_sys.push_back(sysdef("Jet energy corrections", "JESsig", kCorr));
    v_sys.back().shift_index = 1; // JEC Up index in sys_met, etc.

    v_sys.push_back(sysdef("Jet energy resolution", "JER_HH", kSmear));
    v_sys.back().shift_index = 0; // JER index in sys_met, etc.

    v_sys.push_back(sysdef("ISR", "ISR", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_isr["+to_string(i)+"]/w_isr");
    v_sys.push_back(sysdef("Jet ID FS", "JetID_HH", kConst));
    v_sys.back().v_wgts.push_back("0.01");
    v_sys.push_back(sysdef("Pile up", "PUsig", kPU));
    v_sys.push_back(sysdef("Luminosity", "lumi", kConst));
    v_sys.back().v_wgts.push_back("0.026");
  }


  /////////////////////////////  No more changes needed down here to add systematics ///////////////////////
  // prepare the vector of bincuts used to get the yields
  NamedFunc bf_wgt("bf_wgt", [&](const Baby &b){
    float wgt_ = 1;
    if (b.type()==-999999){
      int nh(0), nh_nonbb(0);
      for (unsigned i(0); i<b.mc_id()->size(); i++) {
        if (b.mc_id()->at(i)==25) nh++;
        if (b.mc_mom()->at(i)==25 && abs(b.mc_id()->at(i))!=5) nh_nonbb++;
      }
      nh_nonbb /=2;
      if (incl_nonhh) wgt_ *= (bf*bf/.25*(nh==2) + 2*bf*(1-bf)/.5*(nh==1) + (1-bf)*(1-bf)/.25*(nh==0));
      else wgt_ *= (bf*bf/.25*(nh==2));
      if (!incl_nonbb) wgt_ *= (nh_nonbb==0);
    } else if (incl_nonhh) {
      float mass = b.mgluino();
      double hz_wgt(0), zz_wgt(0);
      if (b.nbdt()==2 && b.nbdm()==2){
        if (b.higd_am()>100 && b.higd_am()<140 && b.higd_dm()<40) {
          hz_wgt = hz_fit["hig_2b"][0]+TMath::Exp(hz_fit["hig_2b"][1]-1e-05*mass*mass);
          zz_wgt = zz_fit["hig_2b"][0]+TMath::Exp(zz_fit["hig_2b"][1]-1e-05*mass*mass);
        } else {
          hz_wgt = hz_fit["sbd_2b"][0]+TMath::Exp(hz_fit["sbd_2b"][1]-1e-05*mass*mass);
          zz_wgt = zz_fit["sbd_2b"][0]+TMath::Exp(zz_fit["sbd_2b"][1]-1e-05*mass*mass);
        }
      } else {
        if (b.higd_am()>100 && b.higd_am()<140 && b.higd_dm()<40) {
          hz_wgt = hz_fit["hig_3b"][0]+TMath::Exp(hz_fit["hig_3b"][1]-1e-05*mass*mass);
          zz_wgt = zz_fit["hig_3b"][0]+TMath::Exp(zz_fit["hig_3b"][1]-1e-05*mass*mass);
        } else {
          hz_wgt = hz_fit["sbd_3b"][0]+TMath::Exp(hz_fit["sbd_3b"][1]-1e-05*mass*mass);
          zz_wgt = zz_fit["sbd_3b"][0]+TMath::Exp(zz_fit["sbd_3b"][1]-1e-05*mass*mass);
        }
      }
      wgt_ *= bf*bf + 2*bf*(1-bf)*hz_wgt + (1-bf)*(1-bf)*zz_wgt;
    } else {
      wgt_ *= bf*bf;
    }
    return wgt_;
  });

  vector<NamedFunc> bcuts;
  sysdef nom = v_sys[0];
  if (nom.tag != "nominal"){
    cerr<<" The first entry in the v_sys vector must be the nominal"<<endl;
    exit(1);
  }
  for (auto &sys: v_sys) {
    sys.ind = bcuts.size(); 
    if (sys.sys_type == kConst){
      continue;
    } else if (sys.sys_type == kWeight) {
      for (auto &bin: v_bins) {
        for (auto &wgt: sys.v_wgts) {
          bcuts.emplace_back(("("+baseline+"&&"+bin.cut+")*"+nom_wgt+"*"+wgt) * bf_wgt);
        }
      }
    } else if (sys.sys_type == kCorr || sys.sys_type == kSmear) {
      for (auto &bin: v_bins) {
        bcuts.emplace_back((nom2sys_bin("("+baseline+"&&"+bin.cut, sys.shift_index)+")*"+nom_wgt) * bf_wgt);
        if (sys.sys_type == kCorr) { //if it is a correction, need to push the 'down' variation as well
          bcuts.emplace_back((nom2sys_bin("("+baseline+"&&"+bin.cut, sys.shift_index+1)+")*"+nom_wgt) * bf_wgt);
        }
      }
    } else if (sys.sys_type == kMetSwap){
      for (auto &bin: v_bins) {
        bcuts.emplace_back((nom2genmet("("+baseline+"&&"+bin.cut)+")*"+nom_wgt) * bf_wgt);
      }
    } else if (sys.sys_type == kPU) {
      for(const auto &bin: v_bins){
        bcuts.emplace_back(("("+baseline+"&&"+bin.cut+"&&"+"npv<=20)*"+nom_wgt) * bf_wgt);
        bcuts.emplace_back(("(npv<=20)*"+nom_wgt) * bf_wgt);
        bcuts.emplace_back(("("+baseline+"&&"+bin.cut+"&&"+"npv>=21)*"+nom_wgt) * bf_wgt);
        bcuts.emplace_back(("(npv>=21)*"+nom_wgt) * bf_wgt);
      }
    }
  }
  
  // get yields from the baby for all the cut strings
  cout<<"Running on: "<<infolder+"/"+infile<<endl;
  Baby_full baby(std::set<std::string>{(infolder+"/"+infile).Data()});
  auto activator = baby.Activate();
  vector<double> yields, w2, entries;
  entries = getYields(baby, baseline, bcuts, yields, w2, luminosity.Atof());

  // --------- Writing datacard -----------------------
  TString outpath = outfolder+"/datacard_SMS-"+TString(model)+"_"+glu_lsp+"_bfH-"+RoundNumber(bf*100, 0);
  if (incl_nonbb) outpath += "_allHigDecays";
  if (incl_nonhh && bf<1.) outpath += "_withZContam";
  outpath += "_"+luminosity.ReplaceAll(".","p")+"ifb";
  if (old_cards) outpath += "_old";
  outpath += ".txt";
  cout<<"open "<<outpath<<endl;
  unsigned wname(25), wdist(5), wbin(15);
  unsigned nbins(v_bins.size());
  unsigned nmet(metbins.size());
  // --------- write header
  ofstream fcard(outpath);
  fcard<<"imax "<<nbins<<"  number of channels\n";
  fcard<<"jmax 1  number of backgrounds\n";
  fcard<<"kmax *  number of nuisance parameters\n";
  fcard<<"shapes * * FAKE\n";
  fcard<<endl<<setw(wname)<<"bin"<<setw(wdist)<<" ";
  for (size_t ibin(0); ibin<nbins; ibin++) fcard<<setw(wbin)<<" "<<setw(wbin)<<v_bins[ibin].tag;
  fcard<<endl<<setw(wname)<<"Observation"<<setw(wdist)<<" ";
  for (size_t ibin(0); ibin<nbins; ibin++) fcard<<setw(wbin)<<" "<<setw(wbin)<<observed[ibin];

  fcard<<endl<<endl<<setw(wname)<<"bin"<<setw(wdist)<<" ";
  for (size_t ibin(0); ibin<nbins; ibin++) fcard<<setw(wbin)<<v_bins[ibin].tag<<setw(wbin)<<v_bins[ibin].tag;
  fcard<<endl<<setw(wname)<<"process"<<setw(wdist)<<" ";
  for (size_t ibin(0); ibin<nbins; ibin++) fcard<<setw(wbin)<<"sig"<<setw(wbin)<<"bkg";
  fcard<<endl<<setw(wname)<<"process"<<setw(wdist)<<" ";
  for (size_t ibin(0); ibin<nbins; ibin++) fcard<<setw(wbin)<<"0"<<setw(wbin)<<"1";
  fcard<<endl<<setw(wname)<<"rate"<<setw(wdist)<<" ";
  if (old_cards) for (size_t ibin(0); ibin<nbins; ibin++) fcard<<setw(wbin)<<Form("%.2f",yields[ibin])<<setw(wbin)<<global_fit[ibin];
  else           for (size_t ibin(0); ibin<nbins; ibin++) fcard<<setw(wbin)<<Form("%.2f",yields[ibin])<<setw(wbin)<<"1";
  fcard<<endl<<endl;
  cout<<"Wrote headers"<<endl;

  //--------- Writing ABCD constraints----------------------------
  unsigned nnb = 2;
  if (!do_3bonly) nnb = 3;

  if (old_cards) {
    for (size_t imet(0); imet<nmet; imet++) {
      double cnb(0); TString cnbstr;
      for (size_t inb(0); inb<nnb; inb++) {
        fcard<<setw(wname)<<"c_"+metbins[imet]+"_"+to_string(inb+2)+"b_HH"<<setw(wdist)<<"lnU";
        if (observed[imet*2*nnb+2*inb]<0.1) cnb = 99999.;
        else cnb = 1+6/sqrt(observed[imet*2*nnb+2*inb]); // <-- based on observed in SBD for this nb category
        cnbstr = Form("%.2f",cnb);
        for (size_t jmet(0); jmet<nmet; jmet++) {    
          for (size_t jnb(0); jnb<nnb; jnb++) {
            if (imet==jmet && inb==jnb) fcard<<setw(wbin)<<"-"<<setw(wbin)<<cnbstr<<setw(wbin)<<"-"<<setw(wbin)<<cnbstr;
            else fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-"<<setw(wbin)<<"-"<<setw(wbin)<<"-";
          }
        }
        fcard<<endl;
      }
      fcard<<setw(wname)<<"c_"+metbins[imet]+"_allnb_HH"<<setw(wdist)<<"lnU";
      for (size_t jmet(0); jmet<nmet; jmet++) {    
        for (size_t jnb(0); jnb<nnb; jnb++) {
          if (imet==jmet) fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-"<<setw(wbin)<<"-"<<setw(wbin)<<cnbstr;
          else fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-"<<setw(wbin)<<"-"<<setw(wbin)<<"-";
        }
      }
      fcard<<endl;
    }
  }

  //--------- Signal statistical uncertainties ----------------------------
  for (size_t ibin(0); ibin<nbins; ibin++) {
    fcard<<setw(wname)<<"sig_stat_"+v_bins[ibin].tag+"_HH"<<setw(wdist)<<"lnN";
    TString sig_stat = Form("%.2f",1.+sqrt(w2[ibin])/yields[ibin]);
    if (yields[ibin]<0.00001) sig_stat = "99999.00";
    for (size_t jbin(0); jbin<nbins; jbin++) {
      if (ibin==jbin) fcard<<setw(wbin)<<sig_stat<<setw(wbin)<<"-";
      else fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-";
    }
    fcard<<endl;
  }
  cout<<"Wrote signal stat. uncertainties"<<endl;

  //--------- MC statistical uncertainty ----------------------------
  vector<double> mcstat_unc = {0.08, 0.15, 0.05, 0.37, 0.19, 0.23, 0.59, 0.75};
  if (nmet<4 || do_3bonly) {
    cout<<"Datacards with uncertainty not supported for partial binning"<<endl; 
  } else {
    for (size_t imet(0); imet<nmet; imet++) {
      for (size_t inb(0); inb<nnb; inb++) {
        if (inb==0) continue;
        fcard<<setw(wname)<<"mcstat_"+metbins[imet]+"_"+to_string(inb+2)+"b_HH"<<setw(wdist)<<"lnN";
        TString uncstr = Form("%.2f",1+mcstat_unc[imet*(nnb-1)+inb-1]);
        for (size_t jmet(0); jmet<nmet; jmet++) {    
          for (size_t jnb(0); jnb<nnb; jnb++) {
            if (imet==jmet && inb==jnb) fcard<<setw(wbin)<<"-"<<setw(wbin)<<uncstr<<setw(wbin)<<"-"<<setw(wbin)<<"-";
            else fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-"<<setw(wbin)<<"-"<<setw(wbin)<<"-";
          }
        }
        fcard<<endl;
      }
    }  
    cout<<"Wrote MC stat. closure uncertainties"<<endl;

    // ------------ Closure uncertainties
    vector<TString> closure_unc_names; vector< vector<double> > closure_unc;
    if (old_cards) {closure_unc_names.push_back("wilks_HH"); closure_unc.push_back({-9999, -9999, -9999, -9999, -9999, -9999, -9999, 10.});}
    closure_unc_names.push_back("ttx_closure_HH"); closure_unc.push_back({0.03, 0.05, 0.03, 0.06, 0.02, 0.04, 0.02, 0.03});
    closure_unc_names.push_back("zll_closure_HH"); closure_unc.push_back({0.01, 0.01, 0.03, 0.01, 0.06, 0.04, 0.08, 0.09});
    closure_unc_names.push_back("qcd_closure_HH"); closure_unc.push_back({0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01});
    closure_unc_names.push_back("bkg_comp_HH"); closure_unc.push_back({-0.03, -0.05,  0.03,  0.01, -0.04,  0.02,  0.06,  0.06});
    for (size_t iunc(0); iunc<closure_unc.size(); iunc++) {
      fcard<<setw(wname)<<closure_unc_names[iunc]<<setw(wdist)<<"lnN";
      for (size_t imet(0); imet<nmet; imet++) {    
        for (size_t inb(0); inb<nnb; inb++) {
          if (inb==0 || closure_unc[iunc][imet*(nnb-1)+inb-1]<-10) fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-"<<setw(wbin)<<"-"<<setw(wbin)<<"-";
          else fcard<<setw(wbin)<<"-"<<setw(wbin)<<Form("%.2f",1+closure_unc[iunc][imet*(nnb-1)+inb-1])<<setw(wbin)<<"-"<<setw(wbin)<<"-";
        }
      }
      fcard<<endl;
    }
    cout<<"Wrote CR-based closure uncertainties"<<endl;

    for (auto &sys: v_sys) {
      if (sys.tag == "nominal") continue;
      fcard<<setw(wname)<<sys.tag<<setw(wdist)<<"lnN";
      for (size_t ibin = 0; ibin<nbins; ++ibin) {
        const double nom_yield(yields[ibin]);
        double up(0.), dn(0.);
        if (sys.sys_type == kConst) {
          up = stod(sys.v_wgts[0].Data());
          dn = -up;
        } else if (sys.sys_type == kWeight) {
          up = yields[sys.ind + 2*ibin]/nom_yield - 1;
          dn = yields[sys.ind + 2*ibin + 1]/nom_yield - 1;
        } else if (sys.sys_type == kSmear) {
          up = yields[sys.ind + ibin]/nom_yield - 1;
          dn = -up;
        } else if (sys.sys_type == kMetSwap) {
          //Use average of met yield and met_tru yield as central value
          dn = yields[sys.ind + ibin]/(0.5*(yields[sys.ind + ibin]+nom_yield)) - 1;
          up = -dn;
        } else if (sys.sys_type == kCorr) {
          up = yields[sys.ind + 2*ibin]/nom_yield - 1;
          dn = yields[sys.ind + 2*ibin + 1]/nom_yield - 1;
        } else if (sys.sys_type == kPU ) {
          double eff_low  = yields[sys.ind+4*ibin+0]/yields[sys.ind+4*ibin+1];
          double eff_high = yields[sys.ind+4*ibin+2]/yields[sys.ind+4*ibin+3];
          double m = (eff_high-eff_low)/(pu_high-pu_low);
          double b = (eff_low*pu_high-eff_high*pu_low)/(pu_high-pu_low);
          double eff_data = 0., eff_mc = 0.;
          for(size_t i = 0; i < v_data_npv.size(); ++i){
            double fx = m*static_cast<double>(i)+b;
            eff_data += fx*h_data_npv.GetBinContent(i+1);
            eff_mc += fx*h_mc_npv.GetBinContent(i+1);
          }
          up = (eff_data-eff_mc)/eff_mc;
          dn = -up;

          //Temporary stand-in until better method available
          if(fake_PU){
            if((v_bins.at(ibin).tag.Contains("lowmet") || v_bins.at(ibin).tag.Contains("medmet"))
               && (v_bins.at(ibin).tag.Contains("lownj") || v_bins.at(ibin).tag.Contains("r1_") || v_bins.at(ibin).tag.Contains("r3_"))){
              up = 0.1;
              dn = -0.1;
            }else{
              up = 0.15;
              dn = -0.15;
            }
          }
        }
        // convert to ra4_stats input and write to file
        // double ln = (up>0 ? 1:-1)*max(up>0 ? up : (1/(1+up)-1), dn>0 ? dn : (1/(1+dn)-1));
        double ln = max(up>0 ? 1+up : 1/(fabs(up)+1), dn>0 ? 1+dn : 1/(fabs(dn)+1));
        if (std::isnan(ln) || std::isinf(ln)) {
          cout <<" Found bad unc. set to 0 -> "<<std::left<<setw(10)<<sys.tag <<std::left<<setw(10)<<v_bins[ibin].tag <<" "<<std::right<<setprecision(0)<<setw(25)<<entries[ibin] <<" "<<setprecision(5)<<setw(15)<<yields[ibin] <<" "<<setprecision(10)<<setw(15)<<w2[ibin] <<endl;  
          ln = 0;
        } 
        if (sys.sys_type == kConst) ln = 1+up;
        fcard<<setw(wbin)<<Form("%.2f",ln)<<setw(wbin)<<"-";
      } // loop over bins
      fcard<<endl;
    } // loop over systematics
    cout<<"Wrote systematics"<<endl;
  }

  if (!old_cards) {
    fcard<<endl;
    for (unsigned imet(0); imet<metbins.size(); imet++) {
      if (metbins[imet]=="met0") {
        fcard<<"rp_hig_3b_met0 rateParam hig_3b_met0 bkg (@0*@1/@2) rp_sbd_3b_met0,rp_hig_2b_met0,rp_sbd_2b_met0"<<endl;
        if (!do_3bonly) fcard<<"rp_hig_4b_met0 rateParam hig_4b_met0 bkg (@0*@1/@2) rp_sbd_4b_met0,rp_hig_2b_met0,rp_sbd_2b_met0"<<endl;
        fcard<<"rp_sbd_2b_met0 rateParam sbd_2b_met0 bkg 1559"<<endl;
        fcard<<"rp_hig_2b_met0 rateParam hig_2b_met0 bkg 658"<<endl;
        fcard<<"rp_sbd_3b_met0 rateParam sbd_3b_met0 bkg 145"<<endl;
        if (!do_3bonly) fcard<<"rp_sbd_4b_met0 rateParam sbd_4b_met0 bkg 45"<<endl<<endl;
      } else if (metbins[imet]=="met1") {
        fcard<<"rp_hig_3b_met1 rateParam hig_3b_met1 bkg (@0*@1/@2) rp_sbd_3b_met1,rp_hig_2b_met1,rp_sbd_2b_met1"<<endl;
        if (!do_3bonly) fcard<<"rp_hig_4b_met1 rateParam hig_4b_met1 bkg (@0*@1/@2) rp_sbd_4b_met1,rp_hig_2b_met1,rp_sbd_2b_met1"<<endl;
        fcard<<"rp_sbd_2b_met1 rateParam sbd_2b_met1 bkg 585"<<endl;
        fcard<<"rp_hig_2b_met1 rateParam hig_2b_met1 bkg 336"<<endl;
        fcard<<"rp_sbd_3b_met1 rateParam sbd_3b_met1 bkg 61"<<endl;
        if (!do_3bonly) fcard<<"rp_sbd_4b_met1 rateParam sbd_4b_met1 bkg 13"<<endl<<endl;
      } else if (metbins[imet]=="met2") {
        fcard<<"rp_hig_3b_met2 rateParam hig_3b_met2 bkg (@0*@1/@2) rp_sbd_3b_met2,rp_hig_2b_met2,rp_sbd_2b_met2"<<endl;
        if (!do_3bonly) fcard<<"rp_hig_4b_met2 rateParam hig_4b_met2 bkg (@0*@1/@2) rp_sbd_4b_met2,rp_hig_2b_met2,rp_sbd_2b_met2"<<endl;
        fcard<<"rp_sbd_2b_met2 rateParam sbd_2b_met2 bkg 74"<<endl;
        fcard<<"rp_hig_2b_met2 rateParam hig_2b_met2 bkg 39"<<endl;
        fcard<<"rp_sbd_3b_met2 rateParam sbd_3b_met2 bkg 4"<<endl;
        if (!do_3bonly) fcard<<"rp_sbd_4b_met2 rateParam sbd_4b_met2 bkg 2"<<endl<<endl;
      } else if (metbins[imet]=="met3") {
        fcard<<"rp_hig_3b_met3 rateParam hig_3b_met3 bkg (@0*@1/@2) rp_sbd_3b_met3,rp_hig_2b_met3,rp_sbd_2b_met3"<<endl;
        if (!do_3bonly) fcard<<"rp_hig_4b_met3 rateParam hig_4b_met3 bkg (@0*@1/@2) rp_sbd_4b_met3,rp_hig_2b_met3,rp_sbd_2b_met3"<<endl;
        fcard<<"rp_sbd_2b_met3 rateParam sbd_2b_met3 bkg 5"<<endl;
        fcard<<"rp_hig_2b_met3 rateParam hig_2b_met3 bkg 5"<<endl;
        fcard<<"rp_sbd_3b_met3 rateParam sbd_3b_met3 bkg 1"<<endl;
        if (!do_3bonly) fcard<<"rp_sbd_4b_met3 rateParam sbd_4b_met3 bkg 0"<<endl;
      }
    }
  } 

  fcard.close();

  cout<<" open "<<outpath<<endl;

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

TString nom2sys_bin(TString ibin, size_t shift_index){
  ibin.ReplaceAll("met", "sys_met["+to_string(shift_index)+"]");
  ibin.ReplaceAll("mt", "sys_mt["+to_string(shift_index)+"]");
  ibin.ReplaceAll("st", "sys_st["+to_string(shift_index)+"]");
  ibin.ReplaceAll("mj14", "sys_mj14["+to_string(shift_index)+"]");
  ibin.ReplaceAll("njets", "sys_njets["+to_string(shift_index)+"]");
  ibin.ReplaceAll("nbm", "sys_nbm["+to_string(shift_index)+"]");
  //replacements for higgsino
  ibin.ReplaceAll("nbdl", "sys_nbdl["+to_string(shift_index)+"]");
  ibin.ReplaceAll("nbdm", "sys_nbdm["+to_string(shift_index)+"]");
  ibin.ReplaceAll("nbdt", "sys_nbdt["+to_string(shift_index)+"]");
  ibin.ReplaceAll("higd_am", "sys_higd_am["+to_string(shift_index)+"]");
  ibin.ReplaceAll("higd_dm", "sys_higd_dm["+to_string(shift_index)+"]");
  ibin.ReplaceAll("higd_drmax", "sys_higd_drmax["+to_string(shift_index)+"]");
  //fix unintended replacement...
  ibin.ReplaceAll("sys_met["+to_string(shift_index)+"]/sys_met["+to_string(shift_index)+"]_calo","met/met_calo");
  return ibin;
}

TString nom2genmet(TString ibin){
  ibin.ReplaceAll("met", "met_tru");
  //fix unintended replacement...
  ibin.ReplaceAll("met_tru/met_tru_calo", "met/met_calo");
  return ibin;

}

vector<double> getYields(Baby_full &baby, const NamedFunc &/*baseline*/, const vector<NamedFunc> &bincuts,
                         vector<double> &yield, vector<double> &w2, double lumi,
                         bool do_trig, const TString &flag){
  for(size_t i = 0; i <v_data_npv.size(); ++i){
    h_data_npv.SetBinContent(i+1, v_data_npv.at(i));
    h_data_npv.SetBinError(i+1, 0.);
    h_mc_npv.SetBinContent(i+1, 0.);
    h_mc_npv.SetBinError(i+1, 0.);
  }
  vector<double> entries = vector<double>(bincuts.size(), 0);
  yield = vector<double>(bincuts.size(), 0);
  w2 = yield;

  long nentries = baby.GetEntries();
  for(long entry = 0; entry < nentries; ++entry){
    baby.GetEntry(entry);
    h_mc_npv.Fill(baby.npv(), baby.weight()*baby.eff_trig());
    if(do_trig){
      if(!baby.pass()) continue;
      if(!baby.trig()->at(4) && !baby.trig()->at(8) && !baby.trig()->at(13) && !baby.trig()->at(33)) continue;
    }
    //if(!baseline.GetScalar(baby)) continue;
    for(size_t ind = 0; ind<bincuts.size(); ++ind){
      float wgt = bincuts.at(ind).GetScalar(baby);
      if(wgt != 0.){
        ++entries.at(ind);

        if(flag=="jer_tail"){
          double jet_res_min = *min_element(baby.jets_pt_res()->begin(), baby.jets_pt_res()->end());
          double jet_res_max = *max_element(baby.jets_pt_res()->begin(), baby.jets_pt_res()->end());
          if((jet_res_min>0&&jet_res_min<0.675) || jet_res_max>1.391)
            wgt *= 1.5;
        }

        yield.at(ind) += wgt;
        w2.at(ind) += wgt*wgt;
      }
    }
  } // Loop over entries
  for(size_t ind = 0; ind<bincuts.size(); ++ind){ 
     yield.at(ind) *= lumi;
     w2.at(ind) *= pow(lumi, 2);
  }
  h_data_npv.Scale(1./h_data_npv.Integral());
  h_mc_npv.Scale(1./h_mc_npv.Integral());
  pu_low = 0.;
  pu_high = 0.;
  double norm = 0.;
  for(size_t npv = 0; npv <= 20 && npv < v_data_npv.size(); ++npv){
    pu_low += npv*h_mc_npv.GetBinContent(npv+1);
    norm += h_mc_npv.GetBinContent(npv+1);
  }
  pu_low /= norm;
  norm = 0.;
  for(size_t npv = 21; npv < v_data_npv.size(); ++npv){
    pu_high += npv*h_mc_npv.GetBinContent(npv+1);
    norm += h_mc_npv.GetBinContent(npv+1);
  }
  pu_high /= norm;
  return entries;
}

void GetOptions(int argc, char *argv[]){
  string blah;
  while(true){
    static struct option long_options[] = {
      {"nosys", required_argument, 0, 'n'},
      {"infolder", required_argument, 0, 'i'},
      {"infile", required_argument, 0, 'f'},
      {"outfolder", required_argument, 0, 'o'},
      {"lumi", required_argument, 0, 'l'},
      {"bf", required_argument, 0, 0},
      {"incl_nonbb", no_argument, 0, 0},
      {"incl_nonhh", no_argument, 0, 0},
      {"old", no_argument, 0, 0},
      {"fake_PU", no_argument, 0, 0},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "ni:f:o:l:b", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'n': nosys = true; break;
    case 'i': infolder = optarg; break;
    case 'f': infile = optarg; break;
    case 'o': outfolder = optarg; break;
    case 'l': luminosity = optarg; break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "fake_PU"){
        fake_PU = true;
      }else if(optname == "incl_nonbb"){
        incl_nonbb = true;
      }else if(optname == "incl_nonhh"){
        incl_nonhh = true;
      }else if(optname == "bf"){
        bf = atof(optarg);
      }else if(optname == "old"){
        old_cards = true;
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default: printf("Bad option! getopt_long returned character code 0%o\n", opt); break;
    }
  }
}
