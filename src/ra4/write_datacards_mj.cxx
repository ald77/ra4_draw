#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>
#include <ctime>
#include <algorithm>
#include <unistd.h> // getopt in Macs
#include <getopt.h>

#include "TSystem.h"
#include "TString.h"
#include "TError.h" // Controls error level reporting

#include "core/named_func.hpp"
#include "core/baby_full.hpp"
#include "core/utilities.hpp"

using namespace std;
namespace {
  bool fake_PU = false;
  TString luminosity = "35.9";
  TString nom_wgt = "weight*eff_trig"; // nominal weight to use
  enum SysType {kConst, kWeight, kSmear, kCorr, kMetSwap, kPU};
  TString syst = "all";

  const vector<double> v_data_npv{6.540e-06, 2.294e-05, 6.322e-05, 8.558e-05, 1.226e-04, 1.642e-04, 1.917e-04, 3.531e-04, 9.657e-04, 2.155e-03, 4.846e-03, 9.862e-03, 1.651e-02, 2.401e-02, 3.217e-02, 4.078e-02, 4.818e-02, 5.324e-02, 5.612e-02, 5.756e-02, 5.841e-02, 5.886e-02, 5.831e-02, 5.649e-02, 5.376e-02, 5.044e-02, 4.667e-02, 4.257e-02, 3.833e-02, 3.406e-02, 2.982e-02, 2.567e-02, 2.169e-02, 1.799e-02, 1.464e-02, 1.170e-02, 9.178e-03, 7.058e-03, 5.306e-03, 3.884e-03, 2.757e-03, 1.890e-03, 1.247e-03, 7.901e-04, 4.795e-04, 2.783e-04, 1.544e-04, 8.181e-05, 4.141e-05, 2.004e-05, 9.307e-06, 4.178e-06, 1.846e-06, 8.350e-07, 4.150e-07, 2.458e-07, 1.779e-07, 1.488e-07, 1.339e-07, 1.238e-07, 1.153e-07, 1.071e-07, 9.899e-08, 9.095e-08, 8.301e-08, 7.527e-08, 6.778e-08, 6.063e-08, 5.387e-08, 4.753e-08, 4.166e-08, 3.627e-08, 3.136e-08, 2.693e-08, 2.297e-08};
  TH1D h_data_npv("h_data_npv", "Data;N_{PV};P(N_{PV})", v_data_npv.size(), -0.5, v_data_npv.size()-0.5);
  TH1D h_mc_npv("h_mc_npv", "MC;N_{PV};P(N_{PV})", v_data_npv.size(), -0.5, v_data_npv.size()-0.5);
  double pu_low = 0.;
  double pu_high = 0.;
  bool do_syst = true;
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
void GetOptions(int argc, char *argv[], TString &infolder, TString &outfolder, TString &infile);
void fillTtbarSys(ofstream &fsys);

vector<double> getYields(Baby_full &baby, const NamedFunc &baseline, const vector<NamedFunc> &bincuts,
                         vector<double> &yield, vector<double> &w2, double lumi,
                         bool do_trig = false, const TString &flag = "");

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  time_t begtime, endtime;
  time(&begtime);
  TString infolder(""), outfolder(""), infile("");
  GetOptions(argc, argv, infolder, outfolder, infile);
  gSystem->mkdir(outfolder, kTRUE);

  // TString infile = "/cms2r0/babymaker/babies/2015_11_27/sms/split_sms/renorm/baby_SMS-T1tttt_mGluino-1500_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15FSPremix-MCRUN2_74_V9_renorm.root";
  string prs = infile.Data();
  int mglu = stoi(prs.substr(prs.find("ino-")+4,prs.find("_mLSP")-prs.find("ino-")-4));
  int mlsp = stoi(prs.substr(prs.find("LSP-")+4,prs.find("_Tune")-prs.find("LSP-")-4));
  cout<<"Working on: mGluino = "<<mglu<<" mLSP = "<<mlsp<<endl;
  string glu_lsp("mGluino-"+to_string(mglu)+"_mLSP-"+to_string(mlsp));
  string model = "T1tttt";
  if(Contains(prs, "T5tttt")) model = "T5tttt";
  if(Contains(prs, "T5tttt-Stop")) model = "T5tttt-Stop";
  if(Contains(prs, "T5tttt-degen")) model = "T5tttt-degen";
  if(Contains(prs, "T2tt")) model = "T2tt";
  if(Contains(prs, "T6ttWW")) model = "T6ttWW";


  //// tables has a vector of the tables you want to print
  TString baseline("st>500 && met>200 && mj14>250 && njets>=6 && nbm>=1 && nleps==1 && nveto==0");
  vector<bindef> v_bins; vector<int> observed;

  observed.push_back(4761); v_bins.push_back(bindef("r1_lowmet_allnb",      "met<=350 && mt<=140 && mj14<=400"));
  observed.push_back(1353); v_bins.push_back(bindef("r2_lowmet_lownj_1b",   "met<=350 && mt<=140 && mj14>400 && nbm==1 && njets<=8"));
  observed.push_back(127);  v_bins.push_back(bindef("r2_lowmet_highnj_1b",  "met<=350 && mt<=140 && mj14>400 && nbm==1 && njets>=9"));
  observed.push_back(910);  v_bins.push_back(bindef("r2_lowmet_lownj_2b",   "met<=350 && mt<=140 && mj14>400 && nbm==2 && njets<=8"));
  observed.push_back(118);  v_bins.push_back(bindef("r2_lowmet_highnj_2b",  "met<=350 && mt<=140 && mj14>400 && nbm==2 && njets>=9"));
  observed.push_back(214);  v_bins.push_back(bindef("r2_lowmet_lownj_3b",   "met<=350 && mt<=140 && mj14>400 && nbm>=3 && njets<=8"));
  observed.push_back(33);   v_bins.push_back(bindef("r2_lowmet_highnj_3b",  "met<=350 && mt<=140 && mj14>400 && nbm>=3 && njets>=9"));
  observed.push_back(247);  v_bins.push_back(bindef("r3_lowmet_allnb",      "met<=350 && mt>140  && mj14<=400"));
  observed.push_back(106);  v_bins.push_back(bindef("r4_lowmet_lownj_1b",   "met<=350 && mt>140  && mj14>400 && nbm==1 && njets<=8"));
  observed.push_back(11);   v_bins.push_back(bindef("r4_lowmet_highnj_1b",  "met<=350 && mt>140  && mj14>400 && nbm==1 && njets>=9"));
  observed.push_back(75);   v_bins.push_back(bindef("r4_lowmet_lownj_2b",   "met<=350 && mt>140  && mj14>400 && nbm==2 && njets<=8"));
  observed.push_back(11);   v_bins.push_back(bindef("r4_lowmet_highnj_2b",  "met<=350 && mt>140  && mj14>400 && nbm==2 && njets>=9"));
  observed.push_back(16);   v_bins.push_back(bindef("r4_lowmet_lownj_3b",   "met<=350 && mt>140  && mj14>400 && nbm>=3 && njets<=8"));
  observed.push_back(2);    v_bins.push_back(bindef("r4_lowmet_highnj_3b",  "met<=350 && mt>140  && mj14>400 && nbm>=3 && njets>=9"));

  observed.push_back(412);  v_bins.push_back(bindef("r1_medmet_allnb",      "met>350&&met<=500 && mt<=140 && mj14<=400"));
  observed.push_back(226);  v_bins.push_back(bindef("r2_medmet_lownj_1b",   "met>350&&met<=500 && mt<=140 && mj14>400 && nbm==1 && njets<=8"));
  observed.push_back(15);   v_bins.push_back(bindef("r2_medmet_highnj_1b",  "met>350&&met<=500 && mt<=140 && mj14>400 && nbm==1 && njets>=9"));
  observed.push_back(155);  v_bins.push_back(bindef("r2_medmet_lownj_2b",   "met>350&&met<=500 && mt<=140 && mj14>400 && nbm==2 && njets<=8"));
  observed.push_back(25);   v_bins.push_back(bindef("r2_medmet_highnj_2b",  "met>350&&met<=500 && mt<=140 && mj14>400 && nbm==2 && njets>=9"));
  observed.push_back(37);   v_bins.push_back(bindef("r2_medmet_lownj_3b",   "met>350&&met<=500 && mt<=140 && mj14>400 && nbm>=3 && njets<=8"));
  observed.push_back(7);    v_bins.push_back(bindef("r2_medmet_highnj_3b",  "met>350&&met<=500 && mt<=140 && mj14>400 && nbm>=3 && njets>=9"));
  observed.push_back(32);   v_bins.push_back(bindef("r3_medmet_allnb",      "met>350&&met<=500 && mt>140  && mj14<=400"));
  observed.push_back(25);   v_bins.push_back(bindef("r4_medmet_lownj_1b",   "met>350&&met<=500 && mt>140  && mj14>400 && nbm==1 && njets<=8"));
  observed.push_back(2);    v_bins.push_back(bindef("r4_medmet_highnj_1b",  "met>350&&met<=500 && mt>140  && mj14>400 && nbm==1 && njets>=9"));
  observed.push_back(10);   v_bins.push_back(bindef("r4_medmet_lownj_2b",   "met>350&&met<=500 && mt>140  && mj14>400 && nbm==2 && njets<=8"));
  observed.push_back(2);    v_bins.push_back(bindef("r4_medmet_highnj_2b",  "met>350&&met<=500 && mt>140  && mj14>400 && nbm==2 && njets>=9"));
  observed.push_back(1);    v_bins.push_back(bindef("r4_medmet_lownj_3b",   "met>350&&met<=500 && mt>140  && mj14>400 && nbm>=3 && njets<=8"));
  observed.push_back(0);    v_bins.push_back(bindef("r4_medmet_highnj_3b",  "met>350&&met<=500 && mt>140  && mj14>400 && nbm>=3 && njets>=9"));
  
  observed.push_back(74);   v_bins.push_back(bindef("r1_highmet_allnb",      "met>500 && mt<=140 && mj14<=400"));
  observed.push_back(71);   v_bins.push_back(bindef("r2_highmet_lownj_1b",   "met>500 && mt<=140 && mj14>400 && nbm==1 && njets<=8"));
  observed.push_back(8);    v_bins.push_back(bindef("r2_highmet_highnj_1b",  "met>500 && mt<=140 && mj14>400 && nbm==1 && njets>=9"));
  observed.push_back(32);   v_bins.push_back(bindef("r2_highmet_lownj_2b",   "met>500 && mt<=140 && mj14>400 && nbm==2 && njets<=8"));
  observed.push_back(5);    v_bins.push_back(bindef("r2_highmet_highnj_2b",  "met>500 && mt<=140 && mj14>400 && nbm==2 && njets>=9"));
  observed.push_back(10);   v_bins.push_back(bindef("r2_highmet_lownj_3b",   "met>500 && mt<=140 && mj14>400 && nbm>=3 && njets<=8"));
  observed.push_back(1);    v_bins.push_back(bindef("r2_highmet_highnj_3b",  "met>500 && mt<=140 && mj14>400 && nbm>=3 && njets>=9"));
  observed.push_back(2);    v_bins.push_back(bindef("r3_highmet_allnb",      "met>500 && mt>140  && mj14<=400"));
  observed.push_back(8);    v_bins.push_back(bindef("r4_highmet_lownj_1b",   "met>500 && mt>140  && mj14>400 && nbm==1 && njets<=8"));
  observed.push_back(2);    v_bins.push_back(bindef("r4_highmet_highnj_1b",  "met>500 && mt>140  && mj14>400 && nbm==1 && njets>=9"));
  observed.push_back(4);    v_bins.push_back(bindef("r4_highmet_lownj_2b",   "met>500 && mt>140  && mj14>400 && nbm==2 && njets<=8"));
  observed.push_back(0);    v_bins.push_back(bindef("r4_highmet_highnj_2b",  "met>500 && mt>140  && mj14>400 && nbm==2 && njets>=9"));
  observed.push_back(1);    v_bins.push_back(bindef("r4_highmet_lownj_3b",   "met>500 && mt>140  && mj14>400 && nbm>=3 && njets<=8"));
  observed.push_back(0);    v_bins.push_back(bindef("r4_highmet_highnj_3b",  "met>500 && mt>140  && mj14>400 && nbm>=3 && njets>=9"));

  unsigned nbins(v_bins.size());


  vector<sysdef> v_sys;
  // order as they will appear in latex table
  // *Nominal must stay in the first spot!!* (will be skipped in table)
  v_sys.push_back(sysdef("Nominal", "nominal", kWeight)); 
  v_sys.back().v_wgts.push_back("1.");

  v_sys.push_back(sysdef("Gen vs reco MET FS", "fs_genmet",kMetSwap));

  if (do_syst) {
    v_sys.push_back(sysdef("Lepton efficiency", "lepeff", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_lep["+to_string(i)+"]/w_lep");
    v_sys.push_back(sysdef("Lepton efficiency FS", "fs_lepeff", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_fs_lep["+to_string(i)+"]/w_fs_lep");
    v_sys.push_back(sysdef("Trigger efficiency", "trig", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_trig["+to_string(i)+"]/eff_trig"); 
    v_sys.push_back(sysdef("B-tag efficiency", "bctag", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_bctag["+to_string(i)+"]/w_btag");
    v_sys.push_back(sysdef("B-tag efficiency FS", "fs_bctag", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_fs_bctag["+to_string(i)+"]/w_btag");
    v_sys.push_back(sysdef("Mistag efficiency", "udsgtag", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_udsgtag["+to_string(i)+"]/w_btag");
    v_sys.push_back(sysdef("Mistag efficiency FS", "fs_udsgtag",kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_fs_udsgtag["+to_string(i)+"]/w_btag");
    v_sys.push_back(sysdef("Jet energy corrections", "jec", kCorr));
    v_sys.back().shift_index = 1; // JEC Up index in sys_met, etc.
    // v_sys.push_back(sysdef("Jet energy resolution", "jer", kSmear));
    // v_sys.back().shift_index = 0; // JER index in sys_met, etc.
    // v_sys.push_back(sysdef("PDFs", "pdf", kWeight));
    // for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_pdf["+to_string(i)+"]");
    // v_sys.push_back(sysdef("RMS PDFs", "rms_pdf", kWeight));
    // for (size_t i = 0; i<100; ++i) v_sys.back().v_wgts.push_back("w_pdf["+to_string(i)+"]");
    v_sys.push_back(sysdef("QCD scales", "murf",kWeight));
    for (size_t i = 0; i<2; ++i) {
      v_sys.back().v_wgts.push_back("sys_mur["+to_string(i)+"]");
      v_sys.back().v_wgts.push_back("sys_muf["+to_string(i)+"]");
      v_sys.back().v_wgts.push_back("sys_murf["+to_string(i)+"]");
    }
    v_sys.push_back(sysdef("ISR", "isr", kWeight));
    for (size_t i = 0; i<2; ++i) v_sys.back().v_wgts.push_back("sys_isr["+to_string(i)+"]/w_isr");
    v_sys.push_back(sysdef("Jet ID FS", "jetid", kConst));
    v_sys.back().v_wgts.push_back("0.01");
    v_sys.push_back(sysdef("Pile up", "pu", kPU));
    v_sys.push_back(sysdef("Luminosity", "lumi", kConst));
    v_sys.back().v_wgts.push_back("0.026");
  } 

  /////////////////////////////  No more changes needed down here to add systematics ///////////////////////
  // prepare the vector of bincuts used to get the yields
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
          bcuts.emplace_back("("+baseline+"&&"+bin.cut+")*"+nom_wgt+"*"+wgt);
        }
      }
    } else if (sys.sys_type == kCorr || sys.sys_type == kSmear) {
      for (auto &bin: v_bins) {
        bcuts.emplace_back(nom2sys_bin("("+baseline+"&&"+bin.cut, sys.shift_index)+")*"+nom_wgt);
        if (sys.sys_type == kCorr) { //if it is a correction, need to push the 'down' variation as well
          bcuts.emplace_back(nom2sys_bin("("+baseline+"&&"+bin.cut, sys.shift_index+1)+")*"+nom_wgt);
        }
      }
    } else if (sys.sys_type == kMetSwap){
      for (auto &bin: v_bins) {
        bcuts.emplace_back(nom2genmet("("+baseline+"&&"+bin.cut)+")*"+nom_wgt);
      }
    } else if (sys.sys_type == kPU) {
      for(const auto &bin: v_bins){
        bcuts.emplace_back("("+baseline+"&&"+bin.cut+"&&"+"ntrupv<=20)*"+nom_wgt);
        bcuts.emplace_back("(ntrupv<=20)*"+nom_wgt);
        bcuts.emplace_back("("+baseline+"&&"+bin.cut+"&&"+"ntrupv>=21)*"+nom_wgt);
        bcuts.emplace_back("(ntrupv>=21)*"+nom_wgt);
      }
    }
  }
  
  // get yields from the baby for all the cut strings
  Baby_full baby(std::set<std::string>{(infolder+"/"+infile).Data()});
  auto activator = baby.Activate();
  vector<double> yields, w2, entries;
  entries = getYields(baby, baseline, bcuts, yields, w2, luminosity.Atof());

  // calculate average of yields with GEN and RECO MET
  vector<float> nom_met_avg_yield, nom_met_avg_w2;
  for (auto &sys: v_sys) {
    if (sys.sys_type == kMetSwap) {
      for (size_t ibin = 0; ibin<nbins; ++ibin) {
        nom_met_avg_yield.push_back(0.5*(yields[sys.ind + ibin]+yields[ibin]));
        nom_met_avg_w2.push_back(0.5*(w2[sys.ind + ibin]+w2[ibin]));
      }
    }
  }

  // --------- Writing datacard -----------------------
  TString outpath = outfolder+"/datacard_SMS-"+TString(model)+"_"+glu_lsp+"_"+luminosity+"ifb.txt";
  if (!do_syst)  outpath.ReplaceAll(".txt","_nosys.txt");
  cout<<"open "<<outpath<<endl;
  unsigned wname(28), wdist(7), wbin(21);
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
  for (size_t ibin(0); ibin<nbins; ibin++) fcard<<setw(wbin)<<Form("%.2f",nom_met_avg_yield[ibin])<<setw(wbin)<<"1";
  fcard<<endl<<endl;
  cout<<"Wrote headers"<<endl;

  //--------- Signal statistical uncertainties ----------------------------
  for (size_t ibin(0); ibin<nbins; ibin++) {
    fcard<<setw(wname)<<"statsig_"+v_bins[ibin].tag<<setw(wdist)<<"lnN";
    TString sig_stat = Form("%.2f",1.+sqrt(nom_met_avg_w2[ibin])/nom_met_avg_yield[ibin]);
    if (yields[ibin]<0.00001) sig_stat = "99999.00";
    for (size_t jbin(0); jbin<nbins; jbin++) {
      if (ibin==jbin) fcard<<setw(wbin)<<sig_stat<<setw(wbin)<<"-";
      else fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-";
    }
    fcard<<endl;
  }
  cout<<"Wrote signal stat. uncertainties"<<endl;

  // ------------ Closure uncertainties
  fcard<<endl<<setw(wname)<<"dilep_lownj"<<setw(wdist)<<"lnN";
  for (size_t ibin(0); ibin<nbins; ibin++) {
    if(v_bins[ibin].tag.Contains("r2_") && v_bins[ibin].tag.Contains("lownj")) fcard<<setw(wbin)<<"-"<<setw(wbin)<<"1.06";
    else fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-";
  }
  fcard<<endl<<setw(wname)<<"dilep_highnj"<<setw(wdist)<<"lnN";
  for (size_t ibin(0); ibin<nbins; ibin++) {
    if(v_bins[ibin].tag.Contains("r2_") && v_bins[ibin].tag.Contains("highnj")) fcard<<setw(wbin)<<"-"<<setw(wbin)<<"1.16";
    else fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-";
  }
  fcard<<endl<<setw(wname)<<"fivejet_lowmet"<<setw(wdist)<<"lnN";
  for (size_t ibin(0); ibin<nbins; ibin++) {
    if(v_bins[ibin].tag.Contains("r2_") && v_bins[ibin].tag.Contains("lowmet")) fcard<<setw(wbin)<<"-"<<setw(wbin)<<"1.15";
    else fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-";
  }
  fcard<<endl<<setw(wname)<<"fivejet_highmet"<<setw(wdist)<<"lnN";
  for (size_t ibin(0); ibin<nbins; ibin++) {
    if(v_bins[ibin].tag.Contains("r2_") && !v_bins[ibin].tag.Contains("lowmet")) fcard<<setw(wbin)<<"-"<<setw(wbin)<<"3.00";
    else fcard<<setw(wbin)<<"-"<<setw(wbin)<<"-";
  }
  cout<<"Wrote CR-based closure uncertainties"<<endl;

  //calculate uncertainties and write results to three files
  ofstream fsys(outpath.ReplaceAll("datacard_","sys_"));
  cout<<"Writing to "<<outpath<<endl;
  fillTtbarSys(fsys);
  // ofstream fsysrms(outpath.ReplaceAll("sys_","sysrms_"));
  // fillTtbarSys(fsysrms);
  ofstream fsysdbg(outpath.ReplaceAll("sys_","sysdbg_"));
  ofstream fsysent(outpath.ReplaceAll("sysdbg_","sysent_"));
  for (auto &sys: v_sys) {
    if (sys.tag != "nominal") {
      if (sys.tag != "rms_pdf") fsys<<"\nSYSTEMATIC "<<sys.tag<<"\n  PROCESSES signal\n";
      // if (sys.tag != "pdf") fsysrms<<"\nSYSTEMATIC "<<sys.tag<<"\n  PROCESSES signal\n";
      fsysdbg<<"\nSYSTEMATIC "<<sys.tag<<"\n  PROCESSES signal\n";
      fcard<<setw(wname)<<sys.tag<<setw(wdist)<<"lnN";
    }
    for (size_t ibin = 0; ibin<nbins; ++ibin) {
      const double nom_yield(yields[ibin]);
      size_t nwgts = sys.v_wgts.size();
      double up(0.), dn(0.); 
      if (sys.sys_type == kConst) {
        up = stod(sys.v_wgts[0].Data());
        dn = -up;
      } else if (sys.sys_type == kWeight) {
        if (sys.tag == "nominal") {
          if (ibin==0) fsysent<<fixed<<"SYSTEMATIC "<<sys.tag<<"\n  PROCESSES signal\n";
          fsysent <<"    " <<std::left<<setw(25)<<v_bins[ibin].tag <<" "<<std::right<<setprecision(0)<<setw(25)<<entries[ibin] <<" "<<setprecision(5)<<setw(15)<<yields[ibin] <<" "<<setprecision(10)<<setw(15)<<w2[ibin] <<endl;
          continue;
        } else if (sys.tag == "rms_pdf") { 
          double sumw2(0), mean(0);
          for (size_t iwgt = 0; iwgt<nwgts; ++iwgt) {
            sumw2 += pow(yields[sys.ind + nwgts*ibin + iwgt],2);
            mean += yields[sys.ind + nwgts*ibin + iwgt];
          }
          mean = mean/nwgts;
          up = sqrt((sumw2-nwgts*pow(mean,2))/(nwgts-1))/nom_yield;  // RMS
          dn = -up;
        } else if (sys.tag == "murf") {
          up = *max_element(yields.begin() + sys.ind + nwgts*ibin, yields.begin() + sys.ind + nwgts*(ibin+1))/nom_yield - 1; //max of all weights mur_up, muf_up and murf_up
          dn = *min_element(yields.begin() + sys.ind + nwgts*ibin, yields.begin() + sys.ind + nwgts*(ibin+1))/nom_yield - 1; //min of all weights mur_up, muf_up and murf_up
        } else {
          up = yields[sys.ind + 2*ibin]/nom_yield - 1;
          dn = yields[sys.ind + 2*ibin + 1]/nom_yield - 1;
        }
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
            dn = 0.1;
          }else{
            up = 0.15;
            dn = 0.15;
          }
        }
      }
      // convert to ra4_stats input and write to file
      double ln = (up>0 ? 1:-1)*max(up>0 ? up : (1/(1+up)-1), dn>0 ? dn : (1/(1+dn)-1));
      if (sys.sys_type == kConst) ln = up;
      if (sys.tag !="rms_pdf") {
        if(sys.tag.Contains("trig")){
          fsys<<"    " <<std::left<<setw(25)<<v_bins[ibin].tag <<" "<<std::right<<setw(10)<<Form("%.3f",ln) <<endl;
        } else {
          fsys<<"    " <<std::left<<setw(25)<<v_bins[ibin].tag <<" "<<std::right<<setw(10)<<Form("%.2f",ln) <<endl;
        }
      }
      // if (sys.tag !="pdf") 
      //   fsysrms<<"    " <<std::left<<setw(25)<<v_bins[ibin].tag <<" "<<std::right<<setw(10)<<Form("%.2f",ln) <<endl;
      if (sys.sys_type == kMetSwap) {
        fsysdbg <<"    " <<std::left<<setw(25)<<v_bins[ibin].tag <<" "<<"mg="<<setw(5)<<mglu <<" "<<"mlsp="<<setw(10)<<mlsp <<" "<<std::right<<setw(10)<<Form("%.2f",up) <<" "<<setw(10)<<Form("%.2f",dn) <<" "<<std::right<<setw(10)<<"Nominal: "<< nom_yield<<" "<<setw(10)<<"met_tru: "<<yields[sys.ind + ibin]<<endl;
      } else {
        fsysdbg <<"    " <<std::left<<setw(25)<<v_bins[ibin].tag <<" "<<"mg="<<setw(5)<<mglu <<" "<<"mlsp="<<setw(10)<<mlsp <<" "<<std::right<<setw(10)<<Form("%.2f",up) <<" "<<setw(10)<<Form("%.2f",dn)<<endl;
      }

      // write systematics to datacard
      ln = max(up>0 ? 1+up : 1/(up+1), dn>0 ? 1+dn : 1/(dn+1));
      if (std::isnan(ln) || std::isinf(ln)) {
        cout <<" Found bad unc. set to 0 -> "<<std::left<<setw(10)<<sys.tag <<std::left<<setw(10)<<v_bins[ibin].tag <<" "<<std::right<<setprecision(0)<<setw(25)<<entries[ibin] <<" "<<setprecision(5)<<setw(15)<<yields[ibin] <<" "<<setprecision(10)<<setw(15)<<w2[ibin] <<endl;  
        ln = 0;
      } 
      if (sys.sys_type == kConst) ln = 1+up;
      fcard<<setw(wbin)<<Form("%.2f",ln)<<setw(wbin)<<"-";
    } // loop over bins
    fcard<<endl;
  } // loop over systematics
  fsys.close();
  // fsysrms.close();
  fsysdbg.close();
  fsysent.close();

  fcard<<"rp_r1_lowmet_allnb      rateParam r1_lowmet_allnb      bkg  4761"<<endl;
  fcard<<"rp_r2_lowmet_lownj_1b   rateParam r2_lowmet_lownj_1b   bkg  1353"<<endl;
  fcard<<"rp_r2_lowmet_highnj_1b  rateParam r2_lowmet_highnj_1b  bkg  127"<<endl;
  fcard<<"rp_r2_lowmet_lownj_2b   rateParam r2_lowmet_lownj_2b   bkg  910"<<endl;
  fcard<<"rp_r2_lowmet_highnj_2b  rateParam r2_lowmet_highnj_2b  bkg  118"<<endl;
  fcard<<"rp_r2_lowmet_lownj_3b   rateParam r2_lowmet_lownj_3b   bkg  214"<<endl;
  fcard<<"rp_r2_lowmet_highnj_3b  rateParam r2_lowmet_highnj_3b  bkg  33"<<endl;
  fcard<<"rp_r3_lowmet_allnb      rateParam r3_lowmet_allnb      bkg  247"<<endl;
  fcard<<"rp_r4_lowmet_lownj_1b   rateParam r4_lowmet_lownj_1b   bkg (@0*@1/@2)*@3 rp_r3_lowmet_allnb,rp_r2_lowmet_lownj_1b,rp_r1_lowmet_allnb,kappa_lowmet_lownj_1b"<<endl;
  fcard<<"rp_r4_lowmet_highnj_1b  rateParam r4_lowmet_highnj_1b  bkg (@0*@1/@2)*@3 rp_r3_lowmet_allnb,rp_r2_lowmet_highnj_1b,rp_r1_lowmet_allnb,kappa_lowmet_highnj_1b"<<endl;
  fcard<<"rp_r4_lowmet_lownj_2b   rateParam r4_lowmet_lownj_2b   bkg (@0*@1/@2)*@3 rp_r3_lowmet_allnb,rp_r2_lowmet_lownj_2b,rp_r1_lowmet_allnb,kappa_lowmet_lownj_2b"<<endl;
  fcard<<"rp_r4_lowmet_highnj_2b  rateParam r4_lowmet_highnj_2b  bkg (@0*@1/@2)*@3 rp_r3_lowmet_allnb,rp_r2_lowmet_highnj_2b,rp_r1_lowmet_allnb,kappa_lowmet_highnj_2b"<<endl;
  fcard<<"rp_r4_lowmet_lownj_3b   rateParam r4_lowmet_lownj_3b   bkg (@0*@1/@2)*@3 rp_r3_lowmet_allnb,rp_r2_lowmet_lownj_3b,rp_r1_lowmet_allnb,kappa_lowmet_lownj_3b"<<endl;
  fcard<<"rp_r4_lowmet_highnj_3b  rateParam r4_lowmet_highnj_3b  bkg (@0*@1/@2)*@3 rp_r3_lowmet_allnb,rp_r2_lowmet_highnj_3b,rp_r1_lowmet_allnb,kappa_lowmet_highnj_3b"<<endl;
  fcard<<"rp_r1_medmet_allnb      rateParam r1_medmet_allnb      bkg  412"<<endl;
  fcard<<"rp_r2_medmet_lownj_1b   rateParam r2_medmet_lownj_1b   bkg  226"<<endl;
  fcard<<"rp_r2_medmet_highnj_1b  rateParam r2_medmet_highnj_1b  bkg  15  [0,100]"<<endl;
  fcard<<"rp_r2_medmet_lownj_2b   rateParam r2_medmet_lownj_2b   bkg  155"<<endl;
  fcard<<"rp_r2_medmet_highnj_2b  rateParam r2_medmet_highnj_2b  bkg  25  [0,100]"<<endl;
  fcard<<"rp_r2_medmet_lownj_3b   rateParam r2_medmet_lownj_3b   bkg  37  [0,100]"<<endl;
  fcard<<"rp_r2_medmet_highnj_3b  rateParam r2_medmet_highnj_3b  bkg  7   [0,100]"<<endl;
  fcard<<"rp_r3_medmet_allnb      rateParam r3_medmet_allnb      bkg  32  [0,100]"<<endl;
  fcard<<"rp_r4_medmet_lownj_1b   rateParam r4_medmet_lownj_1b   bkg (@0*@1/@2)*@3 rp_r3_medmet_allnb,rp_r2_medmet_lownj_1b,rp_r1_medmet_allnb,kappa_medmet_lownj_1b"<<endl;
  fcard<<"rp_r4_medmet_highnj_1b  rateParam r4_medmet_highnj_1b  bkg (@0*@1/@2)*@3 rp_r3_medmet_allnb,rp_r2_medmet_highnj_1b,rp_r1_medmet_allnb,kappa_medmet_highnj_1b"<<endl;
  fcard<<"rp_r4_medmet_lownj_2b   rateParam r4_medmet_lownj_2b   bkg (@0*@1/@2)*@3 rp_r3_medmet_allnb,rp_r2_medmet_lownj_2b,rp_r1_medmet_allnb,kappa_medmet_lownj_2b"<<endl;
  fcard<<"rp_r4_medmet_highnj_2b  rateParam r4_medmet_highnj_2b  bkg (@0*@1/@2)*@3 rp_r3_medmet_allnb,rp_r2_medmet_highnj_2b,rp_r1_medmet_allnb,kappa_medmet_highnj_2b"<<endl;
  fcard<<"rp_r4_medmet_lownj_3b   rateParam r4_medmet_lownj_3b   bkg (@0*@1/@2)*@3 rp_r3_medmet_allnb,rp_r2_medmet_lownj_3b,rp_r1_medmet_allnb,kappa_medmet_lownj_3b"<<endl;
  fcard<<"rp_r4_medmet_highnj_3b  rateParam r4_medmet_highnj_3b  bkg (@0*@1/@2)*@3 rp_r3_medmet_allnb,rp_r2_medmet_highnj_3b,rp_r1_medmet_allnb,kappa_medmet_highnj_3b"<<endl;
  fcard<<"rp_r1_highmet_allnb     rateParam r1_highmet_allnb     bkg 74"<<endl;
  fcard<<"rp_r2_highmet_lownj_1b  rateParam r2_highmet_lownj_1b  bkg 71"<<endl;
  fcard<<"rp_r2_highmet_highnj_1b rateParam r2_highmet_highnj_1b bkg 8  [0,100]"<<endl;
  fcard<<"rp_r2_highmet_lownj_2b  rateParam r2_highmet_lownj_2b  bkg 32 [0,100]"<<endl;
  fcard<<"rp_r2_highmet_highnj_2b rateParam r2_highmet_highnj_2b bkg 5  [0,100]"<<endl;
  fcard<<"rp_r2_highmet_lownj_3b  rateParam r2_highmet_lownj_3b  bkg 10 [0,100]"<<endl;
  fcard<<"rp_r2_highmet_highnj_3b rateParam r2_highmet_highnj_3b bkg 1  [0,100]"<<endl;
  fcard<<"rp_r3_highmet_allnb     rateParam r3_highmet_allnb     bkg 2  [0,100]"<<endl;
  fcard<<"rp_r4_highmet_lownj_1b  rateParam r4_highmet_lownj_1b  bkg (@0*@1/@2)*@3 rp_r3_highmet_allnb,rp_r2_highmet_lownj_1b,rp_r1_highmet_allnb,kappa_highmet_lownj_1b"<<endl;
  fcard<<"rp_r4_highmet_highnj_1b rateParam r4_highmet_highnj_1b bkg (@0*@1/@2)*@3 rp_r3_highmet_allnb,rp_r2_highmet_highnj_1b,rp_r1_highmet_allnb,kappa_highmet_highnj_1b"<<endl;
  fcard<<"rp_r4_highmet_lownj_2b  rateParam r4_highmet_lownj_2b  bkg (@0*@1/@2)*@3 rp_r3_highmet_allnb,rp_r2_highmet_lownj_2b,rp_r1_highmet_allnb,kappa_highmet_lownj_2b"<<endl;
  fcard<<"rp_r4_highmet_highnj_2b rateParam r4_highmet_highnj_2b bkg (@0*@1/@2)*@3 rp_r3_highmet_allnb,rp_r2_highmet_highnj_2b,rp_r1_highmet_allnb,kappa_highmet_highnj_2b"<<endl;
  fcard<<"rp_r4_highmet_lownj_3b  rateParam r4_highmet_lownj_3b  bkg (@0*@1/@2)*@3 rp_r3_highmet_allnb,rp_r2_highmet_lownj_3b,rp_r1_highmet_allnb,kappa_highmet_lownj_3b"<<endl;
  fcard<<"rp_r4_highmet_highnj_3b rateParam r4_highmet_highnj_3b bkg (@0*@1/@2)*@3 rp_r3_highmet_allnb,rp_r2_highmet_highnj_3b,rp_r1_highmet_allnb,kappa_highmet_highnj_3b"<<endl;
  fcard<<"kappa_lowmet_lownj_1b   param 1.21 0.04"<<endl;
  fcard<<"kappa_lowmet_highnj_1b  param 0.98 0.06"<<endl;
  fcard<<"kappa_lowmet_lownj_2b   param 1.17 0.03"<<endl;
  fcard<<"kappa_lowmet_highnj_2b  param 1.24 0.12"<<endl;
  fcard<<"kappa_lowmet_lownj_3b   param 1.49 0.07"<<endl;
  fcard<<"kappa_lowmet_highnj_3b  param 1.37 0.11"<<endl;
  fcard<<"kappa_medmet_lownj_1b   param 0.99 0.05"<<endl;
  fcard<<"kappa_medmet_highnj_1b  param 1.07 0.14"<<endl;
  fcard<<"kappa_medmet_lownj_2b   param 1.13 0.09"<<endl;
  fcard<<"kappa_medmet_highnj_2b  param 0.84 0.13"<<endl;
  fcard<<"kappa_medmet_lownj_3b   param 1.33 0.12"<<endl;
  fcard<<"kappa_medmet_highnj_3b  param 1.16 0.23"<<endl;
  fcard<<"kappa_highmet_lownj_1b  param 0.97 0.11"<<endl;
  fcard<<"kappa_highmet_highnj_1b param 1.04 0.25"<<endl;
  fcard<<"kappa_highmet_lownj_2b  param 1.03 0.13"<<endl;
  fcard<<"kappa_highmet_highnj_2b param 0.97 0.26"<<endl;
  fcard<<"kappa_highmet_lownj_3b  param 1.53 0.30"<<endl;
  fcard<<"kappa_highmet_highnj_3b param 3.03 1.06"<<endl;
  fcard.close();

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
  return ibin;
}

TString nom2genmet(TString ibin){
  ibin.ReplaceAll("met", "met_tru");
  return ibin;

}

void GetOptions(int argc, char *argv[], TString &infolder, TString &outfolder, TString &infile){
  string blah;
  while(true){
    static struct option long_options[] = {
      {"syst", required_argument, 0, 's'},
      {"infolder", required_argument, 0, 'i'},
      {"infile", required_argument, 0, 'f'},
      {"outfolder", required_argument, 0, 'o'},
      {"lumi", required_argument, 0, 'l'},
      {"alt_bin", no_argument, 0, 'b'},
      {"fake_PU", no_argument, 0, 0},
      {"no_syst", no_argument, 0, 0},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s:i:f:o:l:b", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 's': syst = optarg; break;
    case 'i': infolder = optarg; break;
    case 'f': infile = optarg; break;
    case 'o': outfolder = optarg; break;
    case 'l': luminosity = optarg; break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "fake_PU"){
        fake_PU = true;
      }else if(optname == "no_syst"){
        do_syst = false;
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default: printf("Bad option! getopt_long returned character code 0%o\n", opt); break;
    }
  }
}

void fillTtbarSys(ofstream &fsys){
    fsys << "SYSTEMATIC dilep_lownj" << endl;
    fsys << " PROCESSES ttbar,other" << endl;
    fsys << "  r2_lowmet_lownj_1b    0.06" << endl;
    fsys << "  r2_lowmet_lownj_2b    0.06" << endl;
    fsys << "  r2_lowmet_lownj_3b    0.06" << endl;
    fsys << "  r2_medmet_lownj_1b    0.06" << endl;
    fsys << "  r2_medmet_lownj_2b    0.06" << endl;
    fsys << "  r2_medmet_lownj_3b    0.06" << endl;
    fsys << "  r2_highmet_lownj_1b   0.06" << endl;
    fsys << "  r2_highmet_lownj_2b   0.06" << endl;
    fsys << "  r2_highmet_lownj_3b   0.06" << endl;
    fsys << endl;
    fsys << "SYSTEMATIC dilep_highnj" << endl;
    fsys << " PROCESSES ttbar,other" << endl;
    fsys << "  r2_lowmet_highnj_1b   0.16" << endl;
    fsys << "  r2_lowmet_highnj_2b   0.16" << endl;
    fsys << "  r2_lowmet_highnj_3b   0.16" << endl;
    fsys << "  r2_medmet_highnj_1b   0.16" << endl;
    fsys << "  r2_medmet_highnj_2b   0.16" << endl;
    fsys << "  r2_medmet_highnj_3b   0.16" << endl;
    fsys << "  r2_highmet_highnj_1b  0.16" << endl;
    fsys << "  r2_highmet_highnj_2b  0.16" << endl;
    fsys << "  r2_highmet_highnj_3b  0.16" << endl;
    fsys << endl;
    fsys << "SYSTEMATIC fivejet_lowmet" << endl;
    fsys << " PROCESSES ttbar,other" << endl;
    fsys << "  r2_lowmet_lownj_1b    0.15" << endl;
    fsys << "  r2_lowmet_lownj_2b    0.15" << endl;
    fsys << "  r2_lowmet_lownj_3b    0.15" << endl;
    fsys << "  r2_lowmet_highnj_1b   0.15" << endl;
    fsys << "  r2_lowmet_highnj_2b   0.15" << endl;
    fsys << "  r2_lowmet_highnj_3b   0.15" << endl;
    fsys << endl;
    fsys << "SYSTEMATIC fivejet_highmet" << endl;
    fsys << " PROCESSES ttbar,other" << endl;
    fsys << "  r2_medmet_lownj_1b    0.37" << endl;
    fsys << "  r2_medmet_lownj_2b    0.37" << endl;
    fsys << "  r2_medmet_lownj_3b    0.37" << endl;
    fsys << "  r2_medmet_highnj_1b   0.37" << endl;
    fsys << "  r2_medmet_highnj_2b   0.37" << endl;
    fsys << "  r2_medmet_highnj_3b   0.37" << endl;
    fsys << "  r2_highmet_lownj_1b   0.37" << endl;
    fsys << "  r2_highmet_lownj_2b   0.37" << endl;
    fsys << "  r2_highmet_lownj_3b   0.37" << endl;
    fsys << "  r2_highmet_highnj_1b  0.37" << endl;
    fsys << "  r2_highmet_highnj_2b  0.37" << endl;
    fsys << "  r2_highmet_highnj_3b  0.37" << endl;
    fsys << endl;
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
    h_mc_npv.Fill(baby.ntrupv(), baby.weight()*baby.eff_trig());
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
