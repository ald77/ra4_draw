#include "hig/hig_functions.hpp"

#include "TVector2.h"

#include "core/utilities.hpp"
#include "core/config_parser.hpp"

using namespace std;

namespace Higfuncs{

const NamedFunc hig_pt1("hig_pt1",[](const Baby &b) -> NamedFunc::ScalarType{
    float higpt(0);
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (b.mc_id()->at(i)!=25) continue;
      if (b.mc_pt()->at(i)>higpt) higpt = b.mc_pt()->at(i);
    }
    return higpt;
});

const NamedFunc hig_pt2("hig_pt2",[](const Baby &b) -> NamedFunc::ScalarType{
    float higpt1(0), higpt2(0);
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (b.mc_id()->at(i)!=25) continue;
      if (b.mc_pt()->at(i)>higpt1) {
	higpt2 = higpt1;
	higpt1 = b.mc_pt()->at(i);
      } else if (b.mc_pt()->at(i)>higpt2) {
	higpt2 = b.mc_pt()->at(i);
      }
    }
    return higpt2;
});

const NamedFunc ntrub("ntrub",[](const Baby &b) -> NamedFunc::ScalarType{
  int tmp_ntrub(0);
  for (unsigned i(0); i<b.jets_pt()->size(); i++){
    if (!b.jets_h1()->at(i) && !b.jets_h2()->at(i)) continue;
    if (b.jets_hflavor()->at(i)==5) tmp_ntrub++;
  }
  return tmp_ntrub;
});

const NamedFunc hig_bcat("hig_bcat",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.nbt()==2 && b.nbm()==2) return 2;
  else if (b.nbt()>=2 && b.nbm()==3 && b.nbl()==3) return 3;
  else if (b.nbt()>=2 && b.nbm()>=3 && b.nbl()>=4) return 4;
  else return 0;
});

const NamedFunc higd_bcat("higd_bcat",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.nbdt()==2 && b.nbdm()==2) return 2;
  else if (b.nbdt()>=2 && b.nbdm()==3 && b.nbdl()==3) return 3;
  else if (b.nbdt()>=2 && b.nbdm()>=3 && b.nbdl()>=4) return 4;
  else return 0;
});

const NamedFunc higd_bcat_extended("higd_bcat_extended",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.nbdm()==0) return 0;
  else if (b.nbdm()==1) return 1;
  else if (b.nbdt()==2 && b.nbdm()==2) return 2;
  else if (b.nbdt()>=2 && b.nbdm()==3 && b.nbdl()==3) return 3;
  else if (b.nbdt()>=2 && b.nbdm()>=3 && b.nbdl()>=4) return 4;
  else return 99;
});

const NamedFunc higd_bcat_mmmm("higd_bcat_mmmm",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.nbdm()>=2) return min(4,b.nbdm());
  else return 0;
});
const NamedFunc higd_bcat_ttll("higd_bcat_ttll",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.nbdt()==2 && b.nbdl()==2) return 2;
  else if (b.nbdt()>=2 && b.nbdl()==3) return 3;
  else if (b.nbdt()>=2 && b.nbdl()>=4) return 4;
  else return 0;
});
const NamedFunc higd_bcat_tmml("higd_bcat_tmml",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.nbdt()>=1 && b.nbdm()==2) return 2;
  else if (b.nbdt()>=1 && b.nbdm()==3 && b.nbdl()==3) return 3;
  else if (b.nbdt()>=1 && b.nbdm()>=3 && b.nbdl()>=4) return 4;
  else return 0;
});

// apply weights found from the nb and MET data/MC comparisons
const NamedFunc wgt_comp("wgt_comp",[](const Baby &b) -> NamedFunc::ScalarType{
  float wgt = 1;
  if ( (b.type()>=1000 && b.type()<2000) ||  // ttbar
    // (b.type()>=3000 && b.type()<4000) ||     // single top
    (b.type()>=4000 && b.type()<6000) ||     // ttw and ttz
    b.type()==9000  ||                       // ttH
    b.type()==10000  ||                      // ttgamma
    b.type()==11000) {                       // tttt
    //apply normalization factor from MET distribution
    wgt = 1.052;
    if (b.met()<=50)                      wgt*=1.005;
    else if (b.met()>50 && b.met()<=100)  wgt*=1.018;
    else if (b.met()>100 && b.met()<=150) wgt*=0.976;
    else if (b.met()>150 && b.met()<=200) wgt*=0.914;
    else if (b.met()>200 && b.met()<=300) wgt*=0.854;
    else if (b.met()>300)                 wgt*=0.865;
    // nb correction from nb distribution data/mc, inclusive in MET
    if (b.nbdt()==2 && b.nbdm()==2)                     wgt*=0.992;
    else if (b.nbdt()>=2 && b.nbdm()==3 && b.nbdl()==3) wgt*=1.040;
    else if (b.nbdt()>=2 && b.nbdm()>=3 && b.nbdl()>=4) wgt*=1.165;
  } 
  else if((b.type()>=8000 && b.type()<9000) || // zjets
  (b.type()>=2000 && b.type()<3000) ||      // wjets
  (b.type()>=6000 && b.type()<7000)) {      // dyjets
    //apply normalization factor from MET distribution
    wgt = 1.416;
    if (b.met()>150 && b.met()<=200)      wgt*=1.121;
    else if (b.met()>200 && b.met()<=300) wgt*=0.951;
    else if (b.met()>300)                 wgt*=0.722;
    // nb correction from nb distribution data/mc, inclusive in MET
    if (b.nbdt()==2 && b.nbdm()==2)                     wgt*=0.985;
    else if (b.nbdt()>=2 && b.nbdm()==3 && b.nbdl()==3) wgt*=1.176;
    else if (b.nbdt()>=2 && b.nbdm()>=3 && b.nbdl()>=4) wgt*=1.097;
  } 
  else if ( (b.type()>=7000 && b.type()<8000)) { // qcd
  //apply normalization factor from MET distribution
    wgt = 1.700;
    if (b.met()>150 && b.met()<=200)      wgt*=0.927;
    else if (b.met()>200 && b.met()<=300) wgt*=1.199;
    else if (b.met()>300)                 wgt*=1.301;
    // nb correction from nb distribution data/mc, inclusive in MET
    if (b.nbdt()==2 && b.nbdm()==2)                    wgt*=0.982;
    else if (b.nbdt()>=2 && b.nbdm()==3 && b.nbdl()==3) wgt*=1.142;
    else if (b.nbdt()>=2 && b.nbdm()>=3 && b.nbdl()>=4) wgt*=1.069;
  }
  return wgt;
});

// subtract ttbar based on MC prediction reweighted to data in 1l CR
// since ttbar has to be combined in the same process def with data, 
// also apply stitch, json and trigger here
const NamedFunc wgt_subtr_ttx("wgt_subtr_ttx",[](const Baby &b) -> NamedFunc::ScalarType{
  if ( (b.type()>=1000 && b.type()<2000) ||  // ttbar
    // (b.type()>=3000 && b.type()<4000) ||     // single top
    (b.type()>=4000 && b.type()<6000) ||     // ttw and ttz
    b.type()==9000  ||                       // ttH
    b.type()==10000  ||                      // ttgamma
    b.type()==11000) {
    float wgt = -1; 
    // apply lumi 
    wgt*= 35.9;
    // apply ttx weights derived from 1l CR
    wgt = 1.052;
    if (b.met()<=50)                      wgt*=1.005;
    else if (b.met()>50 && b.met()<=100)  wgt*=1.018;
    else if (b.met()>100 && b.met()<=150) wgt*=0.976;
    else if (b.met()>150 && b.met()<=200) wgt*=0.914;
    else if (b.met()>200 && b.met()<=300) wgt*=0.854;
    else if (b.met()>300)                 wgt*=0.865;
    // nb correction from nb distribution data/mc, inclusive in MET
    if (b.nbdt()==2 && b.nbdm()==2)                     wgt*=0.992;
    else if (b.nbdt()>=2 && b.nbdm()==3 && b.nbdl()==3) wgt*=1.040;
    else if (b.nbdt()>=2 && b.nbdm()>=3 && b.nbdl()>=4) wgt*=1.165;
    // apply stitch
    if (b.stitch_met()) return wgt;
    else return 0;
  } else if (b.type()>0 && b.type()<1000){ // apply trigger and json for data
    return trig_hig_decision(b);
  }
  // for all other backgrounds, chill (they are not in the "data" process so no need to apply lumi)
  return 1;
});

// calculate effect of systematics calculated for each background 
// in the data control regions on the total bkg. kappa
const NamedFunc wgt_syst_ttx("wgt_syst_ttx",[](const Baby &b) -> NamedFunc::ScalarType{
  if ( (b.type()>=1000 && b.type()<2000) ||  // ttbar
    // (b.type()>=3000 && b.type()<4000) ||     // single top
    (b.type()>=4000 && b.type()<6000) ||     // ttw and ttz
    b.type()==9000  ||                       // ttH
    b.type()==10000  ||                      // ttgamma
    b.type()==11000) {                       // tttt
    if (b.higd_am()<=100 || (b.higd_am()>140 && b.higd_am()<=200)) {
      if (b.nbdt()>=2 && b.nbdm()==3 && b.nbdl()==3) return 0.03;
      else if (b.nbdt()>=2 && b.nbdm()>=3 && b.nbdl()>=4) return 0.06;
    }
  }
  return 0;
});

const NamedFunc wgt_syst_vjets("wgt_syst_vjets",[](const Baby &b) -> NamedFunc::ScalarType{
  if ( (b.type()>=8000 && b.type()<9000) || // zjets
    (b.type()>=2000 && b.type()<3000) ||    // wjets
    (b.type()>=6000 && b.type()<7000)) {   
    if (b.higd_am()<=100 || (b.higd_am()>140 && b.higd_am()<=200))
      if (b.nbdt()>=2 && b.nbdm()>=3) return 0.19;
  }
  return 0;
});

const NamedFunc wgt_syst_qcd("wgt_syst_qcd",[](const Baby &b) -> NamedFunc::ScalarType{
  if ( (b.type()>=7000 && b.type()<8000)) { // qcd
    if (b.higd_am()<=100 || (b.higd_am()>140 && b.higd_am()<=200))
      if (b.nbdt()>=2 && b.nbdm()>=3) return 0.13;
  }
  return 0;
});

// Definition of analysis trigger
NamedFunc::ScalarType trig_hig_decision(const Baby &b){
    bool mettrig = b.trig()->at(13)||b.trig()->at(33)||b.trig()->at(14)||b.trig()->at(15)
      ||b.trig()->at(30)||b.trig()->at(31);
    bool eltrig = b.trig()->at(22)||b.trig()->at(40)||b.trig()->at(24)||b.trig()->at(41);
    bool mutrig = b.trig()->at(19)||b.trig()->at(55)||b.trig()->at(21);

    if(b.nels()==1 && b.nmus()==0){
      if(mettrig || eltrig) return 1;
      else return -1;
    } else if(b.nels()==0 && b.nmus()==1){
      if(mettrig || mutrig) return 1;
      else return -1;
    } else if(b.nels()==2 && b.nmus()==0){
      if(eltrig) return 1;
      else return -1;
    } else if(b.nels()==0 && b.nmus()==2){
      if(mutrig) return 1;
      else return -1;
    } else if(b.nvleps()==0){
      if(mettrig) return 1;
      else return -1;
    }

    return -1;
}

const NamedFunc trig_hig("trig_hig", [](const Baby &b) -> NamedFunc::ScalarType{
  return trig_hig_decision(b);
  });
  
//// Efficiency of the MET[100||110||120] triggers in all 36.2 ifb
const NamedFunc err_higtrig("eff_higtrig", [](const Baby &b) -> NamedFunc::ScalarType{
    float errup, errdown; // Stat uncertainty. Not used, but for reference
    float uncert = 0., met = b.met(), ht = b.ht();
    errup=0;errdown=0;
    errup+=errdown;

    if(ht>   0 && ht<= 200 && met> 150 && met<= 155) {uncert = 0.072; errup = 0.013; errdown = 0.013;}
    else if(ht> 200 && ht<= 600 && met> 150 && met<= 155) {uncert = 0.071; errup = 0.005; errdown = 0.005;}
    else if(ht> 600 && ht<= 800 && met> 150 && met<= 155) {uncert = 0.075; errup = 0.023; errdown = 0.024;}
    else if(ht> 800 && ht<=1000 && met> 150 && met<= 155) {uncert = 0.083; errup = 0.042; errdown = 0.042;}
    else if(ht>1000 && ht<=9999 && met> 150 && met<= 155) {uncert = 0.089; errup = 0.052; errdown = 0.054;}
    else if(ht>   0 && ht<= 200 && met> 155 && met<= 160) {uncert = 0.072; errup = 0.014; errdown = 0.014;}
    else if(ht> 200 && ht<= 600 && met> 155 && met<= 160) {uncert = 0.071; errup = 0.005; errdown = 0.005;}
    else if(ht> 600 && ht<= 800 && met> 155 && met<= 160) {uncert = 0.074; errup = 0.022; errdown = 0.022;}
    else if(ht> 800 && ht<=1000 && met> 155 && met<= 160) {uncert = 0.083; errup = 0.042; errdown = 0.042;}
    else if(ht>1000 && ht<=9999 && met> 155 && met<= 160) {uncert = 0.091; errup = 0.057; errdown = 0.057;}
    else if(ht>   0 && ht<= 200 && met> 160 && met<= 165) {uncert = 0.047; errup = 0.016; errdown = 0.016;}
    else if(ht> 200 && ht<= 600 && met> 160 && met<= 165) {uncert = 0.044; errup = 0.006; errdown = 0.006;}
    else if(ht> 600 && ht<= 800 && met> 160 && met<= 165) {uncert = 0.050; errup = 0.022; errdown = 0.023;}
    else if(ht> 800 && ht<=1000 && met> 160 && met<= 165) {uncert = 0.060; errup = 0.039; errdown = 0.041;}
    else if(ht>1000 && ht<=9999 && met> 160 && met<= 165) {uncert = 0.066; errup = 0.048; errdown = 0.049;}
    else if(ht>   0 && ht<= 200 && met> 165 && met<= 170) {uncert = 0.047; errup = 0.017; errdown = 0.018;}
    else if(ht> 200 && ht<= 600 && met> 165 && met<= 170) {uncert = 0.044; errup = 0.006; errdown = 0.006;}
    else if(ht> 600 && ht<= 800 && met> 165 && met<= 170) {uncert = 0.050; errup = 0.022; errdown = 0.023;}
    else if(ht> 800 && ht<=1000 && met> 165 && met<= 170) {uncert = 0.062; errup = 0.042; errdown = 0.044;}
    else if(ht>1000 && ht<=9999 && met> 165 && met<= 170) {uncert = 0.077; errup = 0.058; errdown = 0.064;}
    else if(ht>   0 && ht<= 200 && met> 170 && met<= 175) {uncert = 0.048; errup = 0.019; errdown = 0.020;}
    else if(ht> 200 && ht<= 600 && met> 170 && met<= 175) {uncert = 0.044; errup = 0.005; errdown = 0.006;}
    else if(ht> 600 && ht<= 800 && met> 170 && met<= 175) {uncert = 0.050; errup = 0.022; errdown = 0.024;}
    else if(ht> 800 && ht<=1000 && met> 170 && met<= 175) {uncert = 0.063; errup = 0.041; errdown = 0.045;}
    else if(ht>1000 && ht<=9999 && met> 170 && met<= 175) {uncert = 0.075; errup = 0.056; errdown = 0.061;}
    else if(ht>   0 && ht<= 200 && met> 175 && met<= 180) {uncert = 0.049; errup = 0.020; errdown = 0.021;}
    else if(ht> 200 && ht<= 600 && met> 175 && met<= 180) {uncert = 0.044; errup = 0.005; errdown = 0.006;}
    else if(ht> 600 && ht<= 800 && met> 175 && met<= 180) {uncert = 0.050; errup = 0.022; errdown = 0.024;}
    else if(ht> 800 && ht<=1000 && met> 175 && met<= 180) {uncert = 0.062; errup = 0.037; errdown = 0.043;}
    else if(ht>1000 && ht<=9999 && met> 175 && met<= 180) {uncert = 0.076; errup = 0.055; errdown = 0.062;}
    else if(ht>   0 && ht<= 200 && met> 180 && met<= 185) {uncert = 0.049; errup = 0.023; errdown = 0.024;}
    else if(ht> 200 && ht<= 600 && met> 180 && met<= 185) {uncert = 0.042; errup = 0.005; errdown = 0.005;}
    else if(ht> 600 && ht<= 800 && met> 180 && met<= 185) {uncert = 0.048; errup = 0.021; errdown = 0.023;}
    else if(ht> 800 && ht<=1000 && met> 180 && met<= 185) {uncert = 0.064; errup = 0.038; errdown = 0.048;}
    else if(ht>1000 && ht<=9999 && met> 180 && met<= 185) {uncert = 0.069; errup = 0.048; errdown = 0.055;}
    else if(ht>   0 && ht<= 200 && met> 185 && met<= 190) {uncert = 0.049; errup = 0.024; errdown = 0.026;}
    else if(ht> 200 && ht<= 600 && met> 185 && met<= 190) {uncert = 0.042; errup = 0.005; errdown = 0.005;}
    else if(ht> 600 && ht<= 800 && met> 185 && met<= 190) {uncert = 0.048; errup = 0.021; errdown = 0.023;}
    else if(ht> 800 && ht<=1000 && met> 185 && met<= 190) {uncert = 0.065; errup = 0.041; errdown = 0.049;}
    else if(ht>1000 && ht<=9999 && met> 185 && met<= 190) {uncert = 0.069; errup = 0.044; errdown = 0.055;}
    else if(ht>   0 && ht<= 200 && met> 190 && met<= 195) {uncert = 0.051; errup = 0.026; errdown = 0.028;}
    else if(ht> 200 && ht<= 600 && met> 190 && met<= 195) {uncert = 0.042; errup = 0.005; errdown = 0.005;}
    else if(ht> 600 && ht<= 800 && met> 190 && met<= 195) {uncert = 0.048; errup = 0.020; errdown = 0.023;}
    else if(ht> 800 && ht<=1000 && met> 190 && met<= 195) {uncert = 0.062; errup = 0.036; errdown = 0.045;}
    else if(ht>1000 && ht<=9999 && met> 190 && met<= 195) {uncert = 0.073; errup = 0.051; errdown = 0.059;}
    else if(ht>   0 && ht<= 200 && met> 195 && met<= 200) {uncert = 0.055; errup = 0.033; errdown = 0.036;}
    else if(ht> 200 && ht<= 600 && met> 195 && met<= 200) {uncert = 0.042; errup = 0.005; errdown = 0.005;}
    else if(ht> 600 && ht<= 800 && met> 195 && met<= 200) {uncert = 0.046; errup = 0.016; errdown = 0.020;}
    else if(ht> 800 && ht<=1000 && met> 195 && met<= 200) {uncert = 0.059; errup = 0.027; errdown = 0.041;}
    else if(ht>1000 && ht<=9999 && met> 195 && met<= 200) {uncert = 0.072; errup = 0.049; errdown = 0.059;}
    else if(ht>   0 && ht<= 200 && met> 200 && met<= 210) {uncert = 0.035; errup = 0.022; errdown = 0.025;}
    else if(ht> 200 && ht<= 600 && met> 200 && met<= 210) {uncert = 0.024; errup = 0.003; errdown = 0.003;}
    else if(ht> 600 && ht<= 800 && met> 200 && met<= 210) {uncert = 0.028; errup = 0.013; errdown = 0.015;}
    else if(ht> 800 && ht<=1000 && met> 200 && met<= 210) {uncert = 0.036; errup = 0.023; errdown = 0.027;}
    else if(ht>1000 && ht<=9999 && met> 200 && met<= 210) {uncert = 0.049; errup = 0.036; errdown = 0.042;}
    else if(ht>   0 && ht<= 200 && met> 210 && met<= 220) {uncert = 0.040; errup = 0.028; errdown = 0.032;}
    else if(ht> 200 && ht<= 600 && met> 210 && met<= 220) {uncert = 0.024; errup = 0.003; errdown = 0.003;}
    else if(ht> 600 && ht<= 800 && met> 210 && met<= 220) {uncert = 0.027; errup = 0.011; errdown = 0.013;}
    else if(ht> 800 && ht<=1000 && met> 210 && met<= 220) {uncert = 0.039; errup = 0.024; errdown = 0.031;}
    else if(ht>1000 && ht<=9999 && met> 210 && met<= 220) {uncert = 0.036; errup = 0.018; errdown = 0.027;}
    else if(ht>   0 && ht<= 200 && met> 220 && met<= 230) {uncert = 0.044; errup = 0.029; errdown = 0.037;}
    else if(ht> 200 && ht<= 600 && met> 220 && met<= 230) {uncert = 0.024; errup = 0.003; errdown = 0.003;}
    else if(ht> 600 && ht<= 800 && met> 220 && met<= 230) {uncert = 0.026; errup = 0.008; errdown = 0.011;}
    else if(ht> 800 && ht<=1000 && met> 220 && met<= 230) {uncert = 0.035; errup = 0.017; errdown = 0.025;}
    else if(ht>1000 && ht<=9999 && met> 220 && met<= 230) {uncert = 0.036; errup = 0.016; errdown = 0.027;}
    else if(ht>   0 && ht<= 200 && met> 230 && met<= 240) {uncert = 0.055; errup = 0.042; errdown = 0.053;}
    else if(ht> 200 && ht<= 600 && met> 230 && met<= 240) {uncert = 0.015; errup = 0.003; errdown = 0.003;}
    else if(ht> 600 && ht<= 800 && met> 230 && met<= 240) {uncert = 0.020; errup = 0.009; errdown = 0.013;}
    else if(ht> 800 && ht<=1000 && met> 230 && met<= 240) {uncert = 0.027; errup = 0.011; errdown = 0.022;}
    else if(ht>1000 && ht<=9999 && met> 230 && met<= 240) {uncert = 0.046; errup = 0.028; errdown = 0.043;}
    else if(ht>   0 && ht<= 200 && met> 240 && met<= 250) {uncert = 0.067; errup = 0.047; errdown = 0.065;}
    else if(ht> 200 && ht<= 600 && met> 240 && met<= 250) {uncert = 0.015; errup = 0.003; errdown = 0.003;}
    else if(ht> 600 && ht<= 800 && met> 240 && met<= 250) {uncert = 0.018; errup = 0.005; errdown = 0.010;}
    else if(ht> 800 && ht<=1000 && met> 240 && met<= 250) {uncert = 0.030; errup = 0.010; errdown = 0.026;}
    else if(ht>1000 && ht<=9999 && met> 240 && met<= 250) {uncert = 0.047; errup = 0.030; errdown = 0.044;}
    else if(ht>   0 && ht<= 200 && met> 250 && met<= 275) {uncert = 0.049; errup = 0.033; errdown = 0.047;}
    else if(ht> 200 && ht<= 600 && met> 250 && met<= 275) {uncert = 0.012; errup = 0.002; errdown = 0.002;}
    else if(ht> 600 && ht<= 800 && met> 250 && met<= 275) {uncert = 0.013; errup = 0.004; errdown = 0.006;}
    else if(ht> 800 && ht<=1000 && met> 250 && met<= 275) {uncert = 0.020; errup = 0.009; errdown = 0.016;}
    else if(ht>1000 && ht<=9999 && met> 250 && met<= 275) {uncert = 0.026; errup = 0.015; errdown = 0.023;}
    else if(ht>   0 && ht<= 200 && met> 275 && met<= 300) {uncert = 0.096; errup = 0.065; errdown = 0.096;}
    else if(ht> 200 && ht<= 600 && met> 275 && met<= 300) {uncert = 0.012; errup = 0.002; errdown = 0.003;}
    else if(ht> 600 && ht<= 800 && met> 275 && met<= 300) {uncert = 0.015; errup = 0.005; errdown = 0.008;}
    else if(ht> 800 && ht<=1000 && met> 275 && met<= 300) {uncert = 0.027; errup = 0.016; errdown = 0.024;}
    else if(ht>1000 && ht<=9999 && met> 275 && met<= 300) {uncert = 0.023; errup = 0.007; errdown = 0.020;}
    else if(ht>   0 && ht<= 200 && met> 300 && met<=9999) {uncert = 0.089; errup = 0.074; errdown = 0.089;}
    else if(ht> 200 && ht<= 600 && met> 300 && met<=9999) {uncert = 0.006; errup = 0.001; errdown = 0.002;}
    else if(ht> 600 && ht<= 800 && met> 300 && met<=9999) {uncert = 0.007; errup = 0.002; errdown = 0.003;}
    else if(ht> 800 && ht<=1000 && met> 300 && met<=9999) {uncert = 0.007; errup = 0.000; errdown = 0.003;}
    else if(ht>1000 && ht<=9999 && met> 300 && met<=9999) {uncert = 0.010; errup = 0.005; errdown = 0.008;}

    return uncert;
  });
  
  
//// Efficiency of the MET[100||110||120] triggers in all 36.2 ifb
const NamedFunc eff_higtrig("eff_higtrig", [](const Baby &b) -> NamedFunc::ScalarType{
    float errup, errdown; // Not used, but for reference
    float eff = 1., met = b.met(), ht = b.ht();
    errup=0;errdown=0;
    errup+=errdown;
    if(b.type()>0 && b.type()<1000) eff = 1;

    //// Efficiency of the MET[100||110||120] triggers in all 36.2 ifb
    //// "(trig[13]||trig[33]||trig[14]||trig[15]||trig[30]||trig[31])";
    else if(b.nvleps()==0){
      if(b.type()>=7000 && b.type()<8000) { // FAKE MET (QCD)
        if(ht>   0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.242; errup = 0.030; errdown = 0.028;}
        else if(ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff = 0.320; errup = 0.007; errdown = 0.006;}
        else if(ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.359; errup = 0.005; errdown = 0.005;}
        else if(ht> 800 && ht<=1000 && met> 150 && met<= 155) {eff = 0.379; errup = 0.002; errdown = 0.002;}
        else if(ht>1000 && ht<=9999 && met> 150 && met<= 155) {eff = 0.306; errup = 0.002; errdown = 0.002;}
        else if(ht>   0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.217; errup = 0.034; errdown = 0.031;}
        else if(ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff = 0.344; errup = 0.007; errdown = 0.007;}
        else if(ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.394; errup = 0.005; errdown = 0.005;}
        else if(ht> 800 && ht<=1000 && met> 155 && met<= 160) {eff = 0.421; errup = 0.003; errdown = 0.003;}
        else if(ht>1000 && ht<=9999 && met> 155 && met<= 160) {eff = 0.331; errup = 0.002; errdown = 0.002;}
        else if(ht>   0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.211; errup = 0.036; errdown = 0.033;}
        else if(ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff = 0.387; errup = 0.008; errdown = 0.008;}
        else if(ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.434; errup = 0.006; errdown = 0.006;}
        else if(ht> 800 && ht<=1000 && met> 160 && met<= 165) {eff = 0.464; errup = 0.003; errdown = 0.003;}
        else if(ht>1000 && ht<=9999 && met> 160 && met<= 165) {eff = 0.363; errup = 0.002; errdown = 0.002;}
        else if(ht>   0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.263; errup = 0.039; errdown = 0.035;}
        else if(ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff = 0.406; errup = 0.009; errdown = 0.009;}
        else if(ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.469; errup = 0.006; errdown = 0.006;}
        else if(ht> 800 && ht<=1000 && met> 165 && met<= 170) {eff = 0.503; errup = 0.003; errdown = 0.003;}
        else if(ht>1000 && ht<=9999 && met> 165 && met<= 170) {eff = 0.397; errup = 0.002; errdown = 0.002;}
        else if(ht>   0 && ht<= 200 && met> 170 && met<= 175) {eff = 0.328; errup = 0.045; errdown = 0.042;}
        else if(ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff = 0.434; errup = 0.010; errdown = 0.010;}
        else if(ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.507; errup = 0.007; errdown = 0.007;}
        else if(ht> 800 && ht<=1000 && met> 170 && met<= 175) {eff = 0.545; errup = 0.004; errdown = 0.004;}
        else if(ht>1000 && ht<=9999 && met> 170 && met<= 175) {eff = 0.422; errup = 0.003; errdown = 0.003;}
        else if(ht>   0 && ht<= 200 && met> 175 && met<= 180) {eff = 0.256; errup = 0.047; errdown = 0.042;}
        else if(ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff = 0.463; errup = 0.011; errdown = 0.011;}
        else if(ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.520; errup = 0.007; errdown = 0.007;}
        else if(ht> 800 && ht<=1000 && met> 175 && met<= 180) {eff = 0.582; errup = 0.004; errdown = 0.004;}
        else if(ht>1000 && ht<=9999 && met> 175 && met<= 180) {eff = 0.455; errup = 0.003; errdown = 0.003;}
        else if(ht>   0 && ht<= 200 && met> 180 && met<= 185) {eff = 0.234; errup = 0.048; errdown = 0.043;}
        else if(ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff = 0.497; errup = 0.011; errdown = 0.011;}
        else if(ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.558; errup = 0.008; errdown = 0.008;}
        else if(ht> 800 && ht<=1000 && met> 180 && met<= 185) {eff = 0.609; errup = 0.004; errdown = 0.004;}
        else if(ht>1000 && ht<=9999 && met> 180 && met<= 185) {eff = 0.478; errup = 0.003; errdown = 0.003;}
        else if(ht>   0 && ht<= 200 && met> 185 && met<= 190) {eff = 0.330; errup = 0.055; errdown = 0.051;}
        else if(ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff = 0.479; errup = 0.012; errdown = 0.012;}
        else if(ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.576; errup = 0.008; errdown = 0.008;}
        else if(ht> 800 && ht<=1000 && met> 185 && met<= 190) {eff = 0.646; errup = 0.004; errdown = 0.004;}
        else if(ht>1000 && ht<=9999 && met> 185 && met<= 190) {eff = 0.504; errup = 0.003; errdown = 0.003;}
        else if(ht>   0 && ht<= 200 && met> 190 && met<= 195) {eff = 0.406; errup = 0.053; errdown = 0.051;}
        else if(ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff = 0.537; errup = 0.013; errdown = 0.013;}
        else if(ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.592; errup = 0.009; errdown = 0.009;}
        else if(ht> 800 && ht<=1000 && met> 190 && met<= 195) {eff = 0.682; errup = 0.005; errdown = 0.005;}
        else if(ht>1000 && ht<=9999 && met> 190 && met<= 195) {eff = 0.528; errup = 0.003; errdown = 0.003;}
        else if(ht>   0 && ht<= 200 && met> 195 && met<= 200) {eff = 0.530; errup = 0.068; errdown = 0.069;}
        else if(ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff = 0.557; errup = 0.014; errdown = 0.014;}
        else if(ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.609; errup = 0.009; errdown = 0.009;}
        else if(ht> 800 && ht<=1000 && met> 195 && met<= 200) {eff = 0.705; errup = 0.005; errdown = 0.005;}
        else if(ht>1000 && ht<=9999 && met> 195 && met<= 200) {eff = 0.557; errup = 0.004; errdown = 0.004;}
        else if(ht>   0 && ht<= 200 && met> 200 && met<= 210) {eff = 0.407; errup = 0.045; errdown = 0.043;}
        else if(ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff = 0.549; errup = 0.011; errdown = 0.011;}
        else if(ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.626; errup = 0.007; errdown = 0.007;}
        else if(ht> 800 && ht<=1000 && met> 200 && met<= 210) {eff = 0.729; errup = 0.004; errdown = 0.004;}
        else if(ht>1000 && ht<=9999 && met> 200 && met<= 210) {eff = 0.584; errup = 0.003; errdown = 0.003;}
        else if(ht>   0 && ht<= 200 && met> 210 && met<= 220) {eff = 0.462; errup = 0.045; errdown = 0.044;}
        else if(ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff = 0.576; errup = 0.012; errdown = 0.012;}
        else if(ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.660; errup = 0.008; errdown = 0.008;}
        else if(ht> 800 && ht<=1000 && met> 210 && met<= 220) {eff = 0.771; errup = 0.004; errdown = 0.004;}
        else if(ht>1000 && ht<=9999 && met> 210 && met<= 220) {eff = 0.626; errup = 0.003; errdown = 0.003;}
        else if(ht>   0 && ht<= 200 && met> 220 && met<= 230) {eff = 0.505; errup = 0.052; errdown = 0.052;}
        else if(ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff = 0.575; errup = 0.013; errdown = 0.014;}
        else if(ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.667; errup = 0.009; errdown = 0.009;}
        else if(ht> 800 && ht<=1000 && met> 220 && met<= 230) {eff = 0.796; errup = 0.004; errdown = 0.004;}
        else if(ht>1000 && ht<=9999 && met> 220 && met<= 230) {eff = 0.657; errup = 0.003; errdown = 0.003;}
        else if(ht>   0 && ht<= 200 && met> 230 && met<= 240) {eff = 0.462; errup = 0.058; errdown = 0.057;}
        else if(ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff = 0.595; errup = 0.014; errdown = 0.015;}
        else if(ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.705; errup = 0.009; errdown = 0.009;}
        else if(ht> 800 && ht<=1000 && met> 230 && met<= 240) {eff = 0.823; errup = 0.004; errdown = 0.004;}
        else if(ht>1000 && ht<=9999 && met> 230 && met<= 240) {eff = 0.686; errup = 0.004; errdown = 0.004;}
        else if(ht>   0 && ht<= 200 && met> 240 && met<= 250) {eff = 0.500; errup = 0.057; errdown = 0.057;}
        else if(ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff = 0.614; errup = 0.016; errdown = 0.016;}
        else if(ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.722; errup = 0.010; errdown = 0.011;}
        else if(ht> 800 && ht<=1000 && met> 240 && met<= 250) {eff = 0.840; errup = 0.005; errdown = 0.005;}
        else if(ht>1000 && ht<=9999 && met> 240 && met<= 250) {eff = 0.711; errup = 0.004; errdown = 0.004;}
        else if(ht>   0 && ht<= 200 && met> 250 && met<= 275) {eff = 0.535; errup = 0.042; errdown = 0.043;}
        else if(ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff = 0.620; errup = 0.011; errdown = 0.011;}
        else if(ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.736; errup = 0.008; errdown = 0.008;}
        else if(ht> 800 && ht<=1000 && met> 250 && met<= 275) {eff = 0.859; errup = 0.003; errdown = 0.003;}
        else if(ht>1000 && ht<=9999 && met> 250 && met<= 275) {eff = 0.745; errup = 0.003; errdown = 0.003;}
        else if(ht>   0 && ht<= 200 && met> 275 && met<= 300) {eff = 0.569; errup = 0.049; errdown = 0.051;}
        else if(ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff = 0.617; errup = 0.014; errdown = 0.014;}
        else if(ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.737; errup = 0.010; errdown = 0.010;}
        else if(ht> 800 && ht<=1000 && met> 275 && met<= 300) {eff = 0.882; errup = 0.004; errdown = 0.004;}
        else if(ht>1000 && ht<=9999 && met> 275 && met<= 300) {eff = 0.777; errup = 0.003; errdown = 0.003;}
        else if(ht>   0 && ht<= 200 && met> 300 && met<= 350) {eff = 0.558; errup = 0.040; errdown = 0.041;}
        else if(ht> 200 && ht<= 600 && met> 300 && met<= 350) {eff = 0.638; errup = 0.013; errdown = 0.013;}
        else if(ht> 600 && ht<= 800 && met> 300 && met<= 350) {eff = 0.770; errup = 0.010; errdown = 0.010;}
        else if(ht> 800 && ht<=1000 && met> 300 && met<= 350) {eff = 0.902; errup = 0.003; errdown = 0.004;}
        else if(ht>1000 && ht<=9999 && met> 300 && met<= 350) {eff = 0.804; errup = 0.003; errdown = 0.003;}
        else if(ht>   0 && ht<= 200 && met> 350 && met<= 400) {eff = 0.548; errup = 0.056; errdown = 0.057;}
        else if(ht> 200 && ht<= 600 && met> 350 && met<= 400) {eff = 0.642; errup = 0.019; errdown = 0.019;}
        else if(ht> 600 && ht<= 800 && met> 350 && met<= 400) {eff = 0.777; errup = 0.014; errdown = 0.015;}
        else if(ht> 800 && ht<=1000 && met> 350 && met<= 400) {eff = 0.925; errup = 0.005; errdown = 0.005;}
        else if(ht>1000 && ht<=9999 && met> 350 && met<= 400) {eff = 0.838; errup = 0.004; errdown = 0.004;}
        else if(ht>   0 && ht<= 200 && met> 400 && met<= 450) {eff = 0.557; errup = 0.070; errdown = 0.072;}
        else if(ht> 200 && ht<= 600 && met> 400 && met<= 450) {eff = 0.710; errup = 0.024; errdown = 0.025;}
        else if(ht> 600 && ht<= 800 && met> 400 && met<= 450) {eff = 0.864; errup = 0.018; errdown = 0.020;}
        else if(ht> 800 && ht<=1000 && met> 400 && met<= 450) {eff = 0.945; errup = 0.006; errdown = 0.006;}
        else if(ht>1000 && ht<=9999 && met> 400 && met<= 450) {eff = 0.862; errup = 0.005; errdown = 0.005;}
        else if(ht>   0 && ht<= 200 && met> 450 && met<= 500) {eff = 0.735; errup = 0.081; errdown = 0.097;}
        else if(ht> 200 && ht<= 600 && met> 450 && met<= 500) {eff = 0.695; errup = 0.036; errdown = 0.039;}
        else if(ht> 600 && ht<= 800 && met> 450 && met<= 500) {eff = 0.913; errup = 0.022; errdown = 0.028;}
        else if(ht> 800 && ht<=1000 && met> 450 && met<= 500) {eff = 0.970; errup = 0.006; errdown = 0.007;}
        else if(ht>1000 && ht<=9999 && met> 450 && met<= 500) {eff = 0.886; errup = 0.006; errdown = 0.006;}
        else if(ht>   0 && ht<= 200 && met> 500 && met<=9999) {eff = 0.550; errup = 0.088; errdown = 0.091;}
        else if(ht> 200 && ht<= 600 && met> 500 && met<=9999) {eff = 0.731; errup = 0.029; errdown = 0.031;}
        else if(ht> 600 && ht<= 800 && met> 500 && met<=9999) {eff = 0.956; errup = 0.014; errdown = 0.019;}
        else if(ht> 800 && ht<=1000 && met> 500 && met<=9999) {eff = 0.981; errup = 0.004; errdown = 0.005;}
        else if(ht>1000 && ht<=9999 && met> 500 && met<=9999) {eff = 0.907; errup = 0.004; errdown = 0.004;}

        // //////// First 4.3 ifb (FULL STATUS)
        // if(b.met()<=100) {eff = 0.;}
        // if(b.met()> 100 && b.met()<= 125) {eff = 0.107;}
        // if(b.met()> 125 && b.met()<= 150) {eff = 0.245;}
        // if(b.met()> 150 && b.met()<= 175) {eff = 0.416;}
        // if(b.met()> 175 && b.met()<= 200) {eff = 0.556;}
        // if(b.met()> 200 && b.met()<= 225) {eff = 0.659;}
        // if(b.met()> 225 && b.met()<= 250) {eff = 0.727;}
        // if(b.met()> 250 && b.met()<= 275) {eff = 0.773;}
        // if(b.met()> 275 && b.met()<= 300) {eff = 0.791;}
        // if(b.met()> 300 && b.met()<=9999) {eff = 0.819;}
        // //////// First 4.3 ifb (FULL STATUS)

      } else { // TRUE MET
	if(ht>   0 && ht<= 200 && met> 150 && met<= 155) {eff = 0.532; errup = 0.013; errdown = 0.013;}
	else if(ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff = 0.612; errup = 0.005; errdown = 0.005;}
	else if(ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff = 0.589; errup = 0.023; errdown = 0.024;}
	else if(ht> 800 && ht<=1000 && met> 150 && met<= 155) {eff = 0.515; errup = 0.042; errdown = 0.042;}
	else if(ht>1000 && ht<=9999 && met> 150 && met<= 155) {eff = 0.588; errup = 0.052; errdown = 0.054;}
	else if(ht>   0 && ht<= 200 && met> 155 && met<= 160) {eff = 0.591; errup = 0.014; errdown = 0.014;}
	else if(ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff = 0.684; errup = 0.005; errdown = 0.005;}
	else if(ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff = 0.678; errup = 0.022; errdown = 0.022;}
	else if(ht> 800 && ht<=1000 && met> 155 && met<= 160) {eff = 0.537; errup = 0.042; errdown = 0.042;}
	else if(ht>1000 && ht<=9999 && met> 155 && met<= 160) {eff = 0.511; errup = 0.057; errdown = 0.057;}
	else if(ht>   0 && ht<= 200 && met> 160 && met<= 165) {eff = 0.619; errup = 0.016; errdown = 0.016;}
	else if(ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff = 0.727; errup = 0.006; errdown = 0.006;}
	else if(ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff = 0.699; errup = 0.022; errdown = 0.023;}
	else if(ht> 800 && ht<=1000 && met> 160 && met<= 165) {eff = 0.690; errup = 0.039; errdown = 0.041;}
	else if(ht>1000 && ht<=9999 && met> 160 && met<= 165) {eff = 0.568; errup = 0.048; errdown = 0.049;}
	else if(ht>   0 && ht<= 200 && met> 165 && met<= 170) {eff = 0.678; errup = 0.017; errdown = 0.018;}
	else if(ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff = 0.769; errup = 0.006; errdown = 0.006;}
	else if(ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff = 0.732; errup = 0.022; errdown = 0.023;}
	else if(ht> 800 && ht<=1000 && met> 165 && met<= 170) {eff = 0.609; errup = 0.042; errdown = 0.044;}
	else if(ht>1000 && ht<=9999 && met> 165 && met<= 170) {eff = 0.685; errup = 0.058; errdown = 0.064;}
	else if(ht>   0 && ht<= 200 && met> 170 && met<= 175) {eff = 0.670; errup = 0.019; errdown = 0.020;}
	else if(ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff = 0.811; errup = 0.005; errdown = 0.006;}
	else if(ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff = 0.779; errup = 0.022; errdown = 0.024;}
	else if(ht> 800 && ht<=1000 && met> 170 && met<= 175) {eff = 0.736; errup = 0.041; errdown = 0.045;}
	else if(ht>1000 && ht<=9999 && met> 170 && met<= 175) {eff = 0.663; errup = 0.056; errdown = 0.061;}
	else if(ht>   0 && ht<= 200 && met> 175 && met<= 180) {eff = 0.730; errup = 0.020; errdown = 0.021;}
	else if(ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff = 0.838; errup = 0.005; errdown = 0.006;}
	else if(ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff = 0.820; errup = 0.022; errdown = 0.024;}
	else if(ht> 800 && ht<=1000 && met> 175 && met<= 180) {eff = 0.819; errup = 0.037; errdown = 0.043;}
	else if(ht>1000 && ht<=9999 && met> 175 && met<= 180) {eff = 0.736; errup = 0.055; errdown = 0.062;}
	else if(ht>   0 && ht<= 200 && met> 180 && met<= 185) {eff = 0.745; errup = 0.023; errdown = 0.024;}
	else if(ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff = 0.874; errup = 0.005; errdown = 0.005;}
	else if(ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff = 0.848; errup = 0.021; errdown = 0.023;}
	else if(ht> 800 && ht<=1000 && met> 180 && met<= 185) {eff = 0.869; errup = 0.038; errdown = 0.048;}
	else if(ht>1000 && ht<=9999 && met> 180 && met<= 185) {eff = 0.759; errup = 0.048; errdown = 0.055;}
	else if(ht>   0 && ht<= 200 && met> 185 && met<= 190) {eff = 0.777; errup = 0.024; errdown = 0.026;}
	else if(ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff = 0.903; errup = 0.005; errdown = 0.005;}
	else if(ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff = 0.850; errup = 0.021; errdown = 0.023;}
	else if(ht> 800 && ht<=1000 && met> 185 && met<= 190) {eff = 0.839; errup = 0.041; errdown = 0.049;}
	else if(ht>1000 && ht<=9999 && met> 185 && met<= 190) {eff = 0.847; errup = 0.044; errdown = 0.055;}
	else if(ht>   0 && ht<= 200 && met> 190 && met<= 195) {eff = 0.792; errup = 0.026; errdown = 0.028;}
	else if(ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff = 0.907; errup = 0.005; errdown = 0.005;}
	else if(ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff = 0.884; errup = 0.020; errdown = 0.023;}
	else if(ht> 800 && ht<=1000 && met> 190 && met<= 195) {eff = 0.870; errup = 0.036; errdown = 0.045;}
	else if(ht>1000 && ht<=9999 && met> 190 && met<= 195) {eff = 0.781; errup = 0.051; errdown = 0.059;}
	else if(ht>   0 && ht<= 200 && met> 195 && met<= 200) {eff = 0.757; errup = 0.033; errdown = 0.036;}
	else if(ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff = 0.924; errup = 0.005; errdown = 0.005;}
	else if(ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff = 0.921; errup = 0.016; errdown = 0.020;}
	else if(ht> 800 && ht<=1000 && met> 195 && met<= 200) {eff = 0.936; errup = 0.027; errdown = 0.041;}
	else if(ht>1000 && ht<=9999 && met> 195 && met<= 200) {eff = 0.803; errup = 0.049; errdown = 0.059;}
	else if(ht>   0 && ht<= 200 && met> 200 && met<= 210) {eff = 0.841; errup = 0.022; errdown = 0.025;}
	else if(ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff = 0.949; errup = 0.003; errdown = 0.003;}
	else if(ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff = 0.927; errup = 0.013; errdown = 0.015;}
	else if(ht> 800 && ht<=1000 && met> 200 && met<= 210) {eff = 0.894; errup = 0.023; errdown = 0.027;}
	else if(ht>1000 && ht<=9999 && met> 200 && met<= 210) {eff = 0.839; errup = 0.036; errdown = 0.042;}
	else if(ht>   0 && ht<= 200 && met> 210 && met<= 220) {eff = 0.850; errup = 0.028; errdown = 0.032;}
	else if(ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff = 0.966; errup = 0.003; errdown = 0.003;}
	else if(ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff = 0.952; errup = 0.011; errdown = 0.013;}
	else if(ht> 800 && ht<=1000 && met> 210 && met<= 220) {eff = 0.919; errup = 0.024; errdown = 0.031;}
	else if(ht>1000 && ht<=9999 && met> 210 && met<= 220) {eff = 0.959; errup = 0.018; errdown = 0.027;}
	else if(ht>   0 && ht<= 200 && met> 220 && met<= 230) {eff = 0.896; errup = 0.029; errdown = 0.037;}
	else if(ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff = 0.973; errup = 0.003; errdown = 0.003;}
	else if(ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff = 0.979; errup = 0.008; errdown = 0.011;}
	else if(ht> 800 && ht<=1000 && met> 220 && met<= 230) {eff = 0.956; errup = 0.017; errdown = 0.025;}
	else if(ht>1000 && ht<=9999 && met> 220 && met<= 230) {eff = 0.971; errup = 0.016; errdown = 0.027;}
	else if(ht>   0 && ht<= 200 && met> 230 && met<= 240) {eff = 0.844; errup = 0.042; errdown = 0.053;}
	else if(ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff = 0.983; errup = 0.003; errdown = 0.003;}
	else if(ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff = 0.976; errup = 0.009; errdown = 0.013;}
	else if(ht> 800 && ht<=1000 && met> 230 && met<= 240) {eff = 0.983; errup = 0.011; errdown = 0.022;}
	else if(ht>1000 && ht<=9999 && met> 230 && met<= 240) {eff = 0.942; errup = 0.028; errdown = 0.043;}
	else if(ht>   0 && ht<= 200 && met> 240 && met<= 250) {eff = 0.880; errup = 0.047; errdown = 0.065;}
	else if(ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff = 0.985; errup = 0.003; errdown = 0.003;}
	else if(ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff = 0.992; errup = 0.005; errdown = 0.010;}
	else if(ht> 800 && ht<=1000 && met> 240 && met<= 250) {eff = 0.989; errup = 0.010; errdown = 0.026;}
	else if(ht>1000 && ht<=9999 && met> 240 && met<= 250) {eff = 0.931; errup = 0.030; errdown = 0.044;}
	else if(ht>   0 && ht<= 200 && met> 250 && met<= 275) {eff = 0.915; errup = 0.033; errdown = 0.047;}
	else if(ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff = 0.989; errup = 0.002; errdown = 0.002;}
	else if(ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff = 0.992; errup = 0.004; errdown = 0.006;}
	else if(ht> 800 && ht<=1000 && met> 250 && met<= 275) {eff = 0.984; errup = 0.009; errdown = 0.016;}
	else if(ht>1000 && ht<=9999 && met> 250 && met<= 275) {eff = 0.965; errup = 0.015; errdown = 0.023;}
	else if(ht>   0 && ht<= 200 && met> 275 && met<= 300) {eff = 0.862; errup = 0.065; errdown = 0.096;}
	else if(ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff = 0.992; errup = 0.002; errdown = 0.003;}
	else if(ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff = 0.989; errup = 0.005; errdown = 0.008;}
	else if(ht> 800 && ht<=1000 && met> 275 && met<= 300) {eff = 0.963; errup = 0.016; errdown = 0.024;}
	else if(ht>1000 && ht<=9999 && met> 275 && met<= 300) {eff = 0.991; errup = 0.007; errdown = 0.020;}
	else if(ht>   0 && ht<= 200 && met> 300 && met<=9999) {eff = 0.744; errup = 0.074; errdown = 0.089;}
	else if(ht> 200 && ht<= 600 && met> 300 && met<=9999) {eff = 0.994; errup = 0.001; errdown = 0.002;}
	else if(ht> 600 && ht<= 800 && met> 300 && met<=9999) {eff = 0.996; errup = 0.002; errdown = 0.003;}
	else if(ht> 800 && ht<=1000 && met> 300 && met<=9999) {eff = 1.000; errup = 0.000; errdown = 0.003;}
	else if(ht>1000 && ht<=9999 && met> 300 && met<=9999) {eff = 0.987; errup = 0.005; errdown = 0.008;}
      }

      //// MET || Ele27 || Ele105 || Ele115
      //// "(trig[13]||trig[33]||trig[14]||trig[15]||trig[30]||trig[31]||trig[22]||trig[40]||trig[24]||trig[41])"
    } else if(b.nels()==1 && b.nmus()==0){
      vector<float> leps_pt; 
      if (b.leps_pt()->size()>0) leps_pt.push_back(b.leps_pt()->at(0));
      else leps_pt.push_back(0);
      if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met>=0 && met<= 110) {eff=0.160; errup=0.019; errdown = 0.017;}
      else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met>=0 && met<= 110) {eff=0.400; errup=0.024; errdown=0.024;}
      else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met>=0 && met<= 110) {eff=0.728; errup=0.006; errdown=0.006;}
      else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met>=0 && met<= 110) {eff=0.880; errup=0.017; errdown=0.019;}
      else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met>=0 && met<= 110) {eff=0.950; errup=0.003; errdown=0.003;}
      else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met>110 && met<= 120) {eff=0.244; errup=0.024; errdown=0.023;}
      else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met>110 && met<= 120) {eff=0.420; errup=0.027; errdown=0.027;}
      else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met>110 && met<= 120) {eff=0.761; errup=0.007; errdown=0.007;}
      else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met>110 && met<= 120) {eff=0.918; errup=0.015; errdown=0.017;}
      else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met>110 && met<= 120) {eff=0.958; errup=0.003; errdown=0.003;}
      else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met>120 && met<= 130) {eff=0.331; errup=0.030; errdown=0.029;}
      else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met>120 && met<= 130) {eff=0.500; errup=0.031; errdown=0.031;}
      else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met>120 && met<= 130) {eff=0.800; errup=0.007; errdown=0.007;}
      else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met>120 && met<= 130) {eff=0.928; errup=0.015; errdown=0.018;}
      else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met>120 && met<= 130) {eff=0.960; errup=0.003; errdown=0.003;}
      else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met>130 && met<= 140) {eff=0.491; errup=0.031; errdown=0.031;}
      else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met>130 && met<= 140) {eff=0.608; errup=0.033; errdown=0.034;}
      else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met>130 && met<= 140) {eff=0.831; errup=0.007; errdown=0.007;}
      else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met>130 && met<= 140) {eff=0.931; errup=0.016; errdown=0.020;}
      else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met>130 && met<= 140) {eff=0.967; errup=0.003; errdown=0.003;}
      else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met>140 && met<= 150) {eff=0.573; errup=0.033; errdown=0.033;}
      else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met>140 && met<= 150) {eff=0.677; errup=0.035; errdown=0.037;}
      else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met>140 && met<= 150) {eff=0.856; errup=0.007; errdown=0.007;}
      else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met>140 && met<= 150) {eff=0.923; errup=0.018; errdown=0.022;}
      else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met>140 && met<= 150) {eff=0.971; errup=0.003; errdown=0.004;}
      else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met>150 && met<= 160) {eff=0.643; errup=0.037; errdown=0.039;}
      else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met> 150 && met<= 160) {eff=0.738; errup=0.033; errdown=0.036;}
      else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met> 150 && met<= 160) {eff=0.871; errup=0.007; errdown=0.007;}
      else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met> 150 && met<= 160) {eff=0.935; errup=0.016; errdown=0.021;}
      else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met> 150 && met<= 160) {eff=0.982; errup=0.003; errdown=0.003;}
      else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met> 160 && met<= 170) {eff=0.760; errup=0.034; errdown=0.038;}
      else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met> 160 && met<= 170) {eff=0.792; errup=0.034; errdown=0.038;}
      else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met> 160 && met<= 170) {eff=0.910; errup=0.007; errdown=0.007;}
      else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met> 160 && met<= 170) {eff=0.970; errup=0.013; errdown=0.020;}
      else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met> 160 && met<= 170) {eff=0.981; errup=0.003; errdown=0.004;}
      else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met> 170 && met<= 180) {eff=0.829; errup=0.033; errdown=0.038;}
      else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met> 170 && met<= 180) {eff=0.863; errup=0.031; errdown=0.037;}
      else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met> 170 && met<= 180) {eff=0.937; errup=0.006; errdown=0.006;}
      else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met> 170 && met<= 180) {eff=0.988; errup=0.008; errdown=0.016;}
      else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met> 170 && met<= 180) {eff=0.982; errup=0.003; errdown=0.004;}
      else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met> 180 && met<= 190) {eff=0.761; errup=0.041; errdown=0.046;}
      else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met> 180 && met<= 190) {eff=0.863; errup=0.032; errdown=0.038;}
      else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met> 180 && met<= 190) {eff=0.939; errup=0.006; errdown=0.007;}
      else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met> 180 && met<= 190) {eff=0.969; errup=0.015; errdown=0.024;}
      else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met> 180 && met<= 190) {eff=0.984; errup=0.003; errdown=0.004;}
      else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met> 190 && met<= 200) {eff=0.892; errup=0.033; errdown=0.042;}
      else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met> 190 && met<= 200) {eff=0.902; errup=0.030; errdown=0.039;}
      else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met> 190 && met<= 200) {eff=0.956; errup=0.006; errdown=0.006;}
      else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met> 190 && met<= 200) {eff=0.983; errup=0.011; errdown=0.022;}
      else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met> 190 && met<= 200) {eff=0.993; errup=0.002; errdown=0.003;}
      else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met> 200 && met<= 210) {eff=0.950; errup=0.021; errdown=0.032;}
      else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met> 200 && met<= 210) {eff=0.951; errup=0.023; errdown=0.037;}
      else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met> 200 && met<= 210) {eff=0.973; errup=0.005; errdown=0.005;}
      else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met> 200 && met<= 210) {eff=1.000; errup=0.000; errdown=0.018;}
      else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met> 200 && met<= 210) {eff=0.985; errup=0.003; errdown=0.004;}
      else if(leps_pt[0]> 20 && leps_pt[0]<=  25 && met> 210 && met<=9999) {eff=0.974; errup=0.005; errdown=0.006;}
      else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met> 210 && met<=9999) {eff=0.981; errup=0.004; errdown=0.005;}
      else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met> 210 && met<=9999) {eff=0.992; errup=0.001; errdown=0.001;}
      else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met> 210 && met<=9999) {eff=0.997; errup=0.002; errdown=0.003;}
      else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met> 210 && met<=9999) {eff=0.996; errup=0.001; errdown=0.001;}

      //// MET || Mu24 || Mu50
      //// "(trig[13]||trig[33]||trig[14]||trig[15]||trig[30]||trig[31]||trig[19]||trig[55]||trig[21])"
    } else if(b.nels()==0 && b.nmus()==1){
      vector<float> leps_pt; 
      if (b.leps_pt()->size()>0) leps_pt.push_back(b.leps_pt()->at(0));
      else leps_pt.push_back(0);
      if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met>= 0 && met<= 110) {eff=0.271; errup=0.017; errdown = 0.016;}
      else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met>= 0 && met<= 110) {eff=0.725; errup=0.017; errdown=0.018;}
      else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met>= 0 && met<= 110) {eff=0.814; errup=0.008; errdown=0.009;}
      else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met>= 0 && met<= 110) {eff=0.964; errup=0.002; errdown=0.002;}
      else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 110 && met<= 120) {eff=0.363; errup=0.020; errdown=0.020;}
      else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 110 && met<= 120) {eff=0.755; errup=0.018; errdown=0.019;}
      else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 110 && met<= 120) {eff=0.842; errup=0.009; errdown=0.009;}
      else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 110 && met<= 120) {eff=0.969; errup=0.002; errdown=0.002;}
      else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 120 && met<= 130) {eff=0.452; errup=0.022; errdown=0.022;}
      else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 120 && met<= 130) {eff=0.824; errup=0.018; errdown=0.019;}
      else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 120 && met<= 130) {eff=0.869; errup=0.009; errdown=0.009;}
      else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 120 && met<= 130) {eff=0.971; errup=0.002; errdown=0.002;}
      else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 130 && met<= 140) {eff=0.590; errup=0.025; errdown=0.025;}
      else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 130 && met<= 140) {eff=0.875; errup=0.017; errdown=0.019;}
      else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 130 && met<= 140) {eff=0.904; errup=0.008; errdown=0.009;}
      else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 130 && met<= 140) {eff=0.972; errup=0.002; errdown=0.002;}
      else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 140 && met<= 150) {eff=0.660; errup=0.026; errdown=0.027;}
      else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 140 && met<= 150) {eff=0.891; errup=0.017; errdown=0.019;}
      else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 140 && met<= 150) {eff=0.938; errup=0.007; errdown=0.008;}
      else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 140 && met<= 150) {eff=0.980; errup=0.002; errdown=0.002;}
      else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 150 && met<= 160) {eff=0.778; errup=0.024; errdown=0.026;}
      else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 150 && met<= 160) {eff=0.915; errup=0.016; errdown=0.019;}
      else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 150 && met<= 160) {eff=0.940; errup=0.008; errdown=0.009;}
      else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 150 && met<= 160) {eff=0.984; errup=0.002; errdown=0.002;}
      else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 160 && met<= 170) {eff=0.798; errup=0.026; errdown=0.029;}
      else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 160 && met<= 170) {eff=0.946; errup=0.015; errdown=0.020;}
      else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 160 && met<= 170) {eff=0.967; errup=0.006; errdown=0.008;}
      else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 160 && met<= 170) {eff=0.991; errup=0.002; errdown=0.002;}
      else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 170 && met<= 180) {eff=0.885; errup=0.022; errdown=0.025;}
      else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 170 && met<= 180) {eff=0.937; errup=0.016; errdown=0.021;}
      else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 170 && met<= 180) {eff=0.977; errup=0.006; errdown=0.007;}
      else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 170 && met<= 180) {eff=0.987; errup=0.002; errdown=0.002;}
      else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 180 && met<= 190) {eff=0.927; errup=0.019; errdown=0.024;}
      else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 180 && met<= 190) {eff=0.958; errup=0.014; errdown=0.019;}
      else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 180 && met<= 190) {eff=0.974; errup=0.006; errdown=0.008;}
      else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 180 && met<= 190) {eff=0.992; errup=0.002; errdown=0.002;}
      else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 190 && met<= 200) {eff=0.921; errup=0.019; errdown=0.024;}
      else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 190 && met<= 200) {eff=0.965; errup=0.014; errdown=0.020;}
      else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 190 && met<= 200) {eff=0.991; errup=0.004; errdown=0.006;}
      else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 190 && met<= 200) {eff=0.991; errup=0.002; errdown=0.002;}
      else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 200 && met<= 210) {eff=0.926; errup=0.022; errdown=0.028;}
      else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 200 && met<= 210) {eff=0.994; errup=0.005; errdown=0.015;}
      else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 200 && met<= 210) {eff=0.994; errup=0.003; errdown=0.006;}
      else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 200 && met<= 210) {eff=0.994; errup=0.002; errdown=0.002;}
      else if(leps_pt[0]>20 && leps_pt[0]<=  25 && met> 210 && met<=9999) {eff=0.981; errup=0.004; errdown=0.004;}
      else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 210 && met<=9999) {eff=0.994; errup=0.002; errdown=0.003;}
      else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 210 && met<=9999) {eff=0.996; errup=0.001; errdown=0.001;}
      else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 210 && met<=9999) {eff=0.997; errup=0.000; errdown=0.000;}

      //// Ele27 || Ele105 || Ele115
      //// "(trig[22]||trig[40]||trig[24]||trig[41])"
    } else if(b.nels()==2 && b.nmus()==0){
      vector<float> leps_pt; 
      if (b.leps_pt()->size()>0) leps_pt.push_back(b.leps_pt()->at(0));
      else leps_pt.push_back(0);
      if(leps_pt[0]>  40 && leps_pt[0]<=  45) {eff = 0.944; errup = 0.015; errdown = 0.019;}
      else if(leps_pt[0]>  45 && leps_pt[0]<=  50) {eff = 0.910; errup = 0.015; errdown = 0.017;}
      else if(leps_pt[0]>  50 && leps_pt[0]<=  55) {eff = 0.927; errup = 0.013; errdown = 0.015;}
      else if(leps_pt[0]>  55 && leps_pt[0]<=  60) {eff = 0.912; errup = 0.013; errdown = 0.015;}
      else if(leps_pt[0]>  60 && leps_pt[0]<=  65) {eff = 0.941; errup = 0.011; errdown = 0.013;}
      else if(leps_pt[0]>  65 && leps_pt[0]<=  70) {eff = 0.901; errup = 0.014; errdown = 0.016;}
      else if(leps_pt[0]>  70 && leps_pt[0]<=  75) {eff = 0.921; errup = 0.013; errdown = 0.016;}
      else if(leps_pt[0]>  75 && leps_pt[0]<=  80) {eff = 0.947; errup = 0.011; errdown = 0.014;}
      else if(leps_pt[0]>  80 && leps_pt[0]<=  85) {eff = 0.954; errup = 0.011; errdown = 0.013;}
      else if(leps_pt[0]>  85 && leps_pt[0]<=  90) {eff = 0.939; errup = 0.012; errdown = 0.014;}
      else if(leps_pt[0]>  90 && leps_pt[0]<=  95) {eff = 0.940; errup = 0.012; errdown = 0.015;}
      else if(leps_pt[0]>  95 && leps_pt[0]<= 100) {eff = 0.932; errup = 0.014; errdown = 0.017;}
      else if(leps_pt[0]> 100 && leps_pt[0]<= 105) {eff = 0.934; errup = 0.014; errdown = 0.017;}
      else if(leps_pt[0]> 105 && leps_pt[0]<= 110) {eff = 0.965; errup = 0.010; errdown = 0.014;}
      else if(leps_pt[0]> 110 && leps_pt[0]<=9999) {eff = 0.994; errup = 0.001; errdown = 0.001;}

      //// Mu24 || Mu50
      //// "(trig[19]||trig[55]||trig[21])"
    } else if(b.nels()==0 && b.nmus()==2){
      vector<float> leps_pt; 
      if (b.leps_pt()->size()>0) leps_pt.push_back(b.leps_pt()->at(0));
      else leps_pt.push_back(0);
      if(leps_pt[0]>  40 && leps_pt[0]<=  45) {eff = 0.959; errup = 0.010; errdown = 0.012;}
      if(leps_pt[0]>  45 && leps_pt[0]<=  50) {eff = 0.970; errup = 0.006; errdown = 0.007;}
      if(leps_pt[0]>  50 && leps_pt[0]<=9999) {eff = 0.982; errup = 0.000; errdown = 0.000;}
    }
    return eff;
  });

const NamedFunc weight_hig("weight_hig",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.type()>0 && b.type()<1000) return 1;
  else return b.weight()/b.w_btag()*b.w_bhig();
});

const NamedFunc weight_higd("weight_hig_deep",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.type()>0 && b.type()<1000) return 1;
  else return b.weight()/b.w_btag()*b.w_bhig_deep();
});

const NamedFunc mhig("mhig",[](const Baby &b) -> NamedFunc::ScalarType{
  float mass = -999;
  for (unsigned i(0); i<b.mc_mass()->size(); i++){
    if (b.mc_id()->at(i)==1000023) {
      mass = b.mc_mass()->at(i);
      break;
    }
  }
  return mass;
});

const NamedFunc nb_exci("nb_exci",[](const Baby &b) -> NamedFunc::ScalarType{
  int nb=0;
  for (unsigned i(0); i<b.mc_id()->size(); i++){
    if (abs(b.mc_id()->at(i))==5 && b.mc_status()->at(i)==23 
	&& abs(b.mc_mom()->at(i))!=6 && abs(b.mc_mom()->at(i))!=24 && abs(b.mc_mom()->at(i))!=23) 
      nb++;
  }
  return nb;
});

const NamedFunc nb_gs("nb_gs",[](const Baby &b) -> NamedFunc::ScalarType{
  int nb=0;
  for (unsigned i(0); i<b.mc_id()->size(); i++){
    if (abs(b.mc_id()->at(i))==5 && b.mc_status()->at(i)!=23 ) 
      nb++;
  }
  return nb;
});

}
