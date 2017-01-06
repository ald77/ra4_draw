#include "hig/hig_functions.hpp"

#include "TVector2.h"

#include "core/utilities.hpp"
#include "core/config_parser.hpp"

using namespace std;

namespace Higfuncs{

const NamedFunc ntrub("ntrub",[](const Baby &b) -> NamedFunc::ScalarType{
  int tmp_ntrub(0);
  for (unsigned i(0); i<b.jets_pt()->size(); i++){
    if (!b.jets_h1()->at(i) && !b.jets_h2()->at(i)) continue;
    if (b.jets_hflavor()->at(i)==5) tmp_ntrub++;
  }
  return tmp_ntrub;
});

const NamedFunc hig_nb("hig_nb",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.nbt()==2 && b.nbm()==2) return 2;
  else if (b.nbt()>=2 && b.nbm()==3 && b.nbl()==3) return 3;
  else if (b.nbt()>=2 && b.nbm()>=3 && b.nbl()>=4) return 4;
  else return 0;
});

const NamedFunc hig_nb_extended("hig_nb_extended",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.nbm()==0) return 0;
  else if (b.nbm()==1) return 1;
  else if (b.nbt()==2 && b.nbm()==2) return 2;
  else if (b.nbt()>=2 && b.nbm()==3 && b.nbl()==3) return 3;
  else if (b.nbt()>=2 && b.nbm()>=3 && b.nbl()>=4) return 4;
  else return 99;
});

const NamedFunc hig_nb_mmmm("hig_nb_mmmm",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.nbm()>=2) return min(4,b.nbm());
  else return 0;
});
const NamedFunc hig_nb_ttll("hig_nb_ttll",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.nbt()==2 && b.nbl()==2) return 2;
  else if (b.nbt()>=2 && b.nbl()==3) return 3;
  else if (b.nbt()>=2 && b.nbl()>=4) return 4;
  else return 0;
});
const NamedFunc hig_nb_tmml("hig_nb_tmml",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.nbt()>=1 && b.nbm()==2) return 2;
  else if (b.nbt()>=1 && b.nbm()==3 && b.nbl()==3) return 3;
  else if (b.nbt()>=1 && b.nbm()>=3 && b.nbl()>=4) return 4;
  else return 0;
});

// weight used to subtract ttbar based on MC prediction reweighted to data in 1l CR
// since ttbar has to be combined in the same process def with data, also apply stitch here
NamedFunc::ScalarType wgt_subtr_ttx(const Baby &b, string json){
  if ( (b.type()>=1000 && b.type()<2000) ||  // ttbar
    // (b.type()>=3000 && b.type()<4000) ||     // single top
    (b.type()>=4000 && b.type()<6000) ||     // ttw and ttz
    b.type()==9000  ||                       // ttH
    b.type()==10000  ||                      // ttgamma
    b.type()==11000) {
    float wgt = -1; 
    // apply lumi 
    if (json=="json4p0") wgt*= 4.3;
    else wgt*= 36.2;
    // apply weights derived from 1l CR: normalization*ratio(data/mc)
    if (b.nbt()==2 && b.nbm()==2) return wgt*=1.13;
    else if (b.nbt()>=2 && b.nbm()==3 && b.nbl()==3) return wgt*=1.29;
    else if (b.nbt()>=2 && b.nbm()>=3 && b.nbl()>=4) return wgt*=1.41;
    // apply stitch
    if (b.stitch()) return wgt;
    else return 0;
  } else if (b.type()>0 && b.type()<1000){ // apply trigger and json for data
    if (b.trig()->at(13)||b.trig()->at(33)||b.trig()->at(14)||b.trig()->at(15)||b.trig()->at(30)||b.trig()->at(31)
      ||b.trig()->at(22)||b.trig()->at(40)||b.trig()->at(24)||b.trig()->at(41)
      ||b.trig()->at(19)||b.trig()->at(55)||b.trig()->at(21)) {
      if (json=="json4p0") return b.json4p0() ? 1:0;
      else return 1;
    } else {
      return 0;
    }
  }
  // for all other backgrounds, chill
  return 1;
};

// calculate effect of systematics calculated for each background 
// in the data control regions on the total bkg. kappa
const NamedFunc wgt_syst_ttx("wgt_syst_ttx",[](const Baby &b) -> NamedFunc::ScalarType{
  if ( (b.type()>=1000 && b.type()<2000) ||  // ttbar
    // (b.type()>=3000 && b.type()<4000) ||     // single top
    (b.type()>=4000 && b.type()<6000) ||     // ttw and ttz
    b.type()==9000  ||                       // ttH
    b.type()==10000  ||                      // ttgamma
    b.type()==11000) {                       // tttt
    if (b.hig_am()<=100 || (b.hig_am()>140 && b.hig_am()<=200)) {
      if (b.nbt()>=1 && b.nbm()==3 && b.nbl()==3) return 0.19;
      else if (b.nbt()>=1 && b.nbm()>=3 && b.nbl()>=4) return 0.30;
    }
  }
  return 0;
});

const NamedFunc wgt_syst_vjets("wgt_syst_vjets",[](const Baby &b) -> NamedFunc::ScalarType{
  if ( (b.type()>=8000 && b.type()<9000) || // zjets
    (b.type()>=2000 && b.type()<3000) ||    // wjets
    (b.type()>=6000 && b.type()<7000)) {    // dyjets  
    if (b.hig_am()<=100 || (b.hig_am()>140 && b.hig_am()<=200))
      if (b.nbt()>=2 && b.nbm()>=3) return 0.34;
  }
  return 0;
});

const NamedFunc wgt_syst_qcd("wgt_syst_qcd",[](const Baby &b) -> NamedFunc::ScalarType{
  if ( (b.type()>=7000 && b.type()<8000)) { // qcd
    if (b.hig_am()<=100 || (b.hig_am()>140 && b.hig_am()<=200))
      if (b.nbt()>=2 && b.nbm()>=3) return 0.12;
  }
  return 0;
});

// estimate the systematic due to limited knowledge on composition
const NamedFunc wgt_2xhighnb_zjets("wgt_2xhighnb_zjets",[](const Baby &b) -> NamedFunc::ScalarType{
  if (b.type()>=8000 && b.type()<9000)
    if (b.nbt()>=2 && b.nbm()>=3) return 2.;
  return 1;
});

// defintion of analysis trigger
const NamedFunc trig_hig("trig_hig", [](const Baby &b) -> NamedFunc::ScalarType{
	if ( b.trig()->at(13)||b.trig()->at(33)||b.trig()->at(14)||b.trig()->at(15)||b.trig()->at(30)||b.trig()->at(31)
		||b.trig()->at(22)||b.trig()->at(40)||b.trig()->at(24)||b.trig()->at(41)
		||b.trig()->at(19)||b.trig()->at(55)||b.trig()->at(21)) return 1;
	return -1;
});

  //// Efficiency of the MET[100||110||120] triggers in all 36.2 ifb
const NamedFunc eff_higtrig("eff_higtrig", [](const Baby &b) -> NamedFunc::ScalarType{
    float errup, errdown; // Not used, but for reference
    float eff = 1., met = b.met(), ht = b.ht();
    errup=0;errdown=0;
    errup+=errdown;
    if(b.type() < 1000) eff = 1;

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
      	if(ht>   0 && ht<= 200 && met> 150 && met<= 155) {eff=0.456; errup=0.012; errdown=0.012;}
      	else if(ht> 200 && ht<= 600 && met> 150 && met<= 155) {eff=0.525; errup=0.004; errdown=0.004;}
      	else if(ht> 600 && ht<= 800 && met> 150 && met<= 155) {eff=0.459; errup=0.015; errdown=0.015;}
      	else if(ht> 800 && ht<=1000 && met> 150 && met<= 155) {eff=0.455; errup=0.025; errdown=0.024;}
      	else if(ht>1000 && ht<=9999 && met> 150 && met<= 155) {eff=0.396; errup=0.029; errdown=0.028;}
      	else if(ht>   0 && ht<= 200 && met> 155 && met<= 160) {eff=0.529; errup=0.013; errdown=0.013;}
      	else if(ht> 200 && ht<= 600 && met> 155 && met<= 160) {eff=0.584; errup=0.004; errdown=0.004;}
      	else if(ht> 600 && ht<= 800 && met> 155 && met<= 160) {eff=0.540; errup=0.016; errdown=0.016;}
      	else if(ht> 800 && ht<=1000 && met> 155 && met<= 160) {eff=0.440; errup=0.026; errdown=0.026;}
      	else if(ht>1000 && ht<=9999 && met> 155 && met<= 160) {eff=0.490; errup=0.031; errdown=0.031;}
      	else if(ht>   0 && ht<= 200 && met> 160 && met<= 165) {eff=0.552; errup=0.014; errdown=0.014;}
      	else if(ht> 200 && ht<= 600 && met> 160 && met<= 165) {eff=0.640; errup=0.005; errdown=0.005;}
      	else if(ht> 600 && ht<= 800 && met> 160 && met<= 165) {eff=0.598; errup=0.016; errdown=0.016;}
      	else if(ht> 800 && ht<=1000 && met> 160 && met<= 165) {eff=0.496; errup=0.028; errdown=0.027;}
      	else if(ht>1000 && ht<=9999 && met> 160 && met<= 165) {eff=0.441; errup=0.033; errdown=0.033;}
      	else if(ht>   0 && ht<= 200 && met> 165 && met<= 170) {eff=0.594; errup=0.016; errdown=0.016;}
      	else if(ht> 200 && ht<= 600 && met> 165 && met<= 170) {eff=0.686; errup=0.005; errdown=0.005;}
      	else if(ht> 600 && ht<= 800 && met> 165 && met<= 170) {eff=0.611; errup=0.016; errdown=0.016;}
      	else if(ht> 800 && ht<=1000 && met> 165 && met<= 170) {eff=0.613; errup=0.027; errdown=0.028;}
      	else if(ht>1000 && ht<=9999 && met> 165 && met<= 170) {eff=0.548; errup=0.031; errdown=0.032;}
      	else if(ht>   0 && ht<= 200 && met> 170 && met<= 175) {eff=0.635; errup=0.018; errdown=0.018;}
      	else if(ht> 200 && ht<= 600 && met> 170 && met<= 175) {eff=0.730; errup=0.005; errdown=0.005;}
      	else if(ht> 600 && ht<= 800 && met> 170 && met<= 175) {eff=0.665; errup=0.017; errdown=0.017;}
      	else if(ht> 800 && ht<=1000 && met> 170 && met<= 175) {eff=0.621; errup=0.026; errdown=0.027;}
      	else if(ht>1000 && ht<=9999 && met> 170 && met<= 175) {eff=0.511; errup=0.036; errdown=0.036;}
      	else if(ht>   0 && ht<= 200 && met> 175 && met<= 180) {eff=0.630; errup=0.019; errdown=0.019;}
      	else if(ht> 200 && ht<= 600 && met> 175 && met<= 180) {eff=0.770; errup=0.005; errdown=0.005;}
      	else if(ht> 600 && ht<= 800 && met> 175 && met<= 180) {eff=0.670; errup=0.016; errdown=0.017;}
      	else if(ht> 800 && ht<=1000 && met> 175 && met<= 180) {eff=0.683; errup=0.028; errdown=0.030;}
      	else if(ht>1000 && ht<=9999 && met> 175 && met<= 180) {eff=0.574; errup=0.034; errdown=0.035;}
      	else if(ht>   0 && ht<= 200 && met> 180 && met<= 185) {eff=0.715; errup=0.020; errdown=0.021;}
      	else if(ht> 200 && ht<= 600 && met> 180 && met<= 185) {eff=0.800; errup=0.005; errdown=0.005;}
      	else if(ht> 600 && ht<= 800 && met> 180 && met<= 185) {eff=0.731; errup=0.016; errdown=0.016;}
      	else if(ht> 800 && ht<=1000 && met> 180 && met<= 185) {eff=0.713; errup=0.027; errdown=0.028;}
      	else if(ht>1000 && ht<=9999 && met> 180 && met<= 185) {eff=0.648; errup=0.036; errdown=0.038;}
      	else if(ht>   0 && ht<= 200 && met> 185 && met<= 190) {eff=0.701; errup=0.023; errdown=0.024;}
      	else if(ht> 200 && ht<= 600 && met> 185 && met<= 190) {eff=0.831; errup=0.005; errdown=0.005;}
      	else if(ht> 600 && ht<= 800 && met> 185 && met<= 190) {eff=0.766; errup=0.015; errdown=0.016;}
      	else if(ht> 800 && ht<=1000 && met> 185 && met<= 190) {eff=0.742; errup=0.027; errdown=0.029;}
      	else if(ht>1000 && ht<=9999 && met> 185 && met<= 190) {eff=0.692; errup=0.032; errdown=0.034;}
      	else if(ht>   0 && ht<= 200 && met> 190 && met<= 195) {eff=0.767; errup=0.024; errdown=0.026;}
      	else if(ht> 200 && ht<= 600 && met> 190 && met<= 195) {eff=0.857; errup=0.005; errdown=0.005;}
      	else if(ht> 600 && ht<= 800 && met> 190 && met<= 195) {eff=0.833; errup=0.014; errdown=0.015;}
      	else if(ht> 800 && ht<=1000 && met> 190 && met<= 195) {eff=0.760; errup=0.028; errdown=0.030;}
      	else if(ht>1000 && ht<=9999 && met> 190 && met<= 195) {eff=0.734; errup=0.031; errdown=0.033;}
      	else if(ht>   0 && ht<= 200 && met> 195 && met<= 200) {eff=0.699; errup=0.030; errdown=0.032;}
      	else if(ht> 200 && ht<= 600 && met> 195 && met<= 200) {eff=0.868; errup=0.005; errdown=0.005;}
      	else if(ht> 600 && ht<= 800 && met> 195 && met<= 200) {eff=0.776; errup=0.017; errdown=0.017;}
      	else if(ht> 800 && ht<=1000 && met> 195 && met<= 200) {eff=0.774; errup=0.027; errdown=0.029;}
      	else if(ht>1000 && ht<=9999 && met> 195 && met<= 200) {eff=0.713; errup=0.033; errdown=0.035;}
      	else if(ht>   0 && ht<= 200 && met> 200 && met<= 210) {eff=0.774; errup=0.021; errdown=0.023;}
      	else if(ht> 200 && ht<= 600 && met> 200 && met<= 210) {eff=0.897; errup=0.003; errdown=0.003;}
      	else if(ht> 600 && ht<= 800 && met> 200 && met<= 210) {eff=0.865; errup=0.010; errdown=0.011;}
      	else if(ht> 800 && ht<=1000 && met> 200 && met<= 210) {eff=0.823; errup=0.018; errdown=0.019;}
      	else if(ht>1000 && ht<=9999 && met> 200 && met<= 210) {eff=0.726; errup=0.024; errdown=0.026;}
      	else if(ht>   0 && ht<= 200 && met> 210 && met<= 220) {eff=0.799; errup=0.028; errdown=0.030;}
      	else if(ht> 200 && ht<= 600 && met> 210 && met<= 220) {eff=0.926; errup=0.003; errdown=0.003;}
      	else if(ht> 600 && ht<= 800 && met> 210 && met<= 220) {eff=0.888; errup=0.010; errdown=0.010;}
      	else if(ht> 800 && ht<=1000 && met> 210 && met<= 220) {eff=0.864; errup=0.017; errdown=0.019;}
      	else if(ht>1000 && ht<=9999 && met> 210 && met<= 220) {eff=0.793; errup=0.022; errdown=0.024;}
      	else if(ht>   0 && ht<= 200 && met> 220 && met<= 230) {eff=0.828; errup=0.027; errdown=0.031;}
      	else if(ht> 200 && ht<= 600 && met> 220 && met<= 230) {eff=0.940; errup=0.003; errdown=0.003;}
      	else if(ht> 600 && ht<= 800 && met> 220 && met<= 230) {eff=0.898; errup=0.010; errdown=0.011;}
      	else if(ht> 800 && ht<=1000 && met> 220 && met<= 230) {eff=0.906; errup=0.015; errdown=0.017;}
      	else if(ht>1000 && ht<=9999 && met> 220 && met<= 230) {eff=0.845; errup=0.020; errdown=0.023;}
      	else if(ht>   0 && ht<= 200 && met> 230 && met<= 240) {eff=0.824; errup=0.038; errdown=0.044;}
      	else if(ht> 200 && ht<= 600 && met> 230 && met<= 240) {eff=0.956; errup=0.003; errdown=0.003;}
      	else if(ht> 600 && ht<= 800 && met> 230 && met<= 240) {eff=0.924; errup=0.009; errdown=0.010;}
      	else if(ht> 800 && ht<=1000 && met> 230 && met<= 240) {eff=0.919; errup=0.015; errdown=0.017;}
      	else if(ht>1000 && ht<=9999 && met> 230 && met<= 240) {eff=0.900; errup=0.019; errdown=0.022;}
      	else if(ht>   0 && ht<= 200 && met> 240 && met<= 250) {eff=0.756; errup=0.050; errdown=0.057;}
      	else if(ht> 200 && ht<= 600 && met> 240 && met<= 250) {eff=0.967; errup=0.003; errdown=0.003;}
      	else if(ht> 600 && ht<= 800 && met> 240 && met<= 250) {eff=0.961; errup=0.007; errdown=0.008;}
      	else if(ht> 800 && ht<=1000 && met> 240 && met<= 250) {eff=0.940; errup=0.014; errdown=0.018;}
      	else if(ht>1000 && ht<=9999 && met> 240 && met<= 250) {eff=0.883; errup=0.023; errdown=0.027;}
      	else if(ht>   0 && ht<= 200 && met> 250 && met<= 275) {eff=0.857; errup=0.034; errdown=0.041;}
      	else if(ht> 200 && ht<= 600 && met> 250 && met<= 275) {eff=0.973; errup=0.002; errdown=0.002;}
      	else if(ht> 600 && ht<= 800 && met> 250 && met<= 275) {eff=0.968; errup=0.005; errdown=0.005;}
      	else if(ht> 800 && ht<=1000 && met> 250 && met<= 275) {eff=0.969; errup=0.007; errdown=0.009;}
      	else if(ht>1000 && ht<=9999 && met> 250 && met<= 275) {eff=0.893; errup=0.014; errdown=0.015;}
      	else if(ht>   0 && ht<= 200 && met> 275 && met<= 300) {eff=0.877; errup=0.044; errdown=0.060;}
      	else if(ht> 200 && ht<= 600 && met> 275 && met<= 300) {eff=0.985; errup=0.002; errdown=0.002;}
      	else if(ht> 600 && ht<= 800 && met> 275 && met<= 300) {eff=0.969; errup=0.005; errdown=0.006;}
      	else if(ht> 800 && ht<=1000 && met> 275 && met<= 300) {eff=0.958; errup=0.009; errdown=0.011;}
      	else if(ht>1000 && ht<=9999 && met> 275 && met<= 300) {eff=0.943; errup=0.012; errdown=0.014;}
      	else if(ht>   0 && ht<= 200 && met> 300 && met<=9999) {eff=0.804; errup=0.061; errdown=0.076;}
      	else if(ht> 200 && ht<= 600 && met> 300 && met<=9999) {eff=0.990; errup=0.001; errdown=0.001;}
      	else if(ht> 600 && ht<= 800 && met> 300 && met<=9999) {eff=0.995; errup=0.001; errdown=0.001;}
      	else if(ht> 800 && ht<=1000 && met> 300 && met<=9999) {eff=0.991; errup=0.002; errdown=0.003;}
      	else if(ht>1000 && ht<=9999 && met> 300 && met<=9999) {eff=0.971; errup=0.004; errdown=0.004;}

        // //////// First 4.3 ifb (FULL STATUS)
        // if(b.met()<=100) {eff = 0.;}
        // if(b.met()> 100 && b.met()<= 125) {eff = 0.153;}
        // if(b.met()> 125 && b.met()<= 150) {eff = 0.405;}
        // if(b.met()> 150 && b.met()<= 175) {eff = 0.684;}
        // if(b.met()> 175 && b.met()<= 200) {eff = 0.863;}
        // if(b.met()> 200 && b.met()<= 225) {eff = 0.939;}
        // if(b.met()> 225 && b.met()<= 250) {eff = 0.967;}
        // if(b.met()> 250 && b.met()<= 275) {eff = 0.986;}
        // if(b.met()> 275 && b.met()<= 300) {eff = 0.985;}
        // if(b.met()> 300 && b.met()<=9999) {eff = 0.988;}
        // //////// First 4.3 ifb (FULL STATUS)
      }

      //// MET || Ele27 || Ele105 || Ele115
      //// "(trig[13]||trig[33]||trig[14]||trig[15]||trig[30]||trig[31]||trig[22]||trig[40]||trig[24]||trig[41])"
    } else if(b.nels()==1 && b.nmus()==0){
      vector<float> leps_pt; 
      if (b.leps_pt()->size()>0) leps_pt.push_back(b.leps_pt()->at(0));
      else leps_pt.push_back(0);
      if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met> 100 && met<= 110) {eff=0.160; errup=0.019; errdown = 0.017;}
      else if(leps_pt[0]> 25 && leps_pt[0]<=  30 && met>100 && met<= 110) {eff=0.400; errup=0.024; errdown=0.024;}
      else if(leps_pt[0]> 30 && leps_pt[0]<= 110 && met>100 && met<= 110) {eff=0.728; errup=0.006; errdown=0.006;}
      else if(leps_pt[0]>110 && leps_pt[0]<= 120 && met>100 && met<= 110) {eff=0.880; errup=0.017; errdown=0.019;}
      else if(leps_pt[0]>120 && leps_pt[0]<=9999 && met>100 && met<= 110) {eff=0.950; errup=0.003; errdown=0.003;}
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
      if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met> 100 && met<= 110) {eff=0.271; errup=0.017; errdown = 0.016;}
      else if(leps_pt[0]>25 && leps_pt[0]<=  30 && met> 100 && met<= 110) {eff=0.725; errup=0.017; errdown=0.018;}
      else if(leps_pt[0]>30 && leps_pt[0]<=  50 && met> 100 && met<= 110) {eff=0.814; errup=0.008; errdown=0.009;}
      else if(leps_pt[0]>50 && leps_pt[0]<=9999 && met> 100 && met<= 110) {eff=0.964; errup=0.002; errdown=0.002;}
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

}
