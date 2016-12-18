#include "hig/hig_functions.hpp"

#include "TVector2.h"

#include "core/utilities.hpp"
#include "core/config_parser.hpp"

using namespace std;

namespace Functions{

  //// Efficiency of the MET[100||110||120] triggers in all 36.2 ifb
const NamedFunc eff_higtrig("eff_higtrig", [](const Baby &b) -> NamedFunc::ScalarType{
    float errup, errdown; // Not used, but for reference
    float eff = 1., met = b.met();
    vector<float> leps_pt({b.leps_pt()[0]});

    if(b.type() < 1000) eff = 1;

    //// Efficiency of the MET[100||110||120] triggers in all 36.2 ifb
    //// "(trig[13]||trig[33]||trig[14]||trig[15]||trig[30]||trig[31])";
    else if(b.nleps()==0){
      if(b.type()>=7000 && b.type()<8000) { // FAKE MET (QCD)
	if(met> 150 && met<= 155) {eff = 0.328; errup = 0.001; errdown = 0.001;}
	if(met> 155 && met<= 160) {eff = 0.358; errup = 0.001; errdown = 0.001;}
	if(met> 160 && met<= 165) {eff = 0.394; errup = 0.002; errdown = 0.002;}
	if(met> 165 && met<= 170) {eff = 0.428; errup = 0.002; errdown = 0.002;}
	if(met> 170 && met<= 175) {eff = 0.459; errup = 0.002; errdown = 0.002;}
	if(met> 175 && met<= 180) {eff = 0.489; errup = 0.002; errdown = 0.002;}
	if(met> 180 && met<= 185) {eff = 0.516; errup = 0.002; errdown = 0.002;}
	if(met> 185 && met<= 190) {eff = 0.542; errup = 0.002; errdown = 0.002;}
	if(met> 190 && met<= 195) {eff = 0.569; errup = 0.003; errdown = 0.003;}
	if(met> 195 && met<= 200) {eff = 0.594; errup = 0.003; errdown = 0.003;}
	if(met> 200 && met<= 205) {eff = 0.612; errup = 0.003; errdown = 0.003;}
	if(met> 205 && met<= 210) {eff = 0.625; errup = 0.003; errdown = 0.003;}
	if(met> 210 && met<= 215) {eff = 0.649; errup = 0.003; errdown = 0.003;}
	if(met> 215 && met<= 220) {eff = 0.667; errup = 0.003; errdown = 0.003;}
	if(met> 220 && met<= 225) {eff = 0.680; errup = 0.003; errdown = 0.003;}
	if(met> 225 && met<= 230) {eff = 0.687; errup = 0.004; errdown = 0.004;}
	if(met> 230 && met<= 235) {eff = 0.708; errup = 0.004; errdown = 0.004;}
	if(met> 235 && met<= 240) {eff = 0.711; errup = 0.004; errdown = 0.004;}
	if(met> 240 && met<= 245) {eff = 0.730; errup = 0.004; errdown = 0.004;}
	if(met> 245 && met<= 250) {eff = 0.732; errup = 0.004; errdown = 0.004;}
	if(met> 250 && met<= 255) {eff = 0.744; errup = 0.004; errdown = 0.004;}
	if(met> 255 && met<= 260) {eff = 0.751; errup = 0.005; errdown = 0.005;}
	if(met> 260 && met<= 265) {eff = 0.766; errup = 0.005; errdown = 0.005;}
	if(met> 265 && met<= 270) {eff = 0.762; errup = 0.005; errdown = 0.005;}
	if(met> 270 && met<= 275) {eff = 0.772; errup = 0.005; errdown = 0.005;}
	if(met> 275 && met<= 280) {eff = 0.773; errup = 0.005; errdown = 0.005;}
	if(met> 280 && met<= 285) {eff = 0.771; errup = 0.006; errdown = 0.006;}
	if(met> 285 && met<= 290) {eff = 0.777; errup = 0.006; errdown = 0.006;}
	if(met> 290 && met<= 295) {eff = 0.782; errup = 0.006; errdown = 0.006;}
	if(met> 295 && met<= 300) {eff = 0.798; errup = 0.006; errdown = 0.006;}
	if(met> 300 && met<=9999) {eff = 0.815; errup = 0.002; errdown = 0.002;}

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
	if(met> 100 && met<= 105) {eff = 0.074; errup = 0.001; errdown = 0.001;}
	if(met> 105 && met<= 110) {eff = 0.097; errup = 0.001; errdown = 0.001;}
	if(met> 110 && met<= 115) {eff = 0.123; errup = 0.001; errdown = 0.001;}
	if(met> 115 && met<= 120) {eff = 0.157; errup = 0.002; errdown = 0.002;}
	if(met> 120 && met<= 125) {eff = 0.193; errup = 0.002; errdown = 0.002;}
	if(met> 125 && met<= 130) {eff = 0.241; errup = 0.002; errdown = 0.002;}
	if(met> 130 && met<= 135) {eff = 0.292; errup = 0.003; errdown = 0.003;}
	if(met> 135 && met<= 140) {eff = 0.340; errup = 0.003; errdown = 0.003;}
	if(met> 140 && met<= 145) {eff = 0.393; errup = 0.003; errdown = 0.003;}
	if(met> 145 && met<= 150) {eff = 0.453; errup = 0.004; errdown = 0.003;}
	if(met> 150 && met<= 155) {eff = 0.507; errup = 0.004; errdown = 0.004;}
	if(met> 155 && met<= 160) {eff = 0.567; errup = 0.004; errdown = 0.004;}
	if(met> 160 && met<= 165) {eff = 0.619; errup = 0.004; errdown = 0.004;}
	if(met> 165 && met<= 170) {eff = 0.666; errup = 0.004; errdown = 0.004;}
	if(met> 170 && met<= 175) {eff = 0.706; errup = 0.004; errdown = 0.004;}
	if(met> 175 && met<= 180) {eff = 0.742; errup = 0.004; errdown = 0.004;}
	if(met> 180 && met<= 185) {eff = 0.780; errup = 0.004; errdown = 0.004;}
	if(met> 185 && met<= 190) {eff = 0.809; errup = 0.004; errdown = 0.004;}
	if(met> 190 && met<= 195) {eff = 0.841; errup = 0.004; errdown = 0.004;}
	if(met> 195 && met<= 200) {eff = 0.841; errup = 0.004; errdown = 0.005;}
	if(met> 200 && met<= 205) {eff = 0.872; errup = 0.004; errdown = 0.004;}
	if(met> 205 && met<= 210) {eff = 0.888; errup = 0.004; errdown = 0.004;}
	if(met> 210 && met<= 215) {eff = 0.903; errup = 0.004; errdown = 0.004;}
	if(met> 215 && met<= 220) {eff = 0.917; errup = 0.004; errdown = 0.004;}
	if(met> 220 && met<= 225) {eff = 0.920; errup = 0.004; errdown = 0.004;}
	if(met> 225 && met<= 230) {eff = 0.929; errup = 0.004; errdown = 0.004;}
	if(met> 230 && met<= 235) {eff = 0.937; errup = 0.004; errdown = 0.004;}
	if(met> 235 && met<= 240) {eff = 0.947; errup = 0.004; errdown = 0.004;}
	if(met> 240 && met<= 245) {eff = 0.952; errup = 0.004; errdown = 0.004;}
	if(met> 245 && met<= 250) {eff = 0.957; errup = 0.004; errdown = 0.004;}
	if(met> 250 && met<= 255) {eff = 0.955; errup = 0.004; errdown = 0.005;}
	if(met> 255 && met<= 260) {eff = 0.960; errup = 0.004; errdown = 0.005;}
	if(met> 260 && met<= 265) {eff = 0.968; errup = 0.004; errdown = 0.004;}
	if(met> 265 && met<= 270) {eff = 0.971; errup = 0.004; errdown = 0.004;}
	if(met> 270 && met<= 275) {eff = 0.971; errup = 0.004; errdown = 0.005;}
	if(met> 275 && met<= 280) {eff = 0.972; errup = 0.004; errdown = 0.005;}
	if(met> 280 && met<= 285) {eff = 0.973; errup = 0.004; errdown = 0.005;}
	if(met> 285 && met<= 290) {eff = 0.970; errup = 0.005; errdown = 0.006;}
	if(met> 290 && met<= 295) {eff = 0.981; errup = 0.004; errdown = 0.005;}
	if(met> 295 && met<= 300) {eff = 0.976; errup = 0.005; errdown = 0.006;}
	if(met> 300 && met<=9999) {eff = 0.984; errup = 0.001; errdown = 0.001;}

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
      if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met> 150 && met<= 160) {eff = 0.643; errup=0.037; errdown=0.039;}
      if(leps_pt[0]>  25 && leps_pt[0]<=  30 && met> 150 && met<= 160) {eff = 0.738; errup=0.033; errdown=0.036;}
      if(leps_pt[0]>  30 && leps_pt[0]<= 110 && met> 150 && met<= 160) {eff = 0.871; errup=0.007; errdown=0.007;}
      if(leps_pt[0]> 110 && leps_pt[0]<= 120 && met> 150 && met<= 160) {eff = 0.935; errup=0.016; errdown=0.021;}
      if(leps_pt[0]> 120 && leps_pt[0]<=9999 && met> 150 && met<= 160) {eff = 0.982; errup=0.003; errdown=0.003;}
      if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met> 160 && met<= 170) {eff = 0.760; errup=0.034; errdown=0.038;}
      if(leps_pt[0]>  25 && leps_pt[0]<=  30 && met> 160 && met<= 170) {eff = 0.792; errup=0.034; errdown=0.038;}
      if(leps_pt[0]>  30 && leps_pt[0]<= 110 && met> 160 && met<= 170) {eff = 0.910; errup=0.007; errdown=0.007;}
      if(leps_pt[0]> 110 && leps_pt[0]<= 120 && met> 160 && met<= 170) {eff = 0.970; errup=0.013; errdown=0.020;}
      if(leps_pt[0]> 120 && leps_pt[0]<=9999 && met> 160 && met<= 170) {eff = 0.981; errup=0.003; errdown=0.004;}
      if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met> 170 && met<= 180) {eff = 0.829; errup=0.033; errdown=0.038;}
      if(leps_pt[0]>  25 && leps_pt[0]<=  30 && met> 170 && met<= 180) {eff = 0.863; errup=0.031; errdown=0.037;}
      if(leps_pt[0]>  30 && leps_pt[0]<= 110 && met> 170 && met<= 180) {eff = 0.937; errup=0.006; errdown=0.006;}
      if(leps_pt[0]> 110 && leps_pt[0]<= 120 && met> 170 && met<= 180) {eff = 0.988; errup=0.008; errdown=0.016;}
      if(leps_pt[0]> 120 && leps_pt[0]<=9999 && met> 170 && met<= 180) {eff = 0.982; errup=0.003; errdown=0.004;}
      if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met> 180 && met<= 190) {eff = 0.761; errup=0.041; errdown=0.046;}
      if(leps_pt[0]>  25 && leps_pt[0]<=  30 && met> 180 && met<= 190) {eff = 0.863; errup=0.032; errdown=0.038;}
      if(leps_pt[0]>  30 && leps_pt[0]<= 110 && met> 180 && met<= 190) {eff = 0.939; errup=0.006; errdown=0.007;}
      if(leps_pt[0]> 110 && leps_pt[0]<= 120 && met> 180 && met<= 190) {eff = 0.969; errup=0.015; errdown=0.024;}
      if(leps_pt[0]> 120 && leps_pt[0]<=9999 && met> 180 && met<= 190) {eff = 0.984; errup=0.003; errdown=0.004;}
      if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met> 190 && met<= 200) {eff = 0.892; errup=0.033; errdown=0.042;}
      if(leps_pt[0]>  25 && leps_pt[0]<=  30 && met> 190 && met<= 200) {eff = 0.902; errup=0.030; errdown=0.039;}
      if(leps_pt[0]>  30 && leps_pt[0]<= 110 && met> 190 && met<= 200) {eff = 0.956; errup=0.006; errdown=0.006;}
      if(leps_pt[0]> 110 && leps_pt[0]<= 120 && met> 190 && met<= 200) {eff = 0.983; errup=0.011; errdown=0.022;}
      if(leps_pt[0]> 120 && leps_pt[0]<=9999 && met> 190 && met<= 200) {eff = 0.993; errup=0.002; errdown=0.003;}
      if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met> 200 && met<= 210) {eff = 0.950; errup=0.021; errdown=0.032;}
      if(leps_pt[0]>  25 && leps_pt[0]<=  30 && met> 200 && met<= 210) {eff = 0.951; errup=0.023; errdown=0.037;}
      if(leps_pt[0]>  30 && leps_pt[0]<= 110 && met> 200 && met<= 210) {eff = 0.973; errup=0.005; errdown=0.005;}
      if(leps_pt[0]> 110 && leps_pt[0]<= 120 && met> 200 && met<= 210) {eff = 1.000; errup=0.000; errdown=0.018;}
      if(leps_pt[0]> 120 && leps_pt[0]<=9999 && met> 200 && met<= 210) {eff = 0.985; errup=0.003; errdown=0.004;}
      if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met> 210 && met<=9999) {eff = 0.974; errup=0.005; errdown=0.006;}
      if(leps_pt[0]>  25 && leps_pt[0]<=  30 && met> 210 && met<=9999) {eff = 0.981; errup=0.004; errdown=0.005;}
      if(leps_pt[0]>  30 && leps_pt[0]<= 110 && met> 210 && met<=9999) {eff = 0.992; errup=0.001; errdown=0.001;}
      if(leps_pt[0]> 110 && leps_pt[0]<= 120 && met> 210 && met<=9999) {eff = 0.997; errup=0.002; errdown=0.003;}
      if(leps_pt[0]> 120 && leps_pt[0]<=9999 && met> 210 && met<=9999) {eff = 0.996; errup=0.001; errdown=0.001;}

      //// MET || Mu24 || Mu50
      //// "(trig[13]||trig[33]||trig[14]||trig[15]||trig[30]||trig[31]||trig[19]||trig[55]||trig[21])"
    } else if(b.nels()==0 && b.nmus()==1){
      if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met> 150 && met<= 160) {eff = 0.778; errup=0.024; errdown=0.026;}
      if(leps_pt[0]>  25 && leps_pt[0]<=  30 && met> 150 && met<= 160) {eff = 0.915; errup=0.016; errdown=0.019;}
      if(leps_pt[0]>  30 && leps_pt[0]<=  50 && met> 150 && met<= 160) {eff = 0.940; errup=0.008; errdown=0.009;}
      if(leps_pt[0]>  50 && leps_pt[0]<=9999 && met> 150 && met<= 160) {eff = 0.984; errup=0.002; errdown=0.002;}
      if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met> 160 && met<= 170) {eff = 0.798; errup=0.026; errdown=0.029;}
      if(leps_pt[0]>  25 && leps_pt[0]<=  30 && met> 160 && met<= 170) {eff = 0.946; errup=0.015; errdown=0.020;}
      if(leps_pt[0]>  30 && leps_pt[0]<=  50 && met> 160 && met<= 170) {eff = 0.967; errup=0.006; errdown=0.008;}
      if(leps_pt[0]>  50 && leps_pt[0]<=9999 && met> 160 && met<= 170) {eff = 0.991; errup=0.002; errdown=0.002;}
      if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met> 170 && met<= 180) {eff = 0.885; errup=0.022; errdown=0.025;}
      if(leps_pt[0]>  25 && leps_pt[0]<=  30 && met> 170 && met<= 180) {eff = 0.937; errup=0.016; errdown=0.021;}
      if(leps_pt[0]>  30 && leps_pt[0]<=  50 && met> 170 && met<= 180) {eff = 0.977; errup=0.006; errdown=0.007;}
      if(leps_pt[0]>  50 && leps_pt[0]<=9999 && met> 170 && met<= 180) {eff = 0.987; errup=0.002; errdown=0.002;}
      if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met> 180 && met<= 190) {eff = 0.927; errup=0.019; errdown=0.024;}
      if(leps_pt[0]>  25 && leps_pt[0]<=  30 && met> 180 && met<= 190) {eff = 0.958; errup=0.014; errdown=0.019;}
      if(leps_pt[0]>  30 && leps_pt[0]<=  50 && met> 180 && met<= 190) {eff = 0.974; errup=0.006; errdown=0.008;}
      if(leps_pt[0]>  50 && leps_pt[0]<=9999 && met> 180 && met<= 190) {eff = 0.992; errup=0.002; errdown=0.002;}
      if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met> 190 && met<= 200) {eff = 0.921; errup=0.019; errdown=0.024;}
      if(leps_pt[0]>  25 && leps_pt[0]<=  30 && met> 190 && met<= 200) {eff = 0.965; errup=0.014; errdown=0.020;}
      if(leps_pt[0]>  30 && leps_pt[0]<=  50 && met> 190 && met<= 200) {eff = 0.991; errup=0.004; errdown=0.006;}
      if(leps_pt[0]>  50 && leps_pt[0]<=9999 && met> 190 && met<= 200) {eff = 0.991; errup=0.002; errdown=0.002;}
      if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met> 200 && met<= 210) {eff = 0.926; errup=0.022; errdown=0.028;}
      if(leps_pt[0]>  25 && leps_pt[0]<=  30 && met> 200 && met<= 210) {eff = 0.994; errup=0.005; errdown=0.015;}
      if(leps_pt[0]>  30 && leps_pt[0]<=  50 && met> 200 && met<= 210) {eff = 0.994; errup=0.003; errdown=0.006;}
      if(leps_pt[0]>  50 && leps_pt[0]<=9999 && met> 200 && met<= 210) {eff = 0.994; errup=0.002; errdown=0.002;}
      if(leps_pt[0]>  20 && leps_pt[0]<=  25 && met> 210 && met<=9999) {eff = 0.981; errup=0.004; errdown=0.004;}
      if(leps_pt[0]>  25 && leps_pt[0]<=  30 && met> 210 && met<=9999) {eff = 0.994; errup=0.002; errdown=0.003;}
      if(leps_pt[0]>  30 && leps_pt[0]<=  50 && met> 210 && met<=9999) {eff = 0.996; errup=0.001; errdown=0.001;}
      if(leps_pt[0]>  50 && leps_pt[0]<=9999 && met> 210 && met<=9999) {eff = 0.997; errup=0.000; errdown=0.000;}

      //// Ele27 || Ele105 || Ele115
      //// "(trig[22]||trig[40]||trig[24]||trig[41])"
    } else if(b.nels()==2 && b.nmus()==0){
      if(leps_pt[0]>  40 && leps_pt[0]<=  45) {eff = 0.944; errup = 0.015; errdown = 0.019;}
      if(leps_pt[0]>  45 && leps_pt[0]<=  50) {eff = 0.910; errup = 0.015; errdown = 0.017;}
      if(leps_pt[0]>  50 && leps_pt[0]<=  55) {eff = 0.927; errup = 0.013; errdown = 0.015;}
      if(leps_pt[0]>  55 && leps_pt[0]<=  60) {eff = 0.912; errup = 0.013; errdown = 0.015;}
      if(leps_pt[0]>  60 && leps_pt[0]<=  65) {eff = 0.941; errup = 0.011; errdown = 0.013;}
      if(leps_pt[0]>  65 && leps_pt[0]<=  70) {eff = 0.901; errup = 0.014; errdown = 0.016;}
      if(leps_pt[0]>  70 && leps_pt[0]<=  75) {eff = 0.921; errup = 0.013; errdown = 0.016;}
      if(leps_pt[0]>  75 && leps_pt[0]<=  80) {eff = 0.947; errup = 0.011; errdown = 0.014;}
      if(leps_pt[0]>  80 && leps_pt[0]<=  85) {eff = 0.954; errup = 0.011; errdown = 0.013;}
      if(leps_pt[0]>  85 && leps_pt[0]<=  90) {eff = 0.939; errup = 0.012; errdown = 0.014;}
      if(leps_pt[0]>  90 && leps_pt[0]<=  95) {eff = 0.940; errup = 0.012; errdown = 0.015;}
      if(leps_pt[0]>  95 && leps_pt[0]<= 100) {eff = 0.932; errup = 0.014; errdown = 0.017;}
      if(leps_pt[0]> 100 && leps_pt[0]<= 105) {eff = 0.934; errup = 0.014; errdown = 0.017;}
      if(leps_pt[0]> 105 && leps_pt[0]<= 110) {eff = 0.965; errup = 0.010; errdown = 0.014;}
      if(leps_pt[0]> 110 && leps_pt[0]<=9999) {eff = 0.994; errup = 0.001; errdown = 0.001;}

      //// Mu24 || Mu50
      //// "(trig[19]||trig[55]||trig[21])"
    } else if(b.nels()==0 && b.nmus()==2){
      if(leps_pt[0]>  40 && leps_pt[0]<=  45) {eff = 0.959; errup = 0.010; errdown = 0.012;}
      if(leps_pt[0]>  45 && leps_pt[0]<=  50) {eff = 0.970; errup = 0.006; errdown = 0.007;}
      if(leps_pt[0]>  50 && leps_pt[0]<=9999) {eff = 0.982; errup = 0.000; errdown = 0.000;}
    }
    return eff;
  });

}
