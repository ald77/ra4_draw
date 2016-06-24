#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>  // setw

#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "RooStats/NumberCountingUtils.h"

#include "utilities.hpp"
#include "baby.hpp"
#include "process.hpp"
#include "named_func.hpp"
#include "plot_maker.hpp"
#include "palette.hpp"
#include "table.hpp"
using namespace std;

namespace{
  TString method = "met200";
  double lumi(0.815);
  bool do_other(false);
  bool full_lumi(false);
  float syst = 0.001;
  TString blind_s = "$\\spadesuit$";
  bool really_unblind = false;
  bool single_thread = false;
}

template<typename T>
shared_ptr<Process> Proc(const string process_name, Process::Type type,
                         int color, const set<string> &files, const string &cut = "1"){
  return make_shared<Process>(process_name, type, color,
                              unique_ptr<Baby>(new T(files)),
                              cut);
}

void GetOptions(int argc, char *argv[]);
TString Zbi(double Nobs, double Nbkg, double Ebkg){
  double Nsig = Nobs-Nbkg;
  double zbi = RooStats::NumberCountingUtils::BinomialExpZ(Nsig, Nbkg, Ebkg/Nbkg);
  if(Nbkg==0) zbi = RooStats::NumberCountingUtils::BinomialWithTauExpZ(Nsig, Nbkg, 1/Ebkg);
  if(zbi<0) zbi=0;
  TString zbi_s = "$"+RoundNumber(zbi,1)+"\\sigma$";
  if(Nsig<=0) zbi_s = "$-\\sigma$";
  //cout<<"Zbi for Nobs "<<Nobs<<", Nbkg "<<Nbkg<<", Ebkg "<<Ebkg<<" is "<<zbi_s<<endl;
  return zbi_s;
}

int main(int argc, char *argv[]){ 
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  ////// Defining processes
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_standard/");
  string folderdata(bfolder+"/cms2r0/babymaker/babies/2016_06_21/data/skim_standard/");
  if(method.Contains("met150")){
    foldermc = bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_1lht500met150nj5/";
    folderdata = bfolder+"/cms2r0/babymaker/babies/2016_06_14/data/merged_1lht500met150nj5/";
  }
  if(method.Contains("2015")){
    foldermc = bfolder+"/cms2r0/babymaker/babies/2016_04_29/mc/merged_1lht500met200/";
    folderdata = bfolder+"/cms2r0/babymaker/babies/2016_04_29/data/merged_1lht500met200/";
  }

  Palette colors("txt/colors.txt", "default");

  auto tt = Proc<Baby_full>("t#bar{t}", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets*Lept*.root", foldermc+"*_TTJets_HT*.root",
	foldermc+"*_WJetsToLNu*.root",foldermc+"*_ST_*.root",
	foldermc+"*_TTW*.root",foldermc+"*_TTZ*.root",
	foldermc+"*DYJetsToLL*.root",foldermc+"*QCD_HT*.root",
	foldermc+"*_ZJet*.root",foldermc+"*_ttHJetTobb*.root",
	foldermc+"*_TTGJets*.root",foldermc+"*_TTTT*.root",
	foldermc+"*_WH_HToBB*.root",foldermc+"*_ZH_HToBB*.root",
	foldermc+"*_WWTo*.root",foldermc+"*_WZ*.root",foldermc+"*_ZZ_*.root"},
    "stitch");
  auto other = Proc<Baby_full>("Other", Process::Type::background, colors("other"),
    {foldermc+"*_WJetsToLNu*",foldermc+"*_ST_*.root",
	foldermc+"*_TTW*.root",foldermc+"*_TTZ*.root",
	foldermc+"*DYJetsToLL*.root",foldermc+"*QCD_HT*.root",
	foldermc+"*_ZJet*.root",foldermc+"*_ttHJetTobb*.root",
	foldermc+"*_TTGJets*.root",foldermc+"*_TTTT*.root",
	foldermc+"*_WH_HToBB*.root",foldermc+"*_ZH_HToBB*.root",
	foldermc+"*_WWTo*.root",foldermc+"*_WZ*.root",foldermc+"*_ZZ_*.root"});

  string trigs = "(trig[4]||trig[8]||trig[13]||trig[33])";
  if(method.Contains("2015")) trigs = "(trig[4]||trig[8]||trig[28]||trig[14])";
  auto data = Proc<Baby_full>("Data", Process::Type::data, kBlack,
    {folderdata+"*.root"},trigs);

  vector<shared_ptr<Process> > all_procs = {data, tt};
  if (do_other) all_procs.push_back(other);


  ////// Defining cuts
  TString base_s = "mj14>250&&njets>=5&&stitch&&pass&&nonblind";

  vector<TString> njbcuts_stdnob = {"njets>=6&&njets<=8", "njets>=9", "njets>=6&&njets<=8", "njets>=9"}; 
  vector<TString> njbcuts_2l = {"njets>=5&&njets<=7", "njets>=8", "njets>=5&&njets<=7", "njets>=8"}; 
  vector<TString> njbcuts_5j = {"nbm==1&&njets==5", "nbm>=2&&njets==5", "nbm==1&&njets==5", "nbm>=2&&njets==5"}; 
  vector<TString> njbcuts_m1lmet150nb12 = {"nbm==1&&njets>=6&&njets<=8", "nbm==1&&njets>=9", 
					   "nbm>=2&&njets>=6&&njets<=8", "nbm>=2&&njets>=9"}; 
  vector<TString> njbcuts_std = {"nbm==1&&njets>=6&&njets<=8", "nbm==1&&njets>=9", 
				 "nbm==2&&njets>=6&&njets<=8", "nbm==2&&njets>=9", 
				 "nbm>=3&&njets>=6&&njets<=8", "nbm>=3&&njets>=9",
				 "nbm==1&&njets>=6&&njets<=8", "nbm==1&&njets>=9", 
				 "nbm==2&&njets>=6&&njets<=8", "nbm==2&&njets>=9", 
				 "nbm>=3&&njets>=6&&njets<=8", "nbm>=3&&njets>=9"}; 
  vector<TString> njbcuts_nb1 = {"nbm==1&&njets>=6&&njets<=8", "nbm==1&&njets>=9", 
				 "nbm==2&&njets>=6&&njets<=8", "nbm==2&&njets>=9", 
				 "nbm>=3&&njets>=6&&njets<=8", "nbm>=3&&njets>=9"}; 
  vector<TString> njbcuts_met500 = {"nbm==1&&njets>=6&&njets<=8", "nbm==1&&njets>=9", 
				    "nbm==2&&njets>=6&&njets<=8", "nbm==2&&njets>=9", 
				    "nbm>=3&&njets>=6&&njets<=8", "nbm>=3&&njets>=9"}; 

  vector<TString> njbcuts_m1lnob = {"njets>=6&&njets<=8", "njets>=9"}; 
  vector<TString> njbcuts_m2lnob = {"njets>=5&&njets<=7", "njets>=8"}; 
  vector<TString> njbcuts_2lveto_lonj = {"njets>=5", "njets>=5"}; 
  vector<TString> njbcuts_2lveto_hinj = {"njets>=8", "njets>=8"}; 


  size_t ilowmet(2); // njbcuts index up to which metcuts[0] is applied
  vector<TString> metcuts = {"met<=350", "met>350&&met<=500"};

  vector<TString> abcdcuts_std = {"mt<=140&&mj14<=400", 
				  "mt<=140&&mj14>400", 
				  "mt>140&&mj14<=400",          
				  "mt>140&&mj14>400"};
  vector<TString> abcdcuts_2lveto = {"mt<=140&&mj14<=400&&nleps==1&&nveto==0&&nbm>=1&&njets>=6", 
				     "mt<=140&&mj14>400&&nleps==1&&nveto==0&&nbm>=1&&njets>=6", 
				     "mj14<=400&&(nleps==2&&nbm<=2&&njets>=5 || nleps==1&&nveto==1&&nbm>=1&&nbm<=2&&mt>140&&njets>=6)",          
				     "mj14> 400&&(nleps==2&&nbm<=2&&njets>=5 || nleps==1&&nveto==1&&nbm>=1&&nbm<=2&&mt>140&&njets>=6)"};

  vector<TString> abcdcuts_2lveto_el = {"mt<=140&&mj14<=400&&nels==1&&nveto==0&&nbm>=1&&njets>=6", 
					"mt<=140&&mj14>400&&nels==1&&nveto==0&&nbm>=1&&njets>=6", 
					"mj14<=400&&(nels==2&&nbm<=2&&njets>=5 || nels==1&&nveto==1&&nbm>=1&&nbm<=2&&mt>140&&njets>=6)",          
					"mj14> 400&&(nels==2&&nbm<=2&&njets>=5 || nels==1&&nveto==1&&nbm>=1&&nbm<=2&&mt>140&&njets>=6)"};

  vector<TString> abcdcuts_2lveto_mu = {"mt<=140&&mj14<=400&&nmus==1&&nveto==0&&nbm>=1&&njets>=6", 
					"mt<=140&&mj14>400&&nmus==1&&nveto==0&&nbm>=1&&njets>=6", 
					"mj14<=400&&(nmus==2&&nbm<=2&&njets>=5 || nmus==1&&nveto==1&&nbm>=1&&nbm<=2&&mt>140&&njets>=6)",          
					"mj14> 400&&(nmus==2&&nbm<=2&&njets>=5 || nmus==1&&nveto==1&&nbm>=1&&nbm<=2&&mt>140&&njets>=6)"};

  vector<TString> abcdcuts_2lveto_lonj = {"mt<=140&&mj14<=400&&nleps==1&&nveto==0&&nbm>=1&&njets>=6", 
   					  "mt<=140&&mj14>400&&nleps==1&&nveto==0&&nbm>=1&&njets>=6&&njets<=8", 
   					  "mj14<=400&&(nleps==2&&nbm<=2&&njets>=5 || nleps==1&&nveto==1&&nbm>=1&&nbm<=2&&mt>140&&njets>=6)",          
   					  "mj14>400&&(nleps==2&&nbm<=2&&njets>=5&&njets<=7 || nleps==1&&nveto==1&&nbm>=1&&nbm<=2&&mt>140&&njets>=6&&njets<=8)"};


  vector<TString> abcdcuts_2lveto_hinj = {"mt<=140&&mj14<=400&&nleps==1&&nveto==0&&nbm>=1&&njets>=6", 
					  "mt<=140&&mj14>400&&nleps==1&&nveto==0&&nbm>=1&&njets>=9", 
					  "mj14<=400&&(nleps==2&&nbm<=2&&njets>=5 || nleps==1&&nveto==1&&nbm>=1&&nbm<=2&&mt>140&&njets>=6)",          
					  "mj14>400&&(nleps==2&&nbm<=2&&njets>=8 || nleps==1&&nveto==1&&nbm>=1&&nbm<=2&&mt>140&&njets>=9)"};

  vector<TString> abcdcuts_veto = {"mt<=140&&mj14<=400&&njets>=6&&nbm>=1&&nleps==1&&nveto==0", 
				   "mt<=140&&mj14>400&&nbm>=1&&nleps==1&&nveto==0", 
				   "mt>140&&mj14<=400&&njets>=6&&nbm>=1&&nbm<=2&&nleps==1&&nveto==1",          
				   "mt>140&&mj14>400&&njets>=6&&nbm>=1&&nbm<=2&&nleps==1&&nveto==1"};
  vector<TString> abcdcuts_2l = {"mt<=140&&mj14<=400&&njets>=6&&nbm>=1&&nleps==1&&nveto==0", 
				 "mt<=140&&mj14>400&&nbm>=1&&nleps==1&&nveto==0", 
				 "mj14<=400&&nbm<=2&&nleps==2",          
				 "mj14>400&&nbm<=2&&nleps==2"};


  vector<TString> abcdcuts, njbcuts_himt = njbcuts_stdnob;
  vector<TString> njbcuts = njbcuts_stdnob;
  TString region_s = "R", method_s, base_all = "mj14>250&&pass&&nonblind&&stitch&&";
  TString lumi_s = "815 pb$^{-1}$";
  bool unblind = false;

  if(full_lumi){
    base_all = "mj14>250&&pass&&stitch&&";
    if(method.Contains("2015")){
      lumi = 2.3;
      lumi_s = "2.3 fb$^{-1}$";
    } else {
      lumi = 2.6;
      lumi_s = "2.6 fb$^{-1}$";
    }
  }

  if(method=="m2l" || method=="m2l_2015") {
    base_s = base_all+"njets>=5";
    njbcuts_himt = njbcuts_2l;
    abcdcuts = abcdcuts_2l;
    region_s = "D";
    method_s = "D3 and D4 have $2\\ell$";
  } else if(method=="m2lveto" || method=="m2lveto_2015") {
    base_s = base_all+"njets>=5&&nleps>=1";
    njbcuts = njbcuts_2lveto_lonj;
    njbcuts_himt = njbcuts_2lveto_lonj;
    abcdcuts = abcdcuts_2lveto;
    region_s = "D";
    method_s = "D3 and D4 have $\\ell\\ell+\\ell t_{\\rm veto}$";
    ilowmet = 1;
  } else if(method=="m2lveto_el" || method=="m2lveto_el_2015") {
    base_s = base_all+"njets>=5&&nels>=1";
    njbcuts = njbcuts_2lveto_lonj;
    njbcuts_himt = njbcuts_2lveto_lonj;
    abcdcuts = abcdcuts_2lveto_el;
    region_s = "D";
    method_s = "D3 and D4 have $ee+e t_{\\rm veto}$";
    ilowmet = 1;
  } else if(method=="m2lveto_mu" || method=="m2lveto_mu_2015") {
    base_s = base_all+"njets>=5&&nmus>=1";
    njbcuts = njbcuts_2lveto_lonj;
    njbcuts_himt = njbcuts_2lveto_lonj;
    abcdcuts = abcdcuts_2lveto_mu;
    region_s = "D";
    method_s = "D3 and D4 have $\\mu\\mu+\\mu t_{\\rm veto}$";
    ilowmet = 1;
  } else if(method=="m2lveto_lonj") {
    base_s = base_all+"njets>=5";
    njbcuts = njbcuts_2lveto_lonj;
    njbcuts_himt = njbcuts_2lveto_lonj;
    abcdcuts = abcdcuts_2lveto_lonj;
    region_s = "D";
    method_s = "D3 and D4 have $\\ell\\ell+\\ell t_{\\rm veto}$ - low $N_{\\rm jets}$";
    ilowmet = 1;
  } else if(method=="m2lveto_hinj") {
    base_s = base_all+"njets>=5";
    njbcuts = njbcuts_2lveto_hinj;
    njbcuts_himt = njbcuts_2lveto_hinj;
    abcdcuts = abcdcuts_2lveto_hinj;
    region_s = "D";
    method_s = "D3 and D4 have $\\ell\\ell+\\ell t_{\\rm veto}$ - high $N_{\\rm jets}$";
    ilowmet = 1;
  } else if(method=="mveto" || method=="mveto_2015") {
    base_s = base_all+"njets>=6&&nbm>=1&&nleps==1";
    njbcuts_himt = njbcuts_stdnob;
    abcdcuts = abcdcuts_veto;
    region_s = "D";
    method_s = "D3 and D4 have $\\ell t_{\\rm veto}$";
  } else if(method=="m5j") {
    base_s = base_all+"njets==5&&nbm>=1&&nleps==1&&nveto==0";
    njbcuts = njbcuts_5j;
    njbcuts_himt = njbcuts_5j;
    abcdcuts = abcdcuts_std;
    method_s = "$N_{\\rm jets}=5$";
  } else if(method=="m1lmet150") {
    base_s = base_all+"njets>=6&&nbm>=1&&nleps==1&&nveto==0";
    njbcuts = njbcuts_m1lmet150nb12;
    njbcuts_himt = njbcuts_m1lmet150nb12;
    abcdcuts = abcdcuts_std;
    method_s = "$1\\ell, 150<E_{\\rm T}^{\\rm miss}<200$";
    ilowmet = njbcuts.size();
  } else if(method=="mvetomet150" || method=="mvetomet150_2015") {
    base_s = base_all+"njets>=6&&nbm>=1&&nleps==1";
    njbcuts = njbcuts_m1lnob;
    njbcuts_himt = njbcuts_m1lnob;
    abcdcuts = abcdcuts_veto;
    region_s = "D";
    method_s = "$150<E_{\\rm T}^{\\rm miss}<200$, D3 and D4 have $\\ell t_{\\rm veto}$";
  } else if(method=="m2lvetomet150" || method=="m2lvetomet150_2015") {
    base_s = base_all+"njets>=5";
    njbcuts = njbcuts_2lveto_lonj;
    njbcuts_himt = njbcuts_2lveto_lonj;
    abcdcuts = abcdcuts_2lveto;
    region_s = "D";
    method_s = "$150<E_{\\rm T}^{\\rm miss}<200$, D3 and D4 have $\\ell\\ell+\\ell t_{\\rm veto}$";
    ilowmet = 1;
  } else if(method=="m2lmet150" || method=="m2lmet150_2015") {
    base_s = base_all+"njets>=5";
    njbcuts = njbcuts_m1lnob;
    njbcuts_himt = njbcuts_m2lnob;
    abcdcuts = abcdcuts_2l;
    region_s = "D";
    method_s = "$150<E_{\\rm T}^{\\rm miss}<200$, D3 and D4 have $\\ell\\ell+\\ell t_{\\rm veto}$";
  } else if(method=="met200") {
    base_s = base_all+"njets>=6&&nbm>=1&&met>200&&nleps==1&&nveto==0";
    njbcuts = njbcuts_std;
    njbcuts_himt = njbcuts_std;
    abcdcuts = abcdcuts_std;
    method_s = "Signal bins - $1\\ell$, $200<E_{\\rm T}^{\\rm miss}\\leq500$";
    ilowmet = 6;
  } else if(method=="met500") {
    base_s = base_all+"njets>=6&&nbm>=1&&met>500&&nleps==1&&nveto==0";
    njbcuts = njbcuts_met500;
    njbcuts_himt = njbcuts_met500;
    abcdcuts = abcdcuts_std;
    method_s = "Signal bins - $1\\ell$, $E_{\\rm T}^{\\rm miss}>500$";
    metcuts[0] = "met>500";
    ilowmet = njbcuts.size();
  } else if(method=="met200nb1") {
    base_s = base_all+"njets>=6&&nbm>=1&&met>200&&nleps==1&&nveto==0";
    njbcuts = njbcuts_nb1;
    njbcuts_himt = njbcuts_nb1;
    abcdcuts = abcdcuts_std;
    method_s = "Signal bins - $1\\ell$, separated N_{b}=1$, $200<E_{\\rm T}^{\\rm miss}\\leq350$";
    metcuts[0] = "met>200&&met<=350&&nbm==1";
    metcuts[1] = "met>200&&met<=350&&nbm>=2";
    ilowmet = 2;
  } else if(method=="met350nb1") {
    base_s = base_all+"njets>=6&&nbm>=1&&met>350&&met<=500&&nleps==1&&nveto==0";
    njbcuts = njbcuts_nb1;
    njbcuts_himt = njbcuts_nb1;
    abcdcuts = abcdcuts_std;
    method_s = "Signal bins - $1\\ell$, separated N_{b}=1$, $350<E_{\\rm T}^{\\rm miss}\\leq500$";
    metcuts[0] = "met>350&&met<=500&&nbm==1";
    metcuts[1] = "met>350&&met<=500&&nbm>=2";
    ilowmet = 2;
  } else if(method=="met500nb1") {
    base_s = base_all+"njets>=6&&nbm>=1&&met>500&&nleps==1&&nveto==0";
    njbcuts = njbcuts_nb1;
    njbcuts_himt = njbcuts_nb1;
    abcdcuts = abcdcuts_std;
    method_s = "Signal bins - $1\\ell$, separated N_{b}=1$, $E_{\\rm T}^{\\rm miss}>500$";
    metcuts[0] = "met>500&&nbm==1";
    metcuts[1] = "met>500&&nbm>=2";
    ilowmet = 2;
  }else if(method=="agg_himet"){
    base_s = base_all+"met>500&&njets>=6&&nbm>=1";
    njbcuts = vector<TString>{"nbm>=3&&njets>=6"};
    njbcuts_himt = njbcuts;
    abcdcuts = abcdcuts_std;
    method_s = "Agg. Bin: $1\\ell$, MET500, $N_{j}\\geq6$, $N_{b}\\geq3$";
    metcuts = vector<TString>{"met>500"};
    ilowmet = 1;
  }else if(method=="agg_mixed"){
    base_s = base_all+"met>350&&njets>=6&&nbm>=1";
    njbcuts = vector<TString>{"nbm>=2&&njets>=9"};
    njbcuts_himt = njbcuts;
    abcdcuts = abcdcuts_std;
    method_s = "Agg. Bin: $1\\ell$, MET350, $N_{j}\\geq9$, $N_{b}\\geq2$";
    metcuts = vector<TString>{"met>350"};
    ilowmet = 1;
  }else if(method=="agg_himult"){
    base_s = base_all+"met>200&&njets>=6&&nbm>=1";
    njbcuts = vector<TString>{"nbm>=3&&njets>=9"};
    njbcuts_himt = njbcuts;
    abcdcuts = abcdcuts_std;
    method_s = "Agg. Bin: $1\\ell$, MET200, $N_{j}\\geq9$, $N_{b}\\geq3$";
    metcuts = vector<TString>{"met>200"};
    ilowmet = 1;
  }else if(method=="agg_1b"){
    base_s = base_all+"met>500&&njets>=6&&nbm>=1";
    njbcuts = vector<TString>{"nbm>=1&&njets>=9"};
    njbcuts_himt = njbcuts;
    abcdcuts = abcdcuts_std;
    method_s = "Agg. Bin: $1\\ell$, MET500, $N_{j}\\geq9$, $N_{b}\\geq1$";
    metcuts = vector<TString>{"met>500"};
    ilowmet = 1;
  }else {
    cout<<"Method "<<method<<" not available. Exiting"<<endl<<endl; 
    return 0;
  }
  if(method.Contains("2015")) method_s = "2015 data: "+method_s;

  if(really_unblind) unblind = true;

  ////// Finding yields 
  vector<TableRow> bincuts;
  for(size_t ind(0); ind<njbcuts.size(); ind++){
    //cout<<endl<<"New njb"<<endl;
    for(size_t obs(0); obs < abcdcuts.size(); obs++){
      TString totcut(abcdcuts[obs]+"&&"+metcuts[ind>=ilowmet]);
      if(obs == 1) totcut += ("&&"+njbcuts[ind]);
      if(obs == 3) totcut += ("&&"+njbcuts_himt[ind]);
      //cout<<base_s+"&&"+totcut<<endl;
      bincuts.push_back(TableRow((base_s+"&&"+totcut).Data(),(base_s+"&&"+totcut).Data()));
    } // Loop over observables going into kappa
  } // Loop over signal bins
    
  PlotMaker pm;
  pm.Push<Table>("preds",  bincuts, all_procs);
  if(single_thread) pm.multithreaded_ = false;
  pm.MakePlots(lumi);

  Table * yield_table = static_cast<Table*>(pm.Figures().back().get());
  vector<GammaParams> mcyield = yield_table->BackgroundYield(lumi);
  vector<GammaParams> datayield = yield_table->DataYield(1);
  vector<GammaParams> otheryield;
  if (do_other) otheryield = yield_table->Yield(other, lumi); 


  vector<float> pow_pred;
  pow_pred.push_back(-1);  //  mt<=140  mj14<=400   R1
  pow_pred.push_back(1);   //  mt<=140  mj14>400    R2
  pow_pred.push_back(1);   //  mt>140   mj14<=400   D3
  vector<float> pow_tot;
  pow_tot.push_back(-1);   //  mt<=140  mj14<=400   R1
  pow_tot.push_back(1);    //  mt<=140  mj14>400    R2
  pow_tot.push_back(1);    //  mt>140   mj14<=400   D3
  pow_tot.push_back(1);    //  mt<=140  mj14<=400   R1
  pow_tot.push_back(-1);   //  mt<=140  mj14>400    R2
  pow_tot.push_back(-1);   //  mt>140   mj14<=400   D3
  pow_tot.push_back(1);    //  mt>140   mj14>400    D4
  vector<float> pow_k;
  pow_k.push_back(1);    //  mt<=140  mj14<=400   R1
  pow_k.push_back(-1);   //  mt<=140  mj14>400    R2
  pow_k.push_back(-1);   //  mt>140   mj14<=400   D3
  pow_k.push_back(1);    //  mt>140   mj14>400    D4

  float mSigma, pSigma, pred, pred_sys, mSigma_sys, pSigma_sys;
  size_t nabcd(abcdcuts.size()), digits(2);
  vector<vector<float> > preds, kappas, fractions;
  for(size_t ind(0); ind<njbcuts.size(); ind++){
    bool lowjets(ind%2==0);
    vector<vector<float> > entries;
    vector<vector<float> > weights;
    for(size_t obs(0); obs < pow_pred.size(); obs++) {
      size_t index(nabcd*ind+obs);
      entries.push_back(vector<float>());
      weights.push_back(vector<float>());
      entries.back().push_back(datayield[index].Yield());
      weights.back().push_back(1.);
    } // Loop over observables for data
    //pred = calcKappa(entries, weights, pow_pred, mSigma, pSigma);

    vector<vector<float> > kn, kw;
    float k(1.), kup(1.), kdown(1.);    
    fractions.push_back(vector<float>());
    for(size_t obs(0); obs < abcdcuts.size(); obs++){
      size_t index(nabcd*ind+obs);
      if(!do_other){
        k *= pow(mcyield[index].Yield(), pow_tot[3+obs]);
        entries.push_back(vector<float>());
        weights.push_back(vector<float>());
        entries.back().push_back(mcyield[index].NEffective());
        weights.back().push_back(mcyield[index].Weight());

        kn.push_back(vector<float>());
        kw.push_back(vector<float>());
        kn.back().push_back(mcyield[index].NEffective());
        kw.back().push_back(mcyield[index].Weight());
      } else {
        k *= pow(mcyield[index].Yield()+otheryield[index].Yield(), pow_tot[3+obs]);
        float f(otheryield[index].Yield()/(mcyield[index].Yield()+otheryield[index].Yield()));
        fractions.back().push_back(f*100);

        kup *= pow(mcyield[index].Yield()+exp(log(1+syst))*otheryield[index].Yield(), pow_tot[3+obs]);
        kdown *= pow(mcyield[index].Yield()+exp(-log(1+syst))*otheryield[index].Yield(), pow_tot[3+obs]);
      }
    } // Loop over observables for MC
    k = calcKappa(kn, kw, pow_k, kdown, kup);
    cout<<"k = "<<k<<" +"<<kup<<" -"<<kdown<<endl;
    if(!do_other){
      // Calculating predictions without systematics
      pred = calcKappa(entries, weights, pow_tot, mSigma, pSigma);
      if(mSigma<0) mSigma = 0;
      
      // Calculating predictions with systematics
      float totsys = (lowjets?0.51:1.07);
      pred_sys = calcKappa(entries, weights, pow_tot, mSigma_sys, pSigma_sys, false, false, totsys);
      if(mSigma_sys < 0) mSigma_sys = 0;
      preds.push_back(vector<float>({pred, pSigma, mSigma, pred_sys, pSigma_sys, mSigma_sys, k, kup, kdown}));
    }
    if(do_other){
      kup = (kup-k)/k*100;
      kdown = (kdown-k)/k*100;
      kappas.push_back(vector<float>({k, kup, kdown}));
      cout<<"k = "<<RoundNumber(k,2)<<", up "<<setw(5)<<RoundNumber(kup,1)<<"%, down "<<setw(5)
	  <<RoundNumber(kdown,1)<<"%, other fractions ";
      for(size_t oth(0); oth<fractions[ind].size(); oth++) cout <<setw(4)<<RoundNumber(fractions[ind][oth],1)<<"% ";
      cout<<" => cuts "<<bincuts[4*ind+3].label_<<endl;
    }

  } // Loop over signal bins

  ///// Printing table
  TString outname = "txt/table_predictions_lumi0p815_"+method+".tex";
  if(full_lumi) {
    if(method.Contains("2015")) outname.ReplaceAll("lumi0p815", "lumi2p3");
    else outname.ReplaceAll("lumi0p815", "lumi2p6");
  }
  if(unblind) outname.ReplaceAll("lumi", "unblind_lumi");

  if(do_other) outname.ReplaceAll("predictions", "other_sys");
  ofstream out(outname);

  size_t Ncol = 7;
  size_t index;
  out << fixed << setprecision(digits);
  out << "\\documentclass{article}\n";
  out << "\\usepackage{amsmath,graphicx,rotating,geometry}\n";
  out << "\\renewcommand{\\arraystretch}{1.3}\n";
  out << "\\thispagestyle{empty}\n";
  out << "\\begin{document}\n\n";
  out << "\\begin{table}\n";
  out << "\\centering\n";
  //out << "\\resizebox{\\textwidth}{!}{\n";
  if(!do_other){
    out << "\\begin{tabular}[tbp!]{ l ";
    for(size_t ind=0; ind<Ncol-1; ind++) out<<"c";
    out<<"}\\hline\\hline\n";
    out<<" \\multicolumn{"<<Ncol<<"}{c}{"<<method_s<<"} \\\\ \\hline"<<endl;
    out<<"$\\cal{L}=$"<<lumi_s<<" & $\\kappa$ & MC  & Pred. & Obs. & MC/Obs & $Z_{\\rm bi}$ \\\\ \\hline\\hline\n";
    if(method.Contains("met150")) out << " \\multicolumn{"<<Ncol<<"}{c}{$150<\\text{MET}\\leq 200$}  \\\\ \\hline\n";
    else if(method.Contains("met500")) out << " \\multicolumn{"<<Ncol<<"}{c}{$\\text{MET}> 500$}  \\\\ \\hline\n";
    else if(method.Contains("met200nb1")) out << " \\multicolumn{"<<Ncol<<"}{c}{$200<\\text{MET}\\leq 350, N_{b}=1$}  \\\\ \\hline\n";
    else if(method.Contains("met350nb1")) out << " \\multicolumn{"<<Ncol<<"}{c}{$350<\\text{MET}\\leq 500, N_{b}=1$}  \\\\ \\hline\n";
    else if(method.Contains("met200nb1")) out << " \\multicolumn{"<<Ncol<<"}{c}{$\\text{MET}> 500, N_{b}=1$}  \\\\ \\hline\n";
    else out << " \\multicolumn{"<<Ncol<<"}{c}{$200<\\text{MET}\\leq 350$}  \\\\ \\hline\n";

    ///////////////////////////////// Doing R1 ////////////////////////////
    index = 0;
    if(method.Contains("nb1")) out << "R1: $N_b=1,\\text{all }N_j$";
    else  out << "R1: all $N_b,N_j$";
    out<<" &  & "<<mcyield[index].Yield() <<" &  & "
       << setprecision(0) <<datayield[index].Yield()<<setprecision(digits)
       <<" & "<<RoundNumber(mcyield[index].Yield(), digits, datayield[index].Yield())<<" &   \\\\"<<endl;
    ///////////////////////////////// Doing R2 ////////////////////////////
    for(size_t ind(0); ind<ilowmet; ind++){
      index = nabcd*ind+1;
      out<<"R2: "<<cuts2tex(njbcuts[ind])<<" &  & "<<mcyield[index].Yield() <<" &   & "
	 << setprecision(0) <<datayield[index].Yield()<<setprecision(digits)
	 <<" & "<<RoundNumber(mcyield[index].Yield(), digits, datayield[index].Yield())
	 <<" &   \\\\"<<endl;
    }
    //////////////////////////////// Doing R3/D3 ////////////////////////////
    index = 2;
    if(method.Contains("nb1")) out << region_s<<"3: $N_b=1,\\text{all }N_j$";
    else out << region_s<<"3: all $N_b,N_j$";
    out << " &  & "<<mcyield[index].Yield() <<" &  & "
	<< setprecision(0) << datayield[index].Yield() << setprecision(digits)
	<<" & "<<RoundNumber(mcyield[index].Yield(), digits, datayield[index].Yield())<<" &   \\\\"<<endl;
    out << "\\hline"<<endl;
    ///////////////////////////////// Doing R4/D4 ////////////////////////////
    for(size_t ind(0); ind<ilowmet; ind++){
      index = nabcd*ind+3;
      out<<region_s<<"4: "<<cuts2tex(njbcuts_himt[ind])<<" & $"<<preds[ind][6] 
	 << "^{+" << preds[ind][7] <<"}_{-" << preds[ind][8] 
	 <<"}$ & "<<mcyield[index].Yield() <<" & $"<<preds[ind][0] << "^{+" << preds[ind][1] 
	 <<"}_{-" << preds[ind][2] <<"}$ & ";
      if(unblind) out << setprecision(0) <<datayield[index].Yield()<<setprecision(digits)
		      <<" & "<<RoundNumber(mcyield[index].Yield(), digits, datayield[index].Yield())
		      <<" & "<< Zbi(datayield[index].Yield(), preds[ind][0], preds[ind][1])<<" \\\\"<<endl;
      else out << blind_s<<" & "<< blind_s<<" \\\\"<<endl;
    }
    if(ilowmet<njbcuts.size()){
      out<<"\\hline\\hline\n ";
      if(method.Contains("met200nb1"))out<<" \\multicolumn{"<<Ncol<<"}{c}{$200<\\text{MET}\\leq 350, N_{b}\\geq2$} \\\\ \\hline\n";
      else if(method.Contains("met350nb1")) out << " \\multicolumn{"<<Ncol<<"}{c}{$350<\\text{MET}\\leq 500, N_{b}\\geq2$}  \\\\ \\hline\n";
      else if(method.Contains("met500nb1")) out << " \\multicolumn{"<<Ncol<<"}{c}{$\\text{MET}> 500, N_{b}\\geq2$} \\\\ \\hline\n";
      else out << " \\multicolumn{"<<Ncol<<"}{c}{$350<\\text{MET}\\leq 500$}  \\\\ \\hline\n";
      ///////////////////////////////// Doing R1 ////////////////////////////
      index = nabcd*ilowmet;
      if(method.Contains("nb1")) out << "R1: $N_b=1,\\text{all }N_j$";
      else  out << "R1: all $N_b,N_j$";
      out << " &  & "<<mcyield[index].Yield() <<" &  & "
	  << setprecision(0) <<datayield[index].Yield()<<setprecision(digits)
	  <<" & "<<RoundNumber(mcyield[index].Yield(), digits, datayield[index].Yield())<<" &  \\\\"<<endl;
      ///////////////////////////////// Doing R2 ////////////////////////////
      for(size_t ind(ilowmet); ind<njbcuts.size(); ind++){
	index = nabcd*ind+1;
	out<<"R2: "<<cuts2tex(njbcuts[ind])<<" &  & "<<mcyield[index].Yield() <<" &  & "
	   << setprecision(0) <<datayield[index].Yield()<<setprecision(digits)
	   <<" & "<<RoundNumber(mcyield[index].Yield(), digits, datayield[index].Yield())<<" &  \\\\"<<endl;
      }
      
      //////////////////////////////// Doing R3/D3 ////////////////////////////
      index = nabcd*ilowmet+2;
      if(method.Contains("nb1")) out << region_s<<"3: $N_b=1,\\text{all }N_j$";
      else out << region_s<<"3: all $N_b,N_j$";
      out << " &  & "<<mcyield[index].Yield() <<" &  & "
	  << setprecision(0) <<datayield[index].Yield()<<setprecision(digits)
	  <<" & "<<RoundNumber(mcyield[index].Yield(), digits, datayield[index].Yield())<<" &  \\\\"<<endl;
      out << "\\hline"<<endl;
      ///////////////////////////////// Doing R4/D4 ////////////////////////////
      for(size_t ind(ilowmet); ind<njbcuts.size(); ind++){
	index = nabcd*ind+3;
	out<<region_s<<"4: "<<cuts2tex(njbcuts_himt[ind])<<" & $"<<preds[ind][6] 
	   << "^{+" << preds[ind][7] <<"}_{-" << preds[ind][8] 
	   <<"}$ & "<<mcyield[index].Yield() <<" & $"<<preds[ind][0] << "^{+" << preds[ind][1] 
	   <<"}_{-" << preds[ind][2] <<"}$ & ";
	if(unblind) out<< setprecision(0) <<datayield[index].Yield() <<setprecision(digits)
		       <<" & "<<RoundNumber(mcyield[index].Yield(), digits, datayield[index].Yield())
		       <<" & "<< Zbi(datayield[index].Yield(), preds[ind][0], preds[ind][1])<<" \\\\"<<endl;
	else out<< blind_s<<" & "<< blind_s<<" \\\\"<<endl;
      }
    }
  } else {
    out << "\n\\begin{tabular}[tbp!]{ l rrr |";
    for(size_t obs(0); obs < abcdcuts.size(); obs++) out<<"r";
    out << "}\\hline\\hline\n";
    out << " &  & non-$t\\bar{t}\\times"<<RoundNumber(1+syst,1)<<"$ & non-$t\\bar{t}\\times"
	<<RoundNumber(exp(-log(1+syst)),1)<<"$ & \\multicolumn{4}{c}{Fraction of non-$t\\bar{t}$ bkg.} \\\\ \n";
    out << "Bin & $\\kappa$ & $\\Delta\\kappa$ [\\%] & $\\Delta\\kappa$ [\\%]";
    for(size_t obs(0); obs < abcdcuts.size(); obs++) out<<" & $f_{\\text{R"<<obs+1<<"}}$ [\\%]";
    out << " \\\\ \\hline\\hline\n";
    for(size_t ind(0); ind<njbcuts.size(); ind++){
      out<<cuts2tex(njbcuts[ind])<<" & "<<RoundNumber(kappas[ind][0],2) <<" & "<<kappas[ind][1] <<" & "<<kappas[ind][2];
      for(size_t obs(0); obs < abcdcuts.size(); obs++) out<<" & "<<fractions[ind][obs];
      out<<" \\\\"<<endl;
      if(ind==ilowmet-1) out<<"\\hline"<<endl;
    }
  } // if(!do_other)

  out<< "\\hline\\hline\n\\end{tabular}"<<endl;
  //out << "}\n";
  out << "\\end{table}\n\n";
  out << "\\end{document}\n";
  out.close();
  TString pdfname(outname); 
  cout<<endl<<" pdflatex "<<outname<<endl;

  time(&endtime); 
  cout<<endl<<"Calculation took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;

}


void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"method", required_argument, 0, 'm'},
      {"unblind", no_argument, 0, 'u'},
      {"full_lumi", no_argument, 0, 'f'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "m:uf", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'm':
      method = optarg;
      break;
    case 'u':
      really_unblind = true;
      break;
    case 'f':
      full_lumi = true;
      break;
    case 0:
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
