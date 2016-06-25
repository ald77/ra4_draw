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
#include "abcd_method.hpp"

using namespace std;

namespace{
  double lumi(0.815);
  bool single_thread = false;
  bool full_lumi(false);
  bool really_unblind = false;
  bool debug = true;
}

template<typename T>
shared_ptr<Process> Proc(const string process_name, Process::Type type,
                         int color, const set<string> &files, const string &cut = "1"){
  return make_shared<Process>(process_name, type, color,
                              unique_ptr<Baby>(new T(files)),
                              cut);
}

void GetOptions(int argc, char *argv[]);
void findPreds(abcd_method &abcds, vector<GammaParams> &mcyield, vector<GammaParams> &datayield,
	       vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds);
void printDebug(abcd_method &abcds, vector<GammaParams> &mcyield, vector<GammaParams> &datayield,
		vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds);

TString Zbi(double Nobs, double Nbkg, double Ebkg){
  double Nsig = Nobs-Nbkg;
  double zbi = RooStats::NumberCountingUtils::BinomialExpZ(Nsig, Nbkg, Ebkg/Nbkg);
  if(Nbkg==0) zbi = RooStats::NumberCountingUtils::BinomialWithTauExpZ(Nsig, Nbkg, 1/Ebkg);
  if(zbi<0) zbi=0;
  TString zbi_s = "$"+RoundNumber(zbi,1)+"\\sigma$";
  if(Nsig<=0) zbi_s = "-";
  //cout<<"Zbi for Nobs "<<Nobs<<", Nbkg "<<Nbkg<<", Ebkg "<<Ebkg<<" is "<<zbi_s<<endl;
  return zbi_s;
}

int main(int argc, char *argv[]){ 
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_standard/");
  string folderdata(bfolder+"/cms2r0/babymaker/babies/2016_06_21/data/skim_standard/");
  // if(method.Contains("met150")){
  //   foldermc = bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_1lht500met150nj5/";
  //   folderdata = bfolder+"/cms2r0/babymaker/babies/2016_06_21/data/skim_1lmet150/";
  // }
  // if(method.Contains("2015")){
  //   foldermc = bfolder+"/cms2r0/babymaker/babies/2016_04_29/mc/merged_1lht500met200/";
  //   folderdata = bfolder+"/cms2r0/babymaker/babies/2016_04_29/data/merged_1lht500met200/";
  // }

  Palette colors("txt/colors.txt", "default");

  // Cuts in baseline speed up the yield finding
  string baseline = "mj14>250 && nleps>=1 && ht>500 && met>150 && pass && njets>=5";
  if(!full_lumi) baseline += " && nonblind";
  else lumi = 2.6;

  auto tt = Proc<Baby_full>("t#bar{t}", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets*Lept*.root", foldermc+"*_TTJets_HT*.root",
    foldermc+"*_WJetsToLNu*.root",foldermc+"*_ST_*.root",
    foldermc+"*_TTW*.root",foldermc+"*_TTZ*.root",
    foldermc+"*DYJetsToLL*.root",foldermc+"*QCD_HT*.root",
    foldermc+"*_ZJet*.root",foldermc+"*_ttHJetTobb*.root",
    foldermc+"*_TTGJets*.root",foldermc+"*_TTTT*.root",
    foldermc+"*_WH_HToBB*.root",foldermc+"*_ZH_HToBB*.root",
    foldermc+"*_WWTo*.root",foldermc+"*_WZ*.root",foldermc+"*_ZZ_*.root"},
    //{foldermc+"*_ZZ_*.root"},
    baseline+" && stitch");

  string trigs = "(trig[4]||trig[8]||trig[13]||trig[33])";
  //if(method.Contains("2015")) trigs = "(trig[4]||trig[8]||trig[28]||trig[14])";
  auto data = Proc<Baby_full>("Data", Process::Type::data, kBlack,
    {folderdata+"*.root"},baseline+" && "+trigs);

  vector<shared_ptr<Process> > all_procs = {data, tt};


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining basic cuts //////////////////////////////////////////
  // baseline defined above

  ////////////////////////////////////////// MET cuts //////////////////////////
  TString c_vlowmet = "met>150 && met<=200";
  TString c_lowmet  = "met>200 && met<=350";
  TString c_midmet  = "met>350 && met<=500";
  TString c_higmet  = "met>500";


  ////////////////////////////////////////// Nb cuts //////////////////////////
  TString c_lownb = "nbm==1";
  TString c_midnb = "nbm==2";
  TString c_hignb = "nbm>=3";

  ////////////////////////////////////////// Njets cuts ////////////////////////
  TString c_lownj = "njets>=6 && njets<=8";
  TString c_hignj = "njets>=9";
  TString c_nj5   = "njets==5";

  // vector<TString> nj_lohi({c_lownj, c_hignj});
  // vector<TString> nbnj_lohi({c_lownb+"&&"+c_lownj, c_lownb+"&&"+c_hignj, c_hignb+"&&"+c_lownj, c_hignb+"&&"+c_hignj});
  // vector<TString> nj_5({c_nj5});

  ////////////////////////////////////////// ABCD cuts ////////////////////////
  vector<TString> abcdcuts_std = {"mt<=140 && mj14<=400  &&  nj_all_1l", 
				  "mt<=140 && mj14>400   &&  nj_1l", 
				  "mt>140  && mj14<=400  &&  nj_all_1l",          
				  "mt>140  && mj14>400   &&  nj_1l"};

  vector<TString> abcdcuts_veto = {"mt<=140 && mj14<=400 && nleps==1 && nveto==0 && nbm>=1  &&  nj_all_1l", 
				   "mt<=140 && mj14>400  && nleps==1 && nveto==0 && nbm>=1  &&  nj_1l", 
				   "mt>140  && mj14<=400 && nleps==1 && nveto==1 && nbm>=1 && nbm<=2  &&  nj_all_1l",          
				   "mt>140  && mj14>400  && nleps==1 && nveto==1 && nbm>=1 && nbm<=2  &&  nj_1l"};

  vector<TString> abcdcuts_2l = {"mt<=140 && mj14<=400 && nleps==1 && nveto==0 && nbm>=1  &&  nj_all_1l", 
				 "mt<=140 && mj14>400  && nleps==1 && nveto==0 && nbm>=1  &&  nj_1l", 
				 "           mj14<=400 && nleps==2             && nbm<=2  &&  nj_all_2l",          
				 "           mj14>400  && nleps==2             && nbm<=2  &&  nj_2l"};

  vector<TString> abcdcuts_2lveto;
  for(size_t ind=0; ind<2; ind++) abcdcuts_2lveto.push_back(abcdcuts_2l[ind]);
  for(size_t ind=2; ind<abcdcuts_2l.size(); ind++){
    abcdcuts_2lveto.push_back("(("+abcdcuts_2l[ind]+") || ("+abcdcuts_veto[ind]+"))");
    abcdcuts_2lveto.back().ReplaceAll("((  ","((");
  }

  ///////////////////////////////////////// END cut definitons ////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining ABCD methods //////////////////////////////////////////
  vector<abcd_method> abcds;
  vector<TString> abcdcuts, metcuts, bincuts;
  PlotMaker pm;

  ///// Running over these methods
  //vector<TString> methods = {"m2lonly", "signal"};
  vector<TString> methods = {"m2lveto", "signal"};

  for(size_t iabcd=0; iabcd<methods.size(); iabcd++) {
    TString method = methods[iabcd];
    TString basecuts = "", caption = "";
    
    if(method.Contains("2l") || method.Contains("veto")) {
      metcuts = vector<TString>{c_lowmet, c_midmet};
      bincuts = vector<TString>{c_lownj, c_hignj}; // 2l nj cuts automatically lowered in abcd_method
      caption = "D3 and D4 have ";
    } else {
      abcdcuts = abcdcuts_std;
      basecuts = "nleps==1 && nveto==0 && nbm>=1";
      metcuts = vector<TString>{c_lowmet, c_midmet, c_higmet};
      bincuts = vector<TString>{c_lownb+" && "+c_lownj, c_lownb+" && "+c_hignj, 
				c_midnb+" && "+c_lownj, c_midnb+" && "+c_hignj, 
				c_hignb+" && "+c_lownj, c_hignb+" && "+c_hignj}; 
      caption = "Signal search regions";
    }

    if(method.Contains("2lonly")) {
      abcdcuts = abcdcuts_2l;
      caption += "two reconstructed leptons";
    }
    if(method.Contains("2lveto")) {
      abcdcuts = abcdcuts_2lveto;
      caption += "either two reconstructed leptons, or one lepton and one track";
    }

    abcds.push_back(abcd_method(method, metcuts, bincuts, abcdcuts, caption, basecuts));
    //abcds.back().printCuts();

    vector<TableRow> table_cuts;
    for(size_t icut=0; icut < abcds.back().allcuts.size(); icut++)
      table_cuts.push_back(TableRow(abcds.back().allcuts[icut].Data(), abcds.back().allcuts[icut].Data()));

    TString tname = "preds"; tname += iabcd;
    pm.Push<Table>(tname.Data(),  table_cuts, all_procs);   
  } // Loop over ABCD methods

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////// Finding all yields ///////////////////////////////////////////////

  if(single_thread) pm.multithreaded_ = false;
  pm.MakePlots(lumi);


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Calculating and printing //////////////////////////////////////////

  for(size_t imethod=0; imethod<abcds.size(); imethod++) {
    cout<<"Number of figures is "<<pm.Figures().size()<<endl;
    Table * yield_table = static_cast<Table*>(pm.Figures()[imethod].get());
    vector<GammaParams> mcyield = yield_table->BackgroundYield(lumi);
    vector<GammaParams> datayield = yield_table->DataYield(1);

    //// Calculating kappa and Total bkg prediction
    vector<vector<vector<float> > > kappas, preds;
    findPreds(abcds[imethod], mcyield, datayield, kappas, preds);


    //// Print MC/Data yields, cuts applied, kappas, preds
    if(debug) printDebug(abcds[imethod], mcyield, datayield, kappas, preds);

  } // Loop over ABCD methods


  time(&endtime); 
  cout<<endl<<"Calculation took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
} // main
////////////////////////////////////////// End of main //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////




//// Calculating kappa and Total bkg prediction
void findPreds(abcd_method &abcd, vector<GammaParams> &mcyield, vector<GammaParams> &datayield,
	       vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds){
  // Powers for kappa:   ({R1, R2, D3, R4})
  vector<float> pow_kappa({ 1, -1, -1,  1});
  // Powers for TotBkg pred:({R1, R2, D3,  R1, R2, D3, D4})
  vector<float> pow_totpred( {-1,  1,  1,   1, -1, -1,  1});

  float val(1.), valup(1.), valdown(1.);    

  for(size_t iplane=0; iplane < abcd.planecuts.size(); iplane++) {
    kappas.push_back(vector<vector<float> >());
    preds.push_back(vector<vector<float> >());
    for(size_t ibin=0; ibin < abcd.bincuts[iplane].size(); ibin++){
      vector<vector<float> > entries;
      vector<vector<float> > weights;
      //// Pushing data yields for predictions
      for(size_t iabcd=0; iabcd < 3; iabcd++){
	size_t index = abcd.indexBin(iplane, ibin, iabcd);
	entries.push_back(vector<float>());
	weights.push_back(vector<float>());
	entries.back().push_back(datayield[index].Yield());
	weights.back().push_back(1.);
      } // Loop over ABCD cuts

      vector<vector<float> > kentries;
      vector<vector<float> > kweights;
      //// Pushing MC yields for predictions and kappas
      for(size_t iabcd=0; iabcd < 4; iabcd++){
	size_t index = abcd.indexBin(iplane, ibin, iabcd);
	// Yields for predictions
	entries.push_back(vector<float>());
	weights.push_back(vector<float>());
	entries.back().push_back(mcyield[index].NEffective());
	weights.back().push_back(mcyield[index].Weight());
	// Yields for kappas
	kentries.push_back(vector<float>());
	kweights.push_back(vector<float>());
	kentries.back().push_back(mcyield[index].NEffective());
	kweights.back().push_back(mcyield[index].Weight());
      } // Loop over ABCD cuts

	// Throwing toys to find predictions and uncertainties
      val = calcKappa(entries, weights, pow_totpred, valdown, valup);
      if(valdown<0) valdown = 0;
      preds[iplane].push_back(vector<float>({val, valup, valdown}));
      // Throwing toys to find kappas and uncertainties
      val = calcKappa(kentries, kweights, pow_kappa, valdown, valup);
      if(valdown<0) valdown = 0;
      kappas[iplane].push_back(vector<float>({val, valup, valdown}));	  
    } // Loop over bin cuts
  } // Loop over plane cuts

} // findPreds

void printDebug(abcd_method &abcds, vector<GammaParams> &mcyield, vector<GammaParams> &datayield,
		vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds){

  cout<<endl<<endl<<"=================== Printing cuts for method "<<abcds.method<<" ==================="<<endl;  
  for(size_t iplane=0; iplane < abcds.planecuts.size(); iplane++) {
    cout<<endl<<" **** Plane "<<abcds.planecuts[iplane]<<" ***"<<endl;
    for(size_t ibin=0; ibin < abcds.bincuts[iplane].size(); ibin++){
      for(size_t iabcd=0; iabcd < abcds.abcdcuts.size(); iabcd++){
	size_t index = abcds.indexBin(iplane, ibin, iabcd);
	cout<<"MC: "<<setw(7)<<RoundNumber(mcyield[index].Yield(),2)
	    <<"  Data: "<<setw(4)<<RoundNumber(datayield[index].Yield(), 0)
	    <<"  - "<< abcds.allcuts[index]<<endl;
      } // Loop over ABCD cuts
      cout<<"Kappa = "<<RoundNumber(kappas[iplane][ibin][0],2)<<"+"<<RoundNumber(kappas[iplane][ibin][1],2)
	  <<"-"<<RoundNumber(kappas[iplane][ibin][2],2)<<", Prediction = "
	  <<RoundNumber(preds[iplane][ibin][0],2)<<"+"<<RoundNumber(preds[iplane][ibin][1],2)
	  <<"-"<<RoundNumber(preds[iplane][ibin][2],2)<<endl;
      cout<<endl;
    } // Loop over bin cuts
  } // Loop over plane cuts

} // printDebug



void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"unblind", no_argument, 0, 'u'},
      {"full_lumi", no_argument, 0, 'f'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "uf", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
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
