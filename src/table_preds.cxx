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
  bool split_bkg = false;
  bool only_dilepton = false;
  bool do_leptons = false;
  bool do_signal = true;
  bool full_lumi = false;
  bool unblind = false;
  bool debug = false;
  TString skim = "standard";
  float lumi;
}


TString printTable(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
		   vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds);
void findPreds(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
               vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds);
void printDebug(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
                vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds);
TString Zbi(double Nobs, double Nbkg, double Ebkg);
TString cutsToTex(TString cut);

void GetOptions(int argc, char *argv[]);

template<typename T>
shared_ptr<Process> Proc(const string process_name, Process::Type type,
                         int color, const set<string> &files, const string &cut = "1"){
  return make_shared<Process>(process_name, type, color,
                              unique_ptr<Baby>(new T(files)),
                              cut);
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
  if(skim.Contains("met150")){
    foldermc = bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_met150/";
    folderdata = bfolder+"/cms2r0/babymaker/babies/2016_06_21/data/skim_1lmet150/";
  }
  if(skim.Contains("2015")){
    foldermc = bfolder+"/cms2r0/babymaker/babies/2016_04_29/mc/merged_1lht500met200/";
    folderdata = bfolder+"/cms2r0/babymaker/babies/2016_04_29/data/merged_1lht500met200/";
  }

  Palette colors("txt/colors.txt", "default");

  // Cuts in baseline speed up the yield finding
  string baseline = "mj14>250 && nleps>=1 && ht>500 && met>150 && pass && njets>=5";
  if(!full_lumi) {
    baseline += " && nonblind";
    lumi = 0.815;
  } else lumi = 2.6;
  if(skim.Contains("2015")) lumi = 2.3;

  auto proc_bkg = Proc<Baby_full>("All_bkg", Process::Type::background, colors("tt_1l"),
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

  auto proc_t1c = Proc<Baby_full>("T1tttt(C)", Process::Type::signal, colors("t1tttt"),
    {foldermc+"*mGluino-1200_mLSP-800_*root"},
    baseline+" && stitch");
  auto proc_t1nc = Proc<Baby_full>("T1tttt(NC)", Process::Type::signal, colors("t1tttt"),
    {foldermc+"*mGluino-1500_mLSP-100_*root"},
    baseline+" && stitch");
  auto proc_tt1l = Proc<Baby_full>("tt 1lep", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets*Lept*.root", foldermc+"*_TTJets_HT*.root"},
    baseline+" && stitch && ntruleps==1");
  auto proc_tt2l = Proc<Baby_full>("tt 2lep", Process::Type::background, colors("tt_2l"),
    {foldermc+"*_TTJets*Lept*.root", foldermc+"*_TTJets_HT*.root"},
    baseline+" && stitch && ntruleps==2");
  auto proc_other = Proc<Baby_full>("Other", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_WJetsToLNu*.root",foldermc+"*_ST_*.root",
	foldermc+"*_TTW*.root",foldermc+"*_TTZ*.root",
	foldermc+"*DYJetsToLL*.root",foldermc+"*QCD_HT*.root",
	foldermc+"*_ZJet*.root",foldermc+"*_ttHJetTobb*.root",
	foldermc+"*_TTGJets*.root",foldermc+"*_TTTT*.root",
	foldermc+"*_WH_HToBB*.root",foldermc+"*_ZH_HToBB*.root",
	foldermc+"*_WWTo*.root",foldermc+"*_WZ*.root",foldermc+"*_ZZ_*.root"},
    baseline+" && stitch");

  string trigs = "(trig[4]||trig[8]||trig[13]||trig[33])";
  if(skim.Contains("2015")) trigs = "(trig[4]||trig[8]||trig[28]||trig[14])";
  auto proc_data = Proc<Baby_full>("Data", Process::Type::data, kBlack,
    {folderdata+"*.root"},baseline+" && "+trigs);

  vector<shared_ptr<Process> > all_procs = {proc_data, proc_bkg};
  if (do_signal){
    all_procs.push_back(proc_t1nc);
    all_procs.push_back(proc_t1c);
  }
  if (split_bkg){
    all_procs.push_back(proc_tt1l);
    all_procs.push_back(proc_tt2l);
    all_procs.push_back(proc_other);
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining basic cuts //////////////////////////////////////////
  // baseline defined above

  ////// MET cuts
  TString c_vlowmet = "met>150 && met<=200";
  TString c_lowmet  = "met>200 && met<=350";
  TString c_midmet  = "met>350 && met<=500";
  TString c_higmet  = "met>500";

  ////// Nb cuts
  TString c_lownb = "nbm==1";
  TString c_midnb = "nbm==2";
  TString c_hignb = "nbm>=3";

  ////// Njets cuts
  TString c_lownj = "njets>=6 && njets<=8";
  TString c_hignj = "njets>=9";
  TString c_nj5   = "njets==5";

  ////// ABCD cuts
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


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining ABCD methods //////////////////////////////////////////
  vector<abcd_method> abcds;
  vector<TString> abcdcuts, metcuts, bincuts;
  PlotMaker pm;

  ///// Running over these methods
  vector<TString> methods_std = {"m2lveto", "m2lonly", "mvetoonly", "signal", "m5j", 
				 "agg_himet", "agg_mixed", "agg_himult", "agg_1b"};
  vector<TString> methods_met150 = {"m2lvetomet150", "m2lonlymet150", "mvetoonlymet150", "m1lmet150"};
  vector<TString> methods = methods_std;
  if(skim.Contains("met150")) methods = methods_met150;
  if(do_leptons){
    for(auto name: methods){
      name += "_el";
      methods.push_back(name);
      name.ReplaceAll("_el", "_mu");
      methods.push_back(name);
      if(name.Contains("2l")){ 
	name.ReplaceAll("_mu", "_emu");
	methods.push_back(name);
      }
    }
  }

  for(size_t iabcd=0; iabcd<methods.size(); iabcd++) {
    TString method = methods[iabcd];
    TString basecuts = "", caption = "";
    
    //////// General assignments to all methods
    if(method.Contains("2l") || method.Contains("veto")) {
      metcuts = vector<TString>{c_lowmet, c_midmet};
      bincuts = vector<TString>{c_lownj, c_hignj}; // 2l nj cuts automatically lowered in abcd_method
      caption = "Dilepton validation regions. D3 and D4 have ";
    } else {
      if(only_dilepton) continue;
      abcdcuts = abcdcuts_std;
      basecuts = "nleps==1 && nveto==0 && nbm>=1";
    }

    //////// Dilepton methods
    if(method.Contains("2lonly")) {
      abcdcuts = abcdcuts_2l;
      caption += "two reconstructed leptons";
    }
    if(method.Contains("2lveto")) {
      abcdcuts = abcdcuts_2lveto;
      caption += "either two reconstructed leptons, or one lepton and one track";
    }
    if(method.Contains("vetoonly")) {
      abcdcuts = abcdcuts_veto;
      caption += "one lepton and one track";
    }

    //////// Single lepton methods, all use the standard ABCD plane and nleps==1&&nveto==0&&nbm>=1
    if(method.Contains("signal")) {
      metcuts = vector<TString>{c_lowmet, c_midmet, c_higmet};
      bincuts = vector<TString>{c_lownb+" && "+c_lownj, c_lownb+" && "+c_hignj, 
                                c_midnb+" && "+c_lownj, c_midnb+" && "+c_hignj, 
                                c_hignb+" && "+c_lownj, c_hignb+" && "+c_hignj}; 
      caption = "Signal search regions";
    }
    if(method.Contains("m5j")) {
      metcuts = vector<TString>{c_lowmet, c_midmet};
      bincuts = vector<TString>{c_lownb+" && "+c_nj5, c_lownb+" && "+c_nj5, c_hignb+" && "+c_nj5}; 
      caption = "Validation regions with $1\\ell, \\njets=5$";
    }
    ////// Aggregate regions (single lepton). The nbm, njets integration in R1/R3 is done in abcd_method
    if(method.Contains("agg_himet")) {
      metcuts = vector<TString>{"met>500"};
      bincuts = vector<TString>{"nbm>=3&&njets>=6"}; 
      caption = "High-\\met aggregate region with $1\\ell$, $\\met>500\\text{ GeV}$, $\\njets\\geq6$, $\\nb\\geq3$";
    }
    if(method.Contains("agg_mixed")) {
      metcuts = vector<TString>{"met>350"};
      bincuts = vector<TString>{"nbm>=2&&njets>=9"}; 
      caption = "Mixed aggregate region with $1\\ell$, $\\met>350\\text{ GeV}$, $\\njets\\geq9$, $\\nb\\geq2$";
    }
    if(method.Contains("agg_himult")) {
      metcuts = vector<TString>{"met>200"};
      bincuts = vector<TString>{"nbm>=3&&njets>=9"}; 
      caption = "High-multiplicity aggregate region with $1\\ell$, $\\met>200\\text{ GeV}$, $\\njets\\geq9$, $\\nb\\geq3$";
    }
    if(method.Contains("agg_1b")) {
      metcuts = vector<TString>{"met>500"};
      bincuts = vector<TString>{"nbm>=1&&njets>=6"}; 
      caption = "Single b-tag aggregate region with $1\\ell$, $\\met>500\\text{ GeV}$, $\\njets\\geq6$, $\\nb\\geq1$";
    }

    //////// MET150 methods
    if(method.Contains("met150")) {
      metcuts = vector<TString>{c_vlowmet};
    }
    if(method.Contains("m1lmet150")) {
      bincuts = vector<TString>{c_lownb+" && "+c_lownj, c_lownb+" && "+c_hignj, 
                                c_midnb+" && "+c_lownj, c_midnb+" && "+c_hignj}; 
      caption = "Single lepton very low \\met validation region";
    }


    //////// Pushing all cuts to then find the yields
    abcds.push_back(abcd_method(method, metcuts, bincuts, abcdcuts, caption, basecuts));
    if(method.Contains("_el") || method.Contains("_mu") || method.Contains("_emu")) abcds.back().setLeptons();
    if(method.Contains("_el"))  abcds.back().caption += ". Only electrons";
    if(method.Contains("_mu"))  abcds.back().caption += ". Only muons";
    if(method.Contains("_emu")) abcds.back().caption += ". Only $e/\\mu$ pairs in D3 and D4";
    //if(method.Contains("agg_")) abcds.back().int_nbnj = false; // Only changes caption since there is only 1 bin

    vector<TableRow> table_cuts;
    for(size_t icut=0; icut < abcds.back().allcuts.size(); icut++)
      table_cuts.push_back(TableRow(abcds.back().allcuts[icut].Data(), abcds.back().allcuts[icut].Data()));

    TString tname = "preds"; tname += iabcd;
    pm.Push<Table>(tname.Data(),  table_cuts, all_procs, false);   
  } // Loop over ABCD methods

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////// Finding all yields ///////////////////////////////////////////////

  bool single_thread = false;
  if(single_thread) pm.multithreaded_ = false;
  pm.MakePlots(lumi);


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////// Calculating preds/kappas and printing table //////////////////////////////////////
  vector<TString> tablenames;
  for(size_t imethod=0; imethod<abcds.size(); imethod++) {
    Table * yield_table = static_cast<Table*>(pm.Figures()[imethod].get());
    // allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
    // if split_bkg [2/4] Other, [3/5] tt1l, [4/6] tt2l
    vector<vector<GammaParams> > allyields;
    allyields.push_back(yield_table->DataYield(1));
    allyields.push_back(yield_table->Yield(proc_bkg, lumi));
    if(do_signal){
      allyields.push_back(yield_table->Yield(proc_t1nc, lumi));
      allyields.push_back(yield_table->Yield(proc_t1c, lumi));
    }
    if(split_bkg){
      allyields.push_back(yield_table->Yield(proc_other, lumi));
      allyields.push_back(yield_table->Yield(proc_tt1l, lumi));
      allyields.push_back(yield_table->Yield(proc_tt2l, lumi));
    }

    //// Calculating kappa and Total bkg prediction
    vector<vector<vector<float> > > kappas, preds;
    findPreds(abcds[imethod], allyields, kappas, preds);

    //// Makes table MC/Data yields, kappas, preds, Zbi
    TString fullname = printTable(abcds[imethod], allyields, kappas, preds);
    tablenames.push_back(fullname);


    //// Print MC/Data yields, cuts applied, kappas, preds
    if(debug) printDebug(abcds[imethod], allyields, kappas, preds);

  } // Loop over ABCD methods

  //// Printing names of ouput files
  cout<<endl<<"===== Tables to be moved to the AN/PAS/paper:"<<endl;
  for(size_t ind=0; ind<tablenames.size(); ind++){
    TString name=tablenames[ind]; name.ReplaceAll("fulltable","table");
    cout<<" mv "<<name<<"  ${figs_folder}"<<endl;
  }
  cout<<endl<<"===== Tables that can be compiled"<<endl;
  for(size_t ind=0; ind<tablenames.size(); ind++)
    cout<<" pdflatex "<<tablenames[ind]<<"  > /dev/null"<<endl;
  cout<<endl;

  time(&endtime); 
  cout<<endl<<"Finding "<<abcds.size()<<" tables took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
} // main
  
////////////////////////////////////////// End of main //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


//// Prints table with results
// allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
// if split_bkg: [2/4] Other, [3/5] tt1l, [4/6] tt2l
TString printTable(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
		   vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds){
  //// Table general parameters
  int digits = 2;
  TString ump = " & ";
  bool do_zbi = true;
  if(!unblind) do_zbi = false;
  size_t Ncol = 6;
  if(do_zbi) Ncol++;
  if(do_signal) Ncol += 2;
  if(split_bkg) Ncol += 3;
  TString blind_s = "$\\spadesuit$";
  
  //// Setting output file name
  int digits_lumi = 1;
  if(!full_lumi && !skim.Contains("2015")) digits_lumi = 3;
  TString lumi_s = RoundNumber(lumi, digits_lumi); 
  TString outname = "tables/table_pred_lumi"+lumi_s; outname.ReplaceAll(".","p");
  if(skim.Contains("2015")) outname += "_2015";
  if(unblind) outname += "_unblind";
  else outname += "_blind";
  outname += "_"+abcd.method+".tex";
  ofstream out(outname);

  //// Printing main table preamble
  if(abcd.method.Contains("signal") && Ncol>7) out << "\\resizebox{\\textwidth}{!}{\n";
  out << "\\begin{tabular}[tbp!]{ l ";
  for(size_t ind=0; ind<Ncol-1; ind++) out<<"c";
  out<<"}\\hline\\hline\n";
  out<<"${\\cal L}="<<lumi_s<<"$ fb$^{-1}$ & $\\kappa$ ";
  if(do_signal) out << " & T1tttt(NC) & T1tttt(C) ";
  if(split_bkg) out << " & Other & $t\\bar{t}(1\\ell)$ & $t\\bar{t}(2\\ell)$ ";
  out<<"& MC bkg. & Pred. & Obs. & Obs./MC "<<(do_zbi?"& $Z_{\\rm Bi}$":"")<<" \\\\ \\hline\\hline\n";

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////// Printing results////////////////////////////////////////////////
  for(size_t iplane=0; iplane < abcd.planecuts.size(); iplane++) {
    out<<endl<< "\\multicolumn{"<<Ncol<<"}{c}{"<<cutsToTex(abcd.planecuts[iplane])<<"}  \\\\ \\hline\n";
    for(size_t iabcd=0; iabcd < abcd.abcdcuts.size(); iabcd++){
      for(size_t ibin=0; ibin < abcd.bincuts[iplane].size(); ibin++){
        size_t index = abcd.indexBin(iplane, ibin, iabcd);
        if(abcd.int_nbnj && iabcd%2==0 && ibin>0) continue;
	if(iabcd==3 && ibin==0) out << "\\hline" << endl; 
	//// Printing bin name
        out << (iabcd<=1?"R":abcd.rd_letter) << iabcd+1 << ": ";
	if(iabcd%2==0 && abcd.int_nbnj)
	  out << "All "<<(abcd.bincuts[iplane][ibin].Contains("nbm")?"\\nb, ":"")<<"\\njets" ;
	else {
	  if(abcd.method.Contains("2lonly") && iabcd>=2) out<<cutsToTex(abcd.lowerNjets(abcd.bincuts[iplane][ibin]));
	  else if(abcd.method.Contains("2lveto") && iabcd>=2){
	    if(abcd.bincuts[iplane][ibin].Contains("6")) out<<"Low \\njets";
	    else out<<"High \\njets";
	  } else out<<cutsToTex(abcd.bincuts[iplane][ibin]);
	}
	out<<ump;
	//// Printing kappa
        if(iabcd==3) out  << "$"    << RoundNumber(kappas[iplane][ibin][0], digits) 
			 << "^{+"  << RoundNumber(kappas[iplane][ibin][1], digits) 
                         << "}_{-" << RoundNumber(kappas[iplane][ibin][2], digits) <<"}$ ";
	//// Printing signal yields
	if(do_signal) 
	  out<<ump<<RoundNumber(allyields[2][index].Yield(), digits)<< ump <<RoundNumber(allyields[3][index].Yield(), digits);
	//// Printing Other, tt1l, tt2l
	if(split_bkg){
	  size_t offset = (do_signal?2:0);
	  out << ump <<RoundNumber(allyields[offset+2][index].Yield(), digits)
	      << ump <<RoundNumber(allyields[offset+3][index].Yield(), digits)
	      << ump <<RoundNumber(allyields[offset+4][index].Yield(), digits);
	}
	//// Printing MC Bkg yields
	out << ump << RoundNumber(allyields[1][index].Yield(), digits)<< ump;
	//// Printing background predictions
	if(iabcd==3) out << "$"    << RoundNumber(preds[iplane][ibin][0], digits) 
			 << "^{+"  << RoundNumber(preds[iplane][ibin][1], digits) 
                         << "}_{-" << RoundNumber(preds[iplane][ibin][2], digits) <<"}$ ";
	//// Printing observed events in data and Obs/MC ratio
	if(!unblind && iabcd==3) out << ump << blind_s<< ump << blind_s;
	else {
	  out << ump << RoundNumber(allyields[0][index].Yield(), 0);
	  TString ratio_s = "-";
	  double Nobs = allyields[0][index].Yield(), Nmc = allyields[1][index].Yield();
	  double Eobs = sqrt(Nobs), Emc = allyields[1][index].Uncertainty();
	  if(Nobs==0) Eobs=1;
	  if(Emc>0){
	    double ratio = Nobs/Nmc;
	    double Eratio = sqrt(pow(Eobs/Nmc,2) + pow(Nobs*Emc/Nmc/Nmc,2));
	    ratio_s = "$"+RoundNumber(ratio, 2)+"\\pm"+RoundNumber(Eratio,2)+"$";
	  }
	  out << ump << ratio_s;
	}
	//// Printing Zbi significance
	if(do_zbi && iabcd==3) out << ump << Zbi(allyields[0][index].Yield(), preds[iplane][ibin][0], preds[iplane][ibin][1]);
        out << "\\\\ \n";
      } // Loop over bin cuts
    } // Loop over ABCD cuts
    out << "\\hline\\hline\n";
  } // Loop over plane cuts
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //// Printing footer and closing file
  out<< "\\end{tabular}"<<endl;
  if(abcd.method.Contains("signal") && Ncol>7) out << "}\n"; // For resizebox
  out.close();

  //// Copying header and table to the compilable file
  TString fullname = outname; fullname.ReplaceAll("table_","fulltable_");
  ofstream full(fullname);
  ifstream header("txt/header.tex");
  full<<header.rdbuf();
  header.close();
  if(!abcd.method.Contains("signal")) full << "\\usepackage[landscape]{geometry}\n\n";
  full << "\\begin{document}\n\n";
  full << "\\begin{table}\n\\centering\n";
  full << "\\caption{" << abcd.caption <<".}\\vspace{0.1in}\n\\label{tab:"<<abcd.method<<"}\n";
  ifstream outtab(outname);
  full << outtab.rdbuf();
  outtab.close();
  full << "\\end{table}\n\n";
  full << "\\end{document}\n";
  full.close();

  return fullname;
} // printTable

//// Converting ROOT cuts to TeX
TString cutsToTex(TString cut){
  cut.ReplaceAll(" ","");
  cut.ReplaceAll("met>150&&met<=200", "150<met<=200");
  cut.ReplaceAll("met>200&&met<=350", "200<met<=350");
  cut.ReplaceAll("met>350&&met<=500", "350<met<=500");
  cut.ReplaceAll("njets>=5&&njets<=7", "5<=njets<=7");
  cut.ReplaceAll("njets>=6&&njets<=8", "6<=njets<=8");
  cut.ReplaceAll("nbm>=1&&nbm<=2", "1<=nbm<=2");

  cut.ReplaceAll("met","\\met");
  cut.ReplaceAll("njets","\\njets");
  cut.ReplaceAll("nbm","\\nb");
  cut.ReplaceAll("==","=");
  cut.ReplaceAll(">=","\\geq");
  cut.ReplaceAll("<=","\\leq");
  cut.ReplaceAll("&&",", ");

  cut = "$"+cut+"$";
  return cut;
}


//// Estimating significance
TString Zbi(double Nobs, double Nbkg, double Ebkg){
  double Nsig = Nobs-Nbkg;
  double zbi = RooStats::NumberCountingUtils::BinomialExpZ(Nsig, Nbkg, Ebkg/Nbkg);
  if(Nbkg==0) zbi = RooStats::NumberCountingUtils::BinomialWithTauExpZ(Nsig, Nbkg, 1/Ebkg);
  if(zbi<0) zbi=0;
  TString zbi_s = RoundNumber(zbi,1);
  if(zbi_s!="-") zbi_s = "$"+zbi_s+"\\sigma$";
  if(Nsig<=0 || Ebkg<=0) zbi_s = "-";
  //cout<<"Zbi for Nobs "<<Nobs<<", Nbkg "<<Nbkg<<", Ebkg "<<Ebkg<<" is "<<zbi_s<<endl;
  return zbi_s;
}


//// Calculating kappa and Total bkg prediction
// allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
void findPreds(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
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
        entries.back().push_back(allyields[0][index].Yield());
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
        entries.back().push_back(allyields[1][index].NEffective());
        weights.back().push_back(allyields[1][index].Weight());
        // Yields for kappas
        kentries.push_back(vector<float>());
        kweights.push_back(vector<float>());
        kentries.back().push_back(allyields[1][index].NEffective());
        kweights.back().push_back(allyields[1][index].Weight());
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


// allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
void printDebug(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
                vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds){

  cout<<endl<<endl<<"=================== Printing cuts for method "<<abcd.method<<" ==================="<<endl;  
  for(size_t iplane=0; iplane < abcd.planecuts.size(); iplane++) {
    cout<<endl<<" **** Plane "<<abcd.planecuts[iplane]<<" ***"<<endl;
    for(size_t ibin=0; ibin < abcd.bincuts[iplane].size(); ibin++){
      for(size_t iabcd=0; iabcd < abcd.abcdcuts.size(); iabcd++){
        size_t index = abcd.indexBin(iplane, ibin, iabcd);
        cout<<"MC: "<<setw(7)<<RoundNumber(allyields[1][index].Yield(),2)
            <<"  Data: "<<setw(4)<<RoundNumber(allyields[0][index].Yield(), 0)
            <<"  - "<< abcd.allcuts[index]<<endl;
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
      {"skim", required_argument, 0, 's'},
      {"split_bkg", no_argument, 0, 'b'},
      {"no_signal", no_argument, 0, 'n'},
      {"do_leptons", no_argument, 0, 'l'},
      {"unblind", no_argument, 0, 'u'},
      {"full_lumi", no_argument, 0, 'f'},
      {"debug", no_argument, 0, 'd'},
      {"only_dilepton", no_argument, 0, '2'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s:ufdbnl2", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      skim = optarg;
      break;
    case 'b':
      split_bkg = true;
      break;
    case '2':
      only_dilepton = true;
      break;
    case 'l':
      do_leptons = true;
      break;
    case 'n':
      do_signal = false;
      break;
    case 'u':
      unblind = true;
      break;
    case 'd':
      debug = true;
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
