///// plot_kappa_mismeas: Makes kappa plots (and ratio plots) for the various mismeasurement scenarios

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>  // setw

#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "TCanvas.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "RooStats/NumberCountingUtils.h"

#include "utilities.hpp"
#include "baby.hpp"
#include "process.hpp"
#include "named_func.hpp"
#include "plot_maker.hpp"
#include "plot_opt.hpp"
#include "palette.hpp"
#include "table.hpp"
#include "abcd_method.hpp"
#include "styles.hpp"

using namespace std;

namespace{
  bool split_bkg = false;
  bool only_dilepton = false;
  bool do_leptons = false;
  bool do_signal = false;
  bool full_lumi = false;
  bool unblind = false;
  bool debug = false;
  TString skim = "standard";
  TString only_method = "";
  float lumi;
  int Nscen = 2;
}

void plotKappa(abcd_method &abcd, vector<vector<vector<vector<float> > > > &allkappas);
void changeMMCut(TString &cut, vector<TString> &mm_cuts, TString method, TString index_s);

TString printTable(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
                   vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds);
void findPreds(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
               vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds);
void printDebug(abcd_method &abcd, vector<vector<GammaParams> > &allyields, TString baseline,
                vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds);
TString Zbi(double Nobs, double Nbkg, double Ebkg);
TString cutsToTex(TString cut);

void GetOptions(int argc, char *argv[]);

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

  string foldermc(bfolder+"/cms2r0/babymaker/babies/mismeasured/2016_06_14/mc/merged_mm_std_nj5mj250/");
  //foldermc = "/net/cms26/cms26r0/babymaker/babies/mismeasured_v2/2016_06_14/mc/merged_mm_std_nj5mj250/";
  string folderdata(bfolder+"/cms2r0/babymaker/babies/2016_06_26/data/merged_standard/");
  folderdata = foldermc;
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
  //string baseline = "mj14>250 && nleps>=1 && ht>500 && met>150 && pass && njets>=5";
  string baseline = "1";
  if(skim.Contains("2015")) {
    lumi = 2.3;
  }else if(!full_lumi) {
    baseline += " && nonblind";
    lumi = 0.815;
  } else lumi = 2.6;
  if(skim.Contains("mj12")) ReplaceAll(baseline, "mj14","mj");

  lumi = 100; // to easily see digits
  auto proc_bkg = Process::MakeShared<Baby_full>("All_bkg", Process::Type::background, colors("tt_1l"),
    // {foldermc+"*_TTJets*Lept*.root", foldermc+"*_TTJets_HT*.root",
    //  foldermc+"*_WJetsToLNu*.root",foldermc+"*_ST_*.root",
    //  foldermc+"*_TTW*.root",foldermc+"*_TTZ*.root",
    //  foldermc+"*DYJetsToLL*.root",foldermc+"*QCD_HT*.root",
    //  foldermc+"*_ZJet*.root",foldermc+"*_ttHJetTobb*.root",
    //  foldermc+"*_TTGJets*.root",foldermc+"*_TTTT*.root",
    //  foldermc+"*_WH_HToBB*.root",foldermc+"*_ZH_HToBB*.root",
    //  foldermc+"*_WWTo*.root",foldermc+"*_WZ*.root",foldermc+"*_ZZ_*.root"},
    {foldermc+"*TTJets_T*.root"},
    baseline+" && stitch");

  auto proc_t1c = Process::MakeShared<Baby_full>("T1tttt(C)", Process::Type::signal, colors("t1tttt"),
    {foldermc+"*mGluino-1200_mLSP-800_*root"},
    baseline+" && stitch");
  auto proc_t1nc = Process::MakeShared<Baby_full>("T1tttt(NC)", Process::Type::signal, colors("t1tttt"),
    {foldermc+"*mGluino-1500_mLSP-100_*root"},
    baseline+" && stitch");

  auto proc_tt1l = Process::MakeShared<Baby_full>("tt 1lep", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets*SingleLept*.root", foldermc+"*_TTJets_HT*.root"},
    baseline+" && stitch && ntruleps==1");
  auto proc_tt2l = Process::MakeShared<Baby_full>("tt 2lep", Process::Type::background, colors("tt_2l"),
    {foldermc+"*_TTJets*DiLept*.root", foldermc+"*_TTJets_HT*.root"},
    baseline+" && stitch && ntruleps==2");

  // auto proc_tt1l = Process::MakeShared<Baby_full>("tt 1lep", Process::Type::background, colors("tt_1l"),
  //   {foldermc+"*_TTJets*SingleLept*.root"},
  //   baseline+" && ntruleps==1");
  // auto proc_tt2l = Process::MakeShared<Baby_full>("tt 2lep", Process::Type::background, colors("tt_2l"),
  //   {foldermc+"*_TTJets*DiLept*.root"},
  //   baseline+" && ntruleps==2");

  auto proc_other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
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
  auto proc_data = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {folderdata+"*.root"},baseline+" && "+trigs);

  //vector<shared_ptr<Process> > all_procs = {proc_tt1l, proc_tt2l, proc_other};
  vector<shared_ptr<Process> > all_procs = {proc_bkg};

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
  TString c_vlownj = "njets>=4 && njets<=5";
  TString c_lownj = "njets>=6 && njets<=8";
  TString c_hignj = "njets>=9";
  TString c_nj5   = "njets==5";

  ////// ABCD cuts
  vector<TString> abcdcuts_std = {"mt<=140 && mj14<=400  && nj_all_1l",
                                  "mt<=140 && mj14>400   && nj_1l",
                                  "mt>140  && mj14<=400  && nj_all_1l",
                                  "mt>140  && mj14>400   && nj_1l"};

  vector<TString> abcdcuts_veto = {"mt<=140 && mj14<=400 && nleps==1 && nbm>=1  &&  nj_all_1l",
                                   "mt<=140 && mj14>400  && nleps==1 && nbm>=1  &&  nj_1l",
                                   "mt>140  && mj14<=400 && nleps==1 && nbm>=1 && nbm<=2  &&  nj_all_1l",
                                   "mt>140  && mj14>400  && nleps==1 && nbm>=1 && nbm<=2  &&  nj_1l"};

  vector<TString> abcdcuts_2l = {"mt<=140 && mj14<=400 && nleps==1 && nbm>=1  &&  nj_all_1l",
                                 "mt<=140 && mj14>400  && nleps==1 && nbm>=1  &&  nj_1l",
                                 "           mj14<=400 && nleps==2 && nbm<=2  &&  nj_all_2l",
                                 "           mj14>400  && nleps==2 && nbm<=2  &&  nj_2l"};

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
  if(only_method!="") methods = vector<TString>({only_method});
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

  vector<TString> methods_ori = methods;
  methods.clear();
  for(auto name2: methods_ori){
    for(int iscen=0; iscen < Nscen; iscen++){
      //if(iscen!=0 && iscen!=3) continue; // To just do the most extreme ones
      TString name = name2;
      name += "_mm"; name += iscen;
      name += "_lep";
      methods.push_back(name);
      name.ReplaceAll("_lep", "_nolep");
      methods.push_back(name);
    } // Loop over mismeasurement scenarios
  } // Loop over methods

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////// Looping over ABCD methods //////////////////////////////////////////
  for(size_t iabcd=0; iabcd<methods.size(); iabcd++) {
    TString method = methods[iabcd];
    TString basecuts = "1", caption = "";

    //////// General assignments to all methods
    if(method.Contains("2l") || method.Contains("veto")) {
      metcuts = vector<TString>{c_lowmet, c_midmet};
      bincuts = vector<TString>{c_lownj, c_hignj}; // 2l nj cuts automatically lowered in abcd_method
      caption = "Dilepton validation regions (with filters). D3 and D4 have ";
    } else {
      if(only_dilepton) continue;
      abcdcuts = abcdcuts_std;
      basecuts = "nleps==1 && nbm>=1";
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
      bincuts = vector<TString>{c_lownb+" && "+c_nj5, c_midnb+" && "+c_nj5, c_hignb+" && "+c_nj5};
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
      caption.ReplaceAll("regions", "region for very low \\met");
    }
    if(method.Contains("m1lmet150")) {
      bincuts = vector<TString>{c_lownb+" && "+c_lownj, c_lownb+" && "+c_hignj,
                                c_midnb+" && "+c_lownj, c_midnb+" && "+c_hignj};
      caption = "Single lepton validation region for very low \\met";
    }
    if(skim.Contains("2015")) {
      caption += ". Data taken in 2015";
    }

    //////// Changes to cuts due to mismeasurements
    if(method.Contains("_mm")){
      vector<TString> mm_cuts({"met", "nleps", "njets", "nbm", "ht", "mt", "mj14"}); // cuts to change
      basecuts += "&& mj14>250 && nleps>=1 && ht>500 && met>150 && pass && njets>=5"; // Adding the baseline here

      //// Finding scenario index
      TString index_s = method;
      index_s.Remove(0, index_s.Index("_mm")+3);
      index_s.Remove(index_s.First('_'), index_s.Length());

      //// Adjusting cuts
      changeMMCut(basecuts, mm_cuts, method, index_s);
      for(auto &cut : metcuts) changeMMCut(cut, mm_cuts, method, index_s);
      for(auto &cut : bincuts) changeMMCut(cut, mm_cuts, method, index_s);
      for(auto &cut : abcdcuts) changeMMCut(cut, mm_cuts, method, index_s);
    } // If method is mismesurement

    //////// Pushing all cuts to then find the yields
    abcds.push_back(abcd_method(method, metcuts, bincuts, abcdcuts, caption, basecuts));
    if(skim.Contains("mj12")) {
      //      abcds.back().setMj12();
      abcds.back().caption += ". Using $M_J^{1.2}$";
    }
    if(method.Contains("_el") || method.Contains("_mu") || method.Contains("_emu")) abcds.back().setLeptons();
    if(method.Contains("_el"))  abcds.back().caption += ". Only electrons";
    if(method.Contains("_mu"))  abcds.back().caption += ". Only muons";
    if(method.Contains("_emu")) abcds.back().caption += ". Only $e/\\mu$ pairs in D3 and D4";
    //if(method.Contains("agg_")) abcds.back().int_nbnj = false; // Only changes caption since there is only 1 bin

    //cout<<endl<<" ======== Method "<<method<<" =============="<<endl;
    //abcds.back().printCuts();

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
  vector<vector<vector<vector<float> > > > allkappas;
  for(size_t imethod=0; imethod<abcds.size(); imethod++) {
    Table * yield_table = static_cast<Table*>(pm.Figures()[imethod].get());
    // allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
    // if split_bkg [2/4] Other, [3/5] tt1l, [4/6] tt2l
    vector<vector<GammaParams> > allyields;
    allyields.push_back(yield_table->DataYield());
    allyields.push_back(yield_table->BackgroundYield(lumi));
    if(do_signal){
      allyields.push_back(yield_table->Yield(proc_t1nc.get(), lumi));
      allyields.push_back(yield_table->Yield(proc_t1c.get(), lumi));
    }
    if(split_bkg){
      allyields.push_back(yield_table->Yield(proc_other.get(), lumi));
      allyields.push_back(yield_table->Yield(proc_tt1l.get(), lumi));
      allyields.push_back(yield_table->Yield(proc_tt2l.get(), lumi));
    }

    //// Calculating kappa and Total bkg prediction
    vector<vector<vector<float> > > kappas, preds;
    findPreds(abcds[imethod], allyields, kappas, preds);

    allkappas.push_back(kappas);
    plotKappa(abcds[imethod], allkappas);

    // //// Makes table MC/Data yields, kappas, preds, Zbi
    // TString fullname = printTable(abcds[imethod], allyields, kappas, preds);
    // tablenames.push_back(fullname);

    //// Print MC/Data yields, cuts applied, kappas, preds
    if(debug) printDebug(abcds[imethod], allyields, baseline, kappas, preds);

  } // Loop over ABCD methods

  // //// Printing names of ouput files
  // cout<<endl<<"===== Tables to be moved to the AN/PAS/paper:"<<endl;
  // for(size_t ind=0; ind<tablenames.size(); ind++){
  //   TString name=tablenames[ind]; name.ReplaceAll("fulltable","table");
  //   cout<<" mv "<<name<<"  ${tables_folder}"<<endl;
  // }
  // cout<<endl<<"===== Tables that can be compiled"<<endl;
  // for(size_t ind=0; ind<tablenames.size(); ind++)
  //   cout<<" pdflatex "<<tablenames[ind]<<"  > /dev/null"<<endl;
  // cout<<endl;

  time(&endtime);
  cout<<endl<<"Finding "<<abcds.size()<<" tables took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
} // main

////////////////////////////////////////// End of main //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//// Changes standard cuts to mismeasure cuts for a given scenario
void changeMMCut(TString &cut, vector<TString> &mm_cuts, TString method, TString index_s){
  for(auto &mm_cut : mm_cuts) {
    cut.ReplaceAll(mm_cut, "mm_"+mm_cut+"["+index_s+"]");
    if(mm_cut=="mj14"){
      if(method.Contains("_nolep")){
        cut.ReplaceAll("mj14", "mj14_nolep");
        cut.ReplaceAll("_2l", "_1l"); // This avoids lowering number of jets
      } else cut.ReplaceAll("mj14", "mj14_lep");
    }
  } // Loop over mismeasured variables
}

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
  if(skim.Contains("mj12")) outname += "_mj12";
  if(unblind) outname += "_unblind";
  else outname += "_blind";
  outname += "_"+abcd.method+".tex";
  ofstream out(outname);

  //// Printing main table preamble
  if(abcd.method.Contains("signal") && Ncol>7) out << "\\resizebox{\\textwidth}{!}{\n";
  out << "\\begin{tabular}[tbp!]{ l ";
  if(do_signal) out << "|cc";
  if(split_bkg) out << "|ccc";
  out << "|cc|ccc";
  if(do_zbi) out<<"c";
  out<<"}\\hline\\hline\n";
  out<<"${\\cal L}="<<lumi_s<<"$ fb$^{-1}$ ";
  if(do_signal) out << " & T1tttt(NC) & T1tttt(C) ";
  if(split_bkg) out << " & Other & $t\\bar{t}(1\\ell)$ & $t\\bar{t}(2\\ell)$ ";
  out<<"& $\\kappa$ & MC bkg. & Pred. & Obs. & Obs./MC "<<(do_zbi?"& $Z_{\\rm Bi}$":"")<<" \\\\ \\hline\\hline\n";

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
        //// Printing kappa
        out<<ump;
        if(iabcd==3) out  << "$"    << RoundNumber(kappas[iplane][ibin][0], digits)
                          << "^{+"  << RoundNumber(kappas[iplane][ibin][1], digits)
                          << "}_{-" << RoundNumber(kappas[iplane][ibin][2], digits) <<"}$ ";
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
void printDebug(abcd_method &abcd, vector<vector<GammaParams> > &allyields, TString baseline,
                vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds){

  cout<<endl<<endl<<"=================== Printing cuts for method "<<abcd.method<<" ==================="<<endl;
  cout<<"-- Baseline cuts: "<<baseline<<endl;
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

//// Plots kappa if allkappas is size 1 or 2, and ratio of kappa/kappa_nominal if larger
void plotKappa(abcd_method &abcd, vector<vector<vector<vector<float> > > > &allkappas){

  bool label_up = false; //// Putting the MET labels at the bottom

  //// Setting plot style
  PlotOpt opts("txt/plot_styles.txt", "Kappa");
  if(label_up) opts.BottomMargin(0.11);
  setPlotStyle(opts);

  struct kmarker{
    TString cut;
    int color;
    int style;
    vector<float> kappa;
  };
  //// k_ordered has all the kappas group in sets of nb cuts (typically, in bins of njets)
  vector<vector<vector<kmarker> > > k_ordered;
  vector<kmarker> ind_bcuts; // nb cuts actually used in the plot
  vector<float> zz; // Zero length vector for the kmarker constructor
  vector<kmarker> bcuts({{"nbm==1",4,20,zz}, {"nbm==2",2,21,zz}, {"nbm>=3",kGreen+3,22,zz}, {"nbm>=2",kMagenta+2,23,zz}});

  int nbins = 0; // Total number of njets bins (used in the base histo)
  int ind= allkappas.size()-1; // Index of the kappa we'll be plotting, just the last in the vector
  vector<vector<vector<float> > > kappas = allkappas[ind];
  for(size_t iplane=0; iplane < kappas.size(); iplane++) {
    k_ordered.push_back(vector<vector<kmarker> >());
    for(size_t ibin=0; ibin < kappas[iplane].size(); ibin++){
      if(allkappas.size()>2){
        kappas[iplane][ibin][0] /= allkappas[ind%2][iplane][ibin][0];
        kappas[iplane][ibin][1] = 0.;
        kappas[iplane][ibin][2] = 0.;
      }
      TString bincut = abcd.bincuts[iplane][ibin];
      bincut.ReplaceAll(" ","");
      bincut.ReplaceAll("mm_","");
      int index;
      do{
        index = bincut.First('[');
        bincut.Remove(index, bincut.First(']')-index+1);
      }while(index>=0);
      bool found=false;
      for(size_t ib=0; ib<bcuts.size(); ib++){
        if(bincut.Contains(bcuts[ib].cut)){
          //// Storing the number of different nb cuts in ind_bcuts
          bool cutfound=false;
          for(size_t indb=0; indb<ind_bcuts.size(); indb++)
            if(ind_bcuts[indb].color == bcuts[ib].color) cutfound = true;
          if(!cutfound) ind_bcuts.push_back(bcuts[ib]);

          //// Cleaning the nb cut from the bincut
          bincut.ReplaceAll(bcuts[ib].cut+"&&","");
          for(size_t ik=0; ik<k_ordered[iplane].size(); ik++){
            //// Adding point to a given njets cut
            if(bincut==k_ordered[iplane][ik][0].cut){
              k_ordered[iplane][ik].push_back({bincut, bcuts[ib].color, bcuts[ib].style, kappas[iplane][ibin]});
              found = true;
              break;
            } // if same njets cut
          } // Loop over existing ordered kappas
          //// If it doesn't correspond to any njet cut yet, create a new bin
          if(!found) {
            k_ordered[iplane].push_back(vector<kmarker>({{bincut, bcuts[ib].color, bcuts[ib].style, kappas[iplane][ibin]}}));
            found = true;
            nbins++;
          }
        } // if bincut.Contains(bcuts[ib].cut)
      } // Loop over nb cuts

      //// If it doesn't correspond to any nb cut, create a new bin with default (color in [0], blue)
      if(!found) {
        k_ordered[iplane].push_back(vector<kmarker>({{bincut, bcuts[0].color, bcuts[0].style, kappas[iplane][ibin]}}));
        nbins++;
        if(ind_bcuts.size()==0) ind_bcuts.push_back(bcuts[0]);
      }
    } // Loop over bin cuts
  } // Loop over plane cuts

  // //// Printing all values of kappa
  // for(size_t iplane=0; iplane < k_ordered.size(); iplane++) {
  //   cout<<endl<<"New plane"<<endl;
  //   for(size_t ibin=0; ibin < k_ordered[iplane].size(); ibin++){
  //     cout<<"Points for "<<k_ordered[iplane][ibin][0].cut<<endl;
  //     for(size_t ib=0; ib<k_ordered[iplane][ibin].size(); ib++){
  //    cout<<"kappa = "<<k_ordered[iplane][ibin][ib].kappa[0]<<", color "<<k_ordered[iplane][ibin][ib].color<<endl;
  //     }
  //   }
  // }

  //// Plotting kappas
  TCanvas can("can","",1100,700);
  TLine line; line.SetLineWidth(2); line.SetLineStyle(2);
  TLatex label; label.SetTextSize(0.05); label.SetTextFont(132); label.SetTextAlign(23);

  float minx = 0.5, maxx = nbins+0.5, miny = 0, maxy = 2.2;
  if(label_up) maxy = 2.4;
  TH1D histo("histo", abcd.method, nbins, minx, maxx);
  histo.SetMinimum(miny);
  histo.SetMaximum(maxy);
  histo.GetYaxis()->CenterTitle(true);
  histo.GetXaxis()->SetLabelOffset(0.008);
  if(allkappas.size()<=2) histo.SetYTitle("#kappa");
  else histo.SetYTitle("#kappa^{scenario}/#kappa^{nominal}");
  histo.Draw();

  //// Filling vx, vy vectors with kappa coordinates. Each nb cut is stored in a TGraphAsymmetricErrors
  int bin = 0;
  vector<vector<double> > vx(ind_bcuts.size()), vexh(ind_bcuts.size()), vexl(ind_bcuts.size());
  vector<vector<double> > vy(ind_bcuts.size()), veyh(ind_bcuts.size()), veyl(ind_bcuts.size());
  for(size_t iplane=0; iplane < k_ordered.size(); iplane++) {
    for(size_t ibin=0; ibin < k_ordered[iplane].size(); ibin++){
      bin++;
      histo.GetXaxis()->SetBinLabel(bin, cutsToLabel(k_ordered[iplane][ibin][0].cut));
      // xval is the x position of the first marker in the group
      double xval = bin, nbs = k_ordered[iplane][ibin].size(), minxb = 0.15, binw = 0;
      // If there is more than one point in the group, it starts minxb to the left of the center of the bin
      // binw is the distance between points in the njets group
      if(nbs>1) {
        xval -= minxb;
        binw = 2*minxb/(nbs-1);
      }
      for(size_t ib=0; ib<k_ordered[iplane][ibin].size(); ib++){
        //// Finding which TGraph this point goes into by comparing the color of the TGraph and the point
        for(size_t indb=0; indb<ind_bcuts.size(); indb++){
          if(ind_bcuts[indb].color == k_ordered[iplane][ibin][ib].color){
            vx[indb].push_back(xval);
            xval += binw;
            vexl[indb].push_back(0);
            vexh[indb].push_back(0);
            vy[indb].push_back(k_ordered[iplane][ibin][ib].kappa[0]);
            veyl[indb].push_back(k_ordered[iplane][ibin][ib].kappa[1]);
            veyh[indb].push_back(k_ordered[iplane][ibin][ib].kappa[2]);
          }
        } // Loop over nb cuts in ordered TGraphs
      } // Loop over nb cuts in kappa plot
    } // Loop over bin cuts

    // Drawing line separating MET planes
    line.SetLineStyle(2); line.SetLineWidth(2);
    if (iplane<k_ordered.size()-1) line.DrawLine(bin+0.5, miny, bin+0.5, maxy);
    // Drawing MET labels
    if(label_up) label.DrawLatex((2*bin-k_ordered[iplane].size()+1.)/2., maxy-0.1, cutsToLabel(abcd.planecuts[iplane]));
    else label.DrawLatex((2*bin-k_ordered[iplane].size()+1.)/2., -0.25, cutsToLabel(abcd.planecuts[iplane]));
  } // Loop over plane cuts

  //// Drawing legend and TGraphs
  double legX(opts.LeftMargin()+0.026), legY(1-opts.TopMargin()-0.04), legSingle = 0.05;
  if(label_up) legY = 0.8;
  double legW = 0.22, legH = legSingle*(ind_bcuts.size()+1)/2;
  if(ind_bcuts.size()>3) legH = legSingle*((ind_bcuts.size()+1)/2);
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(opts.LegendEntryHeight()); leg.SetFillColor(0);
  leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(42);
  leg.SetNColumns(2);
  TGraphAsymmErrors graph[20]; // There's problems with vectors of TGraphs, so using an array
  for(size_t indb=0; indb<ind_bcuts.size(); indb++){
    graph[indb] = TGraphAsymmErrors(vx[indb].size(), &(vx[indb][0]), &(vy[indb][0]),
                                    &(vexl[indb][0]), &(vexh[indb][0]), &(veyl[indb][0]), &(veyh[indb][0]));
    graph[indb].SetMarkerStyle(ind_bcuts[indb].style); graph[indb].SetMarkerSize(1.65);
    graph[indb].SetMarkerColor(ind_bcuts[indb].color);
    graph[indb].SetLineColor(ind_bcuts[indb].color); graph[indb].SetLineWidth(2);
    graph[indb].Draw("p0 same");
    leg.AddEntry(&graph[indb], cutsToLabel(ind_bcuts[indb].cut), "p");
  } // Loop over TGraphs
  if(ind_bcuts.size()>1) leg.Draw();

  line.SetLineStyle(3); line.SetLineWidth(1);
  line.DrawLine(minx, 1, maxx, 1);

  TString fname="plots/kappa_"+abcd.method+".pdf";
  can.SaveAs(fname);
  cout<<endl<<" open "<<fname<<endl;

}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"method", required_argument, 0, 'm'},  // Method to run on (if you just want one)
      {"skim", required_argument, 0, 's'},    // Which skim to use: standard, met150, 2015 data
      {"split_bkg", no_argument, 0, 'b'},     // Prints Other, tt1l, tt2l contributions
      {"no_signal", no_argument, 0, 'n'},     // Does not print signal columns
      {"do_leptons", no_argument, 0, 'l'},    // Does tables for e/mu/emu as well
      {"unblind", no_argument, 0, 'u'},       // Unblinds R4/D4
      {"full_lumi", no_argument, 0, 'f'},     // Uses all data (does not apply nonblind)
      {"debug", no_argument, 0, 'd'},         // Debug: prints yields and cuts used
      {"only_dilepton", no_argument, 0, '2'}, // Makes tables only for dilepton tests
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "m:s:ufdbnl2", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'm':
      only_method = optarg;
      break;
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
