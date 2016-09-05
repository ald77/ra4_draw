///// table_preds: Generates tables with MC/data yields, bkg predictions
/////              ABCDs are defined in abcd_method, with planecuts (typicaly MET bins),
////               bincuts (typically Nb/Njets bins), and abcdcuts (the cuts for the 4 regions)

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>  // setw
#include <chrono>

#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "RooStats/NumberCountingUtils.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/abcd_method.hpp"
#include "core/styles.hpp"
#include "core/plot_opt.hpp"

using namespace std;

namespace{
  bool only_mc = true;
  bool only_kappa = false;
  bool split_bkg = true;
  bool only_dilepton = false;
  bool do_leptons = false;
  bool do_signal = true;
  bool unblind = true;
  bool debug = false;
  bool do_ht = false;
  TString skim = "standard";
  TString json = "2p6";
  TString only_method = "";
  TString mc_lumi = "40";
  float lumi;
}

TString printTable(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
                   vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds,
		   vector<shared_ptr<Process> > &proc_sigs);
void plotKappa(abcd_method &abcd, vector<vector<vector<float> > >  &kappas);
void findPreds(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
               vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds);
void printDebug(abcd_method &abcd, vector<vector<GammaParams> > &allyields, TString baseline,
                vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds);
TString Zbi(double Nobs, double Nbkg, double Eup_bkg, double Edown_bkg);
TString cutsToTex(TString cut);

void GetOptions(int argc, char *argv[]);

//// Defining st because older ntuples don't have it
const NamedFunc st("st", [](const Baby &b) -> NamedFunc::ScalarType{
    float stvar = b.ht();
    for (const auto &pt: *(b.leps_pt())) stvar += pt; 
    return stvar;
  });

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string ntupletag="";

  //// Capybara
  string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higtight/");
  string foldersig(bfolder+"/cms2r0/babymaker/babies/2016_08_10/TChiHH/merged_higmc_higtight/");
  string folderdata(bfolder+"/cms2r0/babymaker/babies/2016_08_10/data/merged_database_stdnj5/");

  Palette colors("txt/colors.txt", "default");

  // Cuts in baseline speed up the yield finding
  string baseline_s = "hig_drmax<2.2&&ntks==0&&njets>=4&&njets<=5&&!low_dphi&&nvleps==0";
  NamedFunc baseline=baseline_s;

  //// Use this process to make quick plots. Requires being run without split_bkg
  auto proc_bkg = Process::MakeShared<Baby_full>("All_bkg", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets_Tune*"+ntupletag+"*.root"},
    baseline && "stitch && pass");

  vector<string> sigMasses({"225", "300", "400", "700"});
  vector<shared_ptr<Process> > proc_sigs;
  for(size_t ind=0; ind<sigMasses.size(); ind++)
    proc_sigs.push_back(Process::MakeShared<Baby_full>("TChiHH("+sigMasses[ind]+")", Process::Type::signal, 2,
      {foldersig+"*mGluino-"+sigMasses[ind]+"_*"+ntupletag+"*.root"}, baseline && "stitch"));

  auto proc_ttbar = Process::MakeShared<Baby_full>("ttbar", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets*"+ntupletag+"*.root"},
    baseline && "stitch && pass");
  auto proc_singlet = Process::MakeShared<Baby_full>("singlet", Process::Type::background, colors("tt_2l"),
    {foldermc+"*_ST_*"+ntupletag+"*.root"},
    baseline && "stitch && pass");
  
  // Filling all other processes
  vector<string> vnames_other({"_WJetsToLNu", "_TTW", "_TTZ", "DYJetsToLL", 
	"_ZJet", "_ttHJetTobb", "_TTGJets", "_TTTT", 
	"_WH_HToBB", "_ZH_HToBB", "_WWTo", "_WZ", "_ZZ_", "QCD_HT*0_Tune", "QCD_HT*Inf_Tune"});
  set<string> names_other;
  for(auto name : vnames_other)
    names_other.insert(name = foldermc + "*" + name + "*" + ntupletag + "*.root");
  auto proc_other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    names_other, baseline && "stitch && pass");

  //string trigs = "(trig[4]||trig[8]||trig[13]||trig[33])";
  string trigs = "trig_ra4";
  if(skim.Contains("2015")) trigs = "(trig[4]||trig[8]||trig[28]||trig[14])";

  // Setting luminosity
  string jsonCuts = "nonblind";
  if(skim.Contains("2015")) lumi = 2.3;
  else if(json=="0p815"){
    lumi = 0.815;
    jsonCuts = "nonblind";
  } else if(json=="2p6"){
    lumi = 2.6;
    jsonCuts = "json2p6";
  } else if(json=="1p7"){
    lumi = 1.7;
    jsonCuts = "json4p0&&!json2p6";
  } else if(json=="4p3"){
    lumi = 4.3;
    jsonCuts = "json4p0";
  } else if(json=="7p65"){
    lumi = 7.65;
    jsonCuts = "json7p65";
  } else if(json=="12p9"){
    lumi = 12.9;
    jsonCuts = "json12p9";
  }
  if(mc_lumi!="") lumi = mc_lumi.Atof();


  if(only_method.Contains("old")) trigs = "(trig[4]||trig[8])";
  if(!skim.Contains("2015")) trigs += " && "+jsonCuts;

  auto proc_data = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {folderdata+"*"+ntupletag+"*.root"},baseline && trigs && "pass");

  vector<shared_ptr<Process> > all_procs = {proc_ttbar, proc_singlet, proc_other};
  //vector<shared_ptr<Process> > all_procs = {proc_bkg};
  if (do_signal){
    for(size_t ind=0; ind<proc_sigs.size(); ind++)
      all_procs.push_back(proc_sigs[ind]);
  }
  if(!only_mc) all_procs.push_back(proc_data);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining basic cuts //////////////////////////////////////////
  // baseline defined above

  ////// MET cuts
  TString c_lowmet  = "met>100 && met<=200";
  TString c_midmet  = "met>200 && met<=300";
  TString c_higmet  = "met>300";

  ////// Nb cuts
  TString c_2b="nbt==2&&nbm==2";
  TString c_3b="nbt>=2&&nbm==3&&nbl==3";
  TString c_4b="nbt>=2&&nbm>=3&&nbl>=4";

  ////// CR, SR cuts
  TString c_sr="hig_am>100&&hig_am<140&&hig_dm<40";
  TString c_cr="!("+c_sr+")";

  ////// ABCD cuts
  vector<TString> abcdcuts_std = {c_cr + " && 2bcuts",
                                  c_cr + " && nj_1l",
                                  c_sr + " && 2bcuts",
                                  c_sr + " && nj_1l"};


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining ABCD methods //////////////////////////////////////////
  vector<abcd_method> abcds;
  vector<TString> abcdcuts, metcuts, bincuts;
  PlotMaker pm;

  ///// Running over these methods
  vector<TString> methods = {"TTML", "MMMM"};

  if(only_method!="") methods = vector<TString>({only_method});

  for(size_t iabcd=0; iabcd<methods.size(); iabcd++) {
    TString method = methods[iabcd];
    TString basecuts = "", caption = "";

    if(method.Contains("TTML")){
      c_2b="nbt==2&&nbm==2";
      c_3b="nbt>=2&&nbm==3&&nbl==3";
      c_4b="nbt>=2&&nbm>=3&&nbl>=4";
      caption = "TTML method: 2 tight b-tags, 1 medium, 1 loose";
    } // TTML
    if(method.Contains("MMMM")){
      c_2b="nbm==2";
      c_3b="nbm==3";
      c_4b="nbm>=4";
      caption = "MMMM method: all medium b-tags";
    } // MMMM

    //// General assignments to all methods
    abcdcuts = abcdcuts_std;
    for(size_t ind=0; ind<abcdcuts.size(); ind++)
      abcdcuts[ind].ReplaceAll("2bcuts", c_2b);
    metcuts = vector<TString>{c_midmet, c_higmet};
    bincuts = vector<TString>{c_3b, c_4b};


    //////// Pushing all cuts to then find the yields
    abcds.push_back(abcd_method(method, metcuts, bincuts, abcdcuts, caption, basecuts));
    if(method.Contains("noint")) abcds.back().setIntNbNj(false);
 
    vector<TableRow> table_cuts;
    for(size_t icut=0; icut < abcds.back().allcuts.size(); icut++)
      table_cuts.push_back(TableRow(abcds.back().allcuts[icut].Data(), abcds.back().allcuts[icut].Data()));

    TString tname = "preds"; tname += iabcd;
    pm.Push<Table>(tname.Data(),  table_cuts, all_procs, true, false);
  } // Loop over ABCD methods

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////// Finding all yields ///////////////////////////////////////////////

  bool single_thread = false;
  if(single_thread) pm.multithreaded_ = false;
  pm.min_print_ = true;
  pm.MakePlots(lumi);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////// Calculating preds/kappas and printing table //////////////////////////////////////
  vector<TString> tablenames;
  for(size_t imethod=0; imethod<abcds.size(); imethod++) {
    Table * yield_table = static_cast<Table*>(pm.Figures()[imethod].get());
    // allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
    // if split_bkg [2/4] Other, [3/5] tt1l, [4/6] tt2l
    vector<vector<GammaParams> > allyields;
    if(!only_mc) allyields.push_back(yield_table->DataYield());
    else allyields.push_back(yield_table->BackgroundYield(lumi));
    allyields.push_back(yield_table->BackgroundYield(lumi));
    if(do_signal){
      for(size_t ind=0; ind<proc_sigs.size(); ind++)
	allyields.push_back(yield_table->Yield(proc_sigs[ind].get(), lumi));
    }
    if(split_bkg){
      allyields.push_back(yield_table->Yield(proc_other.get(), lumi));
      allyields.push_back(yield_table->Yield(proc_singlet.get(), lumi));
      allyields.push_back(yield_table->Yield(proc_ttbar.get(), lumi));
    }

    //// Calculating kappa and Total bkg prediction
    vector<vector<vector<float> > > kappas, preds;
    findPreds(abcds[imethod], allyields, kappas, preds);

    //// Makes table MC/Data yields, kappas, preds, Zbi
    if(!only_kappa) {
      TString fullname = printTable(abcds[imethod], allyields, kappas, preds, proc_sigs);
      tablenames.push_back(fullname);
    }

    //// Plotting kappa
    plotKappa(abcds[imethod], kappas);

    //// Print MC/Data yields, cuts applied, kappas, preds
    if(debug) printDebug(abcds[imethod], allyields, TString(baseline.Name()), kappas, preds);

  } // Loop over ABCD methods

  if(!only_kappa){
    //// Printing names of ouput files
    cout<<endl<<"===== Tables to be moved to the AN/PAS/paper:"<<endl;
    for(size_t ind=0; ind<tablenames.size(); ind++){
      TString name=tablenames[ind]; name.ReplaceAll("fulltable","table");
      cout<<" mv "<<name<<"  ${tables_folder}"<<endl;
    }
    cout<<endl<<"===== Tables that can be compiled"<<endl;
    for(size_t ind=0; ind<tablenames.size(); ind++)
      cout<<" pdflatex "<<tablenames[ind]<<"  > /dev/null"<<endl;
    cout<<endl;
  }

  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Finding "<<abcds.size()<<" tables took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
} // main

////////////////////////////////////////// End of main //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//// Prints table with results
// allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
// if split_bkg: [2/4] Other, [3/5] tt1l, [4/6] tt2l
TString printTable(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
                   vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds,
		   vector<shared_ptr<Process> > &proc_sigs){
  cout<<endl<<"Printing table (significance estimation can take a bit)"<<endl;

  //// Table general parameters
  int digits = 2;
  TString ump = " & ", blind_s = "$\\spadesuit$";

  size_t Nsig = proc_sigs.size(); // Number of signal points (for now it cannot be changed)
  bool do_zbi = true;
  if(!unblind) do_zbi = false;
  size_t Ncol = 3; // The only colums alwasy printed are the bin names, kappa, and MC bkg.
  if(split_bkg) Ncol += 3;
  if(do_signal) Ncol += Nsig;
  if(only_mc) {
    if(do_zbi) Ncol += Nsig;
  } else {
    Ncol += 3;
    if(do_zbi) Ncol++;
  }

  //// Setting output file name
  int digits_lumi = 1;
  if(lumi < 1) digits_lumi = 3;
  if(lumi >= 15) digits_lumi = 0;
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
  if(split_bkg) out << "|ccc";
  out << "|cc";
  if(!only_mc) {
    out << "|ccc";
    if(do_zbi) out<<"c";
  } else {
    if(do_signal)
      for(size_t ind=0; ind<Nsig; ind++)
	out<<"|c"<<(do_zbi?"c":"");
  }
  out<<"}\\hline\\hline\n";
  out<<"${\\cal L}="<<lumi_s<<"$ fb$^{-1}$ ";
  if(split_bkg) out << " & Other & Single $t$ & $t\\bar{t}$ ";
  out << "& $\\kappa$ & MC bkg.";
  if(!only_mc) out << " & Pred.& Obs. & Obs./MC "<<(do_zbi?"& Signi.":"");
  else if(do_signal) 
    for(size_t ind=0; ind<Nsig; ind++) {
      TString signame = proc_sigs[ind]->name_.c_str();
      if(do_zbi) out << "& \\multicolumn{2}{c"<<(ind<Nsig-1?"|":"")<<"}{" << signame <<"}";
      else  out << "& " << signame;
    }
  out << " \\\\ \\hline\\hline\n";

  vector<TString> binNames({"SBD, 2b", "SBD, xb", "HIG, 2b", "HIG, xb"});
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
	TString binName = binNames[iabcd];
	if(ibin==0) binName.ReplaceAll("xb", "3b");
	else binName.ReplaceAll("xb", "4b");
        out << binName;
        // if(iabcd%2==0 && abcd.int_nbnj)
        //   out << "All "<<(abcd.bincuts[iplane][ibin].Contains("nbm")?"\\nb, ":"")<<"\\njets" ;
        // else {
        //   if(abcd.method.Contains("2lonly") && iabcd>=2) out<<cutsToTex(abcd.lowerNjets(abcd.bincuts[iplane][ibin]));
        //   else if(abcd.method.Contains("2lveto") && iabcd>=2){
        //     if(abcd.bincuts[iplane][ibin].Contains("6")) out<<"Low \\njets";
        //     else out<<"High \\njets";
        //   } else out<<cutsToTex(abcd.bincuts[iplane][ibin]);
        // }
        //// Printing Other, tt1l, tt2l
        if(split_bkg){
          size_t offset = (do_signal?Nsig:0);
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
        if(!only_mc){
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
          if(do_zbi && iabcd==3) out << ump << Zbi(allyields[0][index].Yield(), preds[iplane][ibin][0], 
						   preds[iplane][ibin][1], preds[iplane][ibin][2]);
	  //// Printing signal yields
	  if(do_signal)
	    for(size_t ind=0; ind<Nsig; ind++) 
	      out<<ump<<RoundNumber(allyields[2+ind][index].Yield(), digits);
        } else {// if not only_mc
	  if(do_signal){
	    for(size_t ind=0; ind<Nsig; ind++) {
	      if(ind>0) out<<ump;
	      out<<RoundNumber(allyields[2+ind][index].Yield(), digits);
	      if(do_zbi){
		out << ump;
		if(iabcd==3) 
		  out<<Zbi(allyields[0][index].Yield()+allyields[2+ind][index].Yield(),preds[iplane][ibin][0],
			   preds[iplane][ibin][1], preds[iplane][ibin][2]);
	      } // if do_zbi
	    } // Loop over signals
	  } // if do_signal
	} //if only_mc

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
  cut.ReplaceAll("met>250&&met<=300", "250<met<=300");
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
TString Zbi(double Nobs, double Nbkg, double Eup_bkg, double Edown_bkg){
  TString zbi_s;
  if(false){ // Old, bad Zbi
    double Nsig = Nobs-Nbkg;
    double zbi = RooStats::NumberCountingUtils::BinomialExpZ(Nsig, Nbkg, Eup_bkg/Nbkg);
    if(Nbkg==0) zbi = RooStats::NumberCountingUtils::BinomialWithTauExpZ(Nsig, Nbkg, 1/Eup_bkg);
    if(zbi<0) zbi=0;
    zbi_s = RoundNumber(zbi,1);
    if(zbi_s!="-") zbi_s = "$"+zbi_s+"\\sigma$";
    if(Nsig<=0 || Eup_bkg<=0) zbi_s = "-";
  } else zbi_s = "$"+RoundNumber(Significance(Nobs, Nbkg, Eup_bkg, Edown_bkg),1)+"\\sigma$";
  //cout<<"Zbi for Nobs "<<Nobs<<", Nbkg "<<Nbkg<<", Ebkg "<<Eup_bkg<<" is "<<zbi_s<<endl;
  return zbi_s;
}

//// Makes kappa plots
void plotKappa(abcd_method &abcd, vector<vector<vector<float> > > &kappas){

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
  for(size_t iplane=0; iplane < kappas.size(); iplane++) {
    k_ordered.push_back(vector<vector<kmarker> >());
    for(size_t ibin=0; ibin < kappas[iplane].size(); ibin++){
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

  //// Plotting kappas
  TCanvas can("can","");
  can.SetFillStyle(4000);
  TLine line; line.SetLineWidth(2); line.SetLineStyle(2);
  TLatex label; label.SetTextSize(0.05); label.SetTextFont(42); label.SetTextAlign(23);
  if(k_ordered.size()>3) label.SetTextSize(0.04);


  float minx = 0.5, maxx = nbins+0.5, miny = 0, maxy = 2.4;
  if(label_up) maxy = 2.6;
  TH1D histo("histo", "", nbins, minx, maxx);
  histo.SetMinimum(miny);
  histo.SetMaximum(maxy);
  histo.GetYaxis()->CenterTitle(true);
  histo.GetXaxis()->SetLabelOffset(0.008);
  histo.SetYTitle("#kappa");
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
    else label.DrawLatex((2*bin-k_ordered[iplane].size()+1.)/2., -0.26, cutsToLabel(abcd.planecuts[iplane]));
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
    graph[indb].SetMarkerStyle(ind_bcuts[indb].style); graph[indb].SetMarkerSize(1.1);
    graph[indb].SetMarkerColor(ind_bcuts[indb].color);
    graph[indb].SetLineColor(ind_bcuts[indb].color); graph[indb].SetLineWidth(2);
    graph[indb].Draw("p0 same");
    leg.AddEntry(&graph[indb], cutsToLabel(ind_bcuts[indb].cut), "p");
  } // Loop over TGraphs
  if(ind_bcuts.size()>1) leg.Draw();

  //// Drawing CMS labels and line at 1
  TLatex cmslabel;
  cmslabel.SetTextSize(0.06);
  cmslabel.SetNDC(kTRUE);
  cmslabel.SetTextAlign(11);
  cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}");
  cmslabel.SetTextAlign(31);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.015,"#font[42]{13 TeV}");

  line.SetLineStyle(3); line.SetLineWidth(1);
  line.DrawLine(minx, 1, maxx, 1);

  TString fname="plots/kappa_" + abcd.method;
  fname += ".pdf";
  can.SaveAs(fname);
  cout<<endl<<" open "<<fname<<endl;

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

  int digits = 3;
  cout<<endl<<endl<<"=================== Printing cuts for method "<<abcd.method<<" ==================="<<endl;
  cout<<"-- Baseline cuts: "<<baseline<<endl;
  for(size_t iplane=0; iplane < abcd.planecuts.size(); iplane++) {
    cout<<endl<<" **** Plane "<<abcd.planecuts[iplane]<<" ***"<<endl;
    for(size_t ibin=0; ibin < abcd.bincuts[iplane].size(); ibin++){
      for(size_t iabcd=0; iabcd < abcd.abcdcuts.size(); iabcd++){
        size_t index = abcd.indexBin(iplane, ibin, iabcd);
        cout<<"MC: "<<setw(8)<<RoundNumber(allyields[1][index].Yield(),digits)
            <<"  Data: "<<setw(4)<<RoundNumber(allyields[0][index].Yield(), 0)
            <<"  - "<< abcd.allcuts[index]<<endl;
      } // Loop over ABCD cuts
      cout<<"Kappa = "<<RoundNumber(kappas[iplane][ibin][0],digits)<<"+"<<RoundNumber(kappas[iplane][ibin][1],digits)
          <<"-"<<RoundNumber(kappas[iplane][ibin][2],digits)<<", Prediction = "
          <<RoundNumber(preds[iplane][ibin][0],digits)<<"+"<<RoundNumber(preds[iplane][ibin][1],digits)
          <<"-"<<RoundNumber(preds[iplane][ibin][2],digits)<<endl;
      cout<<endl;
    } // Loop over bin cuts
  } // Loop over plane cuts

} // printDebug

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"method", required_argument, 0, 'm'},  // Method to run on (if you just want one)
      {"lumi", required_argument, 0, 'l'},    // Luminosity to normalize MC with (no data)
      {"skim", required_argument, 0, 's'},    // Which skim to use: standard, 2015 data
      {"json", required_argument, 0, 'j'},    // Which JSON to use: 0p815, 2p6, 4p0, 7p65, 12p9
      {"split_bkg", no_argument, 0, 'b'},     // Prints Other, tt1l, tt2l contributions
      {"no_signal", no_argument, 0, 'n'},     // Does not print signal columns
      {"do_leptons", no_argument, 0, 'p'},    // Does tables for e/mu/emu as well
      {"unblind", no_argument, 0, 'u'},       // Unblinds R4/D4
      {"only_mc", no_argument, 0, 'o'},       // Uses MC as data for the predictions
      {"only_kappa", no_argument, 0, 'k'},    // Only plots kappa (no table)
      {"debug", no_argument, 0, 'd'},         // Debug: prints yields and cuts used
      {"only_dilepton", no_argument, 0, '2'}, // Makes tables only for dilepton tests
      {"ht", no_argument, 0, 0},            // Cuts on ht>500 instead of st>500
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "m:s:j:udbnl:p2ok", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'm':
      only_method = optarg;
      break;
    case 'l':
      mc_lumi = optarg;
      only_mc = true;
      break;
    case 'k':
      only_kappa = true;
      only_mc = true;
      break;
    case 's':
      skim = optarg;
      break;
    case 'j':
      json = optarg;
      break;
    case 'b':
      split_bkg = false;
      break;
    case 'o':
      only_mc = true;
      break;
    case '2':
      only_dilepton = true;
      break;
    case 'p':
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
    case 0:
      optname = long_options[option_index].name;
      if(optname == "ht"){
        do_ht = true;
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
