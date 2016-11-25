///// plot_ratios: plots rMJ and rmT, and combinations of these

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
#include "core/styles.hpp"
#include "core/plot_opt.hpp"
#include "core/functions.hpp"

using namespace std;
using namespace Functions;

namespace{
  bool debug = false;
  TString skim = "both1l";
  float lumi=40.;

  enum Regions {sbd2b, hig2b, sbd3b, hig3b, sbd4b, hig4b};
  enum files {all, ttbar, vjets};
  int ifilesDen = ttbar;
  int ifilesNum = ttbar;
  struct oneplot{
    TString name;
    TString baseline;
    vector<TString> bincuts;
  };
}


void plotRatio(vector<vector<vector<GammaParams> > > &allyields, oneplot &plotdef,
	       vector<vector<vector<int> > > &indices, vector<TString> &leglabels,
	       vector<shared_ptr<Process> > &procs);
void printDebug(vector<vector<TString> > &allcuts, vector<vector<vector<GammaParams> > > &allyields, 
		TString baseline, vector<shared_ptr<Process> > &procs);
void GetOptions(int argc, char *argv[]);

NamedFunc nb_tru("nb_tru",[](const Baby &b) -> NamedFunc::ScalarType{
  int nbtru(0);
  for (unsigned i(0); i<b.jets_pt()->size(); i++){
    if (!b.jets_h1()->at(i) && !b.jets_h2()->at(i)) continue;
    if (b.jets_hflavor()->at(i)==5) nbtru++;
  }
  return nbtru;
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

  string ntupletag="higloose";
  if(skim.Contains("nlep1")) ntupletag="nlep1";
  if(skim.Contains("both")) ntupletag="";
  string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/");

  Palette colors("txt/colors.txt", "default");

  // Cuts in baseline speed up the yield finding
  string zcuts = "nleps==2 && ((mumuv_m>80&&mumuv_m<100) || (elelv_m>80&&elelv_m<100))";
  string baseline="pass && stitch";
  // if(skim.Contains("1l")) baseline += " && nleps==1";
  // if(skim.Contains("zll")) baseline += " && "+zcuts;
  NamedFunc baselinef = baseline;

  set<string> allfiles = {foldermc+"*_TTJets*Lept*"+ntupletag+"*.root", foldermc+"*_TTJets_HT*"+ntupletag+"*.root",
    foldermc+"*_WJetsToLNu*"+ntupletag+"*.root",
    foldermc+"*_ST_*"+ntupletag+"*.root",
    foldermc+"*_TTW*"+ntupletag+"*.root", foldermc+"*_TTZ*"+ntupletag+"*.root",
    foldermc+"*_TTGJets*"+ntupletag+"*.root",foldermc+"*_TTTT*"+ntupletag+"*.root",
    foldermc+"*QCD_HT*Inf_Tune*"+ntupletag+"*.root", 
    foldermc+"*QCD_HT*0_Tune*"+ntupletag+"*.root",foldermc+"*DYJetsToLL*"+ntupletag+"*.root",
    foldermc+"*_ZJet*"+ntupletag+"*.root",foldermc+"*_ttHJetTobb*"+ntupletag+"*.root",
    foldermc+"*_WH_HToBB*"+ntupletag+"*.root",foldermc+"*_ZH_HToBB*"+ntupletag+"*.root",
    foldermc+"*_WWTo*"+ntupletag+"*.root",foldermc+"*_WZ*"+ntupletag+"*.root",
    foldermc+"*_ZZ_*"+ntupletag+"*.root"
  };
  set<string> nonttfiles = {foldermc+"*_WJetsToLNu*"+ntupletag+"*.root",
    foldermc+"*_ST_*"+ntupletag+"*.root",
    foldermc+"*_TTW*"+ntupletag+"*.root", foldermc+"*_TTZ*"+ntupletag+"*.root",
    foldermc+"*_TTGJets*"+ntupletag+"*.root",foldermc+"*_TTTT*"+ntupletag+"*.root",
    foldermc+"*QCD_HT*Inf_Tune*"+ntupletag+"*.root", 
    foldermc+"*QCD_HT*0_Tune*"+ntupletag+"*.root",foldermc+"*DYJetsToLL*"+ntupletag+"*.root",
    foldermc+"*_ZJet*"+ntupletag+"*.root",foldermc+"*_ttHJetTobb*"+ntupletag+"*.root",
    foldermc+"*_WH_HToBB*"+ntupletag+"*.root",foldermc+"*_ZH_HToBB*"+ntupletag+"*.root",
    foldermc+"*_WWTo*"+ntupletag+"*.root",foldermc+"*_WZ*"+ntupletag+"*.root",
    foldermc+"*_ZZ_*"+ntupletag+"*.root"
  };
  set<string> ttfiles = {foldermc+"*_TTJets*Lept*"+ntupletag+"*.root", foldermc+"*_TTJets_HT*"+ntupletag+"*.root",
			 foldermc+"*_TTW*"+ntupletag+"*.root", foldermc+"*_TTZ*"+ntupletag+"*.root",
			 foldermc+"*_TTGJets*"+ntupletag+"*.root",foldermc+"*_ttHJetTobb*"+ntupletag+"*.root",
			 foldermc+"*_TTTT*"+ntupletag+"*.root"};
   set<string> vfiles = {foldermc+"*DYJetsToLL*"+ntupletag+"*.root",foldermc+"*_ZJet*"+ntupletag+"*.root",
   			foldermc+"*_WJetsToLNu*"+ntupletag+"*.root"};
   //set<string> vfiles = {foldermc+"*_ZJet*"+ntupletag+"*.root"};

  //allfiles = set<string>({foldermc+"*_TTJets_Tune*"+ntupletag+"*.root"});
  // allfiles = nonttfiles;

  cout<<endl<<"Doing denominator "<<ifilesDen<<" and numerator "<<ifilesNum<<endl<<endl;

  string dentitle = "All bkg", numtitle="";
  set<string> denfiles = allfiles, numfiles = allfiles;
  if(ifilesDen == ttbar) {
    denfiles = ttfiles;
    dentitle = "t#bar{t}";
  }
  if(ifilesDen == vjets) {
    denfiles = vfiles;
    dentitle = "V+jets";
  }
  if(ifilesNum == ttbar) {
    numfiles = ttfiles;
    numtitle = "t#bar{t}, ";
  }
  if(ifilesNum == vjets) {
    numfiles = vfiles;
    numtitle = "Vjets, ";
  }

  vector<shared_ptr<Process> > procs;
  //// First entry is for low-mT region (kappa denominator)
  procs.push_back(Process::MakeShared<Baby_full>(dentitle, Process::Type::background, 1, denfiles, 
						 baselinef));

  //// Processes for the high-mT region (kappa numerator)
  if(skim.Contains("both")) {
    string leptitle = "N_{lep} = 1",  lepcuts = "nleps==1";
    if(skim.Contains("zll")) {
      leptitle = "Z #rightarrow ll";
      lepcuts = zcuts;
    }
    procs.push_back(Process::MakeShared<Baby_full>
		    (numtitle+"N_{lep} = 0", Process::Type::background, kBlue, 
		     numfiles, baselinef && "nvleps==0"));
    procs.push_back(Process::MakeShared<Baby_full>
		    (numtitle+leptitle, Process::Type::background, kRed+1, 
		     numfiles, baselinef && lepcuts));
  } else {
    // procs.push_back(Process::MakeShared<Baby_full>
    // 		    (numtitle+"#leq1 true B-hadron", Process::Type::background, kPink+2, 
    // 		     numfiles, baselinef && nb_tru<=1));
    // procs.push_back(Process::MakeShared<Baby_full>
    // 		    (numtitle+"2 true B-hadron", Process::Type::background, kOrange-4, 
    // 		     numfiles, baselinef && nb_tru==2));
    // procs.push_back(Process::MakeShared<Baby_full>
    // 		    (numtitle+"3 true B-hadron", Process::Type::background, kTeal-8, 
    // 		     numfiles, baselinef && nb_tru==3));
    // procs.push_back(Process::MakeShared<Baby_full>
    // 		    (numtitle+"#geq4 true B-hadron", Process::Type::background, kAzure-4, 
    // 		     numfiles, baselinef && nb_tru>=4));

    procs.push_back(Process::MakeShared<Baby_full>
		    ("All bkg",Process::Type::background,1,numfiles,
		     baselinef));
    procs.push_back(Process::MakeShared<Baby_full>
		    ("t#bar{t}",Process::Type::background,colors("tt_1l"),ttfiles,
		     baselinef));
  }



/////////////////////////////////////////////////////////////////////////////////////////////////////////
 /////////////////////////////////////////// Defining cuts ///////////////////////////////////////////////
  // baseline defined above
  TString cutlep = "nvleps==0";
  if(skim.Contains("nlep1")) cutlep = "nleps==1&&mt<100";
  if(skim.Contains("both")) cutlep = "1";


  // Makes a plot for each vector in plotcuts
  vector<oneplot> plotcuts({
	{"met",cutlep,{"met>150&&met<=200","met>200&&met<=300", "met>300"}},
	{"met",cutlep+"&&hig_drmax<2.2",{"met>150&&met<=200","met>200&&met<=300", "met>300"}},
	{"met",cutlep+"&&ntks==0 && !low_dphi && hig_drmax<2.2",{"met>150&&met<=200",
	      "met>200&&met<=300", "met>300"}},
	});

  TString c_2b="nbt==2&&nbm==2", c_3b="nbt>=2&&nbm==3&&nbl==3", c_4b="nbt>=2&&nbm>=3&&nbl>=4";
  //TString c_2b="nbm==2", c_3b="nbm==3", c_4b="nbm>=4";
  TString c_hig="hig_am>100&&hig_am<140&&hig_dm<40", c_sbd="!("+c_hig+") && hig_am<200 && hig_dm<40";
  //TString c_hig="hig_am>100&&hig_am<140&&hig_dm<40", c_sbd="!("+c_hig+")";
  TString ump=" && ";

  vector<TString> abcdcuts = {c_2b +ump+ c_sbd, c_2b +ump+ c_hig, c_3b +ump+ c_sbd, c_3b +ump+ c_hig, 
			      c_4b +ump+ c_sbd, c_4b +ump+ c_hig};
  size_t Nabcd = abcdcuts.size();

  PlotMaker pm;
  vector<vector<vector<TString> > > allcuts(plotcuts.size(), vector<vector<TString> > (Nabcd));
  for(size_t iplot=0; iplot<plotcuts.size(); iplot++){
    for(size_t iabcd=0; iabcd<abcdcuts.size(); iabcd++){
      vector<TableRow> table_cuts;
      for(size_t ibin=0; ibin<plotcuts[iplot].bincuts.size(); ibin++){
	TString totcut=plotcuts[iplot].baseline+" && "+plotcuts[iplot].bincuts[ibin]+" && "+abcdcuts[iabcd];
	table_cuts.push_back(TableRow("", totcut.Data()));
	allcuts[iplot][iabcd].push_back(totcut);
      } // Loop over bins
      TString tname = "rmj"; tname += iplot; tname += iabcd;
      pm.Push<Table>(tname.Data(),  table_cuts, procs, false, false);
    } // Loop over abcdcuts
  } // Loop over plots

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////// Finding all yields ///////////////////////////////////////////////

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////// Calculating preds/kappas and printing table //////////////////////////////////////

  for(size_t iplot=0; iplot<plotcuts.size(); iplot++){
    vector<vector<vector<GammaParams> > > allyields(procs.size(), vector<vector<GammaParams> >(Nabcd));
    for(size_t iabcd=0; iabcd<abcdcuts.size(); iabcd++){
      Table * yield_table = static_cast<Table*>(pm.Figures()[iplot*Nabcd+iabcd].get());
      for(size_t ibkg=0; ibkg<procs.size(); ibkg++)
	allyields[ibkg][iabcd] = yield_table->Yield(procs[ibkg].get(), lumi);
    } // Loop over ABCD cuts

    //// Print MC/Data yields, cuts applied, kappas, preds
    if(debug) printDebug(allcuts[iplot], allyields, baseline, procs);

    //// RMJ/RMJ
    vector<vector<vector<int> > > indices;
    vector<TString> leglabels;

    //// 3b kappa
    indices.clear(); leglabels.clear();
    for(int ibkg=1; ibkg<static_cast<int>(procs.size()); ibkg++) {
      int idenom=0;
      if(skim.Contains("both") || skim.Contains("higloose") || skim.Contains("nlep1")) idenom = ibkg;
      indices.push_back(vector<vector<int> >({{ibkg, hig3b, 1}, {ibkg, sbd3b, -1}, 
								  {idenom, hig2b, -1}, {idenom, sbd2b, 1}}));
      leglabels.push_back(procs[ibkg]->name_);
    }
    plotRatio(allyields, plotcuts[iplot], indices, leglabels, procs);

    //// 4b kappa
    indices.clear(); leglabels.clear();
    for(int ibkg=1; ibkg<static_cast<int>(procs.size()); ibkg++) {
      int idenom=0;
      if(skim.Contains("both") || skim.Contains("higloose") || skim.Contains("nlep1")) idenom = ibkg;
      indices.push_back(vector<vector<int> >({{ibkg, hig4b, 1}, {ibkg, sbd4b, -1}, 
								  {idenom, hig2b, -1}, {idenom, sbd2b, 1}}));
      leglabels.push_back(procs[ibkg]->name_);
    }
    plotRatio(allyields, plotcuts[iplot], indices, leglabels, procs);


  } // Loop over plots


  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Finding "<<plotcuts.size()<<" plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
} // main

////////////////////////////////////////// End of main //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////




void plotRatio(vector<vector<vector<GammaParams> > > &allyields, oneplot &plotdef,
	       vector<vector<vector<int> > > &indices, vector<TString> &leglabels, 
	       vector<shared_ptr<Process> > &procs){

  size_t ngraphs = indices.size();
  size_t nbins = allyields[0][0].size();

  //// Finding all ratios for all graphs
  float val(1.), valup(1.), valdown(1.);
  vector<vector<vector<float> > > ratios(ngraphs);
  float maxr=-1., minr=1e6;
  for(size_t igraph=0; igraph<ngraphs; igraph++){
    // Finding powers to calculate ratio
    vector<float> powers;
    for(size_t ipow=0; ipow<indices[igraph].size(); ipow++) powers.push_back(indices[igraph][ipow][2]);

    // Finding ratios for each bin
    for(size_t ibin=0; ibin<nbins; ibin++){
      vector<vector<float> > entries;
      vector<vector<float> > weights;
      for(size_t ind=0; ind<indices[igraph].size(); ind++) {
	size_t ibkg = indices[igraph][ind][0];
	size_t iabcd = indices[igraph][ind][1];
        entries.push_back(vector<float>());
        weights.push_back(vector<float>());
        entries.back().push_back(allyields[ibkg][iabcd][ibin].NEffective());
        weights.back().push_back(allyields[ibkg][iabcd][ibin].Weight());
      } // Loop over indices

      // Throwing toys to find ratios and uncertainties
      val = calcKappa(entries, weights, powers, valdown, valup);
      if(valdown<0) valdown = 0;
      ratios[igraph].push_back(vector<float>({val, valdown, valup}));
      if(maxr < val+valup) maxr = val+valup;
      if(minr > val-valdown) minr = val-valdown;
    } // Loop over bins
  } // Loop over graphs

  //// Finding ytitle
  TString ytitle="Ratio";
  if(indices[0].size()==2){
    size_t ind0=indices[0][0][1], ind1=indices[0][1][1];
    if((ind0==hig2b && ind1==sbd2b)) ytitle = "Hig/Sbd";
  }
  if(indices[0].size()==4){
    size_t ind0=indices[0][0][1], ind1=indices[0][1][1];
    size_t ind2=indices[0][2][1], ind3=indices[0][3][1];
    if((ind0==hig4b&&ind1==sbd4b && ind2==hig2b&&ind3==sbd2b)) {
      ytitle = "[4b Hig/Sbd] / [2b Hig/Sbd("+TString(procs[0]->name_)+")]";
      ytitle = "4b #kappa";
    }
    if((ind0==hig3b&&ind1==sbd3b && ind2==hig2b&&ind3==sbd2b)) {
      ytitle = "[3b Hig/Sbd] / [2b Hig/Sbd("+TString(procs[0]->name_)+")]";
      ytitle = "3b #kappa";
    }
  }
  //// Setting plot style
  PlotOpt opts("txt/plot_styles.txt", "Ratio");
  setPlotStyle(opts);

  //// Plotting kappas
  TCanvas can("can","");
  TLine line; line.SetLineWidth(2); line.SetLineStyle(2);
  TLatex label; label.SetTextSize(0.05); label.SetTextFont(42); label.SetTextAlign(23);

  float minx = 0.5, maxx = nbins+0.5, miny = 0, maxy = maxr*1.2;
  if(maxy>3) maxy = 3;
  if(maxy<2) maxy = 2;
  if(plotdef.baseline=="1") plotdef.baseline = "";
  TH1D histo("histo", CodeToRootTex(plotdef.baseline.Data()).c_str(), nbins, minx, maxx);
  histo.SetMinimum(miny);
  histo.SetMaximum(maxy);
  histo.GetYaxis()->CenterTitle(true);
  histo.GetXaxis()->SetLabelOffset(0.008);
  histo.GetYaxis()->SetTitleSize(0.09);
  histo.GetYaxis()->SetTitleOffset(0.7);
  histo.SetYTitle(ytitle);
  histo.Draw();

  //// Filling vx, vy vectors with kappa coordinates. Each nb cut is stored in a TGraphAsymmetricErrors
  vector<vector<double> > vx(ngraphs), vexh(ngraphs), vexl(ngraphs);
  vector<vector<double> > vy(ngraphs), veyh(ngraphs), veyl(ngraphs);
  for(size_t ibin=0; ibin<nbins; ibin++){
    histo.GetXaxis()->SetBinLabel(ibin+1, CodeToRootTex(plotdef.bincuts[ibin].Data()).c_str());
    // xval is the x position of the first marker in the group
    double xval = ibin+1, minxb = 0.15, binw = 0;
    // If there is more than one point in the group, it starts minxb to the left of the center of the bin
    // binw is the distance between points in the njets group
    if(ngraphs>1) {
      xval -= minxb;
      binw = 2*minxb/(ngraphs-1);
    }
    for(size_t igraph=0; igraph<ngraphs; igraph++){
      vx[igraph].push_back(xval);
      xval += binw;
      vexl[igraph].push_back(0);
      vexh[igraph].push_back(0);
      vy[igraph]  .push_back(ratios[igraph][ibin][0]);
      veyl[igraph].push_back(ratios[igraph][ibin][1]);
      veyh[igraph].push_back(ratios[igraph][ibin][2]);
    } // Loop over TGraphs
  } // Loop over bin cuts

  //// Drawing legend and TGraphs
  double legX(opts.LeftMargin()+0.023), legY(1-opts.TopMargin()-0.03), legSingle = 0.05;
  double legW = 0.19*ngraphs, legH = legSingle;
  int Ncol = ngraphs;
  if(ngraphs>3) {
    legH *= 2;
    Ncol = (ngraphs+1)/2;
    legW = 0.25*Ncol;
  }
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(opts.LegendEntryHeight()*1.2); leg.SetFillColor(0);
  leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(42);
  leg.SetNColumns(Ncol);

  Palette colors("txt/colors.txt", "default");
  vector<int> mcolors, styles;
  for(size_t ibkg=1; ibkg<procs.size(); ibkg++) {
    mcolors.push_back(procs[ibkg]->color_);
    styles.push_back(19+ibkg);
  }
  TGraphAsymmErrors graph[20]; // There's problems with vectors of TGraphs, so using an array
  for(size_t igraph=0; igraph<ngraphs; igraph++){
    graph[igraph] = TGraphAsymmErrors(vx[igraph].size(), &(vx[igraph][0]), &(vy[igraph][0]),
                                    &(vexl[igraph][0]), &(vexh[igraph][0]), &(veyl[igraph][0]), &(veyh[igraph][0]));
    graph[igraph].SetMarkerStyle(styles[igraph]); 
    if(leglabels[igraph].Contains("All")) graph[igraph].SetMarkerStyle(21);
    graph[igraph].SetMarkerSize(1.4);
    graph[igraph].SetMarkerColor(mcolors[igraph]);
    graph[igraph].SetLineColor(mcolors[igraph]); graph[igraph].SetLineWidth(2);
    graph[igraph].Draw("p0 same");
    leg.AddEntry(&graph[igraph], leglabels[igraph], "p");
  } // Loop over TGraphs
  leg.Draw();

  //// Drawing CMS labels and line at 1
  TLatex cmslabel;
  cmslabel.SetTextSize(0.06);
  cmslabel.SetNDC(kTRUE);
  cmslabel.SetTextAlign(11);
  //cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}");
  cmslabel.SetTextAlign(31);
  //cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.015,"#font[42]{13 TeV}");

  line.SetLineStyle(3); line.SetLineWidth(1);
  line.DrawLine(minx, 1, maxx, 1);

  TString fname = "plots/ratio_"+CodeToPlainText(ytitle.Data())+"_"+plotdef.name+"_"
    +CodeToPlainText(plotdef.baseline.Data())+"_NumAllbkg.pdf";
  if(ifilesNum == ttbar) fname.ReplaceAll("NumAllbkg", "Numttbar");
  if(ifilesNum == vjets) fname.ReplaceAll("NumAllbkg", "NumVjets");
  can.SaveAs(fname);
  cout<<endl<<" open "<<fname<<endl;

} // plotRatio



// allyields: [0] All bkg, [1] tt1l, [2] tt2l, [3] other
void printDebug(vector<vector<TString> > &allcuts, vector<vector<vector<GammaParams> > > &allyields, 
		TString baseline, vector<shared_ptr<Process> > &procs){
  int digits = 3;
  cout<<endl<<endl<<"============================ Printing cuts  ============================"<<endl;
  cout<<"-- Baseline cuts: "<<baseline<<endl<<endl;
  for(size_t ibin=0; ibin<allcuts[0].size(); ibin++){
    for(size_t iabcd=0; iabcd<allcuts.size(); iabcd++){
      for(size_t ibkg=0; ibkg<procs.size(); ibkg++)
	cout<<procs[ibkg]->name_<<": "    <<setw(9)<<RoundNumber(allyields[ibkg][iabcd][ibin].Yield(), digits)<<", ";
      cout<<"  - "<< allcuts[iabcd][ibin]<<endl;
    } // Loop over ABCD cuts
    cout<<endl;
  } // Loop over bin cuts

} // printDebug


void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"numerator", required_argument, 0, 'n'},    // Ntuples to use in the numerator
      {"denominator", required_argument, 0, 'd'},    // Ntuples to use in the denominator
      {"lumi", required_argument, 0, 'l'},    // Luminosity to normalize MC with (no data)
      {"skim", required_argument, 0, 's'},    // Which skim to use: standard, met150, 2015 data
      {"debug", no_argument, 0, 'g'},         // Debug: prints yields and cuts used
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s:l:gn:d:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'n':
      ifilesNum = atoi(optarg);
      break;
    case 'd':
      ifilesDen = atoi(optarg);
      break;
    case 'l':
      lumi = atof(optarg);
      break;
    case 's':
      skim = optarg;
      break;
    case 'g':
      debug = true;
      break;
    case 0:
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
