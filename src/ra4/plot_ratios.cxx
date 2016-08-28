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

using namespace std;

namespace{
  bool debug = false;
  TString skim = "both";
  float lumi=40.;

  enum Regions {r1, r2, r3, r4};
  enum Bkgs {bkg, tt1l, tt2l, other};

  struct oneplot{
    TString name;
    TString baseline;
    vector<TString> bincuts;
  };
}

void plotRatio(vector<vector<vector<GammaParams> > > &allyields, oneplot &plotdef,
	       vector<vector<vector<int> > > &indices, vector<TString> &leglabels){

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
      ratios[igraph].push_back(vector<float>({val, valup, valdown}));
      if(maxr < val+valup) maxr = val+valup;
      if(minr > val-valdown) minr = val-valdown;
    } // Loop over bins
  } // Loop over graphs

  //// Finding ytitle
  TString ytitle="Ratio";
  if(indices[0].size()==2){
    size_t ind0=indices[0][0][1], ind1=indices[0][1][1];
    if((ind0==r3&&ind1==r1) || (ind0==r4&&ind1==r2)) ytitle = "R(m_{T})";
    if((ind0==r4&&ind1==r3) || (ind0==r2&&ind1==r1)) ytitle = "R(M_{J})";
  }
  if(indices[0].size()==4){
    size_t ind0=indices[0][0][1], ind1=indices[0][1][1];
    size_t ind2=indices[0][2][1], ind3=indices[0][3][1];
    if((ind0==r4&&ind1==r3&&ind2==r2&&ind3==r1)) ytitle = "R(M_{J}^{high}) / R[M_{J}^{low}(t#bar{t}1l)]";
  }
  //// Setting plot style
  PlotOpt opts("txt/plot_styles.txt", "Ratio");
  setPlotStyle(opts);

  //// Plotting kappas
  TCanvas can("can","");
  TLine line; line.SetLineWidth(2); line.SetLineStyle(2);
  TLatex label; label.SetTextSize(0.05); label.SetTextFont(42); label.SetTextAlign(23);

  float minx = 0.5, maxx = nbins+0.5, miny = 0, maxy = maxr*1.2;
  TH1D histo("histo", "", nbins, minx, maxx);
  histo.SetMinimum(miny);
  histo.SetMaximum(maxy);
  histo.GetYaxis()->CenterTitle(true);
  histo.GetXaxis()->SetLabelOffset(0.008);
  histo.SetYTitle(ytitle);
  histo.Draw();

  //// Filling vx, vy vectors with kappa coordinates. Each nb cut is stored in a TGraphAsymmetricErrors
  vector<vector<double> > vx(ngraphs), vexh(ngraphs), vexl(ngraphs);
  vector<vector<double> > vy(ngraphs), veyh(ngraphs), veyl(ngraphs);
  for(size_t ibin=0; ibin<nbins; ibin++){
    histo.GetXaxis()->SetBinLabel(ibin+1, cutsToLabel(plotdef.bincuts[ibin]));
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
  double legW = 0.19*ngraphs, legH = legSingle*(ngraphs+1)/2;
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(opts.LegendEntryHeight()); leg.SetFillColor(0);
  leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(42);
  leg.SetNColumns(ngraphs);

  Palette colors("txt/colors.txt", "default");
  vector<int> mcolors({4, 2, kGreen+3, kMagenta+2});
  vector<int> styles({20, 21, 22, 23});
  TGraphAsymmErrors graph[20]; // There's problems with vectors of TGraphs, so using an array
  for(size_t igraph=0; igraph<ngraphs; igraph++){
    graph[igraph] = TGraphAsymmErrors(vx[igraph].size(), &(vx[igraph][0]), &(vy[igraph][0]),
                                    &(vexl[igraph][0]), &(vexh[igraph][0]), &(veyl[igraph][0]), &(veyh[igraph][0]));
    graph[igraph].SetMarkerStyle(styles[igraph]); graph[igraph].SetMarkerSize(1.4);
    if(leglabels[igraph]=="All bkg.") mcolors[igraph] = 1;
    if(leglabels[igraph]=="t#bar{t} (2l)") mcolors[igraph] = colors("tt_2l");
    if(leglabels[igraph]=="t#bar{t} (1l)") mcolors[igraph] = colors("tt_1l");
    if(leglabels[igraph]=="Other") mcolors[igraph] = colors("wjets");
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
  cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}");
  cmslabel.SetTextAlign(31);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.015,"#font[42]{13 TeV}");

  line.SetLineStyle(3); line.SetLineWidth(1);
  line.DrawLine(minx, 1, maxx, 1);

  TString fname="plots/ratio_"+cutsToPlain(ytitle)+"_"+plotdef.name+"_"+cutsToPlain(plotdef.baseline)+".pdf";
  can.SaveAs(fname);
  cout<<endl<<" open "<<fname<<endl;

} // plotRatio


void printDebug(vector<vector<TString> > &allcuts, vector<vector<vector<GammaParams> > > &allyields, 
		TString baseline);

void GetOptions(int argc, char *argv[]);

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

  string ntupletag="nleps1met200nj5";
  string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_met150_and_nleps1met200nj5/");
  if(skim.Contains("met150")) ntupletag="_met150";
  if(skim.Contains("both")) ntupletag="";

  Palette colors("txt/colors.txt", "default");

  // Cuts in baseline speed up the yield finding
  string baseline = "pass && stitch && mj14>250 && nleps==1 && nveto==0 && ht>500 && met>150 && njets>=5 && nbm>=1";

  // // For syncing
  // auto proc_tt1l = Process::MakeShared<Baby_full>("tt 1lep", Process::Type::background, colors("tt_1l"),
  //   {foldermc+"*_TTJets*SingleLept*"+ntupletag+"*.root"},
  //   baseline);
  // auto proc_tt2l = Process::MakeShared<Baby_full>("tt 2lep", Process::Type::background, colors("tt_2l"),
  //   {foldermc+"*_TTJets*DiLept*"+ntupletag+"*.root", foldermc+"*_TTJets_HT*"+ntupletag+"*.root"},
  //   baseline);

  auto proc_tt1l = Process::MakeShared<Baby_full>("tt 1lep", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets*SingleLept*"+ntupletag+"*.root", foldermc+"*_TTJets_HT*"+ntupletag+"*.root"},
    baseline+" && stitch && ntruleps==1");
  auto proc_tt2l = Process::MakeShared<Baby_full>("tt 2lep", Process::Type::background, colors("tt_2l"),
    {foldermc+"*_TTJets*DiLept*"+ntupletag+"*.root", foldermc+"*_TTJets_HT*"+ntupletag+"*.root"},
    baseline+" && stitch && ntruleps==2");

  auto proc_other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {foldermc+"*_WJetsToLNu*"+ntupletag+"*.root",foldermc+"*_ST_*"+ntupletag+"*.root",
        foldermc+"*_TTW*"+ntupletag+"*.root",foldermc+"*_TTZ*"+ntupletag+"*.root",
        foldermc+"*DYJetsToLL*"+ntupletag+"*.root",foldermc+"*QCD_HT*"+ntupletag+"*.root",
        foldermc+"*_ZJet*"+ntupletag+"*.root",foldermc+"*_ttHJetTobb*"+ntupletag+"*.root",
        foldermc+"*_TTGJets*"+ntupletag+"*.root",foldermc+"*_TTTT*"+ntupletag+"*.root",
        foldermc+"*_WH_HToBB*"+ntupletag+"*.root",foldermc+"*_ZH_HToBB*"+ntupletag+"*.root",
        foldermc+"*_WWTo*"+ntupletag+"*.root",foldermc+"*_WZ*"+ntupletag+"*.root",foldermc+"*_ZZ_*"+ntupletag+"*.root"},
    baseline+" && stitch");

  // //// Use these processes to make quick plots
  // auto proc_tt1l = Process::MakeShared<Baby_full>("tt 1lep", Process::Type::background, colors("tt_1l"),
  //   {foldermc+"*_TTJets_Tune*.root"},
  //   baseline+" && ntruleps==1");
  // auto proc_tt2l = Process::MakeShared<Baby_full>("tt 2lep", Process::Type::background, colors("tt_2l"),
  //   {foldermc+"*_TTJets_Tune*.root"},
  //   baseline+" && ntruleps==2");

  // auto proc_other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
  //   {foldermc+"*_WWTo*"+ntupletag+"*.root",foldermc+"*_WZ*"+ntupletag+"*.root"},
  //   baseline+" && stitch");

  vector<shared_ptr<Process> > all_procs = {proc_tt1l, proc_tt2l, proc_other};

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////// Defining cuts ///////////////////////////////////////////////
  // baseline defined above

  // Makes a plot for each vector in plotcuts
  vector<oneplot> plotcuts({{"njets", "met>200", {"njets==5", "njets==6", "njets==7", "njets==8", "njets>=9"}},
	{"met", "njets>=6", {"met>150 && met<=200","met>200 && met<=350","met>350 && met<=500","met>500"}}});
  vector<TString> abcdcuts = {"mt<=140 && mj14<=400", "mt<=140 && mj14>400",
			      "mt>140  && mj14<=400", "mt>140  && mj14>400"};
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
      pm.Push<Table>(tname.Data(),  table_cuts, all_procs, false, false);
    } // Loop over abcdcuts
  } // Loop over plots

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////// Finding all yields ///////////////////////////////////////////////

  pm.min_print_ = true;
  pm.MakePlots(lumi);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////// Calculating preds/kappas and printing table //////////////////////////////////////

  for(size_t iplot=0; iplot<plotcuts.size(); iplot++){
    // allyields: [0] All bkg, [1] tt1l, [2] tt2l, [3] other
    vector<vector<vector<GammaParams> > > allyields(4,vector<vector<GammaParams> >(Nabcd));
    for(size_t iabcd=0; iabcd<abcdcuts.size(); iabcd++){
      Table * yield_table = static_cast<Table*>(pm.Figures()[iplot*Nabcd+iabcd].get());
      allyields[bkg]  [iabcd] = yield_table->BackgroundYield(lumi);
      allyields[tt1l] [iabcd] = yield_table->Yield(proc_tt1l.get(), lumi);
      allyields[tt2l] [iabcd] = yield_table->Yield(proc_tt2l.get(), lumi);
      allyields[other][iabcd] = yield_table->Yield(proc_other.get(), lumi);
    } // Loop over ABCD cuts

    //// Print MC/Data yields, cuts applied, kappas, preds
    if(debug) printDebug(allcuts[iplot], allyields, baseline);

    //// RmT
    vector<vector<vector<int> > > indices({
	{{tt2l, r3, 1}, {tt1l, r1, -1}},
	  {{tt2l, r4, 1}, {tt1l, r2, -1}},
	    {{bkg, r3, 1}, {bkg, r1, -1}},
	      {{bkg, r4, 1}, {bkg, r2, -1}}
      });
    vector<TString> leglabels({"tt: M_{J}<400", "tt: M_{J}>400", "bkg: M_{J}<400", "bkg: M_{J}>400"});
    plotRatio(allyields, plotcuts[iplot], indices, leglabels);

    //// RMJ/RMJ
    indices = vector<vector<vector<int> > >({
	{{tt2l, r4, 1}, {tt2l, r3, -1}, {tt1l, r2, -1}, {tt1l, r1, 1}},
	  {{tt1l, r4, 1}, {tt1l, r3, -1}, {tt1l, r2, -1}, {tt1l, r1, 1}},
	    {{other, r4, 1}, {other, r3, -1}, {tt1l, r2, -1}, {tt1l, r1, 1}},
	      {{bkg, r4, 1}, {bkg, r3, -1}, {tt1l, r2, -1}, {tt1l, r1, 1}}
      });
    leglabels = vector<TString>({"t#bar{t} (2l)", "t#bar{t} (1l)", "Other", "All bkg."});
    plotRatio(allyields, plotcuts[iplot], indices, leglabels);


  } // Loop over plots


  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Finding "<<plotcuts.size()<<" plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
} // main

////////////////////////////////////////// End of main //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


// allyields: [0] All bkg, [1] tt1l, [2] tt2l, [3] other
void printDebug(vector<vector<TString> > &allcuts, vector<vector<vector<GammaParams> > > &allyields, 
		TString baseline){
  int digits = 3;
  cout<<endl<<endl<<"============================ Printing cuts  ============================"<<endl;
  cout<<"-- Baseline cuts: "<<baseline<<endl<<endl;
  for(size_t ibin=0; ibin<allcuts[0].size(); ibin++){
    for(size_t iabcd=0; iabcd<allcuts.size(); iabcd++){
      cout<<"MC: "    <<setw(9)<<RoundNumber(allyields[bkg]  [iabcd][ibin].Yield(), digits)
	  <<"  tt1l: "<<setw(9)<<RoundNumber(allyields[tt1l] [iabcd][ibin].Yield(), digits)
	  <<"  tt2l: "<<setw(9)<<RoundNumber(allyields[tt2l] [iabcd][ibin].Yield(), digits)
	  <<" other: "<<setw(9)<<RoundNumber(allyields[other][iabcd][ibin].Yield(), digits)
	  <<"  - "<< allcuts[iabcd][ibin]<<endl;
    } // Loop over ABCD cuts
    cout<<endl;
  } // Loop over bin cuts

} // printDebug


void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"lumi", required_argument, 0, 'l'},    // Luminosity to normalize MC with (no data)
      {"skim", required_argument, 0, 's'},    // Which skim to use: standard, met150, 2015 data
      {"debug", no_argument, 0, 'd'},         // Debug: prints yields and cuts used
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s:l:d", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'l':
      lumi = atof(optarg);
      break;
    case 's':
      skim = optarg;
      break;
    case 'd':
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
