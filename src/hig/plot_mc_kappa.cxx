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
  // options "zll", "qcd", "ttbar", "search"
  string sample = "search";
  float lumi=40.;
  bool do_trim = true;

  enum Regions {sbd2b, hig2b, sbd3b, hig3b, sbd4b, hig4b};
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

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/";
  if (sample=="ttbar") foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higlep1/";
  if (sample=="zll") foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_zisrnjet45/";
  if (sample=="qcd") foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higqcd/";

  set<string> alltags = {"*_TTJets*Lept*.root", "*_TTJets_HT*.root", 
            "*_TTZ*.root", "*_TTW*.root", "*_TTGJets*.root", "*_ttHJetTobb*.root","*_TTTT*.root",
            "*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root",
            "*_ST_*.root",
            "*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root",
            "*_WH_HToBB*.root", "*_ZH_HToBB*.root", "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"};

  set<string> ttxtags = {"*_TTJets*Lept*.root", "*_TTJets_HT*.root", 
            "*_TTZ*.root", "*_TTW*.root", "*_TTGJets*.root", "*_ttHJetTobb*.root","*_TTTT*.root"};

  // Baseline definitions
  string baseline="pass && stitch";
  NamedFunc base_func("pass && stitch");
  if (do_trim) base_func = base_func && "hig_dm<40 && hig_am<200";

  vector<shared_ptr<Process> > procs;
  // All bkg. kappa
  procs.push_back(Process::MakeShared<Baby_full>("All bkg.", Process::Type::background, 
    kBlack, attach_folder(foldermc,alltags), base_func));
  // Sample specific kappa
  if (sample=="zll") 
    procs.push_back(Process::MakeShared<Baby_full>("Z#rightarrow ll", Process::Type::background, 
      kOrange+1, {foldermc+"*DYJetsToLL*.root"}, base_func));
  if (sample=="qcd") 
    procs.push_back(Process::MakeShared<Baby_full>("QCD", Process::Type::background, 
      colors("other"),{foldermc+"*QCD_HT*0_Tune*.root", foldermc+"*QCD_HT*Inf_Tune*.root"}, base_func)); 
  if (sample=="ttbar" || sample=="search") 
    procs.push_back(Process::MakeShared<Baby_full>("t#bar{t}+X", Process::Type::background,
      colors("tt_1l"),attach_folder(foldermc,ttxtags), base_func));

 /////////////////////////////////////////////////////////////////////////////////////////////////////////
 /////////////////////////////////////////// Defining cuts ///////////////////////////////////////////////
  vector<TString> xcuts; // all desired cut combinations
  // zll skim:  ((elelv_m>80&&elelv_m<100)||(mumuv_m>80&&mumuv_m<100))
  // nvleps==2 && nleps>=1 && Max$(leps_pt)>30 && njets>=4&&njets<=5
  if (sample=="zll") {
    xcuts.push_back("nvleps==2 && met<50");
    xcuts.push_back("nvleps==2 && met<50 && hig_drmax<2.2");
  }
  // qcd skim - met>150 && nvleps==0 && (njets==4||njets==5)
  if (sample=="qcd") {
    xcuts.push_back("nvleps==0 && ntks==0 && low_dphi");
    xcuts.push_back("nvleps==0 && ntks==0 && low_dphi && hig_drmax<2.2");
  }
  // ttbar skim - met>100 && nleps==1 && (njets==4||njets==5) && nbm>=2
  if (sample=="ttbar") {
    xcuts.push_back("nleps==1 && mt<100");
    xcuts.push_back("nleps==1 && mt<100 && hig_drmax<2.2");
    xcuts.push_back("nleps==1 && mt<100 && hig_drmax<2.2 && ntks==0 && !low_dphi");
  } 
  // search skim - met>100 && nvleps==0 && (njets==4||njets==5) && nbm>=2
  if (sample=="search") {
    xcuts.push_back("nvleps==0");
    xcuts.push_back("nvleps==0 && ntks==0 && !low_dphi");
    xcuts.push_back("nvleps==0 && ntks==0 && !low_dphi && hig_drmax<2.2");
  }

  vector<TString> metcuts;
  string metdef = "met";
  if (sample=="zll") metdef = "(mumuv_pt*(mumuv_pt>0)+elelv_pt*(elelv_pt>0))";
  // if (sample!="qcd") metcuts.push_back(metdef+">100&&"+metdef+"<=150");
  metcuts.push_back(metdef+">150&&"+metdef+"<=200");
  metcuts.push_back(metdef+">200&&"+metdef+"<=300");
  metcuts.push_back(metdef+">300");

  // Makes a plot for each vector in plotcuts
  vector<oneplot> plotcuts;
  for (auto &ixcut: xcuts) plotcuts.push_back({"met",ixcut,metcuts});

  TString c_ab="nbt==2&&nbm==2";
  TString c_cd="nbt>=2&&nbm==3&&nbl==3";
  TString c_ef="nbt>=2&&nbm>=3&&nbl>=4";
  if (sample=="qcd" || sample=="zll"){
    c_ab="nbm==0";
    c_cd="nbm==1";
    c_ef="nbt==1";  
  }

  TString c_hig="hig_am>100&&hig_am<140&&hig_dm<40";
  TString c_sbd="!("+c_hig+") && hig_am<200 && hig_dm<40";

  vector<TString> abcdcuts = {c_ab+"&&"+c_sbd, c_ab+"&&"+c_hig, 
                              c_cd+"&&"+c_sbd, c_cd+"&&"+c_hig, 
                              c_ef+"&&"+c_sbd, c_ef+"&&"+c_hig};
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
    for(int ibkg=0; ibkg<static_cast<int>(procs.size()); ibkg++) {
      indices.push_back(vector<vector<int> >({{ibkg, hig3b, 1},  {ibkg, sbd3b, -1}, 
                                              {ibkg, hig2b, -1}, {ibkg, sbd2b, 1}}));
      leglabels.push_back(procs[ibkg]->name_);
    }
    plotRatio(allyields, plotcuts[iplot], indices, leglabels, procs);

    //// 4b kappa
    indices.clear(); leglabels.clear();
    for(int ibkg=0; ibkg<static_cast<int>(procs.size()); ibkg++) {
      indices.push_back(vector<vector<int> >({{ibkg, hig4b, 1},  {ibkg, sbd4b, -1}, 
                                              {ibkg, hig2b, -1}, {ibkg, sbd2b, 1}}));
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
      if (sample=="qcd" || sample=="zll") ytitle = "1 tight b #kappa";
      else ytitle = "4b #kappa";
    }
    if((ind0==hig3b&&ind1==sbd3b && ind2==hig2b&&ind3==sbd2b)) {
      if (sample=="qcd" || sample=="zll") ytitle = "1b #kappa";
      else ytitle = "3b #kappa";
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
  for(size_t ibkg=0; ibkg<procs.size(); ibkg++) {
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

  TString fname = "plots/ratio_"+sample+"_"+CodeToPlainText(ytitle.Data())+"_"+plotdef.name+"_"
    +CodeToPlainText(plotdef.baseline.Data())+".pdf";
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
      {"lumi", required_argument, 0, 'l'},    // Luminosity to normalize MC with (no data)
      {"sample", required_argument, 0, 's'},    // Which sample to use: standard, met150, 2015 data
      {"debug", no_argument, 0, 'g'},         // Debug: prints yields and cuts used
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s:l:gn:d:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'l':
      lumi = atof(optarg);
      break;
    case 's':
      sample = optarg;
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
