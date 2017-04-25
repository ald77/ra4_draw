//// plot_ratios: plots rMJ and rmT, and combinations of these

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
#include "TStyle.h"
#include "TBox.h"
#include "TFile.h"

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/styles.hpp"
#include "core/plot_opt.hpp"
#include "hig/hig_functions.hpp"

using namespace std;

namespace{
  bool do_lo = false;
  bool quick = false;
  float lumi = 50.;
  PlotOpt opts("txt/plot_styles.txt", "Eff2D");
}

void GetOptions(int argc, char *argv[]);

void makePlot(TH2D hnum, TH2D hden, TString num_s="", TString den_s="");

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=3000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining baby //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string folder(bfolder+"/cms2r0/babymaker/babies/2017_04_01/mc/merged_bare_dy/");
  if (do_lo) folder = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_bare_dy/";
  string babyname("*DYJetsToLL_*.root");
  if (do_lo) babyname = "*DYJetsToLL_M-50_*.root";
  if (quick) babyname = "*DYJetsToLL_M-50_Tu*.root";
  Baby_full baby(std::set<std::string>{folder+babyname});
  auto activator = baby.Activate();


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  //// Setting plot style
  setPlotStyle(opts);
  gStyle->SetPaintTextFormat("5.2f");

  int ht_nbins(20); float ht_min(0), ht_max(2000);
  int nisr_nbins(4); float nisr_min(-0.5), nisr_max(3.5);
  int njets_nbins(11); float njets_min(-0.5), njets_max(10.5);
  int pt_nbins(20); float pt_min(0), pt_max(1000);

  TString labels ="; Gen H_{T} [GeV]; p_{T}(ll) [GeV]";
  TH2D ht_pt_nj6 = TH2D("ht_pt_nj6", labels, ht_nbins, ht_min, ht_max, pt_nbins, pt_min, pt_max);
  TH2D ht_pt_njl6 = TH2D("ht_pt_njl6", labels, ht_nbins, ht_min, ht_max, pt_nbins, pt_min, pt_max);
  TH2D ht_pt_all = TH2D("ht_pt_all", labels, ht_nbins, ht_min, ht_max, pt_nbins, pt_min, pt_max);
  labels ="; Gen H_{T} [GeV]; ME parton multiplicity";
  TH2D ht_nisr_pt300 = TH2D("ht_nisr_pt300",labels, ht_nbins, ht_min, ht_max, nisr_nbins, nisr_min, nisr_max);
  TH2D ht_nisr_ptl300_ht300 = TH2D("ht_nisr_ptl300_ht300",labels, ht_nbins, ht_min, ht_max, nisr_nbins, nisr_min, nisr_max);
  TH2D ht_nisr_all = TH2D("ht_nisr_all",labels, ht_nbins, ht_min, ht_max, nisr_nbins, nisr_min, nisr_max);
  labels ="; Gen p_{T} (ll) [GeV]; ME parton multiplicity";
  TH2D pt_nisr_ht600 = TH2D("pt_nisr_ht600", labels, pt_nbins, pt_min, pt_max, nisr_nbins, nisr_min, nisr_max);
  TH2D pt_nisr_htl600 = TH2D("pt_nisr_htl600", labels, pt_nbins, pt_min, pt_max, nisr_nbins, nisr_min, nisr_max);
  TH2D pt_nisr_all = TH2D("pt_nisr_all", labels, pt_nbins, pt_min, pt_max, nisr_nbins, nisr_min, nisr_max);
  labels ="; Gen H_{T} [GeV]; Reco. N_{jets}";
  TH2D ht_njets_pt300 = TH2D("ht_njets_pt300",labels, ht_nbins, ht_min, ht_max, njets_nbins, njets_min, njets_max);
  TH2D ht_njets_ptl300_ht300 = TH2D("ht_njets_ptl300_ht300",labels, ht_nbins, ht_min, ht_max, njets_nbins, njets_min, njets_max);
  TH2D ht_njets_all = TH2D("ht_njets_all",labels, ht_nbins, ht_min, ht_max, njets_nbins, njets_min, njets_max);
  labels ="; Gen p_{T} (ll) [GeV]; Reco. N_{jets}";
  TH2D pt_njets_ht600 = TH2D("pt_njets_ht600", labels, pt_nbins, pt_min, pt_max, njets_nbins, njets_min, njets_max);
  TH2D pt_njets_htl600 = TH2D("pt_njets_htl600", labels, pt_nbins, pt_min, pt_max, njets_nbins, njets_min, njets_max);
  TH2D pt_njets_all = TH2D("pt_njets_all", labels, pt_nbins, pt_min, pt_max, njets_nbins, njets_min, njets_max);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////// Loop over entries
  long nentries = baby.GetEntries();
  for(long entry = 0; entry < nentries; ++entry){
    baby.GetEntry(entry);
    if((entry+1)%1000000==0 || entry==nentries-1) {
      cout<<"Doing "<<setw(9)<<AddCommas(entry+1)<<"/"<<AddCommas(nentries)<<"  -> "
          <<setw(5)<<RoundNumber((entry+1)*100,1,nentries)<<"%"<<endl;
    }

    float weight = baby.w_lumi()*lumi;

    if (!quick) {
      if (do_lo && !baby.stitch()) continue;
      else if (baby.type()==6200 && baby.ptll_me()>100) continue;
    }
    
    ht_pt_all.Fill(baby.ht_isr_me(), baby.ptll_me(), weight);
    if (baby.njets()>=6) ht_pt_nj6.Fill(baby.ht_isr_me(), baby.ptll_me(), weight);
    else ht_pt_njl6.Fill(baby.ht_isr_me(), baby.ptll_me(), weight);

    ht_nisr_all.Fill(baby.ht_isr_me(), baby.nisr_me(), weight);
    if (baby.ptll_me()>300) ht_nisr_pt300.Fill(baby.ht_isr_me(), baby.nisr_me(), weight);
    else if (baby.ht_isr_me()>300) ht_nisr_ptl300_ht300.Fill(baby.ht_isr_me(), baby.nisr_me(), weight);

    pt_nisr_all.Fill(baby.ptll_me(), baby.nisr_me(), weight);
    if (baby.ht_isr_me()>600) pt_nisr_ht600.Fill(baby.ptll_me(), baby.nisr_me(), weight);
    // else pt_nisr_htl600.Fill(baby.ptll_me(), baby.nisr_me(), weight);

    ht_njets_all.Fill(baby.ht_isr_me(), baby.njets(), weight);
    if (baby.ptll_me()>300) ht_njets_pt300.Fill(baby.ht_isr_me(), baby.njets(), weight);
    else if (baby.ht_isr_me()>300) ht_njets_ptl300_ht300.Fill(baby.ht_isr_me(), baby.njets(), weight);

    pt_njets_all.Fill(baby.ptll_me(), baby.njets(), weight);
    if (baby.ht_isr_me()>600) pt_njets_ht600.Fill(baby.ptll_me(), baby.njets(), weight);
    else pt_njets_htl600.Fill(baby.ptll_me(), baby.njets(), weight);


  } // Loop over entries
  cout<<endl;

  TH2D hdummy = TH2D("hdummy","hdummy",1,0,1,1,0,1);
  makePlot(ht_pt_nj6, ht_pt_all, "N_{jets} #geq 6", "All events");
  makePlot(ht_nisr_pt300, ht_nisr_all, "p_{T}(ll) > 300 GeV", "All events");
  makePlot(pt_nisr_ht600, pt_nisr_all, "H_{T} > 600 GeV", "All events");
  makePlot(ht_njets_pt300, ht_njets_all, "p_{T}(ll) > 300 GeV", "All events");
  makePlot(pt_njets_ht600, pt_njets_all, "H_{T} > 600 GeV", "All events");
  // 2D distributions
  makePlot(ht_pt_nj6, hdummy, "N_{jets} #geq 6", "");
  makePlot(ht_nisr_pt300, hdummy, "p_{T}(ll) > 300 GeV", "");
  makePlot(pt_nisr_ht600, hdummy, "H_{T} > 600 GeV", "");
  makePlot(ht_njets_pt300, hdummy, "p_{T}(ll) > 300 GeV", "");
  makePlot(pt_njets_ht600, hdummy, "H_{T} > 600 GeV", "");
  //reverse cuts
  makePlot(ht_nisr_ptl300_ht300, hdummy, "p_{T}(ll) < 300 GeV", "");
  makePlot(ht_njets_ptl300_ht300, hdummy, "p_{T}(ll) < 300 GeV", "");

  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Finding plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
} // main

////////////////////////////////////////// End of main //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void makePlot(TH2D hnum, TH2D hden, TString num_s, TString den_s){

  bool dummy = string(hden.GetName())=="hdummy";
  TCanvas can("can","");

  ////// Setting additional histogram style
  hnum.GetXaxis()->SetLabelOffset(0.012);

  hnum.GetZaxis()->SetLabelSize(0.04);
  if (!dummy) {
    hnum.GetZaxis()->SetRangeUser(0.,1.);
    hnum.GetZaxis()->SetTitle("Ratio");
  } else {
    hnum.GetZaxis()->SetTitle("Fraction of events");
    // can.SetLogz();
  }
  hnum.GetZaxis()->SetTitleSize(0.045);
  hnum.GetZaxis()->CenterTitle(true);
  

  ////// Dividing and plot histogram
  if (!dummy) hnum.Divide(&hden);
  else hnum.Scale(1./hnum.Integral(0,hnum.GetNbinsX()+1, 0,hnum.GetNbinsY()+1));
  hnum.Draw("COLZTEXT");


  ////// Drawing line grid
  int nbx = hnum.GetNbinsX(), nby = hnum.GetNbinsY();
  float Xmin = hnum.GetXaxis()->GetBinLowEdge(1), Xmax = hnum.GetXaxis()->GetBinLowEdge(nbx+1);
  float Ymin = hnum.GetYaxis()->GetBinLowEdge(1), Ymax = hnum.GetYaxis()->GetBinLowEdge(nby+1);
  TLine line; line.SetLineWidth(1); line.SetLineStyle(3); line.SetLineColor(kBlack);
  for(int binx=1; binx<nbx; binx++){
    line.DrawLine(hnum.GetXaxis()->GetBinLowEdge(binx+1), Ymin, hnum.GetXaxis()->GetBinLowEdge(binx+1),Ymax);  
  } // Loop over x bins
  for(int biny=1; biny<nby; biny++){
    line.DrawLine(Xmin,  hnum.GetYaxis()->GetBinLowEdge(biny+1), Xmax, hnum.GetYaxis()->GetBinLowEdge(biny+1));  
  } // Loop over y bins
  
  float lMargin = opts.LeftMargin(), rMargin = opts.RightMargin(), tMargin = opts.TopMargin();
  TLatex label; label.SetTextFont(42); label.SetTextSize(0.05);
  label.SetTextAlign(21);

  if (dummy) {
    TString title =  "Selection: " + num_s;
    label.DrawLatexNDC(lMargin+(1-rMargin-lMargin)/2., 1-tMargin+0.025, title);
  } else {
    //// Printing histo title
    TString title = "#font[72]{Den}: "; 
    if(den_s!="") title += den_s;
    else title += hden.GetTitle();
    title += " - #font[72]{Num:} ";
    if(num_s!="") title += num_s;
    else title += hnum.GetTitle();
    label.DrawLatexNDC(lMargin+(1-rMargin-lMargin)/2., 1-tMargin+0.025, title);
  }
  hnum.SetTitle("");

  ////// Saving plot
  TString hname = "plots/ratio_"; 
  if (dummy) hname = "plots/corr_";
  hname += hnum.GetName(); 
  hname.ReplaceAll(".", "p");
  hname += ".pdf";
  can.SaveAs(hname);
  cout<<" open "<<hname<<endl;
}


void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"quick", no_argument, 0, 'q'},      // Plot for high dR between b daughters
      {"lo", required_argument, 0, 0},     // Minimum pT for quark to be in acceptance
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "q", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'q':
      quick = true;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "lo"){
        do_lo = true;
      } else{
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
