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
#include "TStyle.h"

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
  float lumi=35.9;
  struct oneplot{
    NamedFunc den;
    vector<NamedFunc> nums;
  };
}

void GetOptions(int argc, char *argv[]);

void makePlot(TH2D hnum, TH2D hden);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining baby //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string folder(bfolder+"/cms2r0/babymaker/babies/2017_03_17/TChiHH/merged_higmc_unskimmed/");

  Baby_full baby(std::set<std::string>{(folder+"*ino-7*_*root")});
  auto activator = baby.Activate();


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////// Cuts and histograms
  enum Cuts             {all,    nj4,   goodreco,   gr_am};
  vector<TString> scuts({"all", "nj4", "goodreco", "gr_am"});

  //// Setting plot style
  PlotOpt opts("txt/plot_styles.txt", "Eff2D");
  setPlotStyle(opts);

  int nbins = 19;
  float minX = 50, maxX = 1000;//, binW = (maxX-minX)/static_cast<float>(nbins);
  vector<TH2D> histos;
  for(unsigned icut=0; icut<scuts.size(); icut++){
    histos.push_back(TH2D(scuts[icut], scuts[icut], nbins, minX, maxX, nbins, minX, maxX));
  } // Loop over cuts

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////// Loop over entries
  long nentries = baby.GetEntries();
  for(long entry = 0; entry < nentries; ++entry){
    baby.GetEntry(entry);
    if(entry%1000000==0 || entry==nentries-1) {
      cout<<"Doing "<<setw(7)<<entry<<"/"<<nentries<<"  -> "<<setw(5)<<RoundNumber(entry*100,1,nentries)<<"%"<<endl;
    }

    ///////////// Counting b quarks in acceptance and finding Higgs pt
    int nbacc = 0;
    float higpt1(0), higpt2(0);
    for (unsigned imc(0); imc<baby.mc_pt()->size(); imc++){
      if (abs(baby.mc_id()->at(imc))==25) { // It's a Higgs!
	if (baby.mc_pt()->at(imc)>higpt1) {
	  higpt2 = higpt1;
	  higpt1 = baby.mc_pt()->at(imc);
	} else if (baby.mc_pt()->at(imc)>higpt2) {
	  higpt2 = baby.mc_pt()->at(imc);
	}
      } // It's a Higgs!
      if (abs(baby.mc_id()->at(imc))==5 && abs(baby.mc_mom()->at(imc))==25)  // It's a b from a Higgs!
	if(baby.mc_pt()->at(imc)>30 && fabs(baby.mc_eta()->at(imc))<2.4) nbacc++;
    } // Loop over mc particles
    if(nbacc<4) continue;

    ///////////// Calculating other cuts
    bool pass_am = baby.higd_am()>100. && baby.higd_am()<=140.;
    bool pass_goodreco = baby.njets()>=4;
    for(unsigned ijet=0; ijet<baby.jets_pt()->size(); ijet++){
      if(baby.jets_h1d()->at(ijet) || baby.jets_h2d()->at(ijet))
	if(baby.jets_hflavor()->at(ijet)!=5) pass_goodreco = false;
    } // Loop over jets

    ///////////// Calculating yields
    float weight = 1.;
    histos[all].Fill(higpt1, higpt2, weight);
    if(baby.njets()>=4) histos[nj4].Fill(higpt1, higpt2, weight);
    if(pass_goodreco) {
      histos[goodreco].Fill(higpt1, higpt2, weight);
      if(pass_am) histos[gr_am].Fill(higpt1, higpt2, weight);
    } // goodreco
  } // Loop over entries

  cout<<endl;
  makePlot(histos[goodreco], histos[all]);
  makePlot(histos[gr_am], histos[goodreco]);

  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Finding plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
} // main

////////////////////////////////////////// End of main //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void makePlot(TH2D hnum, TH2D hden){
  gStyle->SetPaintTextFormat("5.2f");
  TCanvas can("can","");

  hnum.Divide(&hden);
  hnum.GetXaxis()->SetTitleSize(0.045);
  hnum.GetXaxis()->SetTitleOffset(1.3);
  hnum.GetXaxis()->SetLabelSize(0.04);
  hnum.GetXaxis()->SetLabelOffset(0.015);

  hnum.GetYaxis()->SetTitleSize(0.045);
  hnum.GetYaxis()->SetLabelSize(0.04);
  hnum.GetYaxis()->SetTitleOffset(1.22);

  hnum.GetZaxis()->SetLabelSize(0.04);
  hnum.GetZaxis()->SetRangeUser(0.,1.);
  hnum.GetZaxis()->SetTitleSize(0.045);
  hnum.GetZaxis()->CenterTitle(true);
  hnum.GetXaxis()->SetTitle("Leading Higgs boson p_{T} [GeV]");
  hnum.GetYaxis()->SetTitle("Subleading Higgs boson p_{T} [GeV]");
  hnum.GetZaxis()->SetTitle("Efficiency");
  hnum.Draw("COLZTEXT");

  TString title = "Den: "; title += hden.GetTitle();
  title += "   -   Num: "; title += hnum.GetTitle();
  hnum.SetTitle(title);

  int nbx = hnum.GetNbinsX(), nby = hnum.GetNbinsY();
  float minX = hnum.GetXaxis()->GetBinLowEdge(1), maxX = hnum.GetXaxis()->GetBinLowEdge(nbx+1);
  float minY = hnum.GetYaxis()->GetBinLowEdge(1), maxY = hnum.GetYaxis()->GetBinLowEdge(nby+1);
  TLine line; line.SetLineWidth(1); line.SetLineStyle(3); line.SetLineColor(kBlack);
  for(int binx=1; binx<nbx; binx++){
    line.DrawLine(hnum.GetXaxis()->GetBinLowEdge(binx+1), minY, hnum.GetXaxis()->GetBinLowEdge(binx+1),maxY);  
  } // Loop over x bins
  for(int biny=1; biny<nby; biny++){
    line.DrawLine(minX,  hnum.GetYaxis()->GetBinLowEdge(biny+1), maxX, hnum.GetYaxis()->GetBinLowEdge(biny+1));  
  } // Loop over y bins
  
  TString hname = "plots/eff_"; 
  hname += hnum.GetName(); hname += "_";
  hname += hden.GetName(); hname += ".pdf";
  can.SaveAs(hname);

  cout<<" open "<<hname<<endl;
}


void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"lumi", required_argument, 0, 'l'},    // Luminosity to normalize MC with (no data)
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "l:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'l':
      lumi = atof(optarg);
      break;
    case 0:
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
