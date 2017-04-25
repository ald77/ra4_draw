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
  bool quick = false;
  bool do_hflavor = false;
  bool do_recopt = false;
  bool do_csv = true;
  bool do_deep = true;

  enum Cuts             {loose,   medium,   tight };
  const vector<TString> op_names({"loose", "medium", "tight"});
  const vector<double> op_csv({0.5426, 0.8484, 0.9535});
  const vector<double> op_deep({0.2219, 0.6324, 0.8958});
} 

double sfFullCSV(unsigned op, unsigned jetflav, double pt);
double sfFullDeep(unsigned op, unsigned jetflav, double pt);
TLegend prettylegend(float x0,float y0,float x1,float y1);

void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  gStyle->SetOptStat(0);
  
  GetOptions(argc, argv);

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining baby //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string folder(bfolder+"/cms2r0/babymaker/babies/2017_03_07/mc/merged_btagmc_unskimmed/");
  std::set<std::string> samples = {folder+"*TTJets*"};
  if (quick) samples = {folder+"*T1bbbb_mGluino-1500*"};
  Baby_full baby(samples);
  auto activator = baby.Activate();


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////// Cuts and histograms
  const vector<double> bins = {30, 35, 
    40, 50, 60, 70, 80, 90, 
    100, 125, 150, 175, 200, 225, 250, 275, 
    300, 350, 400, 450, 500, 600, 700, 800, 1000};
  unsigned nbins = bins.size()-1;
  TH1D histos_all = TH1D("denom","denom",nbins, &(bins[0]));
  vector<TH1D> histos_csv, histos_deep;
  for(unsigned icut=0; icut<op_names.size(); icut++){
    histos_csv.push_back(TH1D("CSVv2_algorithm_"+op_names[icut], "CSVv2_algorithm_"+op_names[icut], nbins, &(bins[0])));
    histos_deep.push_back(TH1D("deepCSV_algorithm_"+op_names[icut], "deepCSV_algorithm_"+op_names[icut], nbins, &(bins[0])));
  } // Loop over cuts

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////// Loop over entries
  long nentries = baby.GetEntries();
  for(long entry = 0; entry < nentries; ++entry){
    baby.GetEntry(entry);
    if((entry+1)%1000000==0 || entry==nentries-1) {
      cout<<"Doing "<<setw(9)<<AddCommas(entry+1)<<"/"<<AddCommas(nentries)<<"  -> "
          <<setw(5)<<RoundNumber((entry+1)*100,1,nentries)<<"%"<<endl;
    }

    for (unsigned ijet(0); ijet<baby.jets_pt()->size(); ijet++){
      double pt = baby.jets_pt()->at(ijet);
      if (pt<30 || fabs(baby.jets_eta()->at(ijet))>2.4 || baby.jets_islep()->at(ijet)) continue;

      unsigned flav = floor(fabs(baby.jets_pflavor()->at(ijet))+0.1);
      if (do_hflavor) flav = floor(fabs(baby.jets_hflavor()->at(ijet))+0.1);

      if (flav!=5) continue;

      double genpt = baby.jets_gen_pt()->at(ijet);
      if (do_recopt) genpt = pt;
      if (genpt<0) continue;
      double csv = baby.jets_csv()->at(ijet);
      double csvd = baby.jets_csvd()->at(ijet);

      histos_all.Fill(genpt);
      for (unsigned iop(0); iop<3; iop++){
        if (csv > op_csv[iop])
          histos_csv[iop].Fill(genpt, sfFullCSV(iop, flav, pt));
        if (csvd > op_deep[iop])
          histos_deep[iop].Fill(genpt, sfFullDeep(iop, flav, pt));
      }
    }
  } // Loop over entries

  cout<<endl;
  
  TCanvas can;
  vector<TLegend> legs = {prettylegend(0.15,0.855,0.5, .995), prettylegend(0.5,0.855,0.895,.995)};
  vector<int> colors = {kGreen+2, kAzure+2, kRed+1};
  TH1D hdummy = TH1D("dummy", "", nbins, &(bins[0]));
  TPad grid("grid", "",0,0,1,1);
  grid.Draw();
  grid.cd();
  grid.SetGrid();
  grid.SetFillStyle(4000);
  hdummy.Draw();
  hdummy.GetYaxis()->SetNdivisions(15);
  hdummy.GetYaxis()->SetRangeUser(0.,1.);
  hdummy.GetYaxis()->SetTitle("b-tagging efficiency");
  hdummy.GetXaxis()->SetTitle("Particle-level jet p_{T} [GeV]");
  if (do_recopt) hdummy.GetXaxis()->SetTitle("Reconstructed jet p_{T} [GeV]");
  hdummy.GetXaxis()->SetTitleOffset(1.2);
  grid.SetMargin(.1,.05,.1,.15);

  for (unsigned iop(0); iop<3; iop++){
    if (do_csv) {
      histos_csv[iop].Divide(&histos_all);
      histos_csv[iop].SetLineColor(colors[iop]);
      histos_csv[iop].SetLineWidth(2);
      histos_csv[iop].SetMarkerColor(colors[iop]);
      histos_csv[iop].SetMarkerStyle(20+iop);
      histos_csv[iop].SetMarkerSize(0.8);

      histos_csv[iop].Draw("same");
      legs[0].AddEntry(&histos_csv[iop],"CSVv2, "+op_names[iop],"LP");
    }
    if (do_deep) {
      histos_deep[iop].Divide(&histos_all);
      histos_deep[iop].SetLineColor(colors[iop]);
      histos_deep[iop].SetLineStyle(2);
      histos_deep[iop].SetLineWidth(2);
      histos_deep[iop].SetMarkerColor(colors[iop]);
      histos_deep[iop].SetMarkerStyle(24+iop);
      histos_deep[iop].SetMarkerSize(0.8);

      histos_deep[iop].Draw("same");
      legs[1].AddEntry(&histos_deep[iop],"DeepCSV, "+op_names[iop],"LP");
    }
  }
  
  if (do_csv) legs[0].Draw();
  if (do_deep) legs[1].Draw();
  TString filename = "plots/btag_eff";
  if (do_recopt) filename += "_recopt";
  if (do_hflavor) filename += "_hflav";
  if (do_csv) filename += "_CSVv2";
  if (do_deep) filename += "_DeepCSV";
  filename += ".pdf";
  can.Print(filename);

  filename.ReplaceAll(".pdf", ".root");
  TFile file(filename, "recreate");
  file.cd();
  if (do_csv)
    for (auto &h: histos_csv) h.Write();
  if (do_deep)
    for (auto &h: histos_deep) h.Write();
  file.Close();

  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Finding plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
  cout<<endl<<"open "<<filename.ReplaceAll("root","pdf")<<endl<<endl;
} // main

double sfFullCSV(unsigned op, unsigned jetflav, double pt) {
  if (pt<20) pt = 20.;
  else if (pt>1000) pt = 1000.; 
  if (op==loose){
    if (jetflav >=4) return 0.887973*((1.+(0.0523821*pt))/(1.+(0.0460876*pt)));
    else             return 1.13904+-0.000594946*pt+1.97303e-06*pt*pt+-1.38194e-09*pt*pt*pt;
  } else if (op==medium){            
    if (jetflav >=4) return 0.561694*((1.+(0.31439*pt))/(1.+(0.17756*pt)));
    else             return 1.0589+0.000382569*pt+-2.4252e-07*pt*pt+2.20966e-10*pt*pt*pt;
  } else {
    if (jetflav >=4) return 0.817647*((1.+(0.038703*pt))/(1.+(0.0312388*pt)));
    else             return 0.971945+163.215/(pt*pt)+0.000517836*pt;
  } 
}

double sfFullDeep(unsigned op, unsigned jetflav, double pt) {
  if (pt<20) pt = 20.;
  else if (pt>1000) pt = 1000. ;
  if (op==loose){
    if (jetflav >=4) return 0.733112*((1.+(0.336449*pt))/(1.+(0.246914*pt)));
    else             return 1.06765+0.000317422*pt+-4.61732e-07*pt*pt+2.03608e-10*pt*pt*pt;
  } else if (op==medium){            
    if (jetflav >=4) return 0.637301*((1.+(0.479205*pt))/(1.+(0.311514*pt)));
    else             return 1.03216+0.000504744*pt+2.12276e-08*pt*pt+-2.27663e-10*pt*pt*pt;
  } else {
    if (jetflav >=4) return 0.506673*((1.+(0.464958*pt))/(1.+(0.239689*pt)));
    else             return 1.00762+51.6984/(pt*pt)+0.000370519*pt;
  } 
}

TLegend prettylegend(float x0,float y0,float x1,float y1){
  TLegend leg = TLegend(x0,y0,x1,y1);
  leg.SetTextFont(42);
  leg.SetTextSize(0.04);
  leg.SetFillColor(kWhite);
  leg.SetFillStyle(0);
  leg.SetLineColor(0);
  leg.SetBorderSize(0);
  return leg;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"csv", no_argument, 0, 'c'},
      {"deep", no_argument, 0, 'd'},
      {"quick", no_argument, 0, 'q'},
      {"hflav", no_argument, 0, 0},
      {"reco", no_argument, 0, 0},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "cdq", long_options, &option_index);

    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'c':
      do_csv = true;
      do_deep = false;
      break;
    case 'd':
      do_csv = false;
      do_deep = true;
      break;
    case 'q':
      quick = true;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "hflav"){
        do_hflavor = true;
      }else if(optname == "reco"){
        do_recopt = true;
      }else{
        printf("Bad option! Found option name %s\n", optname.c_str());
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}