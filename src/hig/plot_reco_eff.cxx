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
  float minpt = 30;
  float maxeta = 2.4;

  bool do_aux = false;
  bool low_isr = false;
  bool low_pu = false;
  bool hi_dr = false;
  TString mass = "";
  PlotOpt opts("txt/plot_styles.txt", "Eff2D");
}

void GetOptions(int argc, char *argv[]);

void makePlot(TH2D hnum, TH2D hden, TString num_s="", TString den_s="", TString extra_s="");

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

  Baby_full baby(std::set<std::string>{(folder+"*ino-"+mass.Data()+"*_*root")});
  auto activator = baby.Activate();


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////// Cuts and histograms
  enum Cuts             {all,    acc,   nj4,   goodreco,   gr_am};
  vector<TString> scuts({"all", "acc", "nj4", "goodreco", "gr_am"});

  //// Setting plot style
  setPlotStyle(opts);
  gStyle->SetPaintTextFormat("5.2f");

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
    if((entry+1)%1000000==0 || entry==nentries-1) {
      cout<<"Doing "<<setw(9)<<AddCommas(entry+1)<<"/"<<AddCommas(nentries)<<"  -> "
          <<setw(5)<<RoundNumber((entry+1)*100,1,nentries)<<"%"<<endl;
    }

    float weight = 1.;

    ///////////// Cuts on ISR and PU
    if(low_isr && baby.isr_tru_pt() > 20) continue;
    if(low_pu && baby.npv() > 10) continue;

    ///////////// Counting b quarks in acceptance and finding Higgs pt
    int nbacc = 0;
    float higpt1(0), higpt2(0);
    float dr1=-99, dr2=-99, preveta1=0, prevphi1=0, preveta2=0, prevphi2=0, momidx1=-99, momidx2=-99;
    for (unsigned imc(0); imc<baby.mc_pt()->size(); imc++){
      if (abs(baby.mc_id()->at(imc))==25) { // It's a Higgs!
        if (baby.mc_pt()->at(imc)>higpt1) {
          higpt2 = higpt1;
          higpt1 = baby.mc_pt()->at(imc);
        } else if (baby.mc_pt()->at(imc)>higpt2) {
          higpt2 = baby.mc_pt()->at(imc);
        }
      } // It's a Higgs!
      if (abs(baby.mc_id()->at(imc))==5 && abs(baby.mc_mom()->at(imc))==25){  // It's a b from a Higgs!
        if(baby.mc_pt()->at(imc)>minpt && fabs(baby.mc_eta()->at(imc))<maxeta) nbacc++;
        if(hi_dr){
          float eta = baby.mc_eta()->at(imc), phi = baby.mc_phi()->at(imc);
          int momidx = baby.mc_momidx()->at(imc);
          if(momidx1<0){
            preveta1 = eta;
            prevphi1 = phi;
            momidx1  = momidx;
          } else if(momidx == momidx1) dr1 = deltaR(preveta1, prevphi1, eta, phi);
          if(momidx != momidx1){
            if(momidx2<0){
              preveta2 = eta;
              prevphi2 = phi;
              momidx2 = momidx;
            } else dr2 = deltaR(preveta2, prevphi2, eta, phi);
          }
          
        } // hi_dr
      }// It's a b from a Higgs!
    } // Loop over mc particles

    histos[all].Fill(higpt1, higpt2, weight);
    if(nbacc<4) continue;
    if(hi_dr && (dr1<0.5 || dr2<0.5)) continue;
    if(baby.nleps()>0) continue;

    ///////////// Calculating other cuts
    bool pass_am = baby.higd_am()>100. && baby.higd_am()<=140.;
    bool pass_nj4 = baby.njets()>=4;
    bool pass_goodreco = pass_nj4;
    for(unsigned ijet=0; ijet<baby.jets_pt()->size(); ijet++){
      if(baby.jets_h1d()->at(ijet) || baby.jets_h2d()->at(ijet))
        if(baby.jets_hflavor()->at(ijet)!=5) pass_goodreco = false;
    } // Loop over jets

    ///////////// Calculating yields
    histos[acc].Fill(higpt1, higpt2, weight);
    if(pass_nj4) histos[nj4].Fill(higpt1, higpt2, weight);
    if(pass_goodreco) {
      histos[goodreco].Fill(higpt1, higpt2, weight);
      if(pass_am) histos[gr_am].Fill(higpt1, higpt2, weight);
    } // goodreco
  } // Loop over entries

  cout<<endl;
  TString acc_s = "p_{T} > "+RoundNumber(minpt,0)+" GeV, |#eta| < "+RoundNumber(maxeta,1);
  if(do_aux){
    makePlot(histos[goodreco], histos[acc], "all jets in di-Higgs reconstruction", 
             "all b-quarks from hh#rightarrow4b have "+acc_s, "Numer-matched to B-hadrons");
    makePlot(histos[gr_am], histos[goodreco], "100 < #LTm#GT < 140 GeV",
             "all b-quarks from hh#rightarrow4b have "+acc_s,
             "Denom-+ all jets in di-Higgs reco. matched to B-hadrons");
  } else {
    makePlot(histos[goodreco], histos[acc], "all di-Higgs jets TM-ed", 
             "b-quarks in acceptance");
    makePlot(histos[gr_am], histos[goodreco], "100 < #LTm#GT < 140",
             "all di-Higgs jets TM-ed");
  }

  makePlot(histos[nj4], histos[acc], "N_{jets} #geq 4", "b-quarks in acceptance");
  makePlot(histos[goodreco], histos[nj4], "All b's reco-ed", "N_{jets} #geq 4");

  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Finding plots took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
} // main

////////////////////////////////////////// End of main //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void makePlot(TH2D hnum, TH2D hden, TString num_s, TString den_s, TString extra_s){

  TCanvas can("can","");

  ////// Setting additional histogram style
  hnum.GetXaxis()->SetLabelOffset(0.012);
  hnum.GetXaxis()->SetTitle("Leading Higgs boson p_{T} [GeV]");

  hnum.GetYaxis()->SetTitle("Subleading Higgs boson p_{T} [GeV]");

  hnum.GetZaxis()->SetLabelSize(0.04);
  hnum.GetZaxis()->SetRangeUser(0.,1.);
  hnum.GetZaxis()->SetTitleSize(0.045);
  hnum.GetZaxis()->CenterTitle(true);
  hnum.GetZaxis()->SetTitle("Efficiency");

  ////// Dividing and plot histogram
  hnum.Divide(&hden);
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
  
  ////// Drawing legend
  TBox box;
  Xmin = 85; Xmax = 730; Ymin = 780; Ymax = 975;
  box.SetFillColor(0); box.SetFillStyle(1001);
  box.SetLineColor(1); box.SetLineWidth(2); box.SetLineStyle(1);
  box.DrawBox(Xmin, Ymin, Xmax, Ymax);
  box.SetFillColor(0); box.SetFillStyle(0);
  box.SetLineColor(1); box.SetLineWidth(2); box.SetLineStyle(1);
  box.DrawBox(Xmin, Ymin, Xmax, Ymax);

  float lMargin = opts.LeftMargin(), rMargin = opts.RightMargin(), tMargin = opts.TopMargin();
  TLatex label; label.SetTextFont(42); label.SetTextSize(0.05);
  if(do_aux) {
    TString cmsLogo = "#font[62]{CMS}#scale[0.8]{#font[52]{ Simulation Supplementary}}";
    TString lumiEner = "#font[42]{13 TeV}";
    label.SetTextAlign(11); label.SetTextSize(0.06);
    label.DrawLatexNDC(lMargin+0.005, 1-tMargin+0.015, cmsLogo);
    label.SetTextAlign(31); label.SetTextSize(0.05);
    label.DrawLatexNDC(1-rMargin-0.005, 1-tMargin+0.015, lumiEner);

    //// Printing legend
    label.SetTextAlign(13); label.SetTextSize(0.03);
    den_s = "#font[72]{Den}: "+den_s;
    double xoffset = 15;
    label.DrawLatex(Xmin+xoffset, Ymax-25, den_s);
    double offset = 0;
    if(extra_s.Contains("Denom-")){
      extra_s.ReplaceAll("Denom-","");
      label.DrawLatex(Xmin+xoffset+65, Ymax-80, extra_s);
      offset = 55;
    }
    num_s = "#font[72]{Num}: Den. + "+num_s;
    label.DrawLatex(Xmin+xoffset, Ymax-80-offset, num_s);
    if(extra_s.Contains("Numer-")){
      extra_s.ReplaceAll("Numer-","");
      label.DrawLatex(Xmin+xoffset+70, Ymax-135, extra_s);
    }
  } else {
    //// Printing histo title
    TString title = "#font[72]{Den}: "; 
    if(den_s!="") title += den_s;
    else title += hden.GetTitle();
    title += " - #font[72]{Num:} ";
    if(num_s!="") title += num_s;
    else title += hnum.GetTitle();
    label.SetTextAlign(21);
    label.DrawLatexNDC(lMargin+(1-rMargin-lMargin)/2., 1-tMargin+0.025, title);

    //// Printing legend
    label.SetTextAlign(13); label.SetTextSize(0.032);
    TString acc_s = "All b quarks from hh#rightarrow4b have p_{T} > "+RoundNumber(minpt,0)+" GeV, |#eta| < "
      +RoundNumber(maxeta,1);
    label.DrawLatex(Xmin+20, Ymax-25,acc_s);
    TString cut_s = (low_isr?"ISR p_{T} < 20 GeV":"All ISR");
    cut_s += (low_pu?", N_{PV} < 10":", All PU");
    if(hi_dr) cut_s += ", #DeltaR(b,b)_{h1,h2} > 0.5";
    label.DrawLatex(Xmin+20, Ymax-80,cut_s);
  }
  hnum.SetTitle("");

  ////// Saving plot
  TString hname = "plots/eff_"; 
  hname += hnum.GetName(); hname += "_";
  hname += hden.GetName(); 
  if(minpt!=30) hname += "_minpt"+RoundNumber(minpt,0);
  if(RoundNumber(maxeta,1) != "2.4") hname += "_maxeta"+RoundNumber(maxeta,1);
  if(low_isr) hname += "_lowisr";
  if(low_pu) hname += "_lowpu";
  if(hi_dr) hname += "_hidr";
  if(mass != "") hname += "_mass"+mass;
  if(do_aux) hname += "_aux";
  hname.ReplaceAll(".", "p");
  hname += ".pdf";
  can.SaveAs(hname);
  cout<<" open "<<hname<<endl;

  hname.ReplaceAll(".pdf", ".root");
  hname.ReplaceAll("plots/","");
  TFile file(hname, "recreate");
  file.cd();
  hname.ReplaceAll(".root","");
  hnum.Write(hname);
  file.Close();
  //cout<<"Saved histogram in "<<hname<<endl<<endl;

}


void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"mass", required_argument, 0, 'm'}, // Higgsino mass
      {"low_isr", no_argument, 0, 'i'},    // Plot for low ISR
      {"hi_dr", no_argument, 0, 'r'},      // Plot for high dR between b daughters
      {"aux", no_argument, 0, 0},          // Put AUX labels
      {"low_pu", no_argument, 0, 'p'},     // Plot for low PU
      {"pt", required_argument, 0, 0},     // Minimum pT for quark to be in acceptance
      {"eta", required_argument, 0, 0},    // Maximum eta for quark to be in acceptance
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "m:ipr", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'm':
      mass = optarg;
      break;
    case 'i':
      low_isr = true;
      break;
    case 'r':
      hi_dr = true;
      break;
    case 'p':
      low_pu = true;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "pt"){
        minpt = atof(optarg);
      } else if(optname == "eta"){
        maxeta = atof(optarg);
      } else if(optname == "aux"){
        do_aux = true;
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
