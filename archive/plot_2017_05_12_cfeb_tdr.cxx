#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include <unistd.h>
#include <getopt.h>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"

#include "core/styles.hpp"
#include "core/utilities.hpp"
#include "core/plot_opt.hpp"

using namespace std;

namespace{
  TString filename = "txt/tdr/cfeb_dcfeb_loss.txt";
  PlotOpt opts("txt/plot_styles.txt", "OneCol1D");
}

void GetOptions(int argc, char *argv[]);

void multVector(vector<float> &vals, float factor){
  for(unsigned ind=0; ind<vals.size(); ind++) vals[ind] *= factor;
}

void readFile(TString fname, int ncols, vector<vector<float> > &vals){
  vals = vector<vector<float> >(ncols, vector<float>());
  for(int col=0; col<ncols; col++) vals[col].push_back(0.);

  ifstream infile(fname);
  string line_s;
  float val;
  while(getline(infile, line_s)){
    istringstream iss(line_s);
    for(int col=0; col<ncols; col++) {
      iss >> val;
      if(fabs(val)>1e-10) vals[col].push_back(val);
    } // Loop over columns
  } // Loop over rows
  infile.close();
}

int main(int argc, char *argv[]){
  GetOptions(argc, argv);

  //// Setting plot style
  setPlotStyle(opts);
  gStyle->SetGridStyle(3);


  int linw = 3;
  vector<vector<float> > vals;

  TCanvas can;
  //can.SetGrid(); 
  can.SetFillStyle(4000);
  float minY, maxY;
  TLine line; line.SetLineWidth(2); line.SetLineStyle(2); line.SetLineColor(kOrange+4);
  double legSize = 0.05;
  double legX(0.85), legY(1-opts.TopMargin()-0.05), legSingle = 0.07;
  double legW = 0.12, legH = legSingle*3;
  TString pname;
  TLatex label; label.SetTextSize(0.058); label.SetTextFont(42); label.SetTextAlign(13);
  label.SetTextColor(kOrange+4);

  ///////////////////////////////////////////////////////////////////////
  //////////////////////// CFEBS LOSS MODEL /////////////////////////////
  can.SetLogy(true);

  minY = 0.0005; maxY = 1;
  TH1D histoModel("histoModel", "", 18, 0, 7);
  histoModel.SetMinimum(minY);
  histoModel.SetMaximum(maxY);
  histoModel.GetYaxis()->CenterTitle(true);
  histoModel.GetXaxis()->SetLabelOffset(0.01);
  histoModel.SetXTitle("Instantaneous luminosity [10^{34} s^{-1} cm^{-2}]");
  histoModel.SetYTitle("Event loss fraction");
  histoModel.SetTitle("Event loss for 750 kHz L1A rate, 12.5 #mus latency");
  histoModel.SetTitle("Validation of the statistical model for CFEB losses");
  histoModel.Draw();

  readFile("txt/tdr/cfeb_model.txt", 5, vals);
  multVector(vals[0], 10); multVector(vals[3], 10);

  TGraph gr21b(vals[3].size(), &(vals[3][0]), &(vals[4][0]));
  gr21b.SetLineWidth(linw); gr21b.SetLineColor(kRed); //gr21b.SetLineStyle(2); 
  gr21b.Draw("same"); 
  vector<float> zeroes(vals[0].size());
  TGraphErrors grData(vals[0].size(), &(vals[0][0]), &(vals[1][0]), &(zeroes[0]),  &(vals[2][0]));
  grData.SetMarkerStyle(20); grData.SetMarkerSize(0.8); grData.SetLineColor(kBlack); 
  grData.Draw("same p"); 

  legH = legSingle*2; legY = 0.45; legX = 0.57;
  TLegend legModel(legX-legW, legY-legH, legX, legY);
  legModel.SetTextSize(legSize); legModel.SetFillColor(0); 
  legModel.SetFillStyle(0); legModel.SetBorderSize(0);
  legModel.AddEntry(&grData, "Measured CFEB losses", "p");
  legModel.AddEntry(&gr21b,  "Statistical model", "l");
  legModel.Draw();

  histoModel.Draw("axis same");
  pname = "plots/cfeb_model_loss.pdf";
  can.SaveAs(pname);


  /////////////////////////////////////////////////////////////////
  //////////////////////// CFEBS LOSS /////////////////////////////
  can.SetLogy(false);
  minY = 0; maxY = 0.1;
  TH1D histo("histo", "", 18, 0, 9);
  histo.SetMinimum(minY);
  histo.SetMaximum(maxY);
  histo.GetYaxis()->CenterTitle(true);
  histo.GetXaxis()->SetLabelOffset(0.01);
  histo.SetXTitle("Instantaneous luminosity [10^{34} s^{-1} cm^{-2}]");
  histo.SetYTitle("Event loss fraction");
  histo.SetTitle("Event loss for 750 kHz L1A rate, 12.5 #mus latency");
  histo.SetTitle("CFEB event losses for HL-LHC conditions");
  histo.Draw();

  readFile("txt/tdr/cfeb_loss.txt", 6, vals);
  //multVector(vals[0], 0.1); multVector(vals[2], 0.1); multVector(vals[4], 0.1);

  line.DrawLine(5, minY, 5, maxY);
  label.DrawLatex(5.4, 0.017, "HL-LHC");
  label.DrawLatex(5.4, 0.01, "luminosity");

  //for(unsigned ind=0; ind<vals[0].size(); ind++) cout<<setw(8)<<vals[0][ind]<<", "<<setw(8)<<vals[1][ind]<<endl;
  TGraph gr21(vals[0].size(), &(vals[0][0]), &(vals[1][0]));
  gr21.SetLineWidth(linw); gr21.SetLineColor(kRed); 
  //gr21.Draw("same"); 
  gr21b.Draw("same"); 
  
  TGraph gr31(vals[2].size(), &(vals[2][0]), &(vals[3][0]));
  gr31.SetLineWidth(linw); gr31.SetLineColor(kBlue+1); 
  gr31.Draw("same"); 
  TGraph gr41(vals[4].size(), &(vals[4][0]), &(vals[5][0]));
  gr41.SetLineWidth(linw); gr41.SetLineColor(kGreen+1); 
  gr41.Draw("same"); 

  legH = legSingle*3; legY = 0.7; legX = 0.82;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(legSize); leg.SetFillColor(0); 
  leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.AddEntry(&gr21, "ME2/1", "l");
  leg.AddEntry(&gr31, "ME3/1", "l");
  leg.AddEntry(&gr41, "ME4/1", "l");
  leg.Draw();

  histo.Draw("axis same");
  pname = "plots/cfeb_loss.pdf";
  can.SaveAs(pname);

  ///////////////////////////////////////////////////////////////////////
  //////////////////////// CFEBS VS DCFEB /////////////////////////////
  opts.LoadOptions("txt/plot_styles.txt", "TwoCol1D");
  setPlotStyle(opts);
  gStyle->SetGridStyle(3);
  TCanvas can2; can2.cd();
  can2.SetLogy(false);


  minY = 0; maxY = 0.06;
  TH1D histoDCFEB("histoDCFEB", "", 18, 0, 38);
  histoDCFEB.SetMinimum(minY);
  histoDCFEB.SetMaximum(maxY);
  histoDCFEB.GetYaxis()->CenterTitle(true);
  histoDCFEB.GetXaxis()->SetLabelOffset(0.01);
  histoDCFEB.SetXTitle("Instantaneous luminosity [10^{34} s^{-1} cm^{-2}]");
  histoDCFEB.SetYTitle("Event loss fraction");
  histoDCFEB.SetTitle("(D)CFEB event losses for HL-LHC conditions");
  histoDCFEB.Draw();

  readFile("txt/tdr/cfeb_dcfeb_loss.txt", 4, vals);
  multVector(vals[0], 0.1); multVector(vals[2], 0.1);

  line.DrawLine(5, minY, 5, maxY);
  label.DrawLatex(5.4, 0.013, "HL-LHC");
  label.DrawLatex(5.4, 0.008, "luminosity");

  TGraph grCFEB(vals[0].size(), &(vals[0][0]), &(vals[1][0]));
  grCFEB.SetLineWidth(linw); grCFEB.SetLineColor(kRed); 
  //grCFEB.Draw("same"); 
  gr21b.Draw("same"); 
  TGraph grDCFEB(vals[2].size(), &(vals[2][0]), &(vals[3][0]));
  grDCFEB.SetLineWidth(linw); grDCFEB.SetLineColor(kMagenta+2); 
  grDCFEB.Draw("same"); 

  legH = legSingle*2; legY = 0.8; legX = 0.5;
  TLegend legDCFEB(legX-legW, legY-legH, legX, legY);
  legDCFEB.SetTextSize(legSize); legDCFEB.SetFillColor(0); 
  legDCFEB.SetFillStyle(0); legDCFEB.SetBorderSize(0);
  legDCFEB.AddEntry(&grCFEB,  "ME2/1 CFEB (Phase 1)", "l");
  legDCFEB.AddEntry(&grDCFEB, "ME2/1 DCFEB (Phase 2)", "l");
  legDCFEB.Draw();

  histoDCFEB.Draw("axis same");
  pname = "plots/cfeb_dcfeb_loss.pdf";
  can2.SaveAs(pname);


  ///////////////////////////////////////////////////////////////////////
  //////////////////////// CFEBS VS DCFEB INTRO /////////////////////////////
  minY = 0; maxY = 0.1;
  TH1D histoIntro("histoIntro", "", 18, 0, 33);
  histoIntro.SetMinimum(minY);
  histoIntro.SetMaximum(maxY);
  histoIntro.GetYaxis()->CenterTitle(true);
  histoIntro.GetXaxis()->SetLabelOffset(0.01);
  histoIntro.SetXTitle("Instantaneous luminosity [10^{34} s^{-1} cm^{-2}]");
  histoIntro.SetYTitle("Event loss fraction");
  histoIntro.SetTitle("(D)CFEB event losses for HL-LHC conditions");
  histoIntro.Draw();

  line.DrawLine(5, minY, 5, maxY);
  label.DrawLatex(5.4, 0.017, "HL-LHC");
  label.DrawLatex(5.4, 0.01, "luminosity");

  gr21b.Draw("same"); gr31.Draw("same"); gr41.Draw("same"); 
  grDCFEB.Draw("same"); 

  legH = legSingle*4; legY = 0.8; legX = 0.62;
  TLegend legIntro(legX-legW, legY-legH, legX, legY);
  legIntro.SetTextSize(legSize); legIntro.SetFillColor(0); 
  legIntro.SetFillStyle(0); legIntro.SetBorderSize(0);
  legIntro.AddEntry(&grCFEB, "ME2/1 CFEB (Phase 1)", "l");
  legIntro.AddEntry(&gr31,   "ME3/1 CFEB (Phase 1)", "l");
  legIntro.AddEntry(&gr41,   "ME4/1 CFEB (Phase 1)", "l");
  legIntro.AddEntry(&grDCFEB,  "ME2/1 DCFEB (Phase 2)", "l");
  legIntro.Draw();

  histoIntro.Draw("axis same");
  pname = "plots/csc_cfeb_dcfeb_loss.pdf";
  can2.SaveAs(pname);
}



void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"file", required_argument, 0, 'f'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "f:", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'f':
      filename = optarg;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == ""){
        printf("Bad option! Found option name %s\n", optname.c_str());
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
