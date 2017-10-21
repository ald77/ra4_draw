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
  TString filename = "txt/tdr/tid_me11.txt";
  PlotOpt opts("txt/plot_styles.txt", "TwoCol1D");
}

void GetOptions(int argc, char *argv[]);

void multVector(vector<float> &vals, float factor, float &minval, float &maxval){
  minval = 1e22; maxval = -1e22;
  for(unsigned ind=0; ind<vals.size(); ind++) {
    vals[ind] *= factor;
    if(vals[ind] < minval) minval = vals[ind];
    if(vals[ind] > maxval) maxval = vals[ind];
  }
}

void readFile(TString fname, int ncols, vector<vector<float> > &vals){
  vals = vector<vector<float> >(ncols, vector<float>());

  ifstream infile(fname);
  string line_s;
  float val;
  long nlines = 0;
  while(getline(infile, line_s)){
    istringstream iss(line_s);
    for(int col=0; col<ncols && nlines>0; col++) {
      iss >> val;
      if(fabs(val)>1e-10) vals[col].push_back(val);
    } // Loop over columns
    nlines++;
  } // Loop over rows
  infile.close();
}

  
int main(int argc, char *argv[]){
  GetOptions(argc, argv);

  //// Setting plot style
  setPlotStyle(opts);
  gStyle->SetGridStyle(3);


  vector<vector<float> > vals;

  TCanvas can;
  //can.SetGrid(); 
  can.SetFillStyle(4000);
  float minX, maxX, minY, maxY, delta;
  TLine line; line.SetLineWidth(2); line.SetLineStyle(2); line.SetLineColor(kOrange+4);
  double legSize = 0.05;
  double legX(0.85), legY(1-opts.TopMargin()-0.05), legSingle = 0.07;
  double legW = 0.12, legH = legSingle*3;
  TString pname;
  TLatex label; label.SetTextSize(0.058); label.SetTextFont(42); label.SetTextAlign(13);
  label.SetTextColor(kOrange+4);

  ///////////////////////////////////////////////////////////////////////
  //////////////////////// TID for ME1/1 /////////////////////////////
  can.SetLogy(true);

  readFile("txt/tdr/tid_me11.txt", 3, vals);
  multVector(vals[2], 0.1, minY, maxY);  // Converting Gy to kRad
  for(unsigned ind=0; ind<vals[0].size(); ind++) vals[0][ind] *= vals[2][ind]/100.; // Convertion % error to absolute

  // Plot ranges
  minX = 10*floor(vals[1][0]/10.); maxX = 10*floor(vals[1][vals[1].size()-1]/10.+1);
  delta = maxY/minY; minY = minY/(0.04*delta); maxY = maxY*(0.1*delta);

  cout<<"minX "<<minX<<", maxX "<<maxX<<", minY "<<minY<<", maxY "<<maxY<<endl;
  
  TH1D histoModel("histoModel", "", 18, minX, maxX);
  histoModel.SetMinimum(minY);
  histoModel.SetMaximum(maxY);
  histoModel.GetYaxis()->CenterTitle(true);
  histoModel.GetXaxis()->SetLabelOffset(0.01);
  histoModel.GetYaxis()->SetTitleOffset(0.9);
  histoModel.SetXTitle("r [cm]");
  histoModel.SetYTitle("Dose [kRad]");
  histoModel.SetTitle("TID in ME1/1 for 3000 fb^{-1} of HL-LHC");
  histoModel.Draw();

  vector<float> zeroes(vals[0].size());
  TGraphErrors gr11(vals[0].size(), &(vals[1][0]), &(vals[2][0]), &(zeroes[0]),  &(vals[0][0]));
  gr11.SetMarkerStyle(20); gr11.SetMarkerSize(1); gr11.SetLineColor(kBlack); 
  gr11.Draw("same p"); 

  legH = legSingle*2; legY = 0.4; legX = 0.54;
  TLegend legModel(legX-legW, legY-legH, legX, legY);
  legModel.SetTextSize(legSize); legModel.SetFillColor(0); 
  legModel.SetFillStyle(0); legModel.SetBorderSize(0);
  legModel.AddEntry(&gr11, "Measured CFEB losses", "p");
  // legModel.Draw();

  histoModel.Draw("axis same");
  pname = "plots/tid_me11.pdf";
  can.SaveAs(pname);
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
