//----------------------------------------------------------------------------
// styles - Class to set default plotting styles, read from a text file
//----------------------------------------------------------------------------

#ifndef INT_ROOT
#include "styles.hpp"
#endif

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TROOT.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using std::cout;
using std::endl;

styles::styles(TString group) {
  confFile = "txt/plot_styles.txt";
  TString inames[] = {"nFont", "nDivisions"};
  int *ivalues[] = {&nFont, &nDivisions};
  parseStyleFile("General", 0, 0, 0, inames, ivalues, 2);
  Group = group;
  nPads = 1;

  readGroupStyle();
}

void styles::setDefaultStyle() {
  setGlobalStyle();
  gStyle->SetCanvasDefW(CanvasW);
  gStyle->SetCanvasDefH(CanvasH);
  gStyle->SetTextSize(LegendSize);            // Set global text size
  gStyle->SetTitleFontSize(TitleSize);      // Set top title size
  gStyle->SetTitleSize(LabelSize,"xyz");     // Set the 2 axes title size
  gStyle->SetLabelSize(LabelSize,"xyz");     // Set the 2 axes label size

  gStyle->SetTitleOffset(xTitleOffset,"x");     
  gStyle->SetTitleOffset(yTitleOffset,"y");     
  gStyle->SetTitleOffset(zTitleOffset,"z");     
  gStyle->SetPadRightMargin (PadRightMargin);    
  gStyle->SetPadBottomMargin(PadBottomMargin); 
  gStyle->SetPadTopMargin(PadTopMargin); 
  gStyle->SetPadLeftMargin  (PadLeftMargin); 
  gStyle->SetNdivisions(nDivisions, "xyz");   // 5 primary ticks and 4 secondary ticks

  gStyle->SetTitleFont(nFont,"xyz");          // Set the all 2 axes title font
  gStyle->SetLabelFont(nFont,"xyz");          // Set the all 2 axes label font
  gStyle->SetTextFont(nFont);                // Set global text font
}

void styles::setHistoStyle(TH1 *h) {
  h->SetTitleSize(TitleSize,"xyz");     // Set the 2 axes title size
  h->SetLabelSize(LabelSize,"xyz");     // Set the 2 axes label size

  h->SetTitleOffset(xTitleOffset,"x");     
  h->SetTitleOffset(yTitleOffset,"y");     
  h->SetTitleOffset(zTitleOffset,"z");     
  h->SetNdivisions(nDivisions, "xyz");   // 5 primary ticks and 4 secondary ticks
  h->SetTitleFont(nFont,"xyz");          // Set the all 2 axes title font
  h->SetLabelFont(nFont,"xyz");          // Set the all 2 axes label font
}

// Set default styles globally.   
void styles::setGlobalStyle() {
  gStyle->SetPalette(1);              // Decent colors for 2D plots
  gStyle->SetOptStat(0);              // No Stats box
  gStyle->SetPadTickX(1);             // No ticks at the right
  gStyle->SetPadTickY(1);             // No ticks at the top
}

void styles::setGroup(TString group){
  Group = group;
  readGroupStyle();
}

// Set default style for the specific group 
void styles::readGroupStyle() {
  TString inames[] = {"CanvasW", "CanvasH"};
  TString fnames[] = {"LegendSize", "TitleSize", "LabelSize", "PadRightMargin", "PadTopMargin", "PadBottomMargin",
                      "xTitleOffset", "PadLeftMargin", "yTitleOffset", "zTitleOffset"};
  int   *ivalues[] = {&CanvasW, &CanvasH};
  float *fvalues[] = {&LegendSize,&TitleSize,&LabelSize,&PadRightMargin,&PadTopMargin,&PadBottomMargin,
                      &xTitleOffset,&PadLeftMargin,&yTitleOffset,&zTitleOffset};
  parseStyleFile(Group, fnames, fvalues, 10, inames, ivalues, 2);
}

// Fix for y axes that have too much/little space for the label due to number of digits
void styles::moveYAxisLabel(TH1 *h, float maxi, bool isLog){
  int digits = static_cast<int>((log(maxi)/log(10.)+0.001)+1);
  if(digits<2) digits = 2;
  TString Section = Group; 
  if(isLog)Section += "_Log";
  else {Section += "_Digits_"; Section += digits;  }
  TString fnames[] = {"PadLeftMargin", "yTitleOffset"};
  float *fvalues[] = {&PadLeftMargin, &yTitleOffset};

  parseStyleFile(Section, fnames, fvalues, 2, 0, 0, 0);
  h->SetTitleOffset(yTitleOffset,"y");
  //pad->SetLeftMargin(PadLeftMargin);
}

// Test the global style settings for a generic histogram.  
void styles::testGlobalStyle(bool fixY, float scale) {
  
  readGroupStyle(); setGlobalStyle(); setDefaultStyle();
  
  TH1* h = new TH1F("h", "h", 50, 0, 50);
  TH1* hc[6];
  for (int i=1; i<=50; i++) {
    double value = scale*exp(-0.5*pow(((i-25.)/5.),2));  // Gaussian shape
    h->SetBinContent(i, value);
  }

  TCanvas c;
  if(nPads == 2) c.Divide(2);
  if(nPads == 3) c.Divide(3);
  if(nPads == 4) c.Divide(2,2);
  if(nPads == 6) c.Divide(3,2);
  c.cd(1);
  h->Draw();
  if(fixY) moveYAxisLabel(h,100);
  setTitles(h, "D^{(*)0/+} channels", "xlabel^{2}_{miss} (GeV^{2})", "Events/(10 MeV^{2})");
  float scales[] = {0.1, 10, 0.01};
  for(int pads = 2; pads<=4; pads++){
    if(nPads>=pads){
      c.cd(pads);
      hc[pads-2] = static_cast<TH1F*>(h->Clone());
      hc[pads-2]->Scale(scales[pads-2]); 
      if(fixY) moveYAxisLabel(hc[pads-2],hc[pads-2]->GetMaximum());
      hc[pads-2]->Draw();
      setTitles(hc[pads-2], "D^{(*)0/+} channels", "xlabel^{2}_{miss} (GeV^{2})", "Events/(1000 MeV^{2})");
    }
  }
  TString epsName = "babar_code/styles/Plot_"; epsName += nPads; epsName += "Pads.eps";
  c.Print(epsName);
  
}

void styles::printValues() {
  cout<<"nFont           = " << nFont           << endl;
  cout<<"nPads           = " << nPads           << endl;
  cout<<"nDivisions      = " << nDivisions      << endl;
  cout<<"CanvasW         = " << CanvasW         << endl;   
  cout<<"CanvasH         = " << CanvasH         << endl;   
  cout<<"LegendSize      = " << LegendSize      << endl;   
  cout<<"TitleSize       = " << TitleSize       << endl;
  cout<<"LabelSize       = " << LabelSize       << endl;
  cout<<"PadRightMargin  = " << PadRightMargin  << endl;
  cout<<"PadTopMargin    = " << PadTopMargin    << endl;
  cout<<"PadBottomMargin = " << PadBottomMargin << endl;
  cout<<"xTitleOffset    = " << xTitleOffset    << endl;
  cout<<"PadLeftMargin   = " << PadLeftMargin   << endl;
  cout<<"yTitleOffset    = " << yTitleOffset    << endl;
  cout<<"zTitleOffset    = " << zTitleOffset    << endl;

}

// ----------------------------------------------------------------------
void styles::setMarkers(TH1 *h, float Msize, int Mstyle) {
  h->SetMarkerStyle(Mstyle);
  h->SetMarkerSize(Msize);
}

// ----------------------------------------------------------------------
void styles::setTitles(TH1 *h, TString xTitle, TString yTitle, TString Left, TString Right) {
  if (0==h) {
    cout << " Histogram not defined" << endl;
  } else {
    h->SetXTitle(xTitle); h->SetYTitle(yTitle);
    TLatex label; label.SetNDC(kTRUE);
    label.SetTextAlign(11);
    label.DrawLatex(PadLeftMargin+0.02,1-PadTopMargin+0.02,Left);  
    label.SetTextAlign(31);
    label.DrawLatex(1-PadRightMargin-0.02,1-PadTopMargin+0.02,Right);  
  }
}

// ----------------------------------------------------------------------
void styles::styleHist(TH1 *h, Int_t color, Int_t fillstyle,
                       Int_t symbol, Double_t size, Int_t width) {
  h->SetLineColor(color);   
  h->SetLineWidth(width);
  h->SetMarkerColor(color); 
  h->SetMarkerStyle(symbol);  
  h->SetMarkerSize(size); 
  h->SetStats(kFALSE); 
  h->SetFillStyle(fillstyle); 
  h->SetFillColor(color);
}

// ----------------------------------------------------------------------
void styles::setTitleSizes(TH1 *h,  float size, float lsize, int font,
                           float xoff, float yoff, int divisions) {
  if (0==h) {
    cout << " Histogram not defined" << endl;
  } else {
    h->SetTitleOffset(xoff, "x");      h->SetTitleOffset(yoff, "y");
    h->SetTitleSize(size,   "x");      h->SetTitleSize(size,   "y");
    h->SetLabelSize(lsize,  "x");      h->SetLabelSize(lsize,  "y");
    h->SetLabelFont(font,   "x");      h->SetLabelFont(font,   "y");
    h->GetXaxis()->SetTitleFont(font); h->GetYaxis()->SetTitleFont(font);
    h->SetNdivisions(divisions,"X");   h->SetNdivisions(divisions,   "Y");
  }
}

void styles::parseStyleFile(TString group, TString fnames[], float *fvalues[], int nFloat, 
                            TString inames[], int *ivalues[], int nInt){
  std::ifstream file(confFile);
  TString word, s_value;
  while(file >> word){
    if(word.Contains("[")) {
      word.ReplaceAll("[",""); word.ReplaceAll("]",""); 
      if(word == group){
        while(file >> word){
          if(word.Contains("[")) break;
          file >> s_value; file >> s_value;
          for(int var(0); var < nFloat; var++){
            if(word == fnames[var]){
              *fvalues[var] = s_value.Atof();
              break;
            }
          } // Loop finding requested float variables
          for(int var(0); var < nInt; var++){
            if(word == inames[var]){
              *ivalues[var] = s_value.Atoi();
              break;
            }
          } // Loop finding requested int variables
        } // Loop over group variables
        break;
      }
    }
  } // Loop over all words
}

