//// styles: Utilities to change styles of plots

#include "TStyle.h"

#include "styles.hpp"

void setPlotStyle(PlotOpt opts){
  //// General
  gStyle->SetNdivisions(opts.NDivisions(), "xyz");   
  gStyle->SetOptStat(0);              // No Stats box
  gStyle->SetPadTickX(1);             // Ticks at the right
  gStyle->SetPadTickY(1);             // Ticks at the top

  //// Pad
  gStyle->SetPadRightMargin (opts.RightMargin());    
  gStyle->SetPadBottomMargin(opts.BottomMargin()); 
  gStyle->SetPadTopMargin   (opts.TopMargin()); 
  gStyle->SetPadLeftMargin  (opts.LeftMargin()); 

  //// Histogram
  gStyle->SetTitleOffset(opts.XTitleOffset(),"x");	// Set offset of X title in histogram
  gStyle->SetTitleOffset(opts.YTitleOffset(),"y");      // Set offset of Y title in histogram
  gStyle->SetTextSize(opts.TitleSize());		// Set global text size
  gStyle->SetTitleFontSize(opts.TitleSize());		// Set top title size
  gStyle->SetTitleSize(opts.TitleSize(),"xyz");		// Set the 2 axes title size
  gStyle->SetLabelSize(opts.LabelSize(),"xyz");		// Set the 2 axes label size

}
