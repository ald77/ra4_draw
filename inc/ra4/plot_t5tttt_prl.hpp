// plot_limits_summary: Plots various limit curves on same canvas

#ifndef H_PLOT_LIMITS_SUMMARY
#define H_PLOT_LIMITS_SUMMARY

// System includes
#include <iostream>
#include <vector>

// ROOT includes
#include "TString.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"

class model_limits {
public:
  TString model, title, labMass;
  float legScale;
  std::vector<TString> labels, files, obsnames, expnames;
  std::vector<int> colors;

  void add(TString label, TString file, int color, TString obsname, TString expname);
  model_limits(TString imodel, TString ititle, float ilegScale=1.);
};

void setCanvas(TCanvas &can, float lMargin, float tMargin, float rMargin, float bMargin);
TH2D baseHistogram(float Xmin, float Xmax, float Ymin, float Ymax);
void addLabelsTitle(float lMargin, float tMargin, float rMargin, TString title);
TString altName(const TString &name);
TH2D* getHist2D(TFile &flimit, TString hname, bool allow_name_change = true);
TGraph* getGraph(TFile &flimit, TString gname, bool allow_name_change = true);
TGraph* joinGraphs(TGraph *graph1, TGraph *graph2);
void setGraphStyle(TGraph *graph, int color, int style, int width, double glu_lsp, TString model);
void getModelParams(TString model, float &Xmin, float &Xmax, float &Ymin, float &Ymax, float &glu_lsp);
void reverseGraph(TGraph *graph);
void printExclGlu(TGraph *gobs, TGraph *gexp, std::vector<float> &mLSPs, TString label);
std::vector<float> intersectionLSP(TGraph *graph, std::vector<float> &mLSPs);

#endif
