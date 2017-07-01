#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include <unistd.h>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH2F.h"
#include "TColor.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLatex.h"
#include "TLine.h"

#include "core/styles.hpp"
#include "core/utilities.hpp"
#include "core/plot_opt.hpp"

using namespace std;

namespace{
  PlotOpt opts("txt/plot_styles.txt", "OneCol1D");
  double legSize = 0.05;
  double markerSize = 1.2;
}

void plotGasGain();
void plotDarkRate();
void plotDarkRateProblem();
void plotDarkCurrent();

int main(){

  //// Setting plot style
  setPlotStyle(opts);
  gStyle->SetGridStyle(3);


  plotGasGain();
  plotDarkRate();
  plotDarkRateProblem();
  plotDarkCurrent();

}

void plotGasGain(){

//=========Macro generated from canvas: cME11/ME11
//=========  (Fri Jun  2 14:40:30 2017) by ROOT version6.08/06
   TCanvas *cME11 = new TCanvas("cME11", "ME11");
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   cME11->SetHighLightColor(2);
   cME11->Range(-64.28571,0.5025,364.2857,1.2525);
   cME11->SetFillColor(0);
   cME11->SetBorderMode(0);
   cME11->SetBorderSize(2);
   cME11->SetGridx();
   cME11->SetGridy();
   cME11->SetTickx(1);
   cME11->SetTicky(1);
  
   TH2F *hME1dI1__1 = new TH2F("hME1dI1__1","",300,0,300,1200,0.65,1.2);
   hME1dI1__1->SetFillColor(1);
   hME1dI1__1->SetFillStyle(0);
   hME1dI1__1->SetLineStyle(0);
   hME1dI1__1->SetMarkerStyle(20);
   hME1dI1__1->GetXaxis()->CenterTitle(true);
   hME1dI1__1->GetXaxis()->SetTitle("Accumulated charge [mC/cm]");
   hME1dI1__1->GetYaxis()->CenterTitle(true);
   hME1dI1__1->GetYaxis()->SetTitle("Relative current I/I#kern[0.01]{#lower[-0.1]{_{ref}}}");
   hME1dI1__1->SetTitle("Relative gas gain vs Charge for ME1/1");
   hME1dI1__1->Draw("");
   
   Double_t relI_ME11_L2_fx1001[77] = {
   2.362894,
   5.146731,
   8.677283,
   14.62665,
   24.61029,
   32.05801,
   33.6714,
   38.74127,
   40.40048,
   47.86864,
   50.94259,
   52.68738,
   59.43014,
   62.69608,
   63.43046,
   67.63964,
   70.27712,
   73.80062,
   83.09959,
   113.2732,
   114.1608,
   121.3439,
   124.7864,
   126.9107,
   129.3835,
   135.6164,
   136.4528,
   140.7706,
   142.397,
   142.7824,
   143.5836,
   150.3945,
   152.3572,
   155.2885,
   156.5486,
   159.0907,
   159.9944,
   166.0897,
   171.2794,
   173.1347,
   180.051,
   184.2516,
   187.2563,
   193.5525,
   193.7259,
   194.8045,
   200.1376,
   205.8323,
   207.4543,
   211.6609,
   214.7788,
   222.0638,
   222.5167,
   224.4433,
   227.9533,
   228.5925,
   228.8852,
   229.5719,
   233.1883,
   235.2734,
   241.7202,
   243.4934,
   243.8268,
   250.1803,
   250.8601,
   252.0712,
   257.123,
   259.518,
   260.7698,
   264.4514,
   267.054,
   271.2835,
   272.1236,
   277.3676,
   277.8563,
   279.0993,
   279.9482};
   Double_t relI_ME11_L2_fy1001[77] = {
   1.086147,
   1.086798,
   1.090761,
   1.093724,
   1.090199,
   1.093444,
   1.088389,
   1.092092,
   1.087894,
   1.091505,
   1.085986,
   1.078235,
   1.020346,
   1.074172,
   1.056504,
   1.083106,
   1.075718,
   1.091697,
   1.069545,
   1.080206,
   1.077802,
   1.067187,
   1.073181,
   1.074754,
   1.072584,
   1.078226,
   1.092945,
   1.092982,
   1.100918,
   1.073259,
   1.074688,
   1.092273,
   1.091326,
   1.076357,
   1.088664,
   1.085249,
   1.080634,
   1.086954,
   1.094684,
   1.07907,
   1.087975,
   1.093699,
   1.080183,
   1.093542,
   1.069256,
   1.086597,
   1.069677,
   1.088149,
   1.087544,
   1.0906,
   1.085762,
   1.086903,
   1.086452,
   1.084623,
   1.085892,
   1.081635,
   1.088457,
   1.096183,
   1.086337,
   1.065752,
   1.09051,
   1.083837,
   1.085653,
   1.086185,
   1.092307,
   1.080792,
   1.063676,
   1.066809,
   1.054124,
   1.068858,
   1.071451,
   1.064696,
   1.061853,
   1.080012,
   1.063334,
   1.080454,
   1.07365};
   Double_t relI_ME11_L2_fex1001[77] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t relI_ME11_L2_fey1001[77] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   TGraphErrors *gre = new TGraphErrors(77,relI_ME11_L2_fx1001,relI_ME11_L2_fy1001,relI_ME11_L2_fex1001,relI_ME11_L2_fey1001);
   gre->SetName("relI_ME11_L2");
   gre->SetTitle("relI_ME11_L2");
   gre->SetFillColor(1);

   Int_t ci;      // for color index setting
   //TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#66ccff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(markerSize);
   
   TH1F *Graph_relI_ME1dI1_L21001 = new TH1F("Graph_relI_ME1dI1_L21001","relI_ME11_L2",100,0,307.7068);
   Graph_relI_ME1dI1_L21001->SetMinimum(1.012289);
   Graph_relI_ME1dI1_L21001->SetMaximum(1.108975);
   Graph_relI_ME1dI1_L21001->SetDirectory(0);
   Graph_relI_ME1dI1_L21001->SetStats(0);
   Graph_relI_ME1dI1_L21001->SetFillColor(1);
   Graph_relI_ME1dI1_L21001->SetFillStyle(0);
   Graph_relI_ME1dI1_L21001->SetLineStyle(0);
   Graph_relI_ME1dI1_L21001->SetMarkerStyle(20);
   Graph_relI_ME1dI1_L21001->GetXaxis()->SetNdivisions(505);
   Graph_relI_ME1dI1_L21001->GetXaxis()->SetLabelFont(42);
   Graph_relI_ME1dI1_L21001->GetXaxis()->SetLabelOffset(0.007);
   Graph_relI_ME1dI1_L21001->GetXaxis()->SetLabelSize(0.05);
   Graph_relI_ME1dI1_L21001->GetXaxis()->SetTitleSize(0.06);
   Graph_relI_ME1dI1_L21001->GetXaxis()->SetTitleFont(42);
   Graph_relI_ME1dI1_L21001->GetYaxis()->SetNdivisions(505);
   Graph_relI_ME1dI1_L21001->GetYaxis()->SetLabelFont(42);
   Graph_relI_ME1dI1_L21001->GetYaxis()->SetLabelOffset(0.007);
   Graph_relI_ME1dI1_L21001->GetYaxis()->SetLabelSize(0.05);
   Graph_relI_ME1dI1_L21001->GetYaxis()->SetTitleSize(0.06);
   Graph_relI_ME1dI1_L21001->GetYaxis()->SetTitleOffset(1.25);
   Graph_relI_ME1dI1_L21001->GetYaxis()->SetTitleFont(42);
   Graph_relI_ME1dI1_L21001->GetZaxis()->SetNdivisions(505);
   Graph_relI_ME1dI1_L21001->GetZaxis()->SetLabelFont(42);
   Graph_relI_ME1dI1_L21001->GetZaxis()->SetLabelOffset(0.007);
   Graph_relI_ME1dI1_L21001->GetZaxis()->SetLabelSize(0.05);
   Graph_relI_ME1dI1_L21001->GetZaxis()->SetTitleSize(0.06);
   Graph_relI_ME1dI1_L21001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_relI_ME1dI1_L21001);
   
   gre->Draw("ps");
   
   Double_t relI_ME11_L3_fx1002[77] = {
   2.362894,
   5.146731,
   8.677283,
   14.62665,
   24.61029,
   32.05801,
   33.6714,
   38.74127,
   40.40048,
   47.86864,
   50.94259,
   52.68738,
   59.43014,
   62.69608,
   63.43046,
   67.63964,
   70.27712,
   73.80062,
   83.09959,
   113.2732,
   114.1608,
   121.3439,
   124.7864,
   126.9107,
   129.3835,
   135.6164,
   136.4528,
   140.7706,
   142.397,
   142.7824,
   143.5836,
   150.3945,
   152.3572,
   155.2885,
   156.5486,
   159.0907,
   159.9944,
   166.0897,
   171.2794,
   173.1347,
   180.051,
   184.2516,
   187.2563,
   193.5525,
   193.7259,
   194.8045,
   200.1376,
   205.8323,
   207.4543,
   211.6609,
   214.7788,
   222.0638,
   222.5167,
   224.4433,
   227.9533,
   228.5925,
   228.8852,
   229.5719,
   233.1883,
   235.2734,
   241.7202,
   243.4934,
   243.8268,
   250.1803,
   250.8601,
   252.0712,
   257.123,
   259.518,
   260.7698,
   264.4514,
   267.054,
   271.2835,
   272.1236,
   277.3676,
   277.8563,
   279.0993,
   279.9482};
   Double_t relI_ME11_L3_fy1002[77] = {
   1.068012,
   1.06624,
   1.068679,
   1.073194,
   1.069561,
   1.072523,
   1.061694,
   1.070673,
   1.065652,
   1.069935,
   1.065415,
   1.055926,
   1.004095,
   1.070002,
   1.06786,
   1.079317,
   1.070694,
   1.085861,
   1.069748,
   1.075207,
   1.072331,
   1.067026,
   1.065627,
   1.064376,
   1.053866,
   1.072104,
   1.083457,
   1.090493,
   1.096656,
   1.06745,
   1.066958,
   1.086175,
   1.08025,
   1.070911,
   1.082283,
   1.075336,
   1.068648,
   1.077128,
   1.086805,
   1.070984,
   1.082004,
   1.086697,
   1.078927,
   1.086214,
   1.069417,
   1.082963,
   1.062132,
   1.084125,
   1.080302,
   1.083448,
   1.080113,
   1.080364,
   1.07991,
   1.07805,
   1.07834,
   1.069859,
   1.09138,
   1.090517,
   1.079836,
   1.068228,
   1.083996,
   1.078001,
   1.071197,
   1.084499,
   1.091398,
   1.07825,
   1.065249,
   1.063115,
   1.03801,
   1.069256,
   1.067001,
   1.062708,
   1.061817,
   1.08094,
   1.063916,
   1.079842,
   1.071236};
   Double_t relI_ME11_L3_fex1002[77] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t relI_ME11_L3_fey1002[77] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   gre = new TGraphErrors(77,relI_ME11_L3_fx1002,relI_ME11_L3_fy1002,relI_ME11_L3_fex1002,relI_ME11_L3_fey1002);
   gre->SetName("relI_ME11_L3");
   gre->SetTitle("relI_ME11_L3");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#cc3333");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(21);
   gre->SetMarkerSize(markerSize);
   
   TH1F *Graph_relI_ME1dI1_L31002 = new TH1F("Graph_relI_ME1dI1_L31002","relI_ME11_L3",100,0,307.7068);
   Graph_relI_ME1dI1_L31002->SetMinimum(0.9948386);
   Graph_relI_ME1dI1_L31002->SetMaximum(1.105912);
   Graph_relI_ME1dI1_L31002->SetDirectory(0);
   Graph_relI_ME1dI1_L31002->SetStats(0);
   Graph_relI_ME1dI1_L31002->SetFillColor(1);
   Graph_relI_ME1dI1_L31002->SetFillStyle(0);
   Graph_relI_ME1dI1_L31002->SetLineStyle(0);
   Graph_relI_ME1dI1_L31002->SetMarkerStyle(20);
   Graph_relI_ME1dI1_L31002->GetXaxis()->SetNdivisions(505);
   Graph_relI_ME1dI1_L31002->GetXaxis()->SetLabelFont(42);
   Graph_relI_ME1dI1_L31002->GetXaxis()->SetLabelOffset(0.007);
   Graph_relI_ME1dI1_L31002->GetXaxis()->SetLabelSize(0.05);
   Graph_relI_ME1dI1_L31002->GetXaxis()->SetTitleSize(0.06);
   Graph_relI_ME1dI1_L31002->GetXaxis()->SetTitleFont(42);
   Graph_relI_ME1dI1_L31002->GetYaxis()->SetNdivisions(505);
   Graph_relI_ME1dI1_L31002->GetYaxis()->SetLabelFont(42);
   Graph_relI_ME1dI1_L31002->GetYaxis()->SetLabelOffset(0.007);
   Graph_relI_ME1dI1_L31002->GetYaxis()->SetLabelSize(0.05);
   Graph_relI_ME1dI1_L31002->GetYaxis()->SetTitleSize(0.06);
   Graph_relI_ME1dI1_L31002->GetYaxis()->SetTitleOffset(1.25);
   Graph_relI_ME1dI1_L31002->GetYaxis()->SetTitleFont(42);
   Graph_relI_ME1dI1_L31002->GetZaxis()->SetNdivisions(505);
   Graph_relI_ME1dI1_L31002->GetZaxis()->SetLabelFont(42);
   Graph_relI_ME1dI1_L31002->GetZaxis()->SetLabelOffset(0.007);
   Graph_relI_ME1dI1_L31002->GetZaxis()->SetLabelSize(0.05);
   Graph_relI_ME1dI1_L31002->GetZaxis()->SetTitleSize(0.06);
   Graph_relI_ME1dI1_L31002->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_relI_ME1dI1_L31002);
   
   gre->Draw("ps");
   
   Double_t relI_ME11_L4_fx1003[77] = {
   2.362894,
   5.146731,
   8.677283,
   14.62665,
   24.61029,
   32.05801,
   33.6714,
   38.74127,
   40.40048,
   47.86864,
   50.94259,
   52.68738,
   59.43014,
   62.69608,
   63.43046,
   67.63964,
   70.27712,
   73.80062,
   83.09959,
   113.2732,
   114.1608,
   121.3439,
   124.7864,
   126.9107,
   129.3835,
   135.6164,
   136.4528,
   140.7706,
   142.397,
   142.7824,
   143.5836,
   150.3945,
   152.3572,
   155.2885,
   156.5486,
   159.0907,
   159.9944,
   166.0897,
   171.2794,
   173.1347,
   180.051,
   184.2516,
   187.2563,
   193.5525,
   193.7259,
   194.8045,
   200.1376,
   205.8323,
   207.4543,
   211.6609,
   214.7788,
   222.0638,
   222.5167,
   224.4433,
   227.9533,
   228.5925,
   228.8852,
   229.5719,
   233.1883,
   235.2734,
   241.7202,
   243.4934,
   243.8268,
   250.1803,
   250.8601,
   252.0712,
   257.123,
   259.518,
   260.7698,
   264.4514,
   267.054,
   271.2835,
   272.1236,
   277.3676,
   277.8563,
   279.0993,
   279.9482};
   Double_t relI_ME11_L4_fy1003[77] = {
   1.045887,
   1.041002,
   1.042093,
   1.054449,
   1.049723,
   1.052118,
   1.04342,
   1.05024,
   1.045809,
   1.049272,
   1.048429,
   1.037256,
   1.035005,
   1.032965,
   1.035938,
   1.043329,
   1.032389,
   1.045373,
   1.037967,
   1.045471,
   1.042011,
   1.037874,
   1.032879,
   1.032484,
   1.005717,
   1.039596,
   1.047179,
   1.057304,
   1.060856,
   1.033683,
   1.03447,
   1.056205,
   1.045757,
   1.034377,
   1.048778,
   1.041853,
   1.031022,
   1.041457,
   1.051127,
   1.034632,
   1.04599,
   1.051261,
   1.041651,
   1.046083,
   1.031214,
   1.039793,
   1.024444,
   1.045006,
   1.040544,
   1.045076,
   1.041529,
   1.042133,
   1.040887,
   1.041109,
   1.033112,
   1.028606,
   1.041162,
   1.038042,
   1.040827,
   1.03899,
   1.050638,
   1.042816,
   1.023707,
   1.048157,
   1.052896,
   1.043961,
   1.032716,
   1.030083,
   1.004434,
   1.036416,
   1.026871,
   1.029637,
   1.030381,
   1.047743,
   1.030613,
   1.046184,
   1.037437};
   Double_t relI_ME11_L4_fex1003[77] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t relI_ME11_L4_fey1003[77] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   gre = new TGraphErrors(77,relI_ME11_L4_fx1003,relI_ME11_L4_fy1003,relI_ME11_L4_fex1003,relI_ME11_L4_fey1003);
   gre->SetName("relI_ME11_L4");
   gre->SetTitle("relI_ME11_L4");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#009900");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(22);
   gre->SetMarkerSize(markerSize);
   
   TH1F *Graph_relI_ME1dI1_L41003 = new TH1F("Graph_relI_ME1dI1_L41003","relI_ME11_L4",100,0,307.7068);
   Graph_relI_ME1dI1_L41003->SetMinimum(0.998792);
   Graph_relI_ME1dI1_L41003->SetMaximum(1.066498);
   Graph_relI_ME1dI1_L41003->SetDirectory(0);
   Graph_relI_ME1dI1_L41003->SetStats(0);
   Graph_relI_ME1dI1_L41003->SetFillColor(1);
   Graph_relI_ME1dI1_L41003->SetFillStyle(0);
   Graph_relI_ME1dI1_L41003->SetLineStyle(0);
   Graph_relI_ME1dI1_L41003->SetMarkerStyle(20);
   Graph_relI_ME1dI1_L41003->GetXaxis()->SetNdivisions(505);
   Graph_relI_ME1dI1_L41003->GetXaxis()->SetLabelFont(42);
   Graph_relI_ME1dI1_L41003->GetXaxis()->SetLabelOffset(0.007);
   Graph_relI_ME1dI1_L41003->GetXaxis()->SetLabelSize(0.05);
   Graph_relI_ME1dI1_L41003->GetXaxis()->SetTitleSize(0.06);
   Graph_relI_ME1dI1_L41003->GetXaxis()->SetTitleFont(42);
   Graph_relI_ME1dI1_L41003->GetYaxis()->SetNdivisions(505);
   Graph_relI_ME1dI1_L41003->GetYaxis()->SetLabelFont(42);
   Graph_relI_ME1dI1_L41003->GetYaxis()->SetLabelOffset(0.007);
   Graph_relI_ME1dI1_L41003->GetYaxis()->SetLabelSize(0.05);
   Graph_relI_ME1dI1_L41003->GetYaxis()->SetTitleSize(0.06);
   Graph_relI_ME1dI1_L41003->GetYaxis()->SetTitleOffset(1.25);
   Graph_relI_ME1dI1_L41003->GetYaxis()->SetTitleFont(42);
   Graph_relI_ME1dI1_L41003->GetZaxis()->SetNdivisions(505);
   Graph_relI_ME1dI1_L41003->GetZaxis()->SetLabelFont(42);
   Graph_relI_ME1dI1_L41003->GetZaxis()->SetLabelOffset(0.007);
   Graph_relI_ME1dI1_L41003->GetZaxis()->SetLabelSize(0.05);
   Graph_relI_ME1dI1_L41003->GetZaxis()->SetTitleSize(0.06);
   Graph_relI_ME1dI1_L41003->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_relI_ME1dI1_L41003);
   
   gre->Draw("ps");
   
   Double_t relI_ME11_L5_fx1004[77] = {
   2.362894,
   5.146731,
   8.677283,
   14.62665,
   24.61029,
   32.05801,
   33.6714,
   38.74127,
   40.40048,
   47.86864,
   50.94259,
   52.68738,
   59.43014,
   62.69608,
   63.43046,
   67.63964,
   70.27712,
   73.80062,
   83.09959,
   113.2732,
   114.1608,
   121.3439,
   124.7864,
   126.9107,
   129.3835,
   135.6164,
   136.4528,
   140.7706,
   142.397,
   142.7824,
   143.5836,
   150.3945,
   152.3572,
   155.2885,
   156.5486,
   159.0907,
   159.9944,
   166.0897,
   171.2794,
   173.1347,
   180.051,
   184.2516,
   187.2563,
   193.5525,
   193.7259,
   194.8045,
   200.1376,
   205.8323,
   207.4543,
   211.6609,
   214.7788,
   222.0638,
   222.5167,
   224.4433,
   227.9533,
   228.5925,
   228.8852,
   229.5719,
   233.1883,
   235.2734,
   241.7202,
   243.4934,
   243.8268,
   250.1803,
   250.8601,
   252.0712,
   257.123,
   259.518,
   260.7698,
   264.4514,
   267.054,
   271.2835,
   272.1236,
   277.3676,
   277.8563,
   279.0993,
   279.9482};
   Double_t relI_ME11_L5_fy1004[77] = {
   1.019123,
   1.012465,
   1.011767,
   1.025233,
   1.023127,
   1.024715,
   1.019122,
   1.023119,
   1.022278,
   1.022225,
   1.024598,
   1.014707,
   1.013208,
   1.013388,
   1.014645,
   1.020694,
   1.010907,
   1.016358,
   1.018798,
   1.02759,
   1.023285,
   1.020751,
   1.016554,
   1.014924,
   0.9963414,
   1.02345,
   1.024575,
   1.031543,
   1.033189,
   1.016018,
   1.015577,
   1.026562,
   1.01752,
   1.015243,
   1.023329,
   1.019385,
   1.009833,
   1.02031,
   1.024222,
   1.01658,
   1.021616,
   1.023005,
   1.018187,
   1.021074,
   1.014181,
   1.015068,
   1.00947,
   1.021087,
   1.014246,
   1.019805,
   1.016247,
   1.017017,
   1.016659,
   1.017755,
   1.005311,
   1.011091,
   1.003919,
   1.004948,
   1.016371,
   1.02014,
   1.029232,
   1.021995,
   1.002554,
   1.025053,
   1.030117,
   1.023552,
   1.021518,
   1.014669,
   0.9908729,
   1.017305,
   1.006863,
   1.019734,
   1.021011,
   1.02986,
   1.015924,
   1.025306,
   1.020755};
   Double_t relI_ME11_L5_fex1004[77] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t relI_ME11_L5_fey1004[77] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   gre = new TGraphErrors(77,relI_ME11_L5_fx1004,relI_ME11_L5_fy1004,relI_ME11_L5_fex1004,relI_ME11_L5_fey1004);
   gre->SetName("relI_ME11_L5");
   gre->SetTitle("relI_ME11_L5");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#6600cc");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(23);
   gre->SetMarkerSize(markerSize);
   
   TH1F *Graph_relI_ME1dI1_L51004 = new TH1F("Graph_relI_ME1dI1_L51004","relI_ME11_L5",100,0,307.7068);
   Graph_relI_ME1dI1_L51004->SetMinimum(0.9866412);
   Graph_relI_ME1dI1_L51004->SetMaximum(1.037421);
   Graph_relI_ME1dI1_L51004->SetDirectory(0);
   Graph_relI_ME1dI1_L51004->SetStats(0);
   Graph_relI_ME1dI1_L51004->SetFillColor(1);
   Graph_relI_ME1dI1_L51004->SetFillStyle(0);
   Graph_relI_ME1dI1_L51004->SetLineStyle(0);
   Graph_relI_ME1dI1_L51004->SetMarkerStyle(20);
   Graph_relI_ME1dI1_L51004->GetXaxis()->SetNdivisions(505);
   Graph_relI_ME1dI1_L51004->GetXaxis()->SetLabelFont(42);
   Graph_relI_ME1dI1_L51004->GetXaxis()->SetLabelOffset(0.007);
   Graph_relI_ME1dI1_L51004->GetXaxis()->SetLabelSize(0.05);
   Graph_relI_ME1dI1_L51004->GetXaxis()->SetTitleSize(0.06);
   Graph_relI_ME1dI1_L51004->GetXaxis()->SetTitleFont(42);
   Graph_relI_ME1dI1_L51004->GetYaxis()->SetNdivisions(505);
   Graph_relI_ME1dI1_L51004->GetYaxis()->SetLabelFont(42);
   Graph_relI_ME1dI1_L51004->GetYaxis()->SetLabelOffset(0.007);
   Graph_relI_ME1dI1_L51004->GetYaxis()->SetLabelSize(0.05);
   Graph_relI_ME1dI1_L51004->GetYaxis()->SetTitleSize(0.06);
   Graph_relI_ME1dI1_L51004->GetYaxis()->SetTitleOffset(1.25);
   Graph_relI_ME1dI1_L51004->GetYaxis()->SetTitleFont(42);
   Graph_relI_ME1dI1_L51004->GetZaxis()->SetNdivisions(505);
   Graph_relI_ME1dI1_L51004->GetZaxis()->SetLabelFont(42);
   Graph_relI_ME1dI1_L51004->GetZaxis()->SetLabelOffset(0.007);
   Graph_relI_ME1dI1_L51004->GetZaxis()->SetLabelSize(0.05);
   Graph_relI_ME1dI1_L51004->GetZaxis()->SetTitleSize(0.06);
   Graph_relI_ME1dI1_L51004->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_relI_ME1dI1_L51004);
   
   gre->Draw("ps");
   
   TLegend *leg = new TLegend(0.69, 0.2,0.89,0.5,NULL,"brNDC");
   leg->SetBorderSize(0);
   //leg->SetTextFont(62);
   leg->SetTextSize(legSize);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("relI_ME11_L2","Layer 2","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#66ccff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(markerSize);
   //   entry->SetTextFont(62);
   entry=leg->AddEntry("relI_ME11_L3","Layer 3","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#cc3333");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(markerSize);
   //   entry->SetTextFont(62);
   entry=leg->AddEntry("relI_ME11_L4","Layer 4","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#009900");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(markerSize);
   //   entry->SetTextFont(62);
   entry=leg->AddEntry("relI_ME11_L5","Layer 5","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#6600cc");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(markerSize);
   //   entry->SetTextFont(62);
   leg->Draw();
   cME11->Modified();
   cME11->cd();
   cME11->SetSelected(cME11);

   TString pname = "plots/gifpp_me11_gas_gain.pdf";
   cME11->SaveAs(pname);

}

void plotDarkRate(){
//=========Macro generated from canvas: c1/his1
//=========  (Fri Jun  2 11:19:30 2017) by ROOT version6.08/06
   TCanvas *c1 = new TCanvas("c1", "his1");
   c1->Range(-66.7403,-142.6829,350.3866,954.8781);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetTickx(1);
   c1->SetTicky(1);
   
   Double_t Graph0_fx1[38] = {
   17.42,
   25.23,
   40.4,
   47.87,
   50.94,
   59.38,
   59.39,
   59.43,
   77.94,
   89.81,
   95.44,
   104.24,
   113.27,
   120.9,
   135.62,
   142.4,
   152.13,
   154.06,
   154.07,
   154.42,
   154.58,
   158.07,
   173.13,
   180.46,
   186.99,
   193.55,
   200.14,
   207.45,
   214.78,
   222.52,
   228.59,
   235.27,
   241.72,
   250.86,
   256.42,
   266.13,
   271.28,
   277.37};
   Double_t Graph0_fy1[38] = {
   155.1667,
   149.9667,
   156.1667,
   159.1,
   155.8667,
   160.2,
   158.2,
   151.3333,
   156.6,
   154.0333,
   153.3,
   155.2667,
   153.5667,
   153.9667,
   162.8333,
   152.8333,
   157.8,
   141.6667,
   143.3,
   138.1667,
   145.7333,
   181.9,
   154.9,
   153.3667,
   155.5333,
   147.2667,
   147.7667,
   151.3,
   149.3667,
   148.7,
   149.2667,
   259.9667,
   150.0667,
   179.8,
   152.2,
   154.3,
   149.2333,
   165.4};
   TGraph *graph = new TGraph(38,Graph0_fx1,Graph0_fy1);
   graph->SetName("Graph0");
   graph->SetTitle("; Accumulated charge [mC/cm]; Dark rate [Hz]");
   graph->SetFillColor(1);

   Int_t ci;      // for color index setting
   //   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#cccc00");
   graph->SetLineColor(ci);

   ci = TColor::GetColor("#cccc00");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(22);
   
   TH1F *Graph_Graph_Graph381 = new TH1F("Graph_Graph_Graph381","",100,0,303.365);
   Graph_Graph_Graph381->SetMinimum(0);
   Graph_Graph_Graph381->SetMaximum(900);
   Graph_Graph_Graph381->SetDirectory(0);
   Graph_Graph_Graph381->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph381->SetLineColor(ci);
   Graph_Graph_Graph381->SetLineStyle(0);
   Graph_Graph_Graph381->SetMarkerStyle(20);
   Graph_Graph_Graph381->GetXaxis()->CenterTitle(true);
   Graph_Graph_Graph381->GetXaxis()->SetTitle("Accumulated charge [mC/cm]");
   Graph_Graph_Graph381->GetXaxis()->SetRange(1,99);
   Graph_Graph_Graph381->GetYaxis()->CenterTitle(true);
   Graph_Graph_Graph381->GetYaxis()->SetNdivisions(505);
   Graph_Graph_Graph381->GetYaxis()->SetTitle("Dark rate [Hz]");
   Graph_Graph_Graph381->SetTitle("Dark rate vs Charge for ME1/1");
   graph->SetHistogram(Graph_Graph_Graph381);
   
   graph->Draw("alp");
   
   Double_t Graph1_fx2[38] = {
   17.42,
   25.23,
   40.4,
   47.87,
   50.94,
   59.38,
   59.39,
   59.43,
   77.94,
   89.81,
   95.44,
   104.24,
   113.27,
   120.9,
   135.62,
   142.4,
   152.13,
   154.06,
   154.07,
   154.42,
   154.58,
   158.07,
   173.13,
   180.46,
   186.99,
   193.55,
   200.14,
   207.45,
   214.78,
   222.52,
   228.59,
   235.27,
   241.72,
   250.86,
   256.42,
   266.13,
   271.28,
   277.37};
   Double_t Graph1_fy2[38] = {
   157.5667,
   160.7667,
   157.7333,
   158.9,
   158.7333,
   189.1,
   160.1,
   155.6,
   157.5333,
   159.2333,
   156,
   157.9,
   160.7667,
   160.1667,
   162.4667,
   157.3667,
   159.2667,
   146.7,
   144.7,
   143.0333,
   148.5,
   181.7333,
   155.5,
   154.0667,
   154.1,
   151.1333,
   157.0667,
   154.6333,
   154.4333,
   150.6333,
   151.4333,
   157.2,
   154.0333,
   149.4,
   152.0333,
   155.1667,
   153.5667,
   169.1};
   graph = new TGraph(38,Graph1_fx2,Graph1_fy2);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#66ccff");
   graph->SetLineColor(ci);

   ci = TColor::GetColor("#66ccff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(22);
   
   TH1F *Graph_Graph_Graph392 = new TH1F("Graph_Graph_Graph392","Graph",100,0,303.365);
   Graph_Graph_Graph392->SetMinimum(138.4267);
   Graph_Graph_Graph392->SetMaximum(193.7067);
   Graph_Graph_Graph392->SetDirectory(0);
   Graph_Graph_Graph392->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph392->SetLineColor(ci);
   Graph_Graph_Graph392->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph392->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph392->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph392->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph392->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph392->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph392->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph392->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph392->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph392->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph392->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph392->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph_Graph392);
   
   graph->Draw("lp ");
   
   Double_t Graph2_fx3[38] = {
   17.42,
   25.23,
   40.4,
   47.87,
   50.94,
   59.38,
   59.39,
   59.43,
   77.94,
   89.81,
   95.44,
   104.24,
   113.27,
   120.9,
   135.62,
   142.4,
   152.13,
   154.06,
   154.07,
   154.42,
   154.58,
   158.07,
   173.13,
   180.46,
   186.99,
   193.55,
   200.14,
   207.45,
   214.78,
   222.52,
   228.59,
   235.27,
   241.72,
   250.86,
   256.42,
   266.13,
   271.28,
   277.37};
   Double_t Graph2_fy3[38] = {
   160.7667,
   164,
   160.9,
   160.6,
   161.5333,
   164.3,
   161.3667,
   160.4667,
   165.9333,
   155.1333,
   165,
   157.7,
   167.0333,
   165.5333,
   170.7667,
   158.5,
   166.6667,
   148.1333,
   147.9333,
   143.0333,
   150.4667,
   187.3667,
   162.2,
   157.5333,
   154.1,
   152.5333,
   159.0333,
   159.0667,
   157.8333,
   154.9333,
   155.6667,
   158.1667,
   155.6,
   159.8,
   160.6333,
   158.1,
   154.7333,
   174.8};
   graph = new TGraph(38,Graph2_fx3,Graph2_fy3);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#cc3333");
   graph->SetLineColor(ci);

   ci = TColor::GetColor("#cc3333");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(22);
   
   TH1F *Graph_Graph_Graph403 = new TH1F("Graph_Graph_Graph403","Graph",100,0,303.365);
   Graph_Graph_Graph403->SetMinimum(138.6);
   Graph_Graph_Graph403->SetMaximum(191.8);
   Graph_Graph_Graph403->SetDirectory(0);
   Graph_Graph_Graph403->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph403->SetLineColor(ci);
   Graph_Graph_Graph403->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph403->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph403->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph403->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph403->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph403->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph403->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph403->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph403->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph403->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph403->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph403->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph_Graph403);
   
   graph->Draw("lp ");
   
   Double_t Graph3_fx4[38] = {
   17.42,
   25.23,
   40.4,
   47.87,
   50.94,
   59.38,
   59.39,
   59.43,
   77.94,
   89.81,
   95.44,
   104.24,
   113.27,
   120.9,
   135.62,
   142.4,
   152.13,
   154.06,
   154.07,
   154.42,
   154.58,
   158.07,
   173.13,
   180.46,
   186.99,
   193.55,
   200.14,
   207.45,
   214.78,
   222.52,
   228.59,
   235.27,
   241.72,
   250.86,
   256.42,
   266.13,
   271.28,
   277.37};
   Double_t Graph3_fy4[38] = {
   159.4333,
   164.7667,
   162.7667,
   162.8,
   168.9,
   177.5,
   164,
   158.7667,
   166.6,
   162.1667,
   173.2333,
   170.2333,
   175.0667,
   178.7,
   183.1333,
   175.4,
   160.3333,
   154,
   152.8,
   147.0333,
   151.6333,
   192.1333,
   170.4,
   172.1,
   170.9667,
   158.0667,
   163,
   166.9667,
   168.6,
   178.2,
   171.0667,
   182.9333,
   177.8,
   174.1,
   164.5667,
   167.7,
   170.1,
   182.5333};
   graph = new TGraph(38,Graph3_fx4,Graph3_fy4);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#009900");
   graph->SetLineColor(ci);

   ci = TColor::GetColor("#009900");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(22);
   
   TH1F *Graph_Graph_Graph414 = new TH1F("Graph_Graph_Graph414","Graph",100,0,303.365);
   Graph_Graph_Graph414->SetMinimum(142.5233);
   Graph_Graph_Graph414->SetMaximum(196.6433);
   Graph_Graph_Graph414->SetDirectory(0);
   Graph_Graph_Graph414->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph414->SetLineColor(ci);
   Graph_Graph_Graph414->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph414->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph414->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph414->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph414->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph414->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph414->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph414->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph414->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph414->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph414->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph414->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph_Graph414);
   
   graph->Draw("lp ");
   
   Double_t Graph4_fx5[38] = {
   17.42,
   25.23,
   40.4,
   47.87,
   50.94,
   59.38,
   59.39,
   59.43,
   77.94,
   89.81,
   95.44,
   104.24,
   113.27,
   120.9,
   135.62,
   142.4,
   152.13,
   154.06,
   154.07,
   154.42,
   154.58,
   158.07,
   173.13,
   180.46,
   186.99,
   193.55,
   200.14,
   207.45,
   214.78,
   222.52,
   228.59,
   235.27,
   241.72,
   250.86,
   256.42,
   266.13,
   271.28,
   277.37};
   Double_t Graph4_fy5[38] = {
   166.2333,
   166.3,
   164.8,
   164.9667,
   165.5,
   178.3667,
   188.1333,
   166.0667,
   169.9667,
   170.3667,
   175.6333,
   195.9,
   201.8333,
   189.7,
   181.7333,
   229.3333,
   224.3333,
   162.7667,
   170.7,
   176.6333,
   187.3333,
   203.2667,
   243.6,
   312.7,
   435.7667,
   267.9333,
   259.9,
   263.8333,
   266.6,
   304.8,
   350.5,
   367.4,
   307.4,
   298.6333,
   516.6667,
   265.7,
   393.3667,
   846.3};
   graph = new TGraph(38,Graph4_fx5,Graph4_fy5);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#6600cc");
   graph->SetLineColor(ci);

   ci = TColor::GetColor("#6600cc");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(22);
   
   TH1F *Graph_Graph_Graph425 = new TH1F("Graph_Graph_Graph425","Graph",100,0,303.365);
   Graph_Graph_Graph425->SetMinimum(94.41333);
   Graph_Graph_Graph425->SetMaximum(914.6533);
   Graph_Graph_Graph425->SetDirectory(0);
   Graph_Graph_Graph425->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph425->SetLineColor(ci);
   Graph_Graph_Graph425->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph425->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph425->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph425->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph425->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph425->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph425->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph425->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph425->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph425->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph425->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph425->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph_Graph425);
   
   graph->Draw("lp ");
   
   Double_t Graph5_fx6[38] = {
   17.42,
   25.23,
   40.4,
   47.87,
   50.94,
   59.38,
   59.39,
   59.43,
   77.94,
   89.81,
   95.44,
   104.24,
   113.27,
   120.9,
   135.62,
   142.4,
   152.13,
   154.06,
   154.07,
   154.42,
   154.58,
   158.07,
   173.13,
   180.46,
   186.99,
   193.55,
   200.14,
   207.45,
   214.78,
   222.52,
   228.59,
   235.27,
   241.72,
   250.86,
   256.42,
   266.13,
   271.28,
   277.37};
   Double_t Graph5_fy6[38] = {
   171.2667,
   156.5,
   157.6333,
   158.6667,
   158.1,
   159.8333,
   164.7667,
   151.2,
   160.7667,
   162.6333,
   156.5,
   164.9,
   166.1667,
   157.7333,
   168.6333,
   159.5667,
   162.9333,
   150.5,
   149.1333,
   144.9,
   151.9667,
   182.3,
   154.1333,
   149.6333,
   154.5,
   148.2,
   155.5333,
   152.9,
   157.1333,
   153.4,
   151.5,
   158.3667,
   150.4,
   155.8333,
   154.1333,
   153.7333,
   152.7,
   171.5};
   graph = new TGraph(38,Graph5_fx6,Graph5_fy6);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);

   ci = TColor::GetColor("#cc9933");
   graph->SetLineColor(ci);

   ci = TColor::GetColor("#cc9933");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(22);
   
   TH1F *Graph_Graph_Graph436 = new TH1F("Graph_Graph_Graph436","Graph",100,0,303.365);
   Graph_Graph_Graph436->SetMinimum(141.16);
   Graph_Graph_Graph436->SetMaximum(186.04);
   Graph_Graph_Graph436->SetDirectory(0);
   Graph_Graph_Graph436->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph436->SetLineColor(ci);
   Graph_Graph_Graph436->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph436->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph436->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph436->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph436->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph436->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph436->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph436->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph436->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph436->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph436->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph436->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph_Graph436);
   
   graph->Draw("lp ");
   
   TLegend *leg = new TLegend(0.2,0.46,0.4,0.84,NULL,"brNDC");
   leg->SetBorderSize(0);
   //leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->SetTextSize(legSize);
   TLegendEntry *entry=leg->AddEntry("Graph0","Layer 1","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#cccc00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(markerSize);
   //   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph1","Layer 2","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#66ccff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(markerSize);
   //   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph2","Layer 3","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#cc3333");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(markerSize);
   //   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph3","Layer 4","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#009900");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(markerSize);
   //   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph4","Layer 5","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#6600cc");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(markerSize);
   //   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph5","Layer 6","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#cc9933");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(markerSize);
   //   entry->SetTextFont(62);
   leg->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
   
   TString pname = "plots/gifpp_me11_dark_rate.pdf";
   c1->SaveAs(pname);

}

void plotDarkRateProblem(){
//=========Macro generated from canvas: c1/graph
//=========  (Fri Jun  2 11:27:43 2017) by ROOT version6.08/06
   TCanvas *c1 = new TCanvas("canvasProblem", "graph");
   c1->Range(-71.79638,-0.002397105,357.9707,0.01604216);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetTickx(1);
   c1->SetTicky(1);
   
   Double_t Graph0_fx1[38] = {
   17.42,
   25.23,
   40.4,
   47.87,
   50.94,
   59.38,
   59.39,
   59.43,
   77.94,
   89.81,
   95.44,
   104.24,
   113.27,
   120.9,
   135.62,
   142.4,
   152.13,
   154.06,
   154.07,
   154.42,
   154.58,
   158.07,
   173.13,
   180.46,
   185.52,
   193.55,
   200.14,
   207.45,
   214.78,
   222.52,
   228.59,
   235.27,
   241.72,
   250.86,
   256.42,
   266.13,
   271.28,
   277.37};
   Double_t Graph0_fy1[38] = {
   0.0001062581,
   0.00011865,
   0.0001736865,
   0.0001273858,
   0.000131051,
   0.0002740694,
   0.0001175878,
   0.0003832853,
   0.0003519145,
   0.0004271639,
   0.0006866126,
   0.001546332,
   0.001547412,
   0.001151898,
   0.0005939942,
   0.002803539,
   0.002652112,
   0.0005352006,
   0.001106332,
   0.001398457,
   0.001744453,
   0.000482539,
   0.003494231,
   0.006109568,
   0.009463983,
   0.005038326,
   0.004311398,
   0.004387624,
   0.004459587,
   0.005843096,
   0.007313027,
   0.006674027,
   0.005802442,
   0.005489995,
   0.01158099,
   0.001320907,
   0.008542106,
   0.01375529};
   TGraph *graph = new TGraph(38,Graph0_fx1,Graph0_fy1);
   graph->SetName("Graph0");
   graph->SetTitle("; Accumulated charge [mC/cm]; Singles normalized rate [Hz/cm]");
   graph->SetFillColor(1);
   graph->SetLineColor(2);
   graph->SetMarkerColor(2);
   graph->SetMarkerStyle(22);
   
   TH1F *Graph_Graph_Graph501 = new TH1F("Graph_Graph_Graph501","",100,0,303.365);
   Graph_Graph_Graph501->SetMinimum(0);
   Graph_Graph_Graph501->SetMaximum(0.0151202);
   Graph_Graph_Graph501->SetDirectory(0);
   Graph_Graph_Graph501->SetStats(0);

   Int_t ci;      // for color index setting
   //   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph501->SetLineColor(ci);
   Graph_Graph_Graph501->GetXaxis()->CenterTitle(true);
   Graph_Graph_Graph501->GetXaxis()->SetTitle(" Accumulated charge [mC/cm]");
   Graph_Graph_Graph501->GetXaxis()->SetRange(0,101);
   Graph_Graph_Graph501->GetYaxis()->CenterTitle(true);
   Graph_Graph_Graph501->GetYaxis()->SetTitle("Normalized dark rate [Hz/cm]");
   Graph_Graph_Graph501->GetYaxis()->SetDecimals();
   Graph_Graph_Graph501->GetYaxis()->SetTitleOffset(1.76);
   c1->SetLeftMargin(0.19);
   Graph_Graph_Graph501->SetTitle("Dark rates for problematic wire groups in ME1/1");
   graph->SetHistogram(Graph_Graph_Graph501);
   
   graph->Draw("alp");
   
   Double_t Graph1_fx2[38] = {
   17.42,
   25.23,
   40.4,
   47.87,
   50.94,
   59.38,
   59.39,
   59.43,
   77.94,
   89.81,
   95.44,
   104.24,
   113.27,
   120.9,
   135.62,
   142.4,
   152.13,
   154.06,
   154.07,
   154.42,
   154.58,
   158.07,
   173.13,
   180.46,
   185.52,
   193.55,
   200.14,
   207.45,
   214.78,
   222.52,
   228.59,
   235.27,
   241.72,
   250.86,
   256.42,
   266.13,
   271.28,
   277.37};
   Double_t Graph1_fy2[38] = {
   0.0001541141,
   0.000111315,
   0.0002475807,
   0.0001425542,
   0.0001264287,
   0.0001881214,
   0.0001089787,
   0.0001154841,
   0.0001303098,
   0.0001146759,
   9.179934e-05,
   0.0001139284,
   0.0001246004,
   0.0001184435,
   9.458799e-05,
   0.0001054216,
   9.915431e-05,
   6.634546e-05,
   0.0001099132,
   0.0001113107,
   0.0001189908,
   7.437811e-05,
   0.0001091447,
   0.0001101643,
   7.829553e-05,
   9.603153e-05,
   9.202625e-05,
   0.0001099664,
   8.677653e-05,
   0.0001059153,
   7.610436e-05,
   9.66564e-05,
   0.0001566319,
   9.648693e-05,
   9.161518e-05,
   2.580518e-05,
   8.942546e-05,
   5.592273e-05};
   graph = new TGraph(38,Graph1_fx2,Graph1_fy2);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);
   graph->SetLineColor(3);
   graph->SetMarkerColor(3);
   graph->SetMarkerStyle(23);
   
   TH1F *Graph_Graph_Graph512 = new TH1F("Graph_Graph_Graph512","Graph",100,0,303.365);
   Graph_Graph_Graph512->SetMinimum(3.627628e-06);
   Graph_Graph_Graph512->SetMaximum(0.0002697582);
   Graph_Graph_Graph512->SetDirectory(0);
   Graph_Graph_Graph512->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph512->SetLineColor(ci);
   Graph_Graph_Graph512->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph512->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph512->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph512->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph512->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph512->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph512->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph512->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph512->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph512->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph512->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph512->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph_Graph512);
   
   graph->Draw("lp ");
   
   Double_t Graph2_fx3[38] = {
   17.42,
   25.23,
   40.4,
   47.87,
   50.94,
   59.38,
   59.39,
   59.43,
   77.94,
   89.81,
   95.44,
   104.24,
   113.27,
   120.9,
   135.62,
   142.4,
   152.13,
   154.06,
   154.07,
   154.42,
   154.58,
   158.07,
   173.13,
   180.46,
   185.52,
   193.55,
   200.14,
   207.45,
   214.78,
   222.52,
   228.59,
   235.27,
   241.72,
   250.86,
   256.42,
   266.13,
   271.28,
   277.37};
   Double_t Graph2_fy3[38] = {
   0.0001400133,
   0.0001386195,
   0.0001750282,
   0.0001503649,
   0.0001702995,
   0.0001836784,
   0.0001983086,
   0.0001688356,
   0.0002470049,
   0.0003370222,
   0.0003841339,
   0.0004701146,
   0.0007887728,
   0.0007523964,
   0.0008283081,
   0.0006883217,
   0.000171885,
   0.0002091956,
   0.0002750032,
   0.0002482257,
   0.0002072983,
   0.0001706283,
   0.0004770265,
   0.000497579,
   0.0005714632,
   0.0002013083,
   0.0003488754,
   0.0003848191,
   0.0004967694,
   0.0007025833,
   0.0006385985,
   0.0006740054,
   0.0008401532,
   0.0007893085,
   0.0001130085,
   0.0001702733,
   0.0006066596,
   0.0002966522};
   graph = new TGraph(38,Graph2_fx3,Graph2_fy3);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);
   graph->SetLineColor(2);
   graph->SetMarkerColor(2);
   graph->SetMarkerStyle(26);
   
   TH1F *Graph_Graph_Graph523 = new TH1F("Graph_Graph_Graph523","Graph",100,0,303.365);
   Graph_Graph_Graph523->SetMinimum(4.029407e-05);
   Graph_Graph_Graph523->SetMaximum(0.0009128677);
   Graph_Graph_Graph523->SetDirectory(0);
   Graph_Graph_Graph523->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph523->SetLineColor(ci);
   Graph_Graph_Graph523->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph523->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph523->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph523->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph523->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph523->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph523->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph523->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph523->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph523->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph523->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph523->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph_Graph523);
   
   graph->Draw("lp ");
   
   Double_t Graph3_fx4[38] = {
   17.42,
   25.23,
   40.4,
   47.87,
   50.94,
   59.38,
   59.39,
   59.43,
   77.94,
   89.81,
   95.44,
   104.24,
   113.27,
   120.9,
   135.62,
   142.4,
   152.13,
   154.06,
   154.07,
   154.42,
   154.58,
   158.07,
   173.13,
   180.46,
   185.52,
   193.55,
   200.14,
   207.45,
   214.78,
   222.52,
   228.59,
   235.27,
   241.72,
   250.86,
   256.42,
   266.13,
   271.28,
   277.37};
   Double_t Graph3_fy4[38] = {
   0.000111005,
   0.000150243,
   0.0001209151,
   0.0001101595,
   0.0001503276,
   0.0001375983,
   0.0001368825,
   0.0001013628,
   0.0001378154,
   0.00010311,
   0.0001041907,
   0.0001357722,
   0.0001016949,
   0.0001036579,
   0.000118247,
   0.0001010048,
   0.0001109755,
   8.58902e-05,
   0.0001211338,
   0.0001486304,
   9.997402e-05,
   9.098377e-05,
   0.0001050606,
   0.0001192191,
   0.0001028285,
   9.319794e-05,
   9.062423e-05,
   0.0001402187,
   0.0001287945,
   0.0001042804,
   0.0001100613,
   7.594853e-05,
   0.0001005447,
   8.533331e-05,
   7.290673e-05,
   3.365986e-05,
   9.615451e-05,
   4.458108e-05};
   graph = new TGraph(38,Graph3_fx4,Graph3_fy4);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(1);
   graph->SetLineColor(3);
   graph->SetMarkerColor(3);
   graph->SetMarkerStyle(32);
   
   TH1F *Graph_Graph_Graph534 = new TH1F("Graph_Graph_Graph534","Graph",100,0,303.365);
   Graph_Graph_Graph534->SetMinimum(2.199309e-05);
   Graph_Graph_Graph534->SetMaximum(0.0001619944);
   Graph_Graph_Graph534->SetDirectory(0);
   Graph_Graph_Graph534->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph534->SetLineColor(ci);
   Graph_Graph_Graph534->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph534->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph534->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph534->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph534->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph534->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph534->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph534->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph534->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph534->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph534->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph534->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph_Graph534);
   
   graph->Draw("lp ");
   
   TLegend *leg = new TLegend(0.22,0.66,0.4,0.88,NULL,"brNDC");
   leg->SetBorderSize(0);
   //   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetTextSize(legSize*0.9);
   TLegendEntry *entry=leg->AddEntry("Graph2","Problematic wiregroup in layer 4","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(26);
   entry->SetMarkerSize(1);
   //   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph3","Normal wiregroup in layer 4","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(3);
   entry->SetMarkerStyle(32);
   entry->SetMarkerSize(1);
   //   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph0","Problematic wiregroup in layer 5","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1);
   //   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph1","Normal wiregroup in layer 5","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(3);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(1);
   //   entry->SetTextFont(62);
   leg->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
   
   TString pname = "plots/gifpp_me11_dark_rate_problem.pdf";
   c1->SaveAs(pname);
}

void plotDarkCurrent(){

  TCanvas *c1 = new TCanvas("canvasCurrent", "graph");
  c1->Range(-71.79638,-0.002397105,357.9707,0.01604216);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetTickx(1);
  c1->SetTicky(1);

  int colour[18] = {401, 866, 629, 418, 882, 795};

  float doseV[6] = {118, 187, 198, 220, 280, 289};
  float mean0[12] = {2.3,  0.1,  1.4,  0.1,  1.8,   0.1,  1.9,  0.1,  2.5,  0.1,  2.4,  0.1}; //ME21s1_AllLayer_2016_07_13
  float mean1[12] = {2.15, 0.07, 1.65, 0.05, -0.50, 0.00, 2.09, 0.05, 2.80, 0.13, 2.59, 0.10 };//R_Pi_me21s1_DC_16.11.16
  float mean2[12] = {2.61, 0.05, 2.32, 0.07,  2.43, 0.09, 2.42, 0.07, 2.62, 0.09, 2.59, 0.12}; //R_Pi_me21s1_DC_18.01.17
  float mean3[12] = {2.66, 0.09, 1.73, 0.09,  2.24, 0.10, 2.46, 0.09, 3.38, 0.15, 3.61, 0.16}; //R_Pi_me21s1_DC_8.02.17
  float mean4[12] = {2.72, 0.08, 2.16, 0.07,  2.79, 0.05, 2.97, 0.06, 4.39, 0.09, 3.76, 0.11}; //R_Pi_me21s1_DC_24.04.17
  float mean5[12] = {3.08, 0.08, 2.09, 0.04,  2.63, 0.05, 2.78, 0.07, 3.59, 0.06, 3.45, 0.13}; //R_Pi_me21s1_DC_31.05.17

  TGraphErrors *g[6];

  for(int i=0; i<6; i++){
    g[i] = new TGraphErrors();
    g[i]->SetMarkerColor(colour[i]);
    g[i]->SetLineColor(colour[i]);
    if(i<4) g[i]->SetMarkerStyle(20+i);
    else g[i]->SetMarkerStyle(29+i);
    g[i]->SetMarkerSize(markerSize*1.3);
  };
  g[1]->SetPoint(0, 35, 1.2); //ME21s1_Layer2_2016_03_03 err 0.1
  g[1]->SetPointError(0,  0, 0.1);  
  g[1]->SetPoint(1, 75, 1.1); //ME21s1_Layer2_2016_05_18 err 0.2
  g[1]->SetPointError(1,  0, 0.2);  

  g[2]->SetPoint(0, 52, 1.2);    //ME21s1_Layer3_2016_04_01 err 0.1
  g[2]->SetPointError(0,  0, 0.1);  
  g[2]->SetPoint(1, 75, 1.2);  //ME21s1_Layer3_2016_05_18 err 0.2
  g[2]->SetPointError(1,  0, 0.2);  

  g[3]->SetPoint(0, 75, 1.3);   //ME21s1_Layer4_2016_05_18 err 0.2
  g[3]->SetPointError(0,  0, 0.2);  

  g[4]->SetPoint(0, 75, 1.3);   //ME21s1_Layer5_2016_05_18 err ??? 0.2
  g[4]->SetPointError(0,  0, 0.2);  

  g[5]->SetPoint(0, 75, 1.6);  //ME21s1_Layer6_2016_05_18 err 0.2
  g[5]->SetPointError(0,  0, 0.2);
  
  int start[6]={0,2,2,1,1,1};
  int n=0;
  for(int i=0; i<6; i++){
    g[i]->SetPoint(start[i],      doseV[0], mean0[i*2]);
    g[i]->SetPointError(start[i], 0, mean0[i*2+1]);
    start[i]++; n++;
  };
  
  for(int i=0; i<6; i++){
    if(i!=2){
      g[i]->SetPoint(start[i],      doseV[1], mean1[i*2]);
      g[i]->SetPointError(start[i], 0, mean1[i*2+1]);
      start[i]++; n++;
    };
  };
  
  for(int i=0; i<6; i++){
    g[i]->SetPoint(start[i],      doseV[2], mean2[i*2]);
    g[i]->SetPointError(start[i], 0, mean2[i*2+1]);
    start[i]++; n++;
  };
  for(int i=0; i<6; i++){
    g[i]->SetPoint(start[i],      doseV[3], mean3[i*2]);
    g[i]->SetPointError(start[i], 0, mean3[i*2+1]);
    start[i]++; n++;
  };
  for(int i=0; i<6; i++){
    g[i]->SetPoint(start[i],      doseV[4], mean4[i*2]);
    g[i]->SetPointError(start[i], 0, mean4[i*2+1]);
    start[i]++; n++;
  };
  for(int i=0; i<6; i++){
    g[i]->SetPoint(start[i],      doseV[5], mean5[i*2]);
    g[i]->SetPointError(start[i], 0, mean5[i*2+1]);
    start[i]++; n++;
  };
 
  TH2F * h = new TH2F("h",";Accumulated charge [mC/cm]; Dark current [nA]", 300, 0, 300, 100, 0.0, 8.0);
  h->GetXaxis()->CenterTitle(true);
  h->GetYaxis()->CenterTitle(true);
  h->SetTitle("Dark current vs Charge for ME2/1s1");
  h->GetYaxis()->SetTitleOffset(1.15);
  h->Draw();
  TLegend * l1 = new TLegend(0.22, 0.46, 0.32, 0.85);  
  l1->SetTextSize(legSize);
  l1->SetBorderSize(0);
  l1->SetFillStyle(0);
  for(int i=0; i<6; i++){
    g[i]->Draw("Psames");
    TString la = "Layer "; la+=(i+1);
    l1->AddEntry(g[i], la.Data(),"P");
  }
  l1->Draw();

  TString pname = "plots/gifpp_me21_dark_current.pdf";
  c1->SaveAs(pname);

}
