#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include <unistd.h>
#include <getopt.h>

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
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
  TString lumi = "40";
  TString filename = "txt/limits/limits_TChiHH_lumi"+lumi+".txt";
  TString model = "TChiHH";
}

void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  GetOptions(argc, argv);

  //// Setting plot style
  PlotOpt opts("txt/plot_styles.txt", "Ratio");
  setPlotStyle(opts);
  gStyle->SetGridStyle(3);

  if(filename == "") ERROR("No input file provided");
  ifstream infile(filename);

  vector<double> vmx, vmy, vxsec, vexsec, vobs, vobsup, vobsdown;
  vector<double> vexp, vup, vdown, v2up, v2down, vsigobs, vsigexp, zeroes, ones;
  double maxy=-99., miny=1e99;
  
  string line_s;
  while(getline(infile, line_s)){
    istringstream iss(line_s);
    double pmx, pmy, pxsec, pexsec, pobs, pobsup, pobsdown, pexp, pup, pdown, p2up, p2down, sigobs, sigexp;
    iss >> pmx >> pmy >> pxsec >> pexsec >> pobs >> pobsup >> pobsdown 
	>> pexp >> pup >> pdown >> p2up >> p2down >> sigobs >> sigexp;
    if(pmx<200) continue;
    vmx.push_back(pmx);
    vmy.push_back(pmy);
    vxsec.push_back(pxsec);
    vexsec.push_back(pexsec);
    vobs.push_back(pobs);
    vobsup.push_back(pobsup);
    vobsdown.push_back(pobsdown);
    vexp.push_back(pexp);
    vup.push_back(pup-pexp);
    vdown.push_back(pexp-pdown);
    v2up.push_back(p2up-pexp);
    v2down.push_back(pexp-p2down);
    vsigobs.push_back(sigobs);
    vsigexp.push_back(sigexp);
    zeroes.push_back(0);
    ones.push_back(1);
    if(miny > min(vobs.back(), 1.)) miny = min(vobs.back(), 1.);
    if(maxy < max(vobs.back(), 1.)) maxy = max(vobs.back(), 1.);
  }
  infile.close();

  if(vmx.size() <= 0) ERROR("Need at least 1 model to draw limits");
  if(vmx.size() != vmy.size()
     || vmx.size() != vxsec.size()
     || vmx.size() != vexsec.size()
     || vmx.size() != vobs.size()
     || vmx.size() != vobsup.size()
     || vmx.size() != vobsdown.size()
     || vmx.size() != vexp.size()
     || vmx.size() != vup.size()
     || vmx.size() != vdown.size()
     || vmx.size() != v2up.size()
     || vmx.size() != v2down.size()
     || vmx.size() != vsigobs.size()
     || vmx.size() != vsigexp.size()) ERROR("Error parsing text file. Model point not fully specified");
  
  // Sorting vectors
  vector<size_t> perm = SortPermutation(vmx);
  vmx      = ApplyPermutation(vmx      , perm);
  vmy	   = ApplyPermutation(vmy      , perm);
  vxsec	   = ApplyPermutation(vxsec    , perm);	
  vexsec   = ApplyPermutation(vexsec   , perm);	
  vobs	   = ApplyPermutation(vobs     , perm);	
  vobsup   = ApplyPermutation(vobsup   , perm);	 	
  vobsdown = ApplyPermutation(vobsdown , perm);	
  vexp	   = ApplyPermutation(vexp     , perm);	
  vup	   = ApplyPermutation(vup      , perm);	
  vdown	   = ApplyPermutation(vdown    , perm);	
  v2up	   = ApplyPermutation(v2up     , perm);	
  v2down   = ApplyPermutation(v2down   , perm);	
  vsigobs  = ApplyPermutation(vsigobs  , perm);	 	
  vsigexp  = ApplyPermutation(vsigexp  , perm);	

  TCanvas can;
  //can.SetGrid(); 
  can.SetFillStyle(4000);
  TString chi1n = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  TString chi2n = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.]{#scale[0.85]{_{2}}}";
  TString chi1pm= "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  TString chii= "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0,#pm}}}#kern[-3.]{#scale[0.85]{_{i}}}";
  TString chij= "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0,#mp}}}#kern[-3.]{#scale[0.85]{_{j}}}";
  TString mass_ = "m#kern[0.1]{#lower[-0.12]{_{";
  float minh=200, maxh=1000;
  TH1D histo("histo", "", 18, minh, maxh);
  histo.SetMinimum(0);
  histo.SetMaximum(4.5);
  histo.GetYaxis()->CenterTitle(true);
  histo.GetXaxis()->SetLabelOffset(0.02);
  histo.SetXTitle("Higgsino mass "+mass_+chi1n+"}}} [GeV]");
  histo.SetYTitle("#sigma_{excl}^{95% CL}/#sigma_{theor}");
  histo.Draw();

  TLine line;
  line.SetLineColor(4); line.SetLineStyle(2); line.SetLineWidth(4);
  TLatex cmslabel;
  
  cmslabel.SetNDC(kTRUE);

  TGraphAsymmErrors grexp2(vmx.size(), &(vmx[0]), &(vexp[0]), &(zeroes[0]), &(zeroes[0]), &(v2down[0]), &(v2up[0]));
  grexp2.SetLineColor(1); grexp2.SetFillColor(5); grexp2.SetLineWidth(3); grexp2.SetLineStyle(2);
  grexp2.Draw("e3 same");
  TGraphAsymmErrors grexp1(vmx.size(), &(vmx[0]), &(vexp[0]), &(zeroes[0]), &(zeroes[0]), &(vdown[0]), &(vup[0]));
  grexp1.SetLineColor(1); grexp1.SetFillColor(3); grexp1.SetLineWidth(3); grexp1.SetLineStyle(2);
  grexp1.Draw("e3 same");
  TGraph grexp(vmx.size(), &(vmx[0]), &(vexp[0]));
  grexp.SetLineWidth(3); grexp.SetLineStyle(2);
  grexp.Draw("same"); 
  TGraph grobs(vmx.size(), &(vmx[0]), &(vobs[0]));
  grobs.SetLineWidth(1); 
  grobs.Draw("same"); 

  //// Drawing CMS labels and line at 1
  TString cmsPrel = "#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}";
  TString lumiEner = "#font[42]{"+lumi+" fb^{-1} (13 TeV)}"; lumiEner.ReplaceAll("p",".");
  TString ppChiChi = "pp #rightarrow "+chii+"#kern[0.6]{"+chij+"}  #rightarrow hh#tilde{G}#tilde{G}";
  TString mChis = mass_+chi2n+"}}} #approx "+mass_+chi1pm+"}}} #approx "+mass_+chi1n+"}}}, "+mass_+"#tilde{G}}}} = 1 GeV";
  cmslabel.SetTextAlign(11); cmslabel.SetTextSize(0.06);
  cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015, cmsPrel);
  cmslabel.SetTextAlign(31); cmslabel.SetTextSize(0.056);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.015, lumiEner);
  line.DrawLine(minh, 1, maxh, 1);
  //// Drawing process and masses
  cmslabel.SetTextAlign(31); cmslabel.SetTextSize(0.045);
  cmslabel.SetTextFont(132);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.023, opts.BottomMargin()+0.09, ppChiChi);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.023, opts.BottomMargin()+0.04, mChis);


  double legX(0.54), legY(1-opts.TopMargin()-0.04), legSingle = 0.05;
  double legW = 0.26, legH = legSingle*4;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(0.04); leg.SetFillColor(0); 
  leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.AddEntry(&line, "NLO+NLL", "l");
  leg.AddEntry(&grobs, "\"Observed\"", "l");
  leg.AddEntry(&grexp1, "Expected #pm #sigma");
  leg.AddEntry(&grexp2, "Expected #pm 2#sigma");
  leg.Draw();

  histo.Draw("axis same");
  can.SaveAs("plots/higgsino_limits_lumi"+lumi+".pdf");

  // for(size_t i = 0; i < vxsec.size(); ++i) 
  //   cout<<vmx[i]<<" -> "<<vexp[i]<<"+"<<vup[i]<<"++"<<v2up[i]<<" -"<<vdown[i]<<"--"<<v2down[i]<<endl;

  //////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //// Plotting limits on absolute xsec
  maxy=-99.; miny=1e99;
  for(size_t i = 0; i < vxsec.size(); ++i){
    vxsec[i]   *= 1000; // Converting it to fb
    vexsec[i]  *= vxsec[i];
    vobs[i]    *= vxsec[i];
    vobsup[i]  *= vxsec[i]; 	
    vobsdown[i]*= vxsec[i];
    vexp[i]    *= vxsec[i];
    vup[i]     *= vxsec[i];
    vdown[i]   *= vxsec[i];
    v2up[i]    *= vxsec[i];
    v2down[i]  *= vxsec[i];
    if(miny > min(vexp[i]-v2down[i], vxsec[i])) miny = min(vexp[i]-v2down[i], vxsec[i]);
    if(maxy < max(vexp[i]+v2up[i], vxsec[i])) maxy = max(vexp[i]+v2up[i], vxsec[i]);
    //cout<<vmx[i]<<" -> "<<vobs[i]<<endl;
  }

  histo.GetXaxis()->SetLabelOffset(0.01);
  histo.SetMinimum(miny/2.);
  histo.SetMaximum(maxy*1.2);
  histo.SetYTitle("#sigma_{excl}^{95% CL} #times BF(hh #rightarrow bbbb) [fb]");
  histo.Draw();
  TGraphAsymmErrors gexp2(vmx.size(), &(vmx[0]), &(vexp[0]), &(zeroes[0]), &(zeroes[0]), &(v2down[0]), &(v2up[0]));
  gexp2.SetLineColor(1); gexp2.SetFillColor(5); gexp2.SetLineWidth(3); gexp2.SetLineStyle(2);
  gexp2.Draw("e3 same");
  TGraphAsymmErrors gexp1(vmx.size(), &(vmx[0]), &(vexp[0]), &(zeroes[0]), &(zeroes[0]), &(vdown[0]), &(vup[0]));
  gexp1.SetLineColor(1); gexp1.SetFillColor(3); gexp1.SetLineWidth(3); gexp1.SetLineStyle(2);
  gexp1.Draw("e3 same");
  TGraph gexp(vmx.size(), &(vmx[0]), &(vexp[0]));
  gexp.SetLineWidth(3); gexp.SetLineStyle(2);
  gexp.Draw("same"); 
  TGraph gobs(vmx.size(), &(vmx[0]), &(vobs[0]));
  gobs.SetLineWidth(1); 
  gobs.Draw("same"); 
  TGraph gxsec(vmx.size(), &(vmx[0]), &(vxsec[0]));
  gxsec.SetLineWidth(4); gxsec.SetLineColor(4); gxsec.SetLineStyle(2);
  gxsec.Draw("same");

  can.SetLogy(true);

  legX = 1-opts.RightMargin()-0.02;
  leg.SetX1NDC(legX-legW); leg.SetX2NDC(legX);
  leg.Draw();

  //// Drawing CMS labels
  cmslabel.SetTextAlign(11); cmslabel.SetTextSize(0.06);
  cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015, cmsPrel);
  cmslabel.SetTextAlign(31); cmslabel.SetTextSize(0.056);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.015, lumiEner);

  //// Drawing process and masses
  cmslabel.SetTextAlign(11); cmslabel.SetTextSize(0.045);
  cmslabel.SetTextFont(132);
  cmslabel.DrawLatex(opts.LeftMargin()+0.03, opts.BottomMargin()+0.09, ppChiChi);
  cmslabel.DrawLatex(opts.LeftMargin()+0.03, opts.BottomMargin()+0.04, mChis);


  histo.Draw("axis same");
  can.SaveAs("plots/higgsino_limits_fb_lumi"+lumi+".pdf");

  //////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //// Plotting discovery significance
  can.SetLogy(false);
  histo.GetXaxis()->SetLabelOffset(0.02);
  histo.SetMinimum(0);
  histo.SetMaximum(3.5);
  if(lumi=="40") histo.SetMaximum(4.5);
  histo.SetYTitle(" Expected discovery significance [#sigma]");
  histo.Draw();

  TGraph gsig(vmx.size(), &(vmx[0]), &(vsigexp[0]));
  gsig.SetLineWidth(3);  gsig.SetLineColor(4);
  gsig.Draw("same"); 
  //// Drawing CMS labels
  cmslabel.SetTextAlign(11); cmslabel.SetTextSize(0.06);
  cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015, cmsPrel);
  cmslabel.SetTextAlign(31); cmslabel.SetTextSize(0.056);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.015, lumiEner);
  //// Drawing process and masses
  cmslabel.SetTextAlign(33); cmslabel.SetTextSize(0.045);
  cmslabel.SetTextFont(132);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.023, 1-opts.TopMargin()-0.025, ppChiChi);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.023, 1-opts.TopMargin()-0.09, mChis);

  can.SaveAs("plots/higgsino_significance_lumi"+lumi+".pdf");


}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"model", required_argument, 0, 'm'},
      {"file", required_argument, 0, 'f'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "f:m:", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'm':
      model = optarg;
      break;
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
