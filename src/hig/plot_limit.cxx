#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include <unistd.h>
#include <getopt.h>

#include "TCanvas.h"
#include "TFile.h"
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
  TString lumi = "35p9";
  TString filename = "txt/limits/limits_TChiHH_lumi"+lumi+"_paper.txt";
  TString model = "TChiHH";
  TString datestamp = "";
  bool do_paper = true;

  PlotOpt opts("txt/plot_styles.txt", "Std1D");
  const double hh4b_bf = 0.5824*0.5824;
}

void GetOptions(int argc, char *argv[]);
void DrawCMSLabels(TString type);

int main(int argc, char *argv[]){
  GetOptions(argc, argv);

  //// Setting plot style
  setPlotStyle(opts);
  gStyle->SetGridStyle(3);

  if(filename == "") ERROR("No input file provided");
  ifstream infile(filename);

  vector<double> vmx, vmy, vxsec, vexsec, vobs, vobsup, vobsdown;
  vector<double> vexp, vup, vdown, v2up, v2down, vsigobs, vsigexp, zeroes, ones;
  double maxy=-99., miny=1e99, maxsig = -99., minsig = 1e99;
  vector<double> vxsecup, vxsecdown;
  
  string line_s;
  while(getline(infile, line_s)){
    istringstream iss(line_s);
    double pmx, pmy, pxsec, pexsec, pobs, pobsup, pobsdown, pexp, pup, pdown, p2up, p2down, sigobs, sigexp;
    iss >> pmx >> pmy >> pxsec >> pexsec >> pobs >> pobsup >> pobsdown 
	>> pexp >> pup >> pdown >> p2up >> p2down >> sigobs >> sigexp;
    vmx.push_back(pmx);
    vmy.push_back(pmy);
    vxsec.push_back(pxsec/hh4b_bf);
    vexsec.push_back(pexsec/hh4b_bf);
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
    if(minsig > vsigobs.back()) minsig = vsigobs.back();
    if(maxsig < vsigobs.back()) maxsig = vsigobs.back();

    vxsecup.push_back(1+pexsec);
    vxsecdown.push_back(1-pexsec);
    
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
  vxsecup  = ApplyPermutation(vxsecup  , perm);	
  vxsecdown= ApplyPermutation(vxsecdown, perm);	

  TCanvas can;
  //can.SetGrid(); 
  can.SetFillStyle(4000);
  TString chi1n = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.99]{^{0}}}#kern[-1.3]{#scale[0.99]{_{1}}}";
  TString chi2n = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.99]{^{0}}}#kern[-1.]{#scale[0.99]{_{2}}}";
  TString chi1pm= "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.99]{^{#pm}}}#kern[-1.3]{#scale[0.99]{_{1}}}";
  TString chii= "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.99]{^{0,#pm}}}#kern[-3.]{#scale[0.99]{_{i}}}";
  TString chij= "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.99]{^{0,#mp}}}#kern[-3.]{#scale[0.99]{_{j}}}";
  TString chi10= "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.99]{^{0}}}#kern[-1.3]{#scale[0.99]{_{1}}}";
  TString xsoft= "X#scale[0.85]{_{soft}}";
  TString mass_ = "m#kern[0.1]{#lower[-0.12]{_{";
  float minh=200, maxh=1000, maxXsec = 5e3;
  if(do_paper) {
    minh = 127;
    maxXsec = 1e5;
  }
  TH1D histo("histo", "", 18, minh, maxh);
  histo.SetMinimum(0);
  histo.SetMaximum(4.5);
  histo.GetYaxis()->CenterTitle(true);
  histo.GetXaxis()->SetLabelOffset(0.01);
  histo.SetXTitle("Higgsino mass "+mass_+chi1n+"}}} [GeV]");
  histo.SetYTitle("#sigma_{excl}^{95% CL}/#sigma_{theory}");
  histo.Draw();

  int thcolor = kRed+1, thwidth = 3;
  TLine linXsec;
  linXsec.SetLineColor(thcolor); linXsec.SetLineStyle(1); linXsec.SetLineWidth(thwidth);
  TLatex label;  label.SetNDC(kTRUE);

  //////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //// Plotting limits on relative xsec

  int cyellow = kOrange, cgreen = kGreen+1;
  TGraphAsymmErrors grexp2(vmx.size(), &(vmx[0]), &(vexp[0]), &(zeroes[0]), &(zeroes[0]), &(v2down[0]), &(v2up[0]));
  grexp2.SetLineColor(1); grexp2.SetFillColor(cyellow); grexp2.SetLineWidth(3); grexp2.SetLineStyle(2);
  grexp2.Draw("e3 same");
  TGraphAsymmErrors grexp1(vmx.size(), &(vmx[0]), &(vexp[0]), &(zeroes[0]), &(zeroes[0]), &(vdown[0]), &(vup[0]));
  grexp1.SetLineColor(1); grexp1.SetFillColor(cgreen); grexp1.SetLineWidth(3); grexp1.SetLineStyle(2);
  grexp1.Draw("e3 same");
  TGraph grexp(vmx.size(), &(vmx[0]), &(vexp[0]));
  grexp.SetLineWidth(3); grexp.SetLineStyle(2);
  grexp.Draw("same"); 
  TGraph grobs(vmx.size(), &(vmx[0]), &(vobs[0]));
  grobs.SetLineWidth(3); 
  grobs.Draw("same"); 
  TGraph grxsecup(vmx.size(), &(vmx[0]), &(vxsecup[0]));
  grxsecup.SetLineWidth(1); grxsecup.SetLineStyle(2); grxsecup.SetLineColor(thcolor); 
  grxsecup.Draw("same"); 
  TGraph grxsecdown(vmx.size(), &(vmx[0]), &(vxsecdown[0]));
  grxsecdown.SetLineWidth(1); grxsecdown.SetLineStyle(2); grxsecdown.SetLineColor(thcolor); 
  grxsecdown.Draw("same"); 

  //// Drawing CMS labels and line at 1
  TString ppChiChi = "pp #rightarrow "+chii+"#kern[0.6]{"+chij+"}  #rightarrow "+chi10+"#kern[0.3]{"+chi10+"} + "
    +xsoft+"#rightarrow hh#tilde{G}#tilde{G} + "+xsoft;
  double ppSize = 0.055, ppY = 1-opts.TopMargin()-0.03, ppY2 = 1-opts.TopMargin()-0.11;
  double legSize = 0.044;

  TString mChis = mass_+chi2n+"}}} #approx "+mass_+chi1pm+"}}} #approx "+mass_+chi1n+"}}}, "
    +mass_+"#tilde{G}}}} = 1 GeV";
  DrawCMSLabels("Supplementary");
  linXsec.DrawLine(minh, 1, maxh, 1);

  TLine line;
  int ibox = 0;
  vector<vector<float> > boxes;
  double legX(0.45), legY(1-opts.TopMargin()-0.24), legSingle = 0.053;
  double legW = 0.26, legH = legSingle*5;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetX1NDC(legX-legW); leg.SetX2NDC(legX); // So that GetX1NDC works in getLegendBoxes
  leg.SetY1NDC(legY-legH); leg.SetY2NDC(legY); // So that GetX1NDC works in getLegendBoxes
  leg.SetTextSize(legSize); leg.SetFillColor(0); 
  leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.AddEntry(&linXsec, "NLO+NLL theory #pm s.d.", "l");
  leg.AddEntry(&grobs, " ", "n");
  leg.AddEntry(&grobs, " ", "n");
  leg.AddEntry(&grobs, "Observed", "l");
  leg.AddEntry(&grexp1, "68% expected");
  leg.AddEntry(&grexp2, "95% expected");
  leg.Draw();

  // Drawing theory error lines on legend
  getLegendBoxes(leg, boxes);
  ibox = 0;
  line.SetLineColor(thcolor);line.SetLineWidth(1);line.SetLineStyle(2);
  line.DrawLineNDC(boxes[ibox][0], boxes[ibox][1], boxes[ibox][2], boxes[ibox][1]);
  line.DrawLineNDC(boxes[ibox][0], boxes[ibox][3], boxes[ibox][2], boxes[ibox][3]);

  label.SetTextAlign(12); label.SetTextSize(legSize); label.SetTextFont(42); 
  label.DrawLatex(legX-legW+0.01, legY-legSingle*2, "95% CL upper limits");
  //// Drawing process and masses
  label.SetTextAlign(11); label.SetTextSize(ppSize/1.07);
  label.SetTextFont(132);
  label.DrawLatex(legX-legW+0.01, opts.BottomMargin()+0.70, ppChiChi);
  label.DrawLatex(legX-legW+0.01, opts.BottomMargin()+0.643, mChis);



  histo.Draw("axis same");
  TString basename = "plots/higgsino_limits_lumi"+lumi;
  if(datestamp != "") basename += "_"+datestamp;
  basename += ".pdf";
  TString pname = basename;
  can.SaveAs(pname);

  // // Saving root file
  // pname = "CMS"; if(!do_paper) pname += "-PAS";
  // pname += "-SUS-16-044_AuxFigure_9-b.root";
  // TFile file(pname, "recreate");
  // file.cd();
  // grexp2.Write("ExpLimit_2Sigma");
  // grexp1.Write("ExpLimit_1Sigma");
  // grexp.Write("ExpLimit");
  // grobs.Write("ObsLimit");
  // file.Close();
  // cout<<"Saved graphs in "<<pname<<endl<<endl;

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
    vxsecup[i]  *= vxsec[i];
    vxsecdown[i]  *= vxsec[i];
  }

  histo.GetXaxis()->SetLabelOffset(0.01);
  histo.SetMinimum(miny/2.);
  histo.SetMaximum(maxXsec);
  histo.SetYTitle("#sigma [fb]");
  histo.Draw();
  TGraphAsymmErrors gexp2(vmx.size(), &(vmx[0]), &(vexp[0]), &(zeroes[0]), &(zeroes[0]), &(v2down[0]), &(v2up[0]));
  gexp2.SetLineColor(1); gexp2.SetFillColor(cyellow); gexp2.SetLineWidth(3); gexp2.SetLineStyle(2);
  gexp2.Draw("e3 same");
  TGraphAsymmErrors gexp1(vmx.size(), &(vmx[0]), &(vexp[0]), &(zeroes[0]), &(zeroes[0]), &(vdown[0]), &(vup[0]));
  gexp1.SetLineColor(1); gexp1.SetFillColor(cgreen); gexp1.SetLineWidth(3); gexp1.SetLineStyle(2);
  gexp1.Draw("e3 same");
  TGraph gexp(vmx.size(), &(vmx[0]), &(vexp[0]));
  gexp.SetLineWidth(3); gexp.SetLineStyle(2);
  gexp.Draw("same"); 
  TGraph gobs(vmx.size(), &(vmx[0]), &(vobs[0]));
  gobs.SetLineWidth(3); 
  gobs.Draw("same"); 
  TGraph gxsec(vmx.size(), &(vmx[0]), &(vxsec[0]));
  gxsec.SetLineWidth(thwidth); gxsec.SetLineColor(thcolor); gxsec.SetLineStyle(1);
  gxsec.Draw("same");
  TGraph gxsecup(vmx.size(), &(vmx[0]), &(vxsecup[0]));
  gxsecup.SetLineWidth(1); gxsecup.SetLineStyle(2); gxsecup.SetLineColor(thcolor); 
  gxsecup.Draw("same"); 
  TGraph gxsecdown(vmx.size(), &(vmx[0]), &(vxsecdown[0]));
  gxsecdown.SetLineWidth(1); gxsecdown.SetLineStyle(2); gxsecdown.SetLineColor(thcolor); 
  gxsecdown.Draw("same"); 

  can.SetLogy(true);

  legX = 1-opts.RightMargin()-0.1;
  legY += 0.02;
  leg.SetX1NDC(legX-legW); leg.SetX2NDC(legX);
  leg.SetY1NDC(legY-legH); leg.SetY2NDC(legY);
  leg.Draw();
  label.SetTextAlign(12); label.SetTextSize(legSize); label.SetTextFont(42); 
  label.DrawLatex(legX-legW+0.01, legY-legSingle*2, "95% CL upper limits");

  // Drawing theory error lines on legend
  getLegendBoxes(leg, boxes);
  ibox = 0;
  line.SetLineColor(thcolor);line.SetLineWidth(1);line.SetLineStyle(2);
  line.DrawLineNDC(boxes[ibox][0], boxes[ibox][1], boxes[ibox][2], boxes[ibox][1]);
  line.DrawLineNDC(boxes[ibox][0], boxes[ibox][3], boxes[ibox][2], boxes[ibox][3]);


  //// Drawing CMS labels
  if(do_paper) DrawCMSLabels("");
  else DrawCMSLabels("Preliminary");

  //// Drawing process and masses
  label.SetTextAlign(33); label.SetTextSize(ppSize);
  label.SetTextFont(132);
  label.DrawLatex(1-opts.RightMargin()-0.03, ppY, ppChiChi);
  label.DrawLatex(1-opts.RightMargin()-0.03, ppY2, mChis);

  histo.Draw("axis same");
  pname = basename;
  pname.ReplaceAll("lumi", "fb_lumi");
  can.SaveAs(pname);

  // Saving root file
  pname = "CMS"; if(!do_paper) pname += "-PAS";
  pname += "-SUS-16-044_Figure_9.root";
  TFile file2(pname, "recreate");
  file2.cd();
  gexp2.Write("ExpLimit_2Sigma");
  gexp1.Write("ExpLimit_1Sigma");
  gexp.Write("ExpLimit");
  gobs.Write("ObsLimit");
  gxsec.Write("Xsec");
  gxsecup.Write("XsecUp");
  gxsecdown.Write("XsecDown");
  file2.Close();
  cout<<"Saved graphs in "<<pname<<endl<<endl;

  //////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //// Plotting expected discovery significance

  can.SetLogy(false);
  histo.GetXaxis()->SetLabelOffset(0.01);
  histo.SetMinimum(0);
  histo.SetMaximum(4.5);
  if(lumi=="40") histo.SetMaximum(4.5);
  histo.SetYTitle(" Expected discovery significance [#sigma]");
  histo.Draw();

  TGraph gsig(vmx.size(), &(vmx[0]), &(vsigexp[0]));
  gsig.SetLineWidth(3);  gsig.SetLineColor(4);
  gsig.Draw("same"); 

  //// Drawing CMS labels
  if(do_paper) DrawCMSLabels("Supplementary");
  else DrawCMSLabels("Preliminary");
  //// Drawing process and masses
  label.SetTextAlign(33); label.SetTextSize(ppSize);
  label.SetTextFont(132);
  label.DrawLatex(1-opts.RightMargin()-0.03, ppY, ppChiChi);
  label.DrawLatex(1-opts.RightMargin()-0.03, ppY2, mChis);

  // pname = basename;
  // pname.ReplaceAll("lumi", "exp_significance");
  // can.SaveAs(pname);

  // // Saving root file
  // pname = "CMS"; if(!do_paper) pname += "-PAS";
  // pname += "-SUS-16-044_AuxFigure_3_ExpSignificance.root";
  // TFile file3(pname, "recreate");
  // file3.cd();
  // gsig.Write("ExpSignificance");
  // file3.Close();
  // cout<<"Saved graphs in "<<pname<<endl<<endl;

  //////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //// Plotting observed significance
  can.SetLogy(false);
  histo.GetXaxis()->SetLabelOffset(0.01);
  histo.SetMinimum(-2.5);
  histo.SetMaximum(2.5);
  histo.SetYTitle("Observed signal significance [#sigma]");
  histo.Draw();

  TGraph gobssig(vmx.size(), &(vmx[0]), &(vsigobs[0]));
  gobssig.SetLineWidth(3);  gobssig.SetLineColor(4);
  gobssig.Draw("same"); 

  //// Drawing CMS labels
  if(do_paper) DrawCMSLabels("Supplementary");
  else DrawCMSLabels("Preliminary");

  //// Drawing process and masses
  label.SetTextAlign(33); label.SetTextSize(ppSize);
  label.SetTextFont(132);
  label.DrawLatex(1-opts.RightMargin()-0.03, ppY, ppChiChi);
  label.DrawLatex(1-opts.RightMargin()-0.03, ppY2, mChis);

  // pname = basename;
  // pname.ReplaceAll("lumi", "obs_significance");
  // can.SaveAs(pname);

  // // Saving root file
  // pname = "CMS"; if(!do_paper) pname += "-PAS";
  // pname += "-SUS-16-044_AuxFigure_4_ObsSignificance.root";
  // TFile file4(pname, "recreate");
  // file4.cd();
  // gobssig.Write("ObsSignificance");
  // file4.Close();
  // cout<<"Saved graphs in "<<pname<<endl<<endl;

  //////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //// Plotting observed and expected significance
  can.SetLogy(false);
  histo.GetXaxis()->SetLabelOffset(0.01);
  histo.SetMinimum(-1);
  histo.SetMaximum(5);
  histo.SetYTitle("Signal significance [#sigma]");
  histo.Draw();

  gsig.SetLineStyle(2);  gsig.SetLineColor(kRed+1);  
  gsig.Draw("same"); 
  gobssig.Draw("same"); 

  //// Drawing CMS labels
  if(do_paper) DrawCMSLabels("Supplementary");
  else DrawCMSLabels("Preliminary");

  //// Drawing process and masses
  label.SetTextAlign(33); label.SetTextSize(ppSize);
  label.SetTextFont(132);
  label.DrawLatex(1-opts.RightMargin()-0.03, ppY, ppChiChi);
  label.DrawLatex(1-opts.RightMargin()-0.03, ppY2, mChis);

  //// Drawing legend
  legX = 1- opts.RightMargin()-0.03, legY = ppY2 - 0.2; legSingle = 0.07;
  legW = 0.16; legH = legSingle*2;
  TLegend leg2(legX-legW, legY-legH, legX, legY);
  leg2.SetTextSize(legSize); leg2.SetFillColor(0); 
  leg2.SetFillStyle(0); leg2.SetBorderSize(0);
  leg2.AddEntry(&gsig, "Expected", "l");
  leg2.AddEntry(&gobssig, "Observed", "l");
  leg2.Draw();

  pname = basename;
  pname.ReplaceAll("lumi", "significance");
  can.SaveAs(pname);

  // // Saving root file
  // pname = "CMS"; if(!do_paper) pname += "-PAS";
  // pname += "-SUS-16-044_AuxFigure_3_Significance.root";
  // TFile file5(pname, "recreate");
  // file5.cd();
  // gobssig.Write("ObsSignificance");
  // gsig.Write("ExpSignificance");
  // file5.Close();
  // cout<<"Saved graphs in "<<pname<<endl<<endl;


}

void DrawCMSLabels(TString type){
  TString cmsLogo = "#font[62]{CMS}#scale[0.8]{#font[52]{ "+type+"}}";
  TString lumiEner = "#font[42]{"+lumi+" fb^{-1} (13 TeV)}"; lumiEner.ReplaceAll("p",".");

  TLatex cmslabel;
  cmslabel.SetNDC(kTRUE); cmslabel.SetTextAlign(11); cmslabel.SetTextSize(0.06);
  cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015, cmsLogo);
  cmslabel.SetTextAlign(31); cmslabel.SetTextSize(0.054);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.015, lumiEner);
}




void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"model", required_argument, 0, 'm'},
      {"file", required_argument, 0, 'f'},
      {"datestamp", required_argument, 0, 'd'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "f:m:d:", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 'm':
      model = optarg;
      break;
    case 'f':
      filename = optarg;
      break;
    case 'd':
      datestamp = optarg;
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
