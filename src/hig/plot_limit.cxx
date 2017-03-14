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
  TString lumi = "35p9";
  TString filename = "txt/limits/limits_TChiHH_lumi"+lumi+"_wilk.txt";
  TString model = "TChiHH";
  TString datestamp = "";
}

void higgsinoCrossSection(int hig_mass, float &xsec, float &xsec_unc);
void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  GetOptions(argc, argv);

  //// Setting plot style
  PlotOpt opts("txt/plot_styles.txt", "Std1D");
  setPlotStyle(opts);
  gStyle->SetGridStyle(3);

  if(filename == "") ERROR("No input file provided");
  ifstream infile(filename);

  vector<double> vmx, vmy, vxsec, vexsec, vobs, vobsup, vobsdown;
  vector<double> vexp, vup, vdown, v2up, v2down, vsigobs, vsigexp, zeroes, ones;
  double maxy=-99., miny=1e99;
  vector<double> vxsecup, vxsecdown;
  
  string line_s;
  while(getline(infile, line_s)){
    istringstream iss(line_s);
    double pmx, pmy, pxsec, pexsec, pobs, pobsup, pobsdown, pexp, pup, pdown, p2up, p2down, sigobs, sigexp;
    iss >> pmx >> pmy >> pxsec >> pexsec >> pobs >> pobsup >> pobsdown 
	>> pexp >> pup >> pdown >> p2up >> p2down >> sigobs >> sigexp;
    if(pmx==175) continue;
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
  TString chi1n = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  TString chi2n = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.]{#scale[0.85]{_{2}}}";
  TString chi1pm= "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  TString chii= "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0,#pm}}}#kern[-3.]{#scale[0.85]{_{i}}}";
  TString chij= "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0,#mp}}}#kern[-3.]{#scale[0.85]{_{j}}}";
  TString chi10= "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}";
  TString xsoft= "X#scale[0.85]{_{soft}}";
  TString mass_ = "m#kern[0.1]{#lower[-0.12]{_{";
  float minh=200, maxh=1000;
  TH1D histo("histo", "", 18, minh, maxh);
  histo.SetMinimum(0);
  histo.SetMaximum(4.5);
  histo.GetYaxis()->CenterTitle(true);
  histo.GetXaxis()->SetLabelOffset(0.02);
  histo.SetXTitle("Higgsino mass "+mass_+chi1n+"}}} [GeV]");
  histo.SetYTitle("#sigma_{excl}^{95% CL}/#sigma_{theory}");
  histo.Draw();

  int thcolor = kRed+1, thwidth = 3;
  TLine line;
  line.SetLineColor(thcolor); line.SetLineStyle(1); line.SetLineWidth(thwidth);
  TLatex cmslabel;
  
  cmslabel.SetNDC(kTRUE);

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
  TString cmsPrel = "#font[62]{CMS} #scale[0.8]{#font[52]{Preliminary}}";
  TString cmsSim = "#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}";
  TString lumiEner = "#font[42]{"+lumi+" fb^{-1} (13 TeV)}"; lumiEner.ReplaceAll("p",".");
  TString ppChiChi = "pp #rightarrow "+chii+"#kern[0.6]{"+chij+"}  #rightarrow "+chi10+"#kern[0.3]{"+chi10+"} + "+xsoft+"#rightarrow hh#tilde{G}#tilde{G} + "+xsoft;

  TString mChis = mass_+chi2n+"}}} #approx "+mass_+chi1pm+"}}} #approx "+mass_+chi1n+"}}}, "+mass_+"#tilde{G}}}} = 1 GeV";
  cmslabel.SetTextAlign(11); cmslabel.SetTextSize(0.06);
  cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015, cmsPrel);
  cmslabel.SetTextAlign(31); cmslabel.SetTextSize(0.056);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.015, lumiEner);
  line.DrawLine(minh, 1, maxh, 1);

  double legX(0.5), legY(1-opts.TopMargin()-0.24), legSingle = 0.05;
  double legW = 0.26, legH = legSingle*5;
  TLegend leg(legX-legW, legY-legH, legX, legY);
  leg.SetTextSize(0.04); leg.SetFillColor(0); 
  leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.AddEntry(&line, "NLO+NLL theory #pm s.d.", "l");
  leg.AddEntry(&grobs, " ", "n");
  leg.AddEntry(&grobs, " ", "n");
  leg.AddEntry(&grobs, "Observed", "l");
  leg.AddEntry(&grexp1, "68% expected");
  leg.AddEntry(&grexp2, "95% expected");
  leg.Draw();

  cmslabel.SetTextAlign(12); cmslabel.SetTextSize(0.04); cmslabel.SetTextFont(42); 
  cmslabel.DrawLatex(legX-legW+0.01, legY-legSingle*2, "95% CL upper limits");
  //// Drawing process and masses
  cmslabel.SetTextAlign(11); cmslabel.SetTextSize(0.045);
  cmslabel.SetTextFont(132);
  cmslabel.DrawLatex(legX-legW+0.01, opts.BottomMargin()+0.65, ppChiChi);
  cmslabel.DrawLatex(legX-legW+0.01, opts.BottomMargin()+0.6, mChis);



  histo.Draw("axis same");
  TString pname = "plots/higgsino_limits_lumi"+lumi;
  if(datestamp != "") pname += "_"+datestamp;
  pname += ".pdf";
  can.SaveAs(pname);

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
  histo.SetMaximum(5*1e3);
  histo.SetYTitle("#sigma #times BF(hh #rightarrow bbbb) [fb]");
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

  legX = 1-opts.RightMargin()-0.09;
  legY += 0.05;
  leg.SetX1NDC(legX-legW); leg.SetX2NDC(legX);
  leg.SetY1NDC(legY-legH); leg.SetY2NDC(legY);
  leg.Draw();
  cmslabel.SetTextAlign(12); cmslabel.SetTextSize(0.04); cmslabel.SetTextFont(42); 
  cmslabel.DrawLatex(legX-legW+0.01, legY-legSingle*2, "95% CL upper limits");



  //// Drawing CMS labels
  cmslabel.SetTextAlign(11); cmslabel.SetTextSize(0.06);
  cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015, cmsPrel);
  cmslabel.SetTextAlign(31); cmslabel.SetTextSize(0.056);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.015, lumiEner);

  //// Drawing process and masses
  // cmslabel.SetTextAlign(11); cmslabel.SetTextSize(0.045);
  // cmslabel.SetTextFont(132);
  // cmslabel.DrawLatex(opts.LeftMargin()+0.03, opts.BottomMargin()+0.09, ppChiChi);
  // cmslabel.DrawLatex(opts.LeftMargin()+0.03, opts.BottomMargin()+0.04, mChis);
  cmslabel.SetTextAlign(33); cmslabel.SetTextSize(0.045);
  cmslabel.SetTextFont(132);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.028, 1-opts.TopMargin()-0.03, ppChiChi);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.028, 1-opts.TopMargin()-0.095, mChis);

  histo.Draw("axis same");
  pname.ReplaceAll("lumi", "fb_lumi");
  can.SaveAs(pname);

  //////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //// Plotting discovery significance
  can.SetLogy(false);
  histo.GetXaxis()->SetLabelOffset(0.02);
  histo.SetMinimum(0);
  histo.SetMaximum(4.5);
  if(lumi=="40") histo.SetMaximum(4.5);
  histo.SetYTitle(" Expected discovery significance [#sigma]");
  histo.Draw();

  TGraph gsig(vmx.size(), &(vmx[0]), &(vsigexp[0]));
  gsig.SetLineWidth(3);  gsig.SetLineColor(4);
  gsig.Draw("same"); 
  //// Drawing CMS labels
  cmslabel.SetTextAlign(11); cmslabel.SetTextSize(0.06);
  cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015, cmsSim);
  cmslabel.SetTextAlign(31); cmslabel.SetTextSize(0.056);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.015, lumiEner);
  //// Drawing process and masses
  cmslabel.SetTextAlign(33); cmslabel.SetTextSize(0.045);
  cmslabel.SetTextFont(132);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.023, 1-opts.TopMargin()-0.025, ppChiChi);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.023, 1-opts.TopMargin()-0.09, mChis);

  pname.ReplaceAll("fb_lumi", "significance");
  can.SaveAs(pname);


}

void higgsinoCrossSection(int hig_mass, float &xsec, float &xsec_unc) {
  if(hig_mass == 127) { xsec = .5824*.5824*7602.2/1000; xsec_unc = 0.0393921; return;}
  else if(hig_mass == 150) { xsec = .5824*.5824*3832.31/1000; xsec_unc = 0.0413612; return;}
  else if(hig_mass == 175) { xsec = .5824*.5824*2267.94/1000; xsec_unc = 0.044299;  return;}
  else if(hig_mass == 200) { xsec = .5824*.5824*1335.62/1000; xsec_unc = 0.0474362; return;}
  else if(hig_mass == 225) { xsec = .5824*.5824*860.597/1000; xsec_unc = 0.0504217; return;}
  else if(hig_mass == 250) { xsec = .5824*.5824*577.314/1000; xsec_unc = 0.0532731; return;}
  else if(hig_mass == 275) { xsec = .5824*.5824*400.107/1000; xsec_unc = 0.0560232; return;}
  else if(hig_mass == 300) { xsec = .5824*.5824*284.855/1000; xsec_unc = 0.0586867; return;}
  else if(hig_mass == 325) { xsec = .5824*.5824* 207.36/1000; xsec_unc = 0.0613554; return;}
  else if(hig_mass == 350) { xsec = .5824*.5824*153.841/1000; xsec_unc = 0.0640598; return;}
  else if(hig_mass == 375) { xsec = .5824*.5824*116.006/1000; xsec_unc = 0.066892;  return;}
  else if(hig_mass == 400) { xsec = .5824*.5824*88.7325/1000; xsec_unc = 0.0697517; return;}
  else if(hig_mass == 425) { xsec = .5824*.5824*68.6963/1000; xsec_unc = 0.0723531; return;}
  else if(hig_mass == 450) { xsec = .5824*.5824*53.7702/1000; xsec_unc = 0.0748325; return;}
  else if(hig_mass == 475) { xsec = .5824*.5824*42.4699/1000; xsec_unc = 0.0775146; return;}
  else if(hig_mass == 500) { xsec = .5824*.5824*33.8387/1000; xsec_unc = 0.0802572; return;}
  else if(hig_mass == 525) { xsec = .5824*.5824*27.1867/1000; xsec_unc = 0.0825803; return;}
  else if(hig_mass == 550) { xsec = .5824*.5824*21.9868/1000; xsec_unc = 0.0849278; return;}
  else if(hig_mass == 575) { xsec = .5824*.5824*17.9062/1000; xsec_unc = 0.087561;  return;}
  else if(hig_mass == 600) { xsec = .5824*.5824*14.6677/1000; xsec_unc = 0.0900693; return;}
  else if(hig_mass == 625) { xsec = .5824*.5824* 12.062/1000; xsec_unc = 0.091959;  return;}
  else if(hig_mass == 650) { xsec = .5824*.5824*9.96406/1000; xsec_unc = 0.094065;  return;}
  else if(hig_mass == 675) { xsec = .5824*.5824*8.28246/1000; xsec_unc = 0.0957436; return;}
  else if(hig_mass == 700) { xsec = .5824*.5824*6.89981/1000; xsec_unc = 0.0982894; return;}
  else if(hig_mass == 725) { xsec = .5824*.5824*5.78355/1000; xsec_unc = 0.0999915; return;}
  else if(hig_mass == 750) { xsec = .5824*.5824* 4.8731/1000; xsec_unc = 0.101211;  return;}
  else if(hig_mass == 775) { xsec = .5824*.5824*4.09781/1000; xsec_unc = 0.104646;  return;}
  else if(hig_mass == 800) { xsec = .5824*.5824*3.46143/1000; xsec_unc = 0.107618;  return;}
  else if(hig_mass == 825) { xsec = .5824*.5824* 2.9337/1000; xsec_unc = 0.108353;  return;}
  else if(hig_mass == 850) { xsec = .5824*.5824* 2.4923/1000; xsec_unc = 0.110016;  return;}
  else if(hig_mass == 875) { xsec = .5824*.5824*2.13679/1000; xsec_unc = 0.112636;  return;}
  else if(hig_mass == 900) { xsec = .5824*.5824*1.80616/1000; xsec_unc = 0.1134;    return;}
  else if(hig_mass == 925) { xsec = .5824*.5824*1.55453/1000; xsec_unc = 0.116949;  return;}
  else if(hig_mass == 950) { xsec = .5824*.5824*1.32692/1000; xsec_unc = 0.117027;  return;}
  else if(hig_mass == 975) { xsec = .5824*.5824*1.12975/1000; xsec_unc = 0.121244;  return;}
  else if(hig_mass ==1000) { xsec = .5824*.5824*0.968853/1000; xsec_unc = 0.126209; return;}
  else{ xsec = 0; xsec_unc = 0;}
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
