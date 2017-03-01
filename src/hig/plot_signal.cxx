#include "core/test.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TColor.h"
#include "TVector2.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "hig/hig_functions.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

vector<unsigned> higidx(const Baby &b);

namespace{
  bool single_thread = true;
  double lumi = 35.9;

  vector<string> sigm = {"225","300", "400","700","1000"}; 
  vector<int> sig_colors = {kGreen+1, kTeal-4, kRed, kBlue, kOrange}; // need sigm.size() >= sig_colors.size()
  // extended list
  // vector<string> sigm = {"127","300","400","500","600","700","850","1000"}; 
  // vector<int> sig_colors = {kAzure, kGreen+2, kRed, kViolet-6, kYellow, kMagenta+1, kCyan+1, kOrange-3};
}

int main(int argc, char *argv[]){
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);


  string foldersig = "/cms2r0/babymaker/babies/2017_01_27/TChiHH/merged_higmc_unskimmed/";

  Palette colors("txt/colors.txt", "default");
  vector<shared_ptr<Process> > procs;
  for (unsigned isig(0); isig<sigm.size(); isig++)
    procs.push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
      sig_colors[isig], {foldersig+"*TChiHH_mGluino-"+sigm[isig]+"*.root"}, "1"));

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::off)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info);
  vector<PlotOpt> linplt = {lin_shapes_info};
  vector<PlotOpt> logplt = {log_shapes_info};

  PlotMaker pm;

  NamedFunc lead_higs_pt("lead_higs_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float higpt(0);
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (b.mc_id()->at(i)!=25) continue;
      if (b.mc_pt()->at(i)>higpt) higpt = b.mc_pt()->at(i);
    }
    return higpt;
  });

  NamedFunc rel_lsp_dpt("rel_lsp_dpt",[](const Baby &b) -> NamedFunc::ScalarType{
    vector<unsigned> lspidx;
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (b.mc_id()->at(i)==1000022) lspidx.push_back(i);
      if (lspidx.size()>1) break;
    }
    float max_pt = max(b.mc_pt()->at(lspidx[0]),b.mc_pt()->at(lspidx[1]));
    return fabs(b.mc_pt()->at(lspidx[0])-b.mc_pt()->at(lspidx[1]))/max_pt;
  });

  NamedFunc hig_dphi("hig_dphi",[](const Baby &b) -> NamedFunc::ScalarType{
    vector<unsigned> idx = higidx(b);
    return deltaPhi(b.mc_phi()->at(idx[0]),b.mc_phi()->at(idx[1]));
  });  

  NamedFunc hig_dr("hig_dr",[](const Baby &b) -> NamedFunc::ScalarType{
    vector<unsigned> idx = higidx(b);
    return deltaR(b.mc_eta()->at(idx[0]), b.mc_phi()->at(idx[0]), 
                  b.mc_eta()->at(idx[1]), b.mc_phi()->at(idx[1]));
  });  

  string baseline("njets>=4 && njets<=5 && met>150");

  string c2b = "nbdt==2&&nbdm==2";
  string c3b = "nbdt>=2&&nbdm==3&&nbdl==3";
  string c4b = "nbdt>=2&&nbdm>=3&&nbdl>=4";

  NamedFunc wgt = weight_higd;// * eff_higtrig;

  pm.Push<Hist1D>(Axis(5,0.5,5.5,higd_bcat, "b-tag categories (TTML)"), 
    "met>150" && higd_bcat>0., procs, linplt).Weight(wgt);  
  pm.Push<Hist1D>(Axis(17,150,1000,"met", "E_{T}^{miss} [GeV]", {150, 200, 300, 450}),
    baseline, procs, logplt).Weight(wgt);
  pm.Push<Hist1D>(Axis(20,0,1000,lead_higs_pt, "Leading Higgs p_{T} [GeV]"),
    "met>150", procs, linplt).Weight(wgt);
  pm.Push<Hist1D>(Axis(20,0,4,"higd_drmax", "#DeltaR_{max} [GeV]", {2.2}),
    baseline && higd_bcat>=3 && "higd_am>100 && higd_am<=140", procs, linplt).Weight(wgt);
  pm.Push<Hist1D>(Axis(20,0,4,"higd_drmax", "#DeltaR_{max} [GeV]", {2.2}),
    baseline, procs, linplt).Weight(wgt);
  pm.Push<Hist1D>(Axis(20,0,200,"higd_am", "#LTm#GT [GeV]", {100, 140}),
    baseline && c2b, procs, linplt).Weight(wgt);
  pm.Push<Hist1D>(Axis(20,0,200,"higd_am", "#LTm#GT [GeV]", {100, 140}),
    baseline && c3b, procs, linplt).Weight(wgt);
  pm.Push<Hist1D>(Axis(20,0,200,"higd_am", "#LTm#GT [GeV]", {100, 140}),
    baseline && c4b, procs, linplt).Weight(wgt);

  pm.Push<Hist1D>(Axis(20,0,200,"higd_am", "#LTm#GT [GeV]", {100, 140}),
    baseline && "higd_dm<=40 &&"+c2b, procs, linplt).Weight(wgt);
  pm.Push<Hist1D>(Axis(20,0,200,"higd_am", "#LTm#GT [GeV]", {100, 140}),
    baseline && "higd_dm<=40 &&"+c3b, procs, linplt).Weight(wgt);
  pm.Push<Hist1D>(Axis(20,0,200,"higd_am", "#LTm#GT [GeV]", {100, 140}),
    baseline && "higd_dm<=40 &&"+c4b, procs, linplt).Weight(wgt);

  pm.Push<Hist1D>(Axis(16,0,3.2,hig_dphi, "#Delta#phi(h1, h2) [GeV]"),
    baseline, procs, linplt).Weight(wgt);
  pm.Push<Hist1D>(Axis(20,0,1,rel_lsp_dpt, "#Delta p_{T}(G_{1}, G_{2})/Max(p_{T}(G_{1}), p_{T}(G_{2}))[GeV]"),
    baseline, procs, linplt).Weight(wgt);
  pm.Push<Hist1D>(Axis(20,0,1,rel_lsp_dpt, "#Delta p_{T}(G_{1}, G_{2})/Max(p_{T}(G_{1}), p_{T}(G_{2}))[GeV]"),
    "njets>=4 && njets<=5", procs, linplt).Weight(wgt);

  pm.Push<Hist1D>(Axis(16,0,3.2,"dphi2", "#Delta#Phi_{2}"),
    baseline, procs, linplt).Weight(wgt);
  pm.Push<Hist1D>(Axis(16,0,3.2,"dphi2", "#Delta#Phi_{2}"),
    "njets>=4 && njets<=5", procs, linplt).Weight(wgt);


  if(single_thread) pm.multithreaded_ = false;
  pm.min_print_ = true;
  pm.MakePlots(lumi);
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 's'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s", long_options, &option_index);

    if( opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      single_thread = true;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(false){
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

vector<unsigned> higidx(const Baby &b){
  vector<unsigned> idx;
  for (unsigned i(0); i<b.mc_pt()->size(); i++){
    if (b.mc_id()->at(i)==25) idx.push_back(i);
    if (idx.size()>1) break;
  }
  return idx;
}

// vector<unsigned> bidx(const Baby &b){
//   vector<unsigned> idx;
//   for (unsigned i(0); i<b.mc_pt()->size(); i++){
//     if (b.mc_id()->at(i)==25) higidx.push_back(i);
//     if (higidx.size()>1) break;
//   }
//   return idx;
// }
