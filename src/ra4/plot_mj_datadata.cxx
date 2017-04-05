#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "TError.h"
#include "TColor.h"
#include "TVector2.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist1d.hpp"
#include "core/utilities.hpp"
#include "core/functions.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace {
  double lumi = 35.9;
  bool paper = true;
}

int main(){
  gErrorIgnoreLevel = 6000;
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  Palette colors("txt/colors.txt", "default");

  string ntupletag = "*.root";
  
  string fdata = bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_database_stdnj5/";
  
  string lsp = "{#lower[-0.1]{#tilde{#chi}}#lower[0.2]{#scale[0.95]{^{0}}}#kern[-1.3]{#scale[0.95]{_{1}}}}";
  string t1t_label = "#scale[0.95]{#tilde{g}#kern[0.2]{#tilde{g}}, #tilde{g}#rightarrowt#kern[0.18]{#bar{t}}#kern[0.18]"+lsp;
  string mt = "m#lower[-0.1]{_{T}}";
  string etmiss = "E#lower[-0.1]{_{T}}#kern[-0.25]{#scale[1.15]{#lower[0.15]{^{miss}}}}";
  string nj6 = "N#kern[-0.1]{#scale[1.15]{#lower[-0.15]{_{jets}}}}#kern[0.2]{#geq}#kern[0.2]{6}";
  string nb2 = "N#kern[-0.15]{#scale[1.15]{#lower[-0.15]{_{b}}}}#kern[0.2]{#geq}#kern[0.2]{2}";

  auto data_highmt = Process::MakeShared<Baby_full>("Data, "+mt+" > 140 GeV", Process::Type::data, kBlack,
    {fdata+ntupletag},"pass && trig_ra4 && st>500 && mt>140 && nleps==1 && nveto==0 && njets>=6 && nbm>=1 && met/met_calo<5.0 && pass_ra2_badmu");
  auto data_lowmt = Process::MakeShared<Baby_full>("Data, "+mt+" #leq 140 GeV", Process::Type::background, kBlack,
    {fdata+ntupletag},"pass && trig_ra4 && st>500 && mt<=140 && nleps==1 && nveto==0 && njets>=6 && nbm>=1 && met/met_calo<5.0 && pass_ra2_badmu");
  data_lowmt->SetFillColor(kWhite);
  data_lowmt->SetLineColor(kBlue-7);
  data_lowmt->SetLineWidth(2);

  auto data2lveto = Process::MakeShared<Baby_full>("Data 2l or l+trk", Process::Type::data, kBlue+2,
    {fdata+ntupletag},
    "pass && trig_ra4 && st>500 && ((nleps==2 && njets>=5 && nbm<=2) || (nleps==1 && nveto==1 && njets>=6 && nbm>=1 && mt>140)) && met/met_calo<5.0 && pass_ra2_badmu");

  auto data2l = Process::MakeShared<Baby_full>("Data 2l", Process::Type::data, kMagenta+3,
    {fdata+ntupletag},
    "pass && trig_ra4 && st>500 && (nleps==2 && njets>=5 && nbm<=2) && met/met_calo<5.0 && pass_ra2_badmu");

  auto t1tttt = Process::MakeShared<Baby_full>(t1t_label+" (1800,100)}", Process::Type::signal, colors("t1tttt"),
    {bfolder+"/cms2r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/*SMS-T1tttt_mGluino-1800_mLSP-100_*.root"},"st>500 && mt>140 && nleps==1 && nveto==0 && njets>=6 && nbm>=1 && met/met_calo<5.0 && pass_ra2_badmu");
  t1tttt->SetLineWidth(2);
  auto t1ttttc = Process::MakeShared<Baby_full>(t1t_label+" (1400,1000)}", Process::Type::signal, colors("t1tttt"),
    {bfolder+"/cms2r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/*SMS-T1tttt_mGluino-1400_mLSP-1000_*.root"},"st>500 && mt>140 && nleps==1 && nveto==0 && njets>=6 && nbm>=1 && met/met_calo<5.0 && pass_ra2_badmu");
  t1ttttc->SetLineWidth(2);
  t1ttttc->SetLineStyle(2);


  auto data = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {fdata+ntupletag},"pass && trig_ra4");


  vector<shared_ptr<Process> > data1l_procs = {data_highmt,data_lowmt,t1tttt,t1ttttc};
  vector<shared_ptr<Process> > data2lveto_procs = {data2lveto, data_lowmt};
  vector<shared_ptr<Process> > data2l_procs = {data2l, data_lowmt};

  string style = "Preliminary";
  if(paper) style = "PRLPaper";
  PlotOpt log_lumi("txt/plot_styles.txt", style);
  log_lumi.Title(TitleType::data)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm)
    .RatioMaximum(1.86);
  if(paper){
    log_lumi.Bottom(BottomType::off);
  }else{
    log_lumi=log_lumi.Title(TitleType::preliminary);
  }
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes);
  if(paper){
    log_shapes = log_lumi().Stack(StackType::shapes)
      .Bottom(BottomType::off)
      .ShowBackgroundError(false);
  }
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info);
  vector<PlotOpt> log = {log_lumi_info,log_lumi};
  vector<PlotOpt> lin = {lin_lumi};


  PlotMaker pm;

  //data-to-data
  //  vector<string> metbins = {"met>150 && met<=500", "met>150 && met<=200", "met>200 && met<=350", "met>350 && met<=500", "met>200 && met<=500","met>500","met>200","met>350"};
  vector<string> metbins = { "met>200 && met<=350", "met>350 && met<=500","met>500","met>350"};
  if (paper) metbins = { "met>200 && met<=350","met>350"};
  for (auto &imet: metbins){
    if (paper) {
      string metlabel = CodeToRootTex(imet);
      ReplaceAll(metlabel,"E_{T}^{miss}",etmiss);
      pm.Push<Hist1D>(Axis(19, 50.,1000., "mj14", "M_{J} [GeV]", {250.,400.}),
        imet + "&&nbm>=2", data1l_procs, lin).Tag("data1l2b").RatioTitle("Data, "+mt+" > 140 GeV","Data, "+mt+" #leq 140 GeV")
      .RightLabel({metlabel+" GeV",nj6+", "+nb2}).YAxisZoom(0.85);
    } else {
      pm.Push<Hist1D>(Axis(20, 0.,1000., "mj14", "M_{J} [GeV]",{250.,400.}),
		      imet + "&&nbm==1", data1l_procs, lin).Tag("data1l1b").RatioTitle("Data, "+mt+" > 140 GeV","Data, "+mt+" #leq 140 GeV")
        .RightLabel({CodeToRootTex(imet)+" GeV",CodeToRootTex("njets>=6&&nbm==1")}).YAxisZoom(0.93);
      pm.Push<Hist1D>(Axis(20, 0.,1000., "mj14", "M_{J} [GeV]",{250.,400.}),
		      imet + "&&nbm>=2", data1l_procs, lin).Tag("data1l2b").RatioTitle("Data, "+mt+" > 140 GeV","Data, "+mt+" #leq 140 GeV")
	.RightLabel({CodeToRootTex(imet)+" GeV",CodeToRootTex("njets>=6&&nbm>=2")}).YAxisZoom(0.93);
      pm.Push<Hist1D>(Axis(20, 0.,1000., "mj14", "M_{J} [GeV]", {250.,400.}),
  		    imet,  data2l_procs, lin).Tag("data2l").RatioTitle("Data 2l","Data 1l, "+mt+" #leq 140 GeV");

      pm.Push<Hist1D>(Axis(20, 0.,1000., "mj14", "M_{J} [GeV]", {250.,400.}),
  		    imet, data2lveto_procs, lin).Tag("data2lveto").RatioTitle("Data 2l","Data 1l, "+mt+" #leq 140 GeV");
    }
  } 


  pm.min_print_ = true;
  pm.MakePlots(lumi);

}
