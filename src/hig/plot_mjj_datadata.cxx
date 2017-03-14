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
#include "hig/hig_functions.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace {
  double lumi = 35.9;
}

int main(){
  gErrorIgnoreLevel = 6000;
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms29"; // In laptops, you can't create a /net folder

  Palette colors("txt/colors.txt", "default");
  
  string fdata = bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_higdata_higloose/";
  string fsig = bfolder+"/cms2r0/babymaker/babies/2017_02_26/TChiHH/merged_higmc_higloose/";

  NamedFunc baseline = "pass_ra2_badmu && met/met_calo<5 && nvleps==0 && ntks==0 && !low_dphi && njets>=4 && njets<=5 && nbdt>=2 && met>150 && higd_drmax<=2.2 && higd_am<=200 && higd_dm <= 40";

  auto data_3b = Process::MakeShared<Baby_full>("Data 3b", Process::Type::data, kBlack,
    {fdata+"*.root"},baseline && Higfuncs::trig_hig>0. && "pass && nbdt>=2&&nbdm==3&&nbdl==3");
  auto data_4b = Process::MakeShared<Baby_full>("Data 4b", Process::Type::data, kBlack,
    {fdata+"*.root"},baseline && Higfuncs::trig_hig>0. && "pass && nbdt>=2&&nbdm>=3&&nbdl>=4");
  auto data_2b = Process::MakeShared<Baby_full>("Data 2b", Process::Type::background, kBlack,
    {fdata+"*.root"},baseline && Higfuncs::trig_hig>0. && "pass && nbdt==2&&nbdm==2");
  data_2b->SetFillColor(kWhite);
  data_2b->SetLineColor(kAzure-2);//kBlue-7);
  data_2b->SetLineWidth(2);

  string filters = "pass_ra2_badmu&&pass_goodv&&pass_ecaldeadcell&&pass_hbhe&&pass_hbheiso&&pass_fsmet";
  auto tchi225_3b = Process::MakeShared<Baby_full>("TChiHH(225,1) 3b", Process::Type::signal, kRed,
    {fsig+"*TChiHH*225*.root"}, baseline && filters+"&& nbdt>=2&&nbdm==3&&nbdl==3");
  tchi225_3b->SetLineWidth(2);
  tchi225_3b->SetLineStyle(2);
  auto tchi225_4b = Process::MakeShared<Baby_full>("TChiHH(225,1) 4b", Process::Type::signal, kRed,
    {fsig+"*TChiHH*225*.root"}, baseline && filters+"&& nbdt>=2&&nbdm>=3&&nbdl>=4");
  tchi225_4b->SetLineWidth(2);
  tchi225_4b->SetLineStyle(2);
  auto tchi400_3b = Process::MakeShared<Baby_full>("TChiHH(400,1) 3b", Process::Type::signal, kRed,
    {fsig+"*TChiHH*400*.root"}, baseline && filters+"&& nbdt>=2&&nbdm==3&&nbdl==3");
  tchi400_3b->SetLineWidth(2);
  auto tchi400_4b = Process::MakeShared<Baby_full>("TChiHH(400,1) 4b", Process::Type::signal, kRed,
    {fsig+"*TChiHH*400*.root"}, baseline && filters+"&& nbdt>=2&&nbdm>=3&&nbdl>=4");
  tchi400_4b->SetLineWidth(2);

  vector<shared_ptr<Process> > data3b_procs = {data_3b,tchi225_3b,data_2b,tchi400_3b};
  vector<shared_ptr<Process> > data4b_procs = {data_4b,tchi225_4b,data_2b,tchi400_4b};

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm); 
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  vector<PlotOpt> lin = {lin_lumi};

  NamedFunc wgt = Higfuncs::weight_higd * Higfuncs::eff_higtrig;

  PlotMaker pm;

  //data-to-data
  //  vector<string> metbins = {"met>150 && met<=500", "met>150 && met<=200", "met>200 && met<=350", "met>350 && met<=500", "met>200 && met<=500","met>500","met>200","met>350"};
  vector<string> metbins = {"met>150 && met<=200", "met>200"};
  for (auto &imet: metbins){
    string metlabel = CodeToRootTex(imet)+" GeV"; ReplaceAll(metlabel,"E","p");
    pm.Push<Hist1D>(Axis(20, 0, 200, "higd_am", "#LTm#GT [GeV]",{100.,140.}),
		    imet, data3b_procs, lin).Weight(wgt).Tag("data3b").RatioTitle("Data 3b","Data 2b")
      .RightLabel(metlabel).YAxisZoom(0.9);

    pm.Push<Hist1D>(Axis(20, 0, 200, "higd_am", "#LTm#GT [GeV]", {100.,140.}),
		    imet, data4b_procs, lin).Weight(wgt).Tag("data4b").RatioTitle("Data 4b","Data 2b")
      .RightLabel(metlabel).YAxisZoom(0.9);
  } 


  pm.min_print_ = true;
  pm.MakePlots(lumi);

}
