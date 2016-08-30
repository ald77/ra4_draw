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

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/histo_stack.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  bool single_thread = false;
}

int main(int argc, char *argv[]){
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  double lumi = 20;

  string trig_skim_mc = "/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/";
  string trig_skim_signal = "/cms2r0/babymaker/babies/2016_08_10/TChiHH/merged_higmc_higloose/";

  Palette colors("txt/colors.txt", "default");
  auto qcd = Process::MakeShared<Baby_full>("Inclusive QCD", Process::Type::signal, colors("qcd"),
    {trig_skim_mc+"*_QCD_HT*00_Tune*.root", trig_skim_mc+"*_QCD_HT*Inf_Tune*.root"});
  auto Bqcd = Process::MakeShared<Baby_full>("B-enriched QCD", Process::Type::signal, kAzure,
    {trig_skim_mc+"*_QCD_bEnriched*.root", trig_skim_mc+"*_QCD_HT*BGenFilter*.root"});

  vector<shared_ptr<Process> > full_trig_skim = {qcd, Bqcd};

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
  vector<PlotOpt> all_plot_types = {log_lumi_info, log_shapes_info};

  PlotMaker pm;

  vector<TString> nbcuts;
  nbcuts.push_back(" nbt >= 2 && nbl <= 3");
  nbcuts.push_back(" nbt >= 2 && nbl >= 4");
  TString metskim("nbm>=2&&njets>=4&&njets<=5&&met>100&&nvleps==0");
  TString trkskim("njets>=4&&njets<=5&&met>250&&nvleps==0");
  TString skim("njets>=4&&njets<=5&&met>250&&nvleps==0&&ntks==0");
  TString A("&&");
  TString DeltaR("hig_drmax < 2.2");
  TString AverageM("hig_am > 100 && hig_am < 140");
  TString DeltaM("hig_dm < 40");
  TString LDP("!low_dphi");

  pm.Push<HistoStack>(HistoDef(56,0,2800,"ht_isr_me", "H_{T}^{ISR} [GeV]", metskim, "weight", {9999.}),
			full_trig_skim, all_plot_types);

  pm.Push<Table>("QCD_cutflow", vector<TableRow>{
	  TableRow("$\\text{No b cut}$", "1"),
	  TableRow("$\\text{1M b-tags}$", "nbm>=1"),
	  TableRow("$\\text{2M b-tags}$", "nbm>=2"),
	  TableRow("$\\text{3M b-tags}$", "nbm>=3"),
	  TableRow("$\\text{4M b-tags}$", "nbm>=4"),
	  TableRow("$\\text{5M b-tags}$", "nbm>=5"),
	  TableRow("$\\text{1T b-tags}$", "nbt>=1"),
	  TableRow("$\\text{2T b-tags}$", "nbt>=2"),
	  TableRow("$\\text{3T b-tags}$", "nbt>=3"),
	  TableRow("$\\text{4T b-tags}$", "nbt>=4"),
	  },full_trig_skim, 0);

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
