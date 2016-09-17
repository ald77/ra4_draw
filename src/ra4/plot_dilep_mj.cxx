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
#include "core/hist1d.hpp"
#include "core/utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  bool single_thread = true;
  TString method = "";
}

void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  double lumi = 2.6;

  bool do_rc = true;
  string rcfolder = "";
  if(do_rc) rcfolder = "reclustered/";

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder
  string lowmet_fndata(bfolder+"/cms2r0/babymaker/babies/"+rcfolder+"2016_06_26/data/merged_met150/");
  string himet_fndata(bfolder+"/cms2r0/babymaker/babies/"+rcfolder+"2016_06_26/data/merged_standard/");
  string himet_fodata(bfolder+"/cms2r0/babymaker/babies/"+rcfolder+"2016_04_29/data/merged_1lht500met200/");
  string himet_fmc_74x(bfolder+"/cms2r0/babymaker/babies/"+rcfolder+"2016_04_29/mc/merged_1lht500met200/");
  string lowmet_fmc(bfolder+"/cms2r0/babymaker/babies/"+rcfolder+"2016_06_14/mc/merged_met150/");
  string himet_fmc(bfolder+"/cms2r0/babymaker/babies/"+rcfolder+"2016_06_14/mc/merged_standard/");

  string baseline("nleps>=1 && ht>500 && met>150 && met<500 && pass && njets>=3");

  Palette colors("txt/colors.txt", "default");

  auto ndata = Process::MakeShared<Baby_full>("Data (2016, 2.6 fb^{-1})", Process::Type::data, kBlack,
    {lowmet_fndata+"/*.root", himet_fndata+"/*.root"},
    baseline+"&&json2p6 && (trig[4]||trig[8]||trig[13]||trig[33])");
  auto odata = Process::MakeShared<Baby_full>("Data (2015, 2.3 fb^{-1})", Process::Type::data, kRed+1, {himet_fodata+"/*.root"},
                                              baseline+" && (trig[4]||trig[8]||trig[28]||trig[14])");


  auto tt1l = Process::MakeShared<Baby_full>("tt 1lep", Process::Type::background, colors("tt_1l"),
    {lowmet_fmc+"*_TTJets*Lept*.root", lowmet_fmc+"*_TTJets_HT*.root",
        himet_fmc+"*_TTJets*Lept*.root", himet_fmc+"*_TTJets_HT*.root"},
    baseline+" && stitch && ntruleps==1");
  auto tt2l = Process::MakeShared<Baby_full>("tt 2lep", Process::Type::background, colors("tt_2l"),
    {lowmet_fmc+"*_TTJets*Lept*.root", lowmet_fmc+"*_TTJets_HT*.root",
        himet_fmc+"*_TTJets*Lept*.root", himet_fmc+"*_TTJets_HT*.root"},
    baseline+" && stitch && ntruleps==2");
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {lowmet_fmc+"*_WJetsToLNu*.root",lowmet_fmc+"*_ST_*.root",
        lowmet_fmc+"*_TTW*.root",lowmet_fmc+"*_TTZ*.root",
        lowmet_fmc+"*DYJetsToLL*.root",lowmet_fmc+"*QCD_HT*.root",
        lowmet_fmc+"*_ZJet*.root",lowmet_fmc+"*_ttHJetTobb*.root",
        lowmet_fmc+"*_TTGJets*.root",lowmet_fmc+"*_TTTT*.root",
        lowmet_fmc+"*_WH_HToBB*.root",lowmet_fmc+"*_ZH_HToBB*.root",
        lowmet_fmc+"*_WWTo*.root",lowmet_fmc+"*_WZ*.root",lowmet_fmc+"*_ZZ_*.root",
        himet_fmc+"*_WJetsToLNu*.root",himet_fmc+"*_ST_*.root",
        himet_fmc+"*_TTW*.root",himet_fmc+"*_TTZ*.root",
        himet_fmc+"*DYJetsToLL*.root",himet_fmc+"*QCD_HT*.root",
        himet_fmc+"*_ZJet*.root",himet_fmc+"*_ttHJetTobb*.root",
        himet_fmc+"*_TTGJets*.root",himet_fmc+"*_TTTT*.root",
        himet_fmc+"*_WH_HToBB*.root",himet_fmc+"*_ZH_HToBB*.root",
        himet_fmc+"*_WWTo*.root",himet_fmc+"*_WZ*.root",himet_fmc+"*_ZZ_*.root"},
    baseline+" && stitch");

  auto tt1l_74x = Process::MakeShared<Baby_full>("tt 1lep", Process::Type::background, colors("tt_1l"),
    {himet_fmc_74x+"*_TTJets*Lept*.root", himet_fmc_74x+"*_TTJets_HT*.root"},
    baseline+" && stitch && ntruleps==1");
  auto tt2l_74x = Process::MakeShared<Baby_full>("tt 2lep", Process::Type::background, colors("tt_2l"),
    {himet_fmc_74x+"*_TTJets*Lept*.root", himet_fmc_74x+"*_TTJets_HT*.root"},
    baseline+" && stitch && ntruleps==2");
  auto other_74x = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {himet_fmc_74x+"*_WJetsToLNu*.root",himet_fmc_74x+"*_ST_*.root",
        himet_fmc_74x+"*_TTW*.root",himet_fmc_74x+"*_TTZ*.root",
        himet_fmc_74x+"*DYJetsToLL*.root",himet_fmc_74x+"*QCD_HT*.root",
        himet_fmc_74x+"*_ZJet*.root",himet_fmc_74x+"*_ttHJetTobb*.root",
        himet_fmc_74x+"*_TTGJets*.root",himet_fmc_74x+"*_TTTT*.root",
        himet_fmc_74x+"*_WH_HToBB*.root",himet_fmc_74x+"*_ZH_HToBB*.root",
        himet_fmc_74x+"*_WWTo*.root",himet_fmc_74x+"*_WZ*.root",himet_fmc_74x+"*_ZZ_*.root"},
    baseline+" && stitch");

  vector<shared_ptr<Process> > procs_2015 = {odata, tt1l_74x, tt2l_74x, other_74x};
  vector<shared_ptr<Process> > lowmet_procs = {ndata, tt1l, tt2l, other};

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);

  vector<PlotOpt> plot_types = {lin_lumi_info};

  PlotMaker pm, pm_2015;

  float minx(0), maxx(460);
  int nbins(static_cast<int>((maxx-minx)/10));
  minx = 25; maxx = 1000; nbins = static_cast<int>((maxx-minx)/75);

  string mjname = "M_{J}^{with lep}";
  if(do_rc) mjname = "M_{J}^{no lep}";

  /////////// 1 lepton
  Axis mj_axis(nbins, minx, maxx, "mj14", mjname+" [GeV]", {250., 400.});
  pm_2015.Push<Hist1D>(mj_axis, "met>200&&met<500 && mt<=140 && nleps==1 && nbm>=1 && njets>=3",
                       procs_2015, plot_types);
  if(!do_rc){
    pm_2015.Push<Hist1D>(mj_axis, "met>200&&met<500 && nels==1 && nmus==1 && nbm<=2 && njets>=5",
                         procs_2015, plot_types);
    pm_2015.Push<Hist1D>(mj_axis, "met>200&&met<500 && mt<=140 && nleps==1 && nbm>=1 && njets>=6",
                         procs_2015, plot_types);
  }
  pm.Push<Hist1D>(mj_axis, "met>150&&met<500 && mt<=140 && nleps==1 && nbm>=1 && njets>=3",
                  lowmet_procs, plot_types);
  if(!do_rc){
    pm.Push<Hist1D>(mj_axis, "met>150&&met<500 && nels==1 && nmus==1 && nbm<=2 && njets>=5",
                    lowmet_procs, plot_types);
    pm.Push<Hist1D>(mj_axis, "met>150&&met<500 && mt<=140 && nleps==1 && nbm>=1 && njets>=6",
                    lowmet_procs, plot_types);
  }

  /////////// 2 leptons
  pm_2015.Push<Hist1D>(mj_axis, "met>200&&met<500 && nels==1 && nmus==1 && nbm<=2 && njets>=3",
                       procs_2015, plot_types);
  if(!do_rc){
    pm_2015.Push<Hist1D>(mj_axis, "met>200&&met<500 && mt<=140 && nleps==1 && nbm>=1 && njets>=5",
                         procs_2015, plot_types);
    pm_2015.Push<Hist1D>(mj_axis, "met>200&&met<500 && nels==1 && nmus==1 && nbm<=2 && njets>=6",
                         procs_2015, plot_types);
  }

  pm.Push<Hist1D>(mj_axis, "met>150&&met<500 && nels==1 && nmus==1 && nbm<=2 && njets>=3",
                  lowmet_procs, plot_types);
  if(!do_rc){
    pm.Push<Hist1D>(mj_axis, "met>200&&met<500 && nels==1 && nmus==1 && nbm<=2 && njets==3",
                    lowmet_procs, plot_types);
    pm.Push<Hist1D>(mj_axis, "met>200&&met<500 && nels==1 && nmus==1 && nbm<=2 && njets==4",
                    lowmet_procs, plot_types);
    pm.Push<Hist1D>(mj_axis, "met>150&&met<500 && mt<=140 && nleps==1 && nbm>=1 && njets>=5",
                    lowmet_procs, plot_types);
    pm.Push<Hist1D>(mj_axis, "met>150&&met<500 && nels==1 && nmus==1 && nbm<=2 && njets>=6",
                    lowmet_procs, plot_types);
  }

  if(single_thread) pm.multithreaded_ = false;
  pm.MakePlots(lumi);
  if(single_thread) pm_2015.multithreaded_ = false;
  pm_2015.MakePlots(2.3);

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
