// Plot some kinematic distributions of high MET events to change for
// reconstruction errors and other pathological events

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>
#include <chrono>
#include <string>

#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>

#include "TError.h"

#include "core/utilities.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/event_scan.hpp"
#include "core/styles.hpp"
#include "core/functions.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Functions;

namespace{
  bool single_thread = false;
}



void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  string base_path = "";
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname,"compute-")){
    base_path = "/net/cms2";
  }

  string mc_dir = base_path+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_stdnj5/";
  string data_dir = base_path+"/cms2r0/babymaker/babies/2016_11_08/data/merged_database_met150nj2/";

  Palette colors("txt/colors.txt","default");

  auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {mc_dir+"*_TTJets*Lept*.root", mc_dir+"*_TTJets_HT*.root"}, "ntruleps<=1&&stitch");
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {mc_dir+"*_TTJets*Lept*.root", mc_dir+"*_TTJets_HT*.root"}, "ntruleps>=2&&stitch");
  auto wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {mc_dir+"*_WJetsToLNu*.root"},"stitch");
  auto single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {mc_dir+"*_ST_*.root"});
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {mc_dir+"*_TTWJets*.root", mc_dir+"*_TTZTo*.root"});
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {mc_dir+"*DYJetsToLL*.root", mc_dir+"*_QCD_HT*.root",
        mc_dir+"*_ZJet*.root", mc_dir+"*_WWTo*.root",
        mc_dir+"*ggZH_HToBB*.root", mc_dir+"*ttHJetTobb*.root",
        mc_dir+"*_TTGJets*.root", mc_dir+"*_TTTT_*.root",
        mc_dir+"*_WH_HToBB*.root", mc_dir+"*_WZTo*.root",
        mc_dir+"*_ZH_HToBB*.root", mc_dir+"_ZZ_*.root"});

  auto data = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {data_dir+"*.root"},"pass&&trig_ra4&&(json12p9||mj14<=400.||mt<=140.||njets<=5||nleps!=1)");

  vector<shared_ptr<Process> > processes = {data, tt1l, tt2l, wjets, single_t, ttv, other};

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  vector<PlotOpt> all_opts = {log_lumi, lin_lumi};
  vector<PlotOpt> lin_opts = {lin_lumi};

  NamedFunc baseline = "nleps>=1&&st>500&&met>500.&&njets>=5&&mj14>250.";
  NamedFunc abcd = baseline && "nleps==1&&njets>=6&&nbm>=1"; abcd.Name("abcd");
  NamedFunc abcd_highmt = abcd && "mt>140."; abcd_highmt.Name("abcd_highmt");
  NamedFunc nj5 = baseline && "nleps==1&&njets==5&&nbm>=1&&!(mj14>400&&mt>140)"; nj5.Name("nj5");
  NamedFunc nj5_highmt = nj5 && "mt>140."; nj5_highmt.Name("nj5_highmt");
  NamedFunc dilep = baseline && "nbm<3&&((nleps==2&&nveto==0)||(nleps==1&&nveto==1))"; dilep.Name("dilep");
  NamedFunc dilep_highmt = dilep && "mt>140."; dilep_highmt.Name("dilep_highmt");
  NamedFunc vhigh_met = abcd && "met>2000"; vhigh_met.Name("vhigh_met");

  vector<NamedFunc> all_cuts = {abcd, abcd_highmt, nj5, nj5_highmt, dilep, dilep_highmt, vhigh_met};

  PlotMaker pm;
  for(const auto &cut: all_cuts){
    pm.Push<Hist1D>(Axis(25, 500., 3000., "met", "MET [GeV]"),
                    cut, processes, all_opts);
    pm.Push<Hist1D>(Axis(40, 0., 4000., "st", "S_{T} [GeV]"),
                    cut, processes, all_opts);
    pm.Push<Hist1D>(Axis(41, -0.5, 40.5, "npv", "NPV"),
                    cut, processes, all_opts);
    pm.Push<Hist1D>(Axis(16, -0.5, 15.5, "njets", "N_{jets}"),
                    cut, processes, all_opts);
    pm.Push<Hist1D>(Axis(40, 0., 2000., "jets_pt[0]", "Highest Jet p_{T} [GeV]"),
                    cut, processes, all_opts);
    pm.Push<Hist1D>(Axis(30, 0., 1500., "mt", "m_{T} [GeV]"),
                    cut, processes, all_opts);
    pm.Push<Hist1D>(Axis(30, 0., 1500., "mj14", "M_{J} [GeV]"),
                    cut, processes, all_opts);
    pm.Push<Hist1D>(Axis(30, 0, 1200., "mus_pt", "#mu p_{T} [GeV]"),
                    cut, processes, all_opts);
    pm.Push<Hist1D>(Axis(30, 0, 1200., "els_pt", "#mu p_{T} [GeV]"),
                    cut, processes, all_opts);
    pm.Push<Hist1D>(Axis(30, 0., acos(-1.), min_dphi_met_jet, "Min #Delta#phi(jet,MET)"),
                    cut, processes, all_opts);
    pm.Push<Hist1D>(Axis(30, 0., acos(-1.), max_dphi_met_jet, "Max #Delta#phi(jet,MET)"),
                    cut, processes, all_opts);
    pm.Push<Hist2D>(Axis(20, -acos(-1.), acos(-1.), "jets_phi", "Jet #phi"),
                    Axis(20, -3., 3., "jets_eta", "Jet #eta"),
                    cut, processes, lin_opts);
    pm.Push<Hist2D>(Axis(20, 500., 3000., "met", "MET [GeV]"),
                    Axis(50, 0., acos(-1.), min_dphi_met_jet, "Min #Delta#phi(jet,MET)"),
                    cut, processes, lin_opts);
    pm.Push<Hist2D>(Axis(20, 500., 3000., "met", "MET [GeV]"),
                    Axis(50, 0., acos(-1.), max_dphi_met_jet, "Max #Delta#phi(jet,MET)"),
                    cut, processes, lin_opts);
    pm.Push<EventScan>(cut.Name(), cut, vector<NamedFunc>{
        "type", "run", "lumiblock", "event",
          "nmus", "nels", "nveto", "njets", "mj14", "mt",
          "met", "met_phi",
          min_dphi_met_jet, max_dphi_met_jet,
          "jets_pt", "jets_eta", "jets_phi", "jets_m", "jets_islep", "jets_csv",
          "mus_pt", "mus_eta", "mus_phi",
          "els_pt", "els_sceta", "els_phi",
          }, processes, 10);
                                         
  }

  pm.multithreaded_ = !single_thread;
  pm.MakePlots(12.9);
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"single_thread", no_argument, 0, 's'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "l:su", long_options, &option_index);
    if(opt == -1) break;

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
	exit(1);
      }
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
