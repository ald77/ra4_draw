#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include <unistd.h>
#include <getopt.h>

#include "TError.h"
#include "TColor.h"

#include "baby.hpp"
#include "process.hpp"
#include "named_func.hpp"
#include "plot_maker.hpp"
#include "plot_opt.hpp"
#include "palette.hpp"
#include "table.hpp"
#include "histo_stack.hpp"

using namespace std;
using namespace PlotOptTypes;

int main(){
  gErrorIgnoreLevel = 6000;

  double lumi = 2.7;

  string mc_folder = "/net/cms29/cms29r0/cawest/skims/ht1200/";
  string sig_folder = "/net/cms9/cms9r0/rohan/babies/2016_07_13/T1tbs/split/renorm/";

  Palette colors("txt/colors.txt", "default");

  // Background samples
  auto tt = Process::MakeShared<Baby_full>("t#bar{t}", Process::Type::background, colors("tt_1l"),
    {mc_folder+"*TTJets*Lept*"},
    "ntruleps>=1");

  auto qcd = Process::MakeShared<Baby_full>("QCD", Process::Type::background, colors("tt_1l"),
    {mc_folder+"*QCD_HT*",mc_folder+"*TTJets_TuneCUETP8M1_13TeV-madgraphMLM*"},
    "ntruleps==0");

  auto wjets = Process::MakeShared<Baby_full>("W+Jets", Process::Type::background, colors("tt_1l"),
    {mc_folder+"*WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8*"});

  auto singlet = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("tt_1l"),
    {mc_folder+"*ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1*",mc_folder+"*ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1*",
        mc_folder+"*ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1*", mc_folder+"*ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1*",
        mc_folder+"*ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1*"});

  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("tt_1l"),
    {mc_folder+"*DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8*",mc_folder+"*TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8*",
        mc_folder+"*TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8*",mc_folder+"*TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8*",
        mc_folder+"*TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8*",mc_folder+"*ttHJetTobb_M125_13TeV_amcatnloFXFX_madspin_pythia8*",
        mc_folder+"*TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8*",mc_folder+"*WJetsToQQ_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8*",
        mc_folder+"*ZJetsToQQ_HT600toInf_13TeV-madgraph*"});

  //Signal Samples
  auto tbs750 = Process::MakeShared<Baby_full>("m_{glu}=750$ GeV$", Process::Type::signal, colors("tt_1l"),
    {sig_folder+"*RPV_mGluino-750*"});

  auto tbs1200 = Process::MakeShared<Baby_full>("m_{glu}=1200$ GeV$", Process::Type::signal, colors("tt_1l"),
    {sig_folder+"*RPV_mGluino-1200*"});

  auto tbs1500 = Process::MakeShared<Baby_full>("m_{glu}=1500$ GeV$", Process::Type::signal, colors("tt_1l"),
    {sig_folder+"*RPV_mGluino-1500*"});

  vector<shared_ptr<Process> > samples = {other, singlet, wjets, qcd, tt, tbs750, tbs1200, tbs1500};

  PlotMaker pm;
  pm.Push<Table>("rpv_cutflow", vector<TableRow>{
      TableRow("$N_{leps}=0$"),
        TableRow("$H_{T}>1500$","nleps==0&&ht>1500",0,0,"weight*w_pu_rpv/eff_trig"),
        TableRow("$N_{b}\\geq1$","nleps==0&&ht>1500&&nbm>=1",0,0,"weight*w_pu_rpv/eff_trig"),
        TableRow("$N_{jets}\\geq4$","nleps==0&&ht>1500&&nbm>=1&&njets>=4",0,0,"weight*w_pu_rpv/eff_trig"),
        TableRow("$M_{J}>500$","nleps==0&&ht>1500&&nbm>=1&&njets>=4&&mj>500",0,1,"weight*w_pu_rpv/eff_trig"),
        TableRow("$M_{J}>800$","nleps==0&&ht>1500&&nbm>=1&&njets>=4&&mj>800",0,0,"weight*w_pu_rpv/eff_trig"),
        TableRow("$N_{jets}\\geq10$","nleps==0&&ht>1500&&nbm>=1&&njets>=10&&mj>800",0,0,"weight*w_pu_rpv/eff_trig"),
        TableRow("$N_{b}\\geq3$","nleps==0&&ht>1500&&nbm>=3&&njets>=8&&mj>800",0,1,"weight*w_pu_rpv/eff_trig"),
        TableRow("$N_{leps}=1$"),
        TableRow("$H_{T}>1200$","nleps==1&&ht>1200",1,0,"weight*w_pu_rpv/eff_trig"),
        TableRow("$N_{b}\\geq1$","nleps==1&&ht>1200&&nbm>=1",0,0,"weight*w_pu_rpv/eff_trig"),
        TableRow("$N_{jets}\\geq4$","nleps==1&&ht>1200&&nbm>=1&&njets>=4",0,0,"weight*w_pu_rpv/eff_trig"),
        TableRow("$M_{J}>500$","nleps==1&&ht>1200&&nbm>=1&&njets>=4&&mj>500",0,1,"weight*w_pu_rpv/eff_trig"),
        TableRow("$M_{J}>800$","nleps==1&&ht>1200&&nbm>=1&&njets>=4&&mj>800",0,0,"weight*w_pu_rpv/eff_trig"),
        TableRow("$N_{jets}\\geq8$","nleps==1&&ht>1200&&nbm>=1&&njets>=8&&mj>800",0,0,"weight*w_pu_rpv/eff_trig"),
        TableRow("$N_{b}\\geq3$","nleps==1&&ht>1200&&nbm>=3&&njets>=8&&mj>800",0,1,"weight*w_pu_rpv/eff_trig"),
        },samples,false,true);

  pm.MakePlots(lumi);
}
