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


template<typename T>
shared_ptr<Process> Proc(const string process_name, Process::Type type,
                         int color, const set<string> &files, const string &cut = "1"){
  return make_shared<Process>(process_name, type, color,
                              unique_ptr<Baby>(new T(files)),
                              cut);
}

int main(){
  gErrorIgnoreLevel = 6000;

  double lumi = 2.7;

  string mc_folder = "/net/cms29/cms29r0/cawest/skims/ht1200/";
  string sig_folder = "/net/cms9/cms9r0/rohan/babies/2016_07_13/T1tbs/split/renorm/";

  Palette colors("txt/colors.txt", "default");

  // Background samples
  auto tt = Proc<Baby_full>("t#bar{t}", Process::Type::background, colors("tt_1l"),
    {mc_folder+"*TTJets*Lept*"},
    "ntruleps>=1");

  auto qcd = Proc<Baby_full>("QCD", Process::Type::background, colors("tt_1l"),
    {mc_folder+"*QCD_HT*",mc_folder+"*TTJets_TuneCUETP8M1_13TeV-madgraphMLM*"},
    "ntruleps==0");

  auto wjets = Proc<Baby_full>("W+Jets", Process::Type::background, colors("tt_1l"),
    {mc_folder+"*WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8*"});

  auto singlet = Proc<Baby_full>("Single t", Process::Type::background, colors("tt_1l"),
    {mc_folder+"*ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1*",mc_folder+"*ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1*",
	mc_folder+"*ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1*", mc_folder+"*ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1*",
	mc_folder+"*ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1*"});
  
  auto other = Proc<Baby_full>("Other", Process::Type::background, colors("tt_1l"),
    {mc_folder+"*DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8*",mc_folder+"*TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8*",
	mc_folder+"*TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8*",mc_folder+"*TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8*",
	mc_folder+"*TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8*",mc_folder+"*ttHJetTobb_M125_13TeV_amcatnloFXFX_madspin_pythia8*",
	mc_folder+"*TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8*",mc_folder+"*WJetsToQQ_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8*",
	mc_folder+"*ZJetsToQQ_HT600toInf_13TeV-madgraph*"});
  
  //Signal Samples
  auto tbs1000 = Proc<Baby_full>("m_{glu}=1000$ GeV$", Process::Type::signal, colors("tt_1l"),
    {sig_folder+"*RPV_mGluino-1000*"});
  
  auto tbs1400 = Proc<Baby_full>("m_{glu}=1400$ GeV$", Process::Type::signal, colors("tt_1l"),
    {sig_folder+"*RPV_mGluino-1400*"});
  
  vector<shared_ptr<Process> > samples = {other, singlet, wjets, qcd, tt, tbs1000, tbs1400};
  
  PlotMaker pm;
  pm.Push<Table>("rpv_regions", vector<TableRow>{
      //0-lepton
      TableRow("$N_{leps}=0$, $H_{T}>1500$, $N_{b}\\geq1$"),
	TableRow("Baseline"),
	TableRow("$M_{J}>500$, $N_{jets}\\geq4$", "nleps==0&&ht>1500&&nbm>=1&&mj>500&&njets>=4",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("Control Regions"),
	TableRow("$500<M_{J}\\leq800$, $4\\leq N_{jets}\\leq5$", "nleps==0&&ht>1500&&nbm>=1&&mj>500&&mj<=800&&njets>=4&&njets<=5",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("$500<M_{J}\\leq800$, $6\\leq N_{jets}\\leq7$", "nleps==0&&ht>1500&&nbm>=1&&mj>500&&mj<=800&&njets>=6&&njets<=7",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("$M_{J}\\geq800$, $4\\leq N_{jets}\\leq5$", "nleps==0&&ht>1500&&nbm>=1&&mj>800&&njets>=4&&njets<=5",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("$M_{J}\\geq800$, $6\\leq N_{jets}\\leq7$", "nleps==0&&ht>1500&&nbm>=1&&mj>800&&njets>=6&&njets<=7",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("Signal Regions"),
	TableRow("$500<M_{J}\\leq800$, $8\\leq N_{jets}\\leq9$", "nleps==0&&ht>1500&&nbm>=1&&mj>500&&mj<=800&&njets>=8&&njets<=9",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("$500<M_{J}\\leq800$, $N_{jets}\\geq10$", "nleps==0&&ht>1500&&nbm>=1&&mj>500&&mj<=800&&njets>=10",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("$M_{J}\\geq800$, $8\\leq N_{jets}\\leq9$", "nleps==0&&ht>1500&&nbm>=1&&mj>800&&njets>=8&&njets<=9",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("$M_{J}>800$, $N_{jets}\\geq10$", "nleps==0&&ht>1500&&nbm>=1&&mj>800&&njets>=10",0,1,"weight*w_pu_rpv/eff_trig"),
	TableRow("$\\quad N_{b}=1$", "nleps==0&&ht>1500&&nbm==1&&mj>800&&njets>=10",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("$\\quad N_{b}=2$", "nleps==0&&ht>1500&&nbm==2&&mj>800&&njets>=10",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("$\\quad N_{b}=3$", "nleps==0&&ht>1500&&nbm==3&&mj>800&&njets>=10",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("$\\quad N_{b}\\geq4$", "nleps==0&&ht>1500&&nbm>=4&&mj>800&&njets>=10",0,1,"weight*w_pu_rpv/eff_trig"),

	//1-lepton
	TableRow("$N_{leps}=1$, $H_{T}>1200$, $N_{b}\\geq1$"),
	TableRow("Baseline"),
	TableRow("$M_{J}>500$, $N_{jets}\\geq4$", "nleps==1&&ht>1200&&nbm>=1&&mj>500&&njets>=4",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("Control Regions"),
	TableRow("$500<M_{J}\\leq800$, $4\\leq N_{jets}\\leq5$", "nleps==1&&ht>1200&&nbm>=1&&mj>500&&mj<=800&&njets>=4&&njets<=5",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("$M_{J}\\geq800$, $4\\leq N_{jets}\\leq5$", "nleps==1&&ht>1200&&nbm>=1&&mj>800&&njets>=4&&njets<=5",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("Signal Regions"),
	TableRow("$500<M_{J}\\leq800$, $6\\leq N_{jets}\\leq7$", "nleps==1&&ht>1200&&nbm>=1&&mj>500&&mj<=800&&njets>=6&&njets<=7",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("$500<M_{J}\\leq800$, $N_{jets}\\geq8$", "nleps==1&&ht>1200&&nbm>=1&&mj>500&&mj<=800&&njets>=8",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("$M_{J}\\geq800$, $6\\leq N_{jets}\\leq7$", "nleps==1&&ht>1200&&nbm>=1&&mj>800&&njets>=6&&njets<=7",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("$M_{J}>800$, $N_{jets}\\geq8$", "nleps==1&&ht>1200&&nbm>=1&&mj>800&&njets>=8",0,1,"weight*w_pu_rpv/eff_trig"),
	TableRow("$\\quad N_{b}=1$", "nleps==1&&ht>1200&&nbm==1&&mj>800&&njets>=8",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("$\\quad N_{b}=2$", "nleps==1&&ht>1200&&nbm==2&&mj>800&&njets>=8",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("$\\quad N_{b}=3$", "nleps==1&&ht>1200&&nbm==3&&mj>800&&njets>=8",0,0,"weight*w_pu_rpv/eff_trig"),
	TableRow("$\\quad N_{b}\\geq4$", "nleps==1&&ht>1200&&nbm>=4&&mj>800&&njets>=8",0,0,"weight*w_pu_rpv/eff_trig"),
	},samples,false);
  
  pm.MakePlots(lumi);
}


