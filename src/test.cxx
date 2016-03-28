#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "TError.h"

#include "baby.hpp"
#include "process.hpp"
#include "named_func.hpp"
#include "plot_maker.hpp"

using namespace std;

template<typename T>
shared_ptr<Process> Proc(const string process_name, Process::Type type,
                         int color, const set<string> &files){
  return make_shared<Process>(process_name, type, color,
                              unique_ptr<Baby>(new T(files)));
}

int main(){
  gErrorIgnoreLevel = 6000;

  auto qcd = Proc<Baby_basic>("QCD", Process::Type::background, kGreen+2,
    {"~/ntuples/2015_09_28_ana/skim/*_QCD*.root"});
  auto ttjets = Proc<Baby_basic>("t#bar{t}", Process::Type::background, kBlue+2,
    {"~/ntuples/2015_09_28_ana/skim/*_TTJets*.root"});
  auto ttjets2 = Proc<Baby_basic>("t#bar{t}2", Process::Type::background, kBlue+1,
    {"~/ntuples/2015_09_28_ana/skim2/*_TTJets*.root"});
  auto sms_nc = Proc<Baby_basic>("T1tttt(1500,100)", Process::Type::signal, kRed,
    {"~/ntuples/2015_09_28_ana/skim/*_SMS-1500*100*.root"});
  auto sms_c = Proc<Baby_basic>("T1tttt(1200,800)", Process::Type::signal, kRed+2,
    {"~/ntuples/2015_09_28_ana/skim/*_SMS-1200*800*.root"});

  PlotMaker pm;
  for(int i = 0; i < 1; ++i){
    pm.AddPlot({qcd, ttjets, ttjets2, sms_nc, sms_c},
               HistoDef(20, 0., 1000., "ht", "H_{T} [GeV]"));
    pm.AddPlot({qcd, ttjets, ttjets2, sms_nc, sms_c},
               HistoDef(20, 0., 1000., "ht+met", "H_{T} [GeV]"));
    pm.AddPlot({qcd, ttjets, ttjets2, sms_nc, sms_c},
               HistoDef(20, 0., 1000., "njets", "H_{T} [GeV]"));
    pm.AddPlot({qcd, ttjets, ttjets2, sms_nc, sms_c},
               HistoDef(20, 0., 1000., "jets_pt", "H_{T} [GeV]"));
    pm.AddPlot({qcd, ttjets, ttjets2, sms_nc, sms_c},
               HistoDef(20, 0., 1000., "jets_pt[0]", "H_{T} [GeV]"));
  }
  pm.MakePlots();
}
