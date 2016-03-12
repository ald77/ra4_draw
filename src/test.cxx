#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "TError.h"

#include "baby.hpp"
#include "process.hpp"
#include "named_func.hpp"
#include "plot_maker.hpp"
#include "styles.hpp"

using namespace std;

template<typename T>
unique_ptr<Baby> BabyPtr(const set<string> &names){
  return unique_ptr<Baby>(new T(names));
}

int main(){
  gErrorIgnoreLevel = 6000;
  styles style("RA4");
  style.setDefaultStyle();
  
  auto qcd = make_shared<Process>("QCD", Process::Type::background, kGreen+2,
                                  BabyPtr<Baby_basic>({"~/ntuples/2015_09_28_ana/skim/*_QCD*.root"}));
  auto ttjets = make_shared<Process>("t#bar{t}", Process::Type::background, kBlue+2,
                                     BabyPtr<Baby_basic>({"~/ntuples/2015_09_28_ana/skim/*_TTJets*.root"}));
  auto ttjets2 = make_shared<Process>("t#bar{t}2", Process::Type::background, kBlue+1,
                                      BabyPtr<Baby_basic>({"~/ntuples/2015_09_28_ana/skim2/*_TTJets*.root"}));
  auto sms_nc = make_shared<Process>("T1tttt(1500,100)", Process::Type::signal, kRed,
                                     BabyPtr<Baby_basic>({"~/ntuples/2015_09_28_ana/skim/*_SMS-1500*100*.root"}));
  auto sms_c = make_shared<Process>("T1tttt(1200,800)", Process::Type::signal, kRed+2,
                                    BabyPtr<Baby_basic>({"~/ntuples/2015_09_28_ana/skim/*_SMS-1200*800*.root"}));

  PlotMaker pm;
  pm.AddPlot({qcd, ttjets, ttjets2, sms_nc, sms_c},
             HistoDef(20, 0., 1000., "ht", "H_{T} [GeV]"));
  pm.MakePlots();
}
