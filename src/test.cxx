#include <cmath>

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <functional>
#include <stdexcept>
#include <memory>

#include "TError.h"

#include "baby.hpp"
#include "process.hpp"
#include "named_func.hpp"
#include "plot_maker.hpp"
#include "styles.hpp"
#include "utilities.hpp"

#include "thread_pool.hpp"

using namespace std;

template<typename T>
unique_ptr<Baby> BabyPtr(const set<string> &names){
  return unique_ptr<Baby>(new T(names));
}

int dumb = 0;
mutex m;

int Dumb(){
  this_thread::sleep_for(chrono::seconds(2));
  lock_guard<mutex> lock(m);
  return ++dumb;
}

int main(){
  if(false){
    ThreadPool tp;
    vector<future<int> > v(10);
    for(size_t i = 0; i < v.size(); ++i){
      v.at(i) = tp.Push(Dumb);
    }
    for(size_t i = 0; i < v.size(); ++i){
      cout << "Waiting for " << i << endl;
      double x = v.at(i).get();
      cout << "Found x = " << x << " for i = " << i << endl;
    }
  }
  //return 0;

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
             HistoDef(20, 0., 1000., FUNC(b.ht()), "H_{T} [GeV]"));
  pm.AddPlot({qcd, ttjets, ttjets2, sms_nc, sms_c},
             HistoDef(20, 0., 1000., FUNC(b.met()+b.jets_csv()->size()+b.jets_eta()->size()+b.jets_m()->size()+b.jets_phi()->size()+b.jets_pt()->size()+b.jets_islep()->size()), "MET [GeV]", FUNC(b.ht()>500)));
  pm.MakePlots();
}
