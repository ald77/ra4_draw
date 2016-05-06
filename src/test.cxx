#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "TError.h"

#include "baby.hpp"
#include "process.hpp"
#include "named_func.hpp"
#include "plot_maker.hpp"
#include "plot_opt.hpp"

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

  auto bkg1 = Proc<Baby_basic>("Background 1", Process::Type::background, kBlue,
    {"~/ntuples/2015_09_28_ana/skim/*_QCD*.root"}, "nbm>2");
  auto bkg2 = Proc<Baby_basic>("Background 2", Process::Type::background, kGreen,
    {"~/ntuples/2015_09_28_ana/skim/*_QCD*.root"}, "nbm>2");
  auto sig = Proc<Baby_basic>("Signal", Process::Type::signal, kRed,
    {"~/ntuples/2015_09_28_ana/skim/*_QCD*.root"}, "nbm>2");
  auto data = Proc<Baby_basic>("Data", Process::Type::data, kBlack,
    {"~/ntuples/2015_09_28_ana/skim/*_QCD*.root"}, "nbm>2");

  PlotOpt opt("txt/plot_styles.txt", "CMSPaper");
  
  PlotMaker pm;
  for(int i = 0; i < 1; ++i){
    pm.AddPlot({bkg1, bkg2, sig, data},
               HistoDef(20, 0., 2000., "ht", "H_{T}", "GeV"), opt().Bottom(BottomType::ratio));
    pm.AddPlot({bkg1, bkg2, sig, data},
               HistoDef(20, 0., 2000., "ht+met", "H_{T}+MET", "GeV"), opt().Bottom(BottomType::diff));
    pm.AddPlot({bkg1, bkg2, sig, data},
               HistoDef(11, -0.5, 10.5, "njets", "N_{jets}"), opt);
    pm.AddPlot({bkg1, bkg2, sig, data},
               HistoDef(10, 0., 1000., "jets_pt", "Jet p_{T}", "GeV"), opt);
    pm.AddPlot({bkg1, bkg2, sig, data},
               HistoDef(10, 0., 1000., "jets_pt[0]", "Lead Jet p_{T}", "GeV"), opt);
  }
  pm.MakePlots();
}
