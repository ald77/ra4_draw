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
#include "TVector2.h"

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/event_scan.hpp"
#include "core/utilities.hpp"
#include "hig/hig_functions.hpp"

using namespace std;
using namespace PlotOptTypes;
using namespace Higfuncs;

vector<unsigned> higidx(const Baby &b);

namespace{
  bool single_thread = true;
  double lumi = 35.9;
}

int main(int argc, char *argv[]){
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  const NamedFunc isgood("isgood",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned i(0); i<b.mc_mass()->size(); i++){
      if (b.mc_id()->at(i)==1000023 || b.mc_id()->at(i)==1000025) {
        int mom = abs(b.mc_mom()->at(i));
        if (mom==35 || mom==36 || mom==37) 
          return false;
      }
    }
    return true;
  });

  const NamedFunc isgood_nob("isgood_nob",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned i(0); i<b.mc_mass()->size(); i++){
      if (b.mc_id()->at(i)==1000023 || b.mc_id()->at(i)==1000025) {
        int mom = abs(b.mc_mom()->at(i));
        if (mom==5) 
          return false;
      }
    }
    return true;
  });

  string foldersig = "/cms2r0/babymaker/babies/2017_02_26/TChiHH/merged_higmc_unskimmed/";
  string folderbad = "/cms2r0/babymaker/babies/2017_02_26/TChiHH/bad/";

  Palette colors("txt/colors.txt", "default");
  vector<vector<shared_ptr<Process> > > procs;
  procs.push_back(vector<shared_ptr<Process> >());
  procs.back().push_back(Process::MakeShared<Baby_full>("m=200: filter 35,36,5", Process::Type::signal, 
    kMagenta, {foldersig+"*TChiHH_mGluino-200*.root"}, isgood_nob==true));
  procs.back().back()->SetLineStyle(2);
  procs.back().push_back(Process::MakeShared<Baby_full>("m=200: only mom=5", Process::Type::signal, 
    kBlue, {foldersig+"*TChiHH_mGluino-200*.root"}, isgood_nob==false));
  procs.back().back()->SetLineStyle(2);

  procs.back().push_back(Process::MakeShared<Baby_full>("m=200: New GEN", Process::Type::signal, 
    kGray+2, {folderbad+"genbaby_tchi200new.root"}, "1"));
  procs.back().back()->SetLineStyle(2);
  procs.back().push_back(Process::MakeShared<Baby_full>("m=200: filter 35,36", Process::Type::signal, 
    kRed+1, {foldersig+"*TChiHH_mGluino-200*.root"}, "1"));

  procs.push_back(vector<shared_ptr<Process> >());
  procs.back().push_back(Process::MakeShared<Baby_full>("m=225: filter 5", Process::Type::signal, 
    kCyan, {foldersig+"*TChiHH_mGluino-225*.root"}, isgood_nob==true));
  procs.back().back()->SetLineStyle(2);
  procs.back().push_back(Process::MakeShared<Baby_full>("m=225: filter 5", Process::Type::signal, 
    kBlue, {foldersig+"*TChiHH_mGluino-225*.root"}, isgood_nob==false));
  procs.back().back()->SetLineStyle(2);
  procs.back().push_back(Process::MakeShared<Baby_full>("m=225: New GEN", Process::Type::signal, 
    kGray+2, {folderbad+"genbaby_tchi225new.root"}, "1"));
  procs.back().back()->SetLineStyle(2);
  procs.back().push_back(Process::MakeShared<Baby_full>("m=225: No filter", Process::Type::signal, 
    kGreen+1, {foldersig+"*TChiHH_mGluino-225*.root"}, "1"));



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
  vector<PlotOpt> linplt = {lin_shapes_info};
  vector<PlotOpt> logplt = {log_shapes_info};

  PlotMaker pm;

  NamedFunc lead_higs_pt("lead_higs_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float higpt(0);
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (b.mc_id()->at(i)!=25) continue;
      if (b.mc_pt()->at(i)>higpt) higpt = b.mc_pt()->at(i);
    }
    return higpt;
  });

  NamedFunc nbacc("nbacc",[](const Baby &b) -> NamedFunc::ScalarType{
    int tmp_nbacc(0);
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (abs(b.mc_id()->at(i))!=5 || b.mc_mom()->at(i)!=25) continue;
      if (b.mc_pt()->at(i)>30 && abs(b.mc_eta()->at(i))<2.5) tmp_nbacc++;
    }
    return tmp_nbacc;
  });

  NamedFunc bpt("bpt",[](const Baby &b) -> NamedFunc::VectorType{
    vector<double> tmp_bpt(0);
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (abs(b.mc_id()->at(i))!=5 || b.mc_mom()->at(i)!=25) continue;
      if (b.mc_pt()->at(i)>30 && abs(b.mc_eta()->at(i))<2.5) {
        tmp_bpt.push_back(b.mc_pt()->at(i));
      }
    }
    return tmp_bpt;
  });

  NamedFunc lead_higsino_pt("lead_higsino_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    float higpt(0);
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (b.mc_id()->at(i)!=1000025 && b.mc_id()->at(i)!=1000023) continue;
      if (b.mc_pt()->at(i)>higpt) higpt = b.mc_pt()->at(i);
    }
    return higpt;
  });

  NamedFunc n3mom("n3mom",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (b.mc_id()->at(i)!=1000025) continue;
      return abs(b.mc_mom()->at(i));
    }
    return 0;
  });
 
  NamedFunc rel_lsp_dpt("rel_lsp_dpt",[](const Baby &b) -> NamedFunc::ScalarType{
    vector<unsigned> lspidx;
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (b.mc_id()->at(i)==1000022) lspidx.push_back(i);
      if (lspidx.size()>1) break;
    }
    float max_pt = max(b.mc_pt()->at(lspidx[0]),b.mc_pt()->at(lspidx[1]));
    return fabs(b.mc_pt()->at(lspidx[0])-b.mc_pt()->at(lspidx[1]))/max_pt;
  });

  NamedFunc hig_dphi("hig_dphi",[](const Baby &b) -> NamedFunc::ScalarType{
    vector<unsigned> idx = higidx(b);
    return deltaPhi(b.mc_phi()->at(idx[0]),b.mc_phi()->at(idx[1]));
  });  

  NamedFunc hig_dr("hig_dr",[](const Baby &b) -> NamedFunc::ScalarType{
    vector<unsigned> idx = higidx(b);
    return deltaR(b.mc_eta()->at(idx[0]), b.mc_phi()->at(idx[0]), 
                  b.mc_eta()->at(idx[1]), b.mc_phi()->at(idx[1]));
  });  

  NamedFunc trumet("trumet",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.mlsp()==1) return b.met_tru();
    else return b.gen_met();
  });  

  NamedFunc truht("truht",[](const Baby &b) -> NamedFunc::ScalarType{
    if (b.mlsp()==1) return b.ht_tru();
    else return b.genht();
  });  

  NamedFunc kill_bprod("kill_bprod",[](const Baby &b) -> NamedFunc::ScalarType{
    for (unsigned i(0); i<b.mc_pt()->size(); i++){
      if (b.mc_id()->at(i)==1000023 || b.mc_id()->at(i)==1000025)
        if (abs(b.mc_mom()->at(i))==5)
          return 0;
    }
    return 1;
  });  

  string chi20 = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.]{#scale[0.85]{_{2}}}";
  string chi30 = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.]{#scale[0.85]{_{3}}}";

  string baseline("met>150");

  string c2b = "nbdt==2&&nbdm==2";
  string c3b = "nbdt>=2&&nbdm==3&&nbdl==3";
  string c4b = "nbdt>=2&&nbdm>=3&&nbdl>=4";

  NamedFunc wgt = "1"; //weight_higd;// * eff_higtrig;

  vector<string> lbl = {"m200","m225"};
  for (unsigned ipr(0); ipr<procs.size(); ipr++) {
    pm.Push<Hist1D>(Axis(25,0,25, n3mom, "Higgsino Mother PDG Id"),
      "1", procs[ipr], linplt).Tag(lbl[ipr]).Weight(wgt);
    pm.Push<Hist1D>(Axis(20,0,600,lead_higsino_pt, "Leading Higgsino p_{T} [GeV]"),
      "1", procs[ipr], linplt).Tag(lbl[ipr]).Weight(wgt);

    pm.Push<Hist1D>(Axis(16,0,400,trumet, "Gen E_{T}^{miss} [GeV]", {150, 200, 300, 450}),
      "1", procs[ipr], linplt).Tag(lbl[ipr]).Weight(wgt);

    pm.Push<Hist1D>(Axis(16,0,400,bpt, "p_{T} of b-quarks [GeV]"),
      "1", procs[ipr], linplt).Tag(lbl[ipr]).Weight(wgt);

    pm.Push<Hist1D>(Axis(6,-0.5,5.5,nbacc, "Number of b-quarks within acceptance"),
      "1", procs[ipr], linplt).Tag(lbl[ipr]).Weight(wgt);
  }

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

vector<unsigned> higidx(const Baby &b){
  vector<unsigned> idx;
  for (unsigned i(0); i<b.mc_pt()->size(); i++){
    if (b.mc_id()->at(i)==25) idx.push_back(i);
    if (idx.size()>1) break;
  }
  return idx;
}

// vector<unsigned> bidx(const Baby &b){
//   vector<unsigned> idx;
//   for (unsigned i(0); i<b.mc_pt()->size(); i++){
//     if (b.mc_id()->at(i)==25) higidx.push_back(i);
//     if (higidx.size()>1) break;
//   }
//   return idx;
// }
