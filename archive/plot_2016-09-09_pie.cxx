///// table_preds: Makes piecharts

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

#include "TError.h" // Controls error level reporting

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/plot_opt.hpp"

namespace{
  //bool do_met150 = true;
}

using namespace std;

int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_stdnj5/");
  //if(do_met150) foldermc = (bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_met150/");
  Palette colors("txt/colors.txt", "default");

  // Cuts in baseline speed up the yield finding
  string base1l = "mj14>250 && nleps==1 && nveto==0 && st>500 && met>200 && pass && weight<1 && njets>=6 && nbm>=1"; // Excluding one QCD event
  string base2l = "mj14>250 && ((nleps==1 && nveto==1 && njets>=6 && nbm>=1 && mt>140) || (nleps==2  && njets>=5 && nbm<=2)) && st>500 && met>200 && met<500 && pass && weight<1"; // Excluding one QCD event
  
  string ntupletag = "*metG200*.root";  
  auto proc_tt1l = Process::MakeShared<Baby_full>("t#bar{t} (l)", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets*SingleLept"+ntupletag, foldermc+"*_TTJets_HT"+ntupletag},
    "stitch && ntruleps==1");
  auto proc_tt2l = Process::MakeShared<Baby_full>("t#bar{t} (ll)", Process::Type::background, colors("tt_2l"),
    {foldermc+"*_TTJets*DiLept"+ntupletag, foldermc+"*_TTJets_HT"+ntupletag},
    "stitch && ntruleps==2 ");
  auto proc_wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {foldermc+"*_WJetsToLNu"+ntupletag}, "stitch");
  auto proc_single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {foldermc+"*_ST_"+ntupletag}, "1");
  auto proc_ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {foldermc+"*_TTWJets"+ntupletag, foldermc+"*_TTZ"+ntupletag}, "1");
  auto proc_other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {foldermc+"*DYJetsToLL"+ntupletag,foldermc+"*QCD_HT*0_Tune"+ntupletag,foldermc+"*QCD_HT*Inf_Tune"+ntupletag,
        foldermc+"*_ZJet"+ntupletag,foldermc+"*_ttHJetTobb"+ntupletag,
        foldermc+"*_TTGJets"+ntupletag,foldermc+"*_TTTT"+ntupletag,
        foldermc+"*_WH_HToBB"+ntupletag,foldermc+"*_ZH_HToBB"+ntupletag,
        foldermc+"*_WWTo"+ntupletag,foldermc+"*_WZ"+ntupletag,foldermc+"*_ZZ_"+ntupletag},
    "stitch");

  vector<shared_ptr<Process> > all_procs = {proc_tt1l, proc_tt2l, proc_wjets, proc_single_t, proc_ttv, proc_other};

  vector<TString> cuts;
  cuts.push_back(base1l + "&& mt<=140 && mj14<=400");
  cuts.push_back(base1l + "&& mt<=140 && mj14>400");
  cuts.push_back(base1l + "&& mt>140  && mj14<=400");
  cuts.push_back(base1l + "&& mt>140  && mj14>400");
  cuts.push_back(base2l + "&& mj14<=400");
  cuts.push_back(base2l + "&& mj14>400");

  vector<TableRow> table_cuts;
  for(size_t icut=0; icut<cuts.size(); icut++){
    table_cuts.push_back(TableRow(cuts2tex(cuts[icut]).Data(), cuts[icut].Data()));
  }

  PlotMaker pm;
  pm.Push<Table>("chart",  table_cuts, all_procs, true, true, true);
  pm.min_print_ = true;
  pm.MakePlots(40.);

  time(&endtime);
  cout<<endl<<"Making "<<table_cuts.size()<<" piecharts took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
