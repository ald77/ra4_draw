// Cheatsheet:
//    To obtain plots with no mismeasurement and kappa's with expected data stats:
//         ./run/hig/plot_systematics.exe --mm mc_as_data -l 36.2
//    There are 4 possibilities for the skim, requested with option -s. These are: search, zll, qcd, ttbar
//    Option -t plots the kappas with a tighter selection, see basecuts below, e.g.
//         ./run/hig/plot_systematics.exe --mm mc_as_data -t -s zll -l 36.2

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>  // setw
#include <chrono>
#include <string>

#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "RooStats/NumberCountingUtils.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/abcd_method.hpp"
#include "core/styles.hpp"
#include "core/plot_opt.hpp"
#include "core/config_parser.hpp"
#include "core/functions.hpp"
#include "hig/hig_functions.hpp"

using namespace std;

namespace{
  bool actualZbi = false;
  bool only_mc = false;
  bool only_kappa = false;
  bool split_bkg = false;
  bool only_dilepton = false;
  bool do_leptons = false;
  bool do_signal = true;
  bool unblind = false;
  bool debug = false;
  bool do_ht = false;
  bool do_correction = false;
  bool do_loose = true;
  bool do_trim = true;
  bool do_highnb = false;
  bool do_midnb = false;
  bool do_onemet = false;
  TString skim = "search";
  TString json = "2p6";
  TString only_method = "";
  TString mc_lumi = "";
  string sys_wgts_file = "txt/sys_weights.cfg";
  string mm_scen = "";
  float lumi=36.2;
  bool quick_test = false;
  // for office use only
  vector<TString> syst_names;
  vector<vector<float>> syst_values;
}

string GetScenario(const string &method);

TString printTable(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
                   vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds, 
                   vector<vector<float> > yieldsPlane, vector<shared_ptr<Process> > &proc_sigs);
void plotKappa(abcd_method &abcd, vector<vector<vector<float> > >  &kappas, 
               vector<vector<vector<float> > >  &kappas_mm, vector<vector<vector<float> > >  &kmcdat);
vector<vector<float> > findPreds(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
                                 vector<vector<vector<float> > > &kappas, 
                                 vector<vector<vector<float> > > &kappas_mm, 
                                 vector<vector<vector<float> > > &kmcdat, 
                                 vector<vector<vector<float> > > &datapreds,
                                 vector<vector<vector<float> > > &preds);
void printDebug(abcd_method &abcd, vector<vector<GammaParams> > &allyields, TString baseline,
                vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &kappas_mm, 
                vector<vector<vector<float> > > &preds);
TString Zbi(double Nobs, double Nbkg, double Eup_bkg, double Edown_bkg);

void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  chrono::high_resolution_clock::time_point begTime;
  begTime = chrono::high_resolution_clock::now();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  if(mm_scen == ""){
    cout << " ======== Doing all mis-measurement scenarios ======== \n" << endl;
    only_mc = true;
  }else if(mm_scen == "data"){
    cout << " ======== Comparing MC and actual DATA ======== \n" << endl;
    only_mc = false;
  }else if(mm_scen == "no_mismeasurement" || mm_scen == "off" || mm_scen == "mc_as_data"){
    cout << " ======== No mismeasurement applied ======== \n" << endl;
    // if(mm_scen == "mc_as_data") only_mc = true;
  }else{
    cout << " ======== Doing mis-measurement scenario " << mm_scen << " ======== \n" << endl;
    only_mc = true;
  }


  vector<string> scenarios = ConfigParser::GetOptSets(sys_wgts_file);
  //NamedFunc w = "weight*eff_trig";
  NamedFunc w = "weight" * Higfuncs::eff_higtrig;
  map<string, NamedFunc> weights, corrections;
  auto central = Functions::Variation::central;
  weights.emplace("no_mismeasurement", w);
  corrections.emplace("no_mismeasurement", 1.);
  if(mm_scen == ""){
    for(const auto &scen: scenarios){
      weights.emplace(scen, w*Functions::MismeasurementWeight(sys_wgts_file, scen));
      corrections.emplace(scen, Functions::MismeasurementCorrection(sys_wgts_file, scen, central));
    }
  }else if(mm_scen == "totunc"){
    scenarios = vector<string>();
    do_correction = true;
    // scenarios.push_back("syst_ttx_up");
    // weights.emplace("syst_ttx_up", w*(Higfuncs::wgt_comp)*(1+Higfuncs::wgt_syst_ttx));
    // corrections.emplace("syst_ttx_up", Higfuncs::wgt_comp);
    // scenarios.push_back("syst_ttx_dn");
    // weights.emplace("syst_ttx_dn", w*(Higfuncs::wgt_comp)*(1-Higfuncs::wgt_syst_ttx));
    // corrections.emplace("syst_ttx_dn", Higfuncs::wgt_comp);

    // scenarios.push_back("syst_vjets_up");
    // weights.emplace("syst_vjets_up", w*(Higfuncs::wgt_comp)*(1+Higfuncs::wgt_syst_vjets));
    // corrections.emplace("syst_vjets_up", Higfuncs::wgt_comp);
    // scenarios.push_back("syst_vjets_dn");
    // weights.emplace("syst_vjets_dn", w*(Higfuncs::wgt_comp)*(1-Higfuncs::wgt_syst_vjets));
    // corrections.emplace("syst_vjets_dn", Higfuncs::wgt_comp);

    // scenarios.push_back("syst_qcd_up");
    // weights.emplace("syst_qcd_up", w*(Higfuncs::wgt_comp)*(1+Higfuncs::wgt_syst_qcd));
    // corrections.emplace("syst_qcd_up", Higfuncs::wgt_comp);
    // scenarios.push_back("syst_qcd_dn");
    // weights.emplace("syst_qcd_dn", w*(Higfuncs::wgt_comp)*(1-Higfuncs::wgt_syst_qcd));
    // corrections.emplace("syst_qcd_dn", Higfuncs::wgt_comp);

    // scenarios.push_back("syst_comp");
    // weights.emplace("syst_comp", w*(Higfuncs::wgt_comp)); 
    // corrections.emplace("syst_comp", 1.);

    scenarios.push_back("syst_mcstat");
    weights.emplace("syst_mcstat", w);
    corrections.emplace("syst_mcstat", 1.);

    // scenarios.push_back("syst_bctag");
    // weights.emplace("syst_bctag", w*"sys_bctag[0]");
    // corrections.emplace("syst_bctag", 1.);

    // scenarios.push_back("syst_udsgtag");
    // weights.emplace("syst_udsgtag", w*"sys_udsgtag[0]");
    // corrections.emplace("syst_udsgtag", 1.);
  }else if(mm_scen == "data"){
    scenarios = vector<string>{mm_scen};
  }else if(mm_scen != "no_mismeasurement"){
    scenarios = vector<string>{mm_scen};
    weights.emplace(mm_scen, w*Functions::MismeasurementWeight(sys_wgts_file, mm_scen));
    corrections.emplace(mm_scen, Functions::MismeasurementCorrection(sys_wgts_file, mm_scen, central));
  }

  //// Capybara
  string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higloose/");
  if (skim=="ttbar") foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higlep1met0/";
  if (skim=="zll") foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higlep2/";
  if (skim=="qcd") foldermc = bfolder+"/cms2r0/babymaker/babies/2016_08_10/mc/merged_higmc_higqcd/";
  string folderdata(bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_higloose/");
  if (skim=="ttbar") folderdata = bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_higlep1met0/";
  if (skim=="zll") folderdata = bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_higlep2/";
  if (skim=="qcd") folderdata = bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/merged_higdata_higqcd/";
  string foldersig(bfolder+"/cms2r0/babymaker/babies/2016_08_10/TChiHH/merged_higmc_higloose/");

  Palette colors("txt/colors.txt", "default");

  // Cuts in baseline speed up the yield finding
  string baseline_s = "met/met_calo<5 && njets>=4 && njets<=5";
  NamedFunc baseline=baseline_s;

  map<string, set<string>> mctags; 
  mctags["ttx"]     = set<string>({"*_TTJets*Lept*.root", "*_TTJets_HT*.root", "*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*_ttHJetTobb*.root","*_TTTT*.root"});
  mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
  mctags["other"]   = set<string>({"*_ST_*.root",
                                   "*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root",
                                   "*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                   "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  
  vector<string> sigMasses({"225", "300", "400", "700"});
  vector<shared_ptr<Process> > proc_sigs;
  for(size_t ind=0; ind<sigMasses.size(); ind++)
    proc_sigs.push_back(Process::MakeShared<Baby_full>("TChiHH("+sigMasses[ind]+")", Process::Type::signal, 2,
      {foldersig+"*mGluino-"+sigMasses[ind]+"_*.root"}, baseline && "stitch && pass_ra2_badmu"));

  auto proc_tt1l = Process::MakeShared<Baby_full>("tt+X", Process::Type::background, colors("tt_1l"),
    attach_folder(foldermc, mctags["ttx"]), baseline && "stitch && pass && pass_ra2_badmu");
  auto proc_tt2l = Process::MakeShared<Baby_full>("V+jets", Process::Type::background, kOrange+1,
    attach_folder(foldermc, mctags["vjets"]), baseline && "stitch && pass && pass_ra2_badmu");
  auto proc_other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    attach_folder(foldermc, mctags["other"]), baseline && "stitch && pass && pass_ra2_badmu");

  //// All MC files to make pseudodata
  set<string> names_allmc;
  for (auto &iset: mctags) names_allmc.insert(iset.second.begin(), iset.second.end());

  // Setting luminosity
  string jsonCuts = "nonblind";
  if(skim.Contains("2015")) lumi = 2.3;
  else if(json=="0p869"){
    lumi = 0.869;
    jsonCuts = "nonblind";
  } else if(json=="12p9"){
    lumi = 12.9;
    jsonCuts = "json12p9";
  } else if(json=="full"){
    lumi = 36.2;
    jsonCuts = "1";
  } 
  if(mc_lumi!="") lumi = mc_lumi.Atof();

  set<string> names_data({folderdata+"*.root"});
  if(only_mc){
    names_data = attach_folder(foldermc, names_allmc);
    if(quick_test) names_data = set<string>({foldermc+"*_TTJets_SingleLeptFromT_*.root"});
  }

  NamedFunc base_data = baseline && Higfuncs::trig_hig>0. && jsonCuts+"&& pass && pass_ra2_badmu";
  if (only_mc) base_data = baseline && "stitch && pass && pass_ra2_badmu";
   if(mm_scen == "data")
     cout<<"Data files are "<<*(names_data.begin())<<" with cuts "<<baseline<<"&&"<< Higfuncs::trig_hig << "&&pass&&pass_ra2_badmu"<<endl<<endl;
  auto proc_data = Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    names_data,base_data);

  //// Use this process to make quick plots. Requires being run without split_bkg
  auto proc_bkg = Process::MakeShared<Baby_full>("All_bkg", Process::Type::background, colors("tt_1l"),
    {foldermc+"*_TTJets_SingleLeptFromT_*.root"}, baseline && " pass && pass_ra2_badmu");

  vector<shared_ptr<Process> > all_procs;
  if(!quick_test) all_procs = vector<shared_ptr<Process> >{proc_tt1l, proc_tt2l, proc_other};
  else {
    all_procs = vector<shared_ptr<Process> >{proc_bkg};
    split_bkg = false;
  }
  if (do_signal){
    for(size_t ind=0; ind<proc_sigs.size(); ind++)
      all_procs.push_back(proc_sigs[ind]);
  }
  all_procs.push_back(proc_data);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining basic cuts //////////////////////////////////////////
  // baseline defined above
  vector<TString> abcdcuts, metcuts, nbcuts;

  ////// MET cuts
  string metdef = "met";
  if (skim=="zll") metdef = "(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))";
  if (skim=="ttbar" || skim=="zll"){
    metcuts.push_back(metdef+">0&&"+metdef+"<=75");
    metcuts.push_back(metdef+">75&&"+metdef+"<=150");
  }
  metcuts.push_back(metdef+">150&&"+metdef+"<=200");
  metcuts.push_back(metdef+">200&&"+metdef+"<=300");
  metcuts.push_back(metdef+">300&&"+metdef+"<=450");
  metcuts.push_back(metdef+">450");
  if (skim=="qcd" || skim=="search") {
    if (do_onemet) metcuts.push_back(metdef+">150");
  } else if (skim=="ttbar" || skim=="zll"){
    if (do_onemet) metcuts.push_back(metdef+">0");
  }
  

  ////// Nb cuts
  if (skim=="ttbar" || skim=="search" || do_highnb){
    nbcuts.push_back("nbt==2&&nbm==2");
    nbcuts.push_back("nbt>=2&&nbm==3&&nbl==3");
    nbcuts.push_back("nbt>=2&&nbm>=3&&nbl>=4");
  } else if (do_midnb){
    nbcuts.push_back("nbm==1");
    nbcuts.push_back("nbt==2&&nbm==2");
  } else {
    nbcuts.push_back("nbm==0");
    nbcuts.push_back("nbm==1");
  }

  ////// CR, SR cuts
  TString c_sr="hig_am>100&&hig_am<140&&hig_dm<40";
  TString c_cr="!("+c_sr+") && hig_am<200 && hig_dm<40";//"(hig_am<=100 || hig_am>=140) && hig_am<200 && hig_dm<40";
  if (!do_trim) c_cr="hig_am<=100 || hig_am>=140 || hig_dm>=40";

  ////// One loose and one tight selection option for each region
  TString basecuts("njets>=4 && njets<=5"); 
  // zll skim:  ((elel_m>80&&elel_m<100)||(mumu_m>80&&mumu_m<100))
  // nleps==2 && Max$(leps_pt)>40 
  if (skim=="zll") {
    if (do_loose) basecuts += " && nleps==2 && ((elel_m>80&&elel_m<100)||(mumu_m>80&&mumu_m<100)) && met<50";
    else basecuts += " && nleps==2 && ((elel_m>80&&elel_m<100)||(mumu_m>80&&mumu_m<100)) && met<50 && hig_drmax<2.2";
  }
  // qcd skim - met>150 && nvleps==0 && (njets==4||njets==5)
  if (skim=="qcd") {
    if (do_loose) basecuts += " && nvleps==0 && ntks==0 && low_dphi";
    else basecuts += " && nvleps==0 && ntks==0 && low_dphi && hig_drmax<2.2";
  }
  // ttbar skim - met>100 && nleps==1 && (njets==4||njets==5) && nbm>=2
  if (skim=="ttbar") {
    if (do_loose) basecuts += " && nleps==1 && mt<100";
    else basecuts += " && nleps==1 && mt<100 && hig_drmax<2.2";
  } 
  // search skim - met>100 && nvleps==0 && (njets==4||njets==5) && nbm>=2
  if (skim=="search") {
    if (do_loose) basecuts += " && nvleps==0 && ntks==0 && !low_dphi";
    else basecuts += " && nvleps==0 && ntks==0 && !low_dphi && hig_drmax<2.2";
  }

  ////// ABCD cuts
  vector<TString> abcdcuts_std = {c_cr + " && 2bcuts",
                                  c_cr + " && nj_1l",
                                  c_sr + " && 2bcuts",
                                  c_sr + " && nj_1l"};

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining ABCD methods //////////////////////////////////////////
  PlotMaker pm;

  ///// Running over these methods
  vector<TString> methods = {"TTML"};

  if(only_method!="") methods = vector<TString>({only_method});

  vector<TString> methods_tmp = methods;
  methods.clear();
  for(const auto &scenario: scenarios){
    for(const auto &method: methods_tmp){
      methods.push_back(method+"_scen_"+scenario.c_str());
    }
  }

  TString njets = "N#lower[-0.1]{_{jets}}";
  TString nleps = "N#lower[-0.1]{_{leps}}";
  TString nbs = "N#lower[-0.1]{_{b}}";
  TString llpt = "(mumu_pt*(mumu_pt>0)+elel_pt*(elel_pt>0))";
  vector<abcd_method> abcds;
  TString caption = "", abcd_title;
  for(size_t iabcd=0; iabcd<methods.size(); iabcd++) {
    TString method = methods[iabcd];
    mm_scen = GetScenario(method.Data());

    if(skim.Contains("zll")){
      caption = "$N_{\\rm leps}=2$ CR";
      abcd_title = "Control region: "+nleps+" = 2";
   } else if(skim.Contains("qcd")){
      caption = "Low $\\Delta\\phi$ CR";
      abcd_title = "Control region: Low #Delta#phi";
    } else if(skim.Contains("ttbar")){
      caption = "$N_{\\rm leps}=1$ CR";
      abcd_title = "Control region: "+nleps+" = 1";
    } else {
      caption = "Search bins";
      abcd_title = "Search bins";
    }
    if (do_loose) abcd_title += " (no #DeltaR#lower[-0.1]{_{max}} cut)";

    if(method.Contains("MMMM")){
      caption = "MMMM method: all medium b-tags";
      abcd_title = "MMMM";
      nbcuts = {"nbm==2", "nbm==3", "nbm>=4"};
    } // MMMM

    //// General assignments to all methods
    abcdcuts = abcdcuts_std;
    for(size_t ind=0; ind<abcdcuts.size(); ind++)
      abcdcuts[ind].ReplaceAll("2bcuts", nbcuts[0]);
    vector<TString> bincuts = vector<TString>(nbcuts.begin()+1, nbcuts.end());

    //////// Pushing all cuts to then find the yields
    if (Contains(mm_scen,"syst_qcd") || Contains(mm_scen,"syst_comp")){ 
      // loosen selection for propagating qcd systematics by: removing delta phi and using only 3b
      // nbcut filled twice to avoid complications with printing table
      vector<TString> bincuts_tmp = vector<TString>({"nbt>=2&&nbm>=3", "nbt>=2&&nbm>=3"});
      TString basecuts_tmp = basecuts; basecuts_tmp.ReplaceAll("!low_dphi","!(dphi1<0.3 || dphi2<0.3)");
      abcds.push_back(abcd_method(method, metcuts, bincuts_tmp, abcdcuts, caption, basecuts_tmp, abcd_title));
    } else {
      abcds.push_back(abcd_method(method, metcuts, bincuts, abcdcuts, caption, basecuts, abcd_title));
    }
    if(method.Contains("noint")) abcds.back().setIntNbNj(false);

    vector<TableRow> table_cuts, table_cuts_mm;
    NamedFunc correction = do_correction ? corrections.at(mm_scen) : NamedFunc(1.);
    for(size_t icut=0; icut < abcds.back().allcuts.size(); icut++){
      table_cuts.push_back(TableRow(abcds.back().allcuts[icut].Data(), abcds.back().allcuts[icut].Data(),
                                    0,0,weights.at("no_mismeasurement")*correction));
      if(only_mc) table_cuts_mm.push_back(TableRow(abcds.back().allcuts[icut].Data(), abcds.back().allcuts[icut].Data(),
                                                   0,0,weights.at(mm_scen)));
    }
    TString tname = "preds"; tname += iabcd;
    pm.Push<Table>(tname.Data(),  table_cuts, all_procs, true, false);
    tname += mm_scen;
    if(only_mc) pm.Push<Table>(tname.Data(),  table_cuts_mm, all_procs, true, false);
  } // Loop over ABCD methods

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////// Finding all yields ///////////////////////////////////////////////

  bool single_thread = false;
  if(single_thread) pm.multithreaded_ = false;
  pm.min_print_ = true;
  pm.MakePlots(lumi);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////// Calculating preds/kappas and printing table //////////////////////////////////////
  vector<TString> tablenames;
  for(size_t imethod=0; imethod<abcds.size(); imethod++) {
    mm_scen = GetScenario(methods.at(imethod).Data());
    // allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
    // if split_bkg [2/4] Other, [3/5] tt1l, [4/6] tt2l
    vector<vector<GammaParams> > allyields;
    Table * yield_table;
    if(only_mc){
      yield_table = static_cast<Table*>(pm.Figures()[imethod*2].get());
      Table * yield_table_mm = static_cast<Table*>(pm.Figures()[imethod*2+1].get());
      allyields.push_back(yield_table_mm->Yield(proc_data.get(), lumi));
    } else {
      yield_table = static_cast<Table*>(pm.Figures()[imethod].get());
      allyields.push_back(yield_table->DataYield());
    }
    allyields.push_back(yield_table->BackgroundYield(lumi));
    if(do_signal){
      for(size_t ind=0; ind<proc_sigs.size(); ind++)
        allyields.push_back(yield_table->Yield(proc_sigs[ind].get(), lumi));
    }
    if(split_bkg){
      allyields.push_back(yield_table->Yield(proc_other.get(), lumi));
      allyields.push_back(yield_table->Yield(proc_tt1l.get(), lumi));
      allyields.push_back(yield_table->Yield(proc_tt2l.get(), lumi));
    }

    if (debug) {
      for (unsigned i=0; i<allyields.size(); i++)
        for (unsigned j=0; j<allyields[i].size(); j++)
          cout<<"allyield["<<i<<"]["<<j<<"]"<<allyields[i][j]<<endl;
      }

    //// Calculating kappa and Total bkg prediction
    vector<vector<vector<float> > > kappas, kappas_mm, kmcdat, datapreds, preds;
    vector<vector<float> > yieldsPlane = findPreds(abcds[imethod], allyields, kappas, kappas_mm, kmcdat, datapreds, preds);

    //// Print MC/Data yields, cuts applied, kappas, preds
    if(debug) printDebug(abcds[imethod], allyields, TString(baseline.Name()), kappas, kappas_mm, preds);

    //// Makes table MC/Data yields, kappas, preds, Zbi
    if(!only_kappa) {
      TString fullname = printTable(abcds[imethod], allyields, kappas, preds, yieldsPlane, proc_sigs);
      tablenames.push_back(fullname);
    }
    
    // reserve vectors to be filled in plotKappa()
    if (mm_scen=="syst_mcstat"){
      syst_names.push_back("syst_mcstat_up");   syst_values.push_back(vector<float>()); 
      syst_names.push_back("syst_mcstat_dn");   syst_values.push_back(vector<float>()); 
    } else {
      syst_names.push_back(mm_scen);
      syst_values.push_back(vector<float>());
    }

    //// Plotting kappa
    plotKappa(abcds[imethod], kappas, kappas_mm, kmcdat);

    // piggy back on one of the scenarios to get the expected data stat unc.
    if (mm_scen=="syst_mcstat"){
      syst_names.push_back("syst_datastat_up"); syst_values.push_back(vector<float>());
      syst_names.push_back("syst_datastat_dn"); syst_values.push_back(vector<float>());
      unsigned nsys = syst_values.size();
      for (auto &iplane: datapreds) {
        for (auto &ibin: iplane) {
          if (ibin[0]==0) ibin[0] = 1e-6; 
          syst_values[nsys-1].push_back(ibin[1]/ibin[0]);
          syst_values[nsys-2].push_back(ibin[2]/ibin[0]);
        }
      }
    }
  } // Loop over ABCD methods

  // print AN systematics table
  if (Contains(mm_scen,"syst")) {
    cout<<endl<<"===== AN systematics table rows:"<<endl;
    for (unsigned isys(0); isys<syst_names.size(); isys++) {
      if (syst_names[isys].Contains("_ttx")) cout<<setw(20)<<"$t\\bar{t}$+X closure";
      else if (syst_names[isys].Contains("_vjets")) cout<<setw(20)<<"V+jets closure";
      else if (syst_names[isys].Contains("_qcd")) cout<<setw(20)<<"QCD closure";
      else if (syst_names[isys].Contains("_comp")) cout<<setw(20)<<"Bkg. composition";
      else if (syst_names[isys].Contains("_bctag")) cout<<setw(20)<<"B-tag SFs";
      else if (syst_names[isys].Contains("_udsgtag")) cout<<setw(20)<<"Mistag SFs";
      else if (syst_names[isys].Contains("_mcstat")) cout<<setw(20)<<"MC stat.";
      else if (syst_names[isys].Contains("_datastat")) cout<<setw(20)<<"Data stat.";
      else cout<<setw(20)<<syst_names[isys];
      if (syst_names[isys].EndsWith("_up")) {
        if (isys>=syst_names.size()) {
          cout<<"Variation up without matching variation down "<<syst_names[isys]<<endl;
        } else {
          for (unsigned ibin(0); ibin<syst_values[isys].size(); ibin++) {
            float up(syst_values[isys][ibin]), dn(syst_values[isys+1][ibin]);
            float isysval = (up>0 ? 1 : -1)*max(up>0 ? up: 1/(1+up)-1, dn>0 ? dn: 1/(1+dn)-1);
            if (fabs(isysval)>0.01) {
              // the closure propagation come out negative because of the way we introduce the shift
              if (syst_names[isys].Contains("syst_qcd") || syst_names[isys].Contains("syst_vjets")|| syst_names[isys].Contains("syst_"))
                isysval *=-1;
              cout<<" & "<<setw(7)<<RoundNumber(isysval*100,0);
            } else {
              cout<<" & "<<setw(7)<<"$<$ 1";
            }
          }
          isys++;
        }
        cout<<" \\\\"<<endl;
      } else { 
        for (unsigned ibin(0); ibin<syst_values[isys].size(); ibin++) {
          float isysval = syst_values[isys][ibin];
          if (fabs(isysval)>0.01) {
            if (isysval<0) isysval = 1-1/(1+isysval);
            cout<<" & "<<setw(7)<<RoundNumber(isysval*100,0);
          } else { 
            cout<<" & "<<setw(7)<<"$<$ 1";
          }
        }
        cout<<" \\\\"<<endl;
      }
    }
  }

  if(!only_kappa){
    //// Printing names of ouput files
    cout<<endl<<"===== Tables to be moved to the AN/PAS/paper:"<<endl;
    for(size_t ind=0; ind<tablenames.size(); ind++){
      TString name=tablenames[ind]; name.ReplaceAll("fulltable","table");
      cout<<" mv "<<name<<"  ${tables_folder}"<<endl;
    }
    cout<<endl<<"===== Tables that can be compiled"<<endl;
    for(size_t ind=0; ind<tablenames.size(); ind++)
      cout<<" pdflatex "<<tablenames[ind]<<"  > /dev/null"<<endl;
    cout<<endl;
  }

  double seconds = (chrono::duration<double>(chrono::high_resolution_clock::now() - begTime)).count();
  TString hhmmss = HoursMinSec(seconds);
  cout<<endl<<"Finding "<<abcds.size()<<" tables took "<<round(seconds)<<" seconds ("<<hhmmss<<")"<<endl<<endl;
} // main

////////////////////////////////////////// End of main //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

string GetScenario(const string &method){
  string key = "_scen_";
  return method.substr(method.rfind(key)+key.size());
}

//// Prints table with results
// allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
// if split_bkg: [2/4] Other, [3/5] tt1l, [4/6] tt2l
TString printTable(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
                   vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &preds, 
                   vector<vector<float> > yieldsPlane, vector<shared_ptr<Process> > &proc_sigs){
  //cout<<endl<<"Printing table (significance estimation can take a bit)"<<endl;
  //// Table general parameters
  int digits = 2;
  TString ump = " & ";

  size_t Nsig = proc_sigs.size(); // Number of signal points (for now it cannot be changed)
  bool do_zbi = true;
  if(!unblind) do_zbi = false;
  size_t Ncol = 6;
  if(do_signal) Ncol += Nsig;
  if(split_bkg) Ncol += 3;
  if(only_mc) {
    Ncol -= 3;
    if(do_signal && do_zbi) Ncol += Nsig;
  } else {
    if(do_zbi) Ncol++;
  }
  TString blind_s = "$\\spadesuit$";

  //// Setting output file name
  int digits_lumi = 1;
  if(lumi < 1) digits_lumi = 3;
  if(lumi-floor(lumi)==0) digits_lumi = 0;
  TString lumi_s = RoundNumber(lumi, digits_lumi);
  TString outname = "tables/table_pred_"+skim+"_tight_lumi"+lumi_s; outname.ReplaceAll(".","p");
  if (do_loose) outname = "tables/table_pred_"+skim+"_loose_lumi"+lumi_s; outname.ReplaceAll(".","p");
  if (do_highnb) outname +="_highnb";
  if (do_midnb) outname +="_midnb";
  if(unblind) outname += "_unblind";
  else outname += "_blind";
  if(do_ht) outname += "_ht500";
  outname += "_"+abcd.method+".tex";
  ofstream out(outname);

  //// Printing main table preamble
  if(abcd.method.Contains("signal") && Ncol>7) out << "\\resizebox{\\textwidth}{!}{\n";
  out << "\\begin{tabular}[tbp!]{ l ";
  if(split_bkg) out << "|ccc";
  out << "|cc";
  if(!only_mc) {
    out << "|ccc";
    if(do_zbi) out<<"c";
  } else {
    if(do_signal)
      for(size_t ind=0; ind<Nsig; ind++)
        out<<"|c"<<(do_zbi?"c":"");
  }
  out<<"}\\hline\\hline\n";
  out<<"${\\cal L}="<<lumi_s<<"$ fb$^{-1}$ ";
  if(split_bkg) out << " & Other & Single $t$ & $t\\bar{t}$ ";
  out << "& $\\kappa$ & MC bkg.";
  if(!only_mc) out << " & Pred.& Obs. & Obs./MC "<<(do_zbi?"& Signi.":"");
  else if(do_signal) 
    for(size_t ind=0; ind<Nsig; ind++) {
      TString signame = proc_sigs[ind]->name_.c_str();
      if(do_zbi) out << "& \\multicolumn{2}{c"<<(ind<Nsig-1?"|":"")<<"}{" << signame <<"}";
      else  out << "& " << signame;
    }
  out << " \\\\ \\hline\\hline\n";

  vector<TString> binNames({"SBD, 2b", "SBD, xb", "HIG, 2b", "HIG, xb"});
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////// Printing results////////////////////////////////////////////////
  for(size_t iplane=0; iplane < abcd.planecuts.size(); iplane++) {
    out<<endl<< "\\multicolumn{"<<Ncol<<"}{c}{$"<<CodeToLatex(abcd.planecuts[iplane].Data())
       <<"$ (Obs/MC = $"<<RoundNumber(yieldsPlane[iplane][0],2)<<"\\pm"<<RoundNumber(yieldsPlane[iplane][1],2)
       <<"$)}  \\\\ \\hline\n";
    for(size_t iabcd=0; iabcd < abcd.abcdcuts.size(); iabcd++){
      for(size_t ibin=0; ibin < abcd.bincuts[iplane].size(); ibin++){
        size_t index = abcd.indexBin(iplane, ibin, iabcd);
        if(abcd.int_nbnj && iabcd%2==0 && ibin>0) continue;
        if(iabcd==3 && ibin==0) out << "\\hline" << endl;
        //// Printing bin name
        TString binName = binNames[iabcd];
        if(ibin==0) binName.ReplaceAll("xb", "3b");
        else binName.ReplaceAll("xb", "4b");
        out << binName;
        //// Printing Other, tt1l, tt2l
        if(split_bkg){
          size_t offset = (do_signal?Nsig:0);
          out << ump <<RoundNumber(allyields[offset+2][index].Yield(), digits)
              << ump <<RoundNumber(allyields[offset+3][index].Yield(), digits)
              << ump <<RoundNumber(allyields[offset+4][index].Yield(), digits);
        }
        //// Printing kappa
        out<<ump;
        if(iabcd==3) out  << "$"    << RoundNumber(kappas[iplane][ibin][0], digits)
                          << "^{+"  << RoundNumber(kappas[iplane][ibin][1], digits)
                          << "}_{-" << RoundNumber(kappas[iplane][ibin][2], digits) <<"}$ ";
        //// Printing MC Bkg yields
        out << ump << RoundNumber(allyields[1][index].Yield(), digits);
        // //// Printing background predictions
        // if(iabcd==3) out << "$"    << RoundNumber(preds[iplane][ibin][0], digits)
        //                  << "^{+"  << RoundNumber(preds[iplane][ibin][1], digits)
        //                  << "}_{-" << RoundNumber(preds[iplane][ibin][2], digits) <<"}$ ";
        if(!only_mc){
          //// Printing observed events in data and Obs/MC ratio
          if(!unblind && iabcd==3) out << ump << blind_s<< ump << blind_s;
          else {
            out << ump << RoundNumber(allyields[0][index].Yield(), 0);
            TString ratio_s = "-";
            double Nobs = allyields[0][index].Yield(), Nmc = allyields[1][index].Yield();
            double Eobs = sqrt(Nobs), Emc = allyields[1][index].Uncertainty();
            if(Nobs==0) Eobs=1;
            if(Emc>0){
              double ratio = Nobs/Nmc;
              double Eratio = sqrt(pow(Eobs/Nmc,2) + pow(Nobs*Emc/Nmc/Nmc,2));
              ratio_s = "$"+RoundNumber(ratio, 2)+"\\pm"+RoundNumber(Eratio,2)+"$";
            }
            out << ump << ratio_s;
          }
          //// Printing Zbi significance
          if(do_zbi && iabcd==3) out << ump << Zbi(allyields[0][index].Yield(), preds[iplane][ibin][0], 
                                                   preds[iplane][ibin][1], preds[iplane][ibin][2]);
          //// Printing signal yields
          if(do_signal)
            for(size_t ind=0; ind<Nsig; ind++) 
              out<<ump<<RoundNumber(allyields[2+ind][index].Yield(), digits);
        } else {// if not only_mc
          if(do_signal){
            for(size_t ind=0; ind<Nsig; ind++) {
              out<<ump<<RoundNumber(allyields[2+ind][index].Yield(), digits);
              if(do_zbi){
                out << ump;
                if(iabcd==3) 
                  out<<Zbi(allyields[0][index].Yield()+allyields[2+ind][index].Yield(),preds[iplane][ibin][0],
                           preds[iplane][ibin][1], preds[iplane][ibin][2]);
              } // if do_zbi
            } // Loop over signals
          } // if do_signal
        }
        out << "\\\\ \n";
      } // Loop over bin cuts
    } // Loop over ABCD cuts
    out << "\\hline\\hline\n";
  } // Loop over plane cuts
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //// Printing footer and closing file
  out<< "\\end{tabular}"<<endl;
  if(abcd.method.Contains("signal") && Ncol>7) out << "}\n"; // For resizebox
  out.close();

  //// Copying header and table to the compilable file
  TString fullname = outname; fullname.ReplaceAll("table_","fulltable_");
  ofstream full(fullname);
  ifstream header("txt/header.tex");
  full<<header.rdbuf();
  header.close();
  if(!abcd.method.Contains("signal")) full << "\\usepackage[landscape]{geometry}\n\n";
  full << "\\begin{document}\n\n";
  full << "\\begin{table}\n\\centering\n";
  full << "\\caption{" << abcd.caption <<".}\\vspace{0.1in}\n\\label{tab:"<<abcd.method<<"}\n";
  ifstream outtab(outname);
  full << outtab.rdbuf();
  outtab.close();
  full << "\\end{table}\n\n";
  full << "\\end{document}\n";
  full.close();

  return fullname;
} // printTable
//// Estimating significance
TString Zbi(double Nobs, double Nbkg, double Eup_bkg, double Edown_bkg){
  TString zbi_s;
  if(actualZbi){ // Old, bad Zbi
    double Nsig = Nobs-Nbkg;
    double zbi = RooStats::NumberCountingUtils::BinomialExpZ(Nsig, Nbkg, Eup_bkg/Nbkg);
    if(Nbkg==0) zbi = RooStats::NumberCountingUtils::BinomialWithTauExpZ(Nsig, Nbkg, 1/Eup_bkg);
    if(zbi<0) zbi=0;
    zbi_s = RoundNumber(zbi,1);
    if(zbi_s!="-") zbi_s = "$"+zbi_s+"\\sigma$";
    if(Nsig<=0 || Eup_bkg<=0) zbi_s = "-";
  } else zbi_s = "$"+RoundNumber(Significance(Nobs, Nbkg, Eup_bkg, Edown_bkg),1)+"\\sigma$";
  //cout<<"Zbi for Nobs "<<Nobs<<", Nbkg "<<Nbkg<<", Ebkg "<<Eup_bkg<<" is "<<zbi_s<<endl;
  return zbi_s;
}

//// Makes kappa plots
void plotKappa(abcd_method &abcd, vector<vector<vector<float> > > &kappas, 
               vector<vector<vector<float> > > &kappas_mm, vector<vector<vector<float> > > &kmcdat){

  bool label_up = false; //// Putting the MET labels at the bottom
  double markerSize = 1.1;

  //// Setting plot style
  PlotOpt opts("txt/plot_styles.txt", "Kappa");
  if(label_up) opts.BottomMargin(0.11);
  if(kappas.size() >= 1) { // Used to be 4
    opts.CanvasWidth(1600);
    markerSize = 1.5;
  }
  setPlotStyle(opts);

  struct kmarker{
    TString cut;
    int color;
    int style;
    vector<float> kappa;
  };
  //// k_ordered has all the kappas grouped in sets of nb cuts (typically, in bins of njets)
  vector<vector<vector<kmarker> > > k_ordered, kmd_ordered, k_ordered_mm;
  vector<kmarker> ind_bcuts; // nb cuts actually used in the plot
  vector<float> zz; // Zero length vector for the kmarker constructor
  // vector<kmarker> bcuts({{"nbm==1",2,21,zz}, {"nbm==2",4,20,zz}, {"nbm>=3",kGreen+3,22,zz}, 
  //                                                               {"nbm==0",kMagenta+2,23,zz}, 
  //                                                               {"nbl==0",kMagenta+2,23,zz}});
  vector<kmarker> bcuts({{"none",2,21,zz}});
   
  int cSignal = kBlue;
  float maxy = 2.4, fYaxis = 1.3;
  int nbins = 0; // Total number of njets bins (used in the base histo)
  for(size_t iplane=0; iplane < kappas.size(); iplane++) {
    k_ordered.push_back(vector<vector<kmarker> >());
    kmd_ordered.push_back(vector<vector<kmarker> >());
    k_ordered_mm.push_back(vector<vector<kmarker> >());
    for(size_t ibin=0; ibin < kappas[iplane].size(); ibin++){
      if(maxy < fYaxis*(kappas[iplane][ibin][0]+kappas[iplane][ibin][1])) 
        maxy = fYaxis*(kappas[iplane][ibin][0]+kappas[iplane][ibin][1]);
      if(maxy < fYaxis*(kmcdat[iplane][ibin][0]+kmcdat[iplane][ibin][1])) 
        maxy = fYaxis*(kmcdat[iplane][ibin][0]+kmcdat[iplane][ibin][1]);
      if(maxy < fYaxis*(kappas_mm[iplane][ibin][0])) 
        maxy = fYaxis*(kappas_mm[iplane][ibin][0]);
      TString bincut = abcd.bincuts[iplane][ibin];
      bincut.ReplaceAll(" ","");
      bincut.ReplaceAll("mm_","");
      int index;
      do{
        index = bincut.First('[');
        bincut.Remove(index, bincut.First(']')-index+1);
      }while(index>=0);
      bool found=false;
      for(size_t ib=0; ib<bcuts.size(); ib++){
        if(bincut.Contains(bcuts[ib].cut)){
          //// Storing the number of different nb cuts in ind_bcuts
          bool cutfound=false;
          for(size_t indb=0; indb<ind_bcuts.size(); indb++)
            if(ind_bcuts[indb].color == bcuts[ib].color) cutfound = true;
          if(!cutfound) ind_bcuts.push_back(bcuts[ib]);

          //// Cleaning the nb cut from the bincut
          bincut.ReplaceAll(bcuts[ib].cut+"&&","");
          for(size_t ik=0; ik<k_ordered[iplane].size(); ik++){
            //// Adding point to a given njets cut
            if(bincut==k_ordered[iplane][ik][0].cut){
              k_ordered[iplane][ik].push_back({bincut, bcuts[ib].color, bcuts[ib].style, kappas[iplane][ibin]});
              kmd_ordered[iplane][ik].push_back({bincut, bcuts[ib].color, bcuts[ib].style, kmcdat[iplane][ibin]});
              k_ordered_mm[iplane][ik].push_back({bincut, 1, bcuts[ib].style, kappas_mm[iplane][ibin]});
              found = true;
              break;
            } // if same njets cut
          } // Loop over existing ordered kappas
          //// If it doesn't correspond to any njet cut yet, create a new bin
          if(!found) {
            k_ordered[iplane].push_back(vector<kmarker>({{bincut, bcuts[ib].color, bcuts[ib].style, kappas[iplane][ibin]}}));
            kmd_ordered[iplane].push_back(vector<kmarker>({{bincut, bcuts[ib].color, bcuts[ib].style, kmcdat[iplane][ibin]}}));
            k_ordered_mm[iplane].push_back(vector<kmarker>({{bincut, 1, bcuts[ib].style, kappas_mm[iplane][ibin]}}));
            found = true;
            nbins++;
          }
        } // if bincut.Contains(bcuts[ib].cut)
      } // Loop over nb cuts

      //// If it doesn't correspond to any nb cut, create a new bin with default (color in [0], blue)
      if(!found) {
        k_ordered[iplane].push_back(vector<kmarker>({{bincut, bcuts[0].color, bcuts[0].style, kappas[iplane][ibin]}}));
        kmd_ordered[iplane].push_back(vector<kmarker>({{bincut, bcuts[0].color, bcuts[0].style, kmcdat[iplane][ibin]}}));
        k_ordered_mm[iplane].push_back(vector<kmarker>({{bincut, 1, bcuts[0].style, kappas_mm[iplane][ibin]}}));
        nbins++;
        if(ind_bcuts.size()==0) ind_bcuts.push_back(bcuts[0]);
      }
    } // Loop over bin cuts
  } // Loop over plane cuts

  //// Plotting kappas
  TCanvas can("can","");
  can.SetFillStyle(4000);
  TLine line; line.SetLineWidth(2); line.SetLineStyle(2);
  TLatex label; label.SetTextSize(0.05); label.SetTextFont(42); label.SetTextAlign(23);
  if(k_ordered.size()>3) label.SetTextSize(0.04);
  TLatex klab; klab.SetTextFont(42); klab.SetTextAlign(23);


  float minx = 0.5, maxx = nbins+0.5, miny = 0;
  if(label_up) maxy = 2.6;
  if(maxy > 4) maxy = 4;
  if(mm_scen=="syst_mcstat") maxy = 3;
  TH1D histo("histo", "", nbins, minx, maxx);
  histo.SetMinimum(miny);
  histo.SetMaximum(maxy);
  histo.GetYaxis()->CenterTitle(true);
  histo.GetXaxis()->SetLabelOffset(0.008);
  TString ytitle = "#kappa";
  if(mm_scen!="data" && mm_scen!="syst_mcstat") ytitle += " (Scen. = "+mm_scen+")";
  histo.SetTitleOffset(0.45,"y");
  histo.SetTitleSize(0.06,"y");
  histo.SetYTitle(ytitle);
  histo.Draw();

  //// Filling vx, vy vectors with kappa coordinates. Each nb cut is stored in a TGraphAsymmetricErrors
  int bin = 0;
  vector<vector<double> > vx(ind_bcuts.size()), vexh(ind_bcuts.size()), vexl(ind_bcuts.size());
  vector<vector<double> > vy(ind_bcuts.size()), veyh(ind_bcuts.size()), veyl(ind_bcuts.size());
  vector<vector<double> > vx_mm(ind_bcuts.size()), vexh_mm(ind_bcuts.size()), vexl_mm(ind_bcuts.size());
  vector<vector<double> > vy_mm(ind_bcuts.size()), veyh_mm(ind_bcuts.size()), veyl_mm(ind_bcuts.size());
  vector<vector<double> > vx_kmd(ind_bcuts.size()), vexh_kmd(ind_bcuts.size()), vexl_kmd(ind_bcuts.size());
  vector<vector<double> > vy_kmd(ind_bcuts.size()), veyh_kmd(ind_bcuts.size()), veyl_kmd(ind_bcuts.size());
  for(size_t iplane=0; iplane < k_ordered.size(); iplane++) {
    for(size_t ibin=0; ibin < k_ordered[iplane].size(); ibin++){
      bin++;
      histo.GetXaxis()->SetBinLabel(bin, CodeToRootTex(k_ordered[iplane][ibin][0].cut.Data()).c_str());
      // xval is the x position of the first marker in the group
      double xval = bin, nbs = k_ordered[iplane][ibin].size(), minxb = 0.15, binw = 0;
      // If there is more than one point in the group, it starts minxb to the left of the center of the bin
      // binw is the distance between points in the njets group
      if(nbs>1) {
        xval -= minxb;
        binw = 2*minxb/(nbs-1);
      }
      for(size_t ib=0; ib<k_ordered[iplane][ibin].size(); ib++){
        //// Finding which TGraph this point goes into by comparing the color of the TGraph and the point
        for(size_t indb=0; indb<ind_bcuts.size(); indb++){
          if(ind_bcuts[indb].color == k_ordered[iplane][ibin][ib].color){
            vx[indb].push_back(xval);
            vexl[indb].push_back(0);
            vexh[indb].push_back(0);
            vy[indb].push_back(k_ordered[iplane][ibin][ib].kappa[0]);
            veyh[indb].push_back(k_ordered[iplane][ibin][ib].kappa[1]);
            veyl[indb].push_back(k_ordered[iplane][ibin][ib].kappa[2]);

            //// MC kappas with data uncertainties
            vx_kmd[indb].push_back(xval);
            vexl_kmd[indb].push_back(0);
            vexh_kmd[indb].push_back(0);
            vy_kmd[indb].push_back(kmd_ordered[iplane][ibin][ib].kappa[0]);
            float ekmdUp = sqrt(pow(k_ordered[iplane][ibin][ib].kappa[1],2) +
                                pow(kmd_ordered[iplane][ibin][ib].kappa[1],2));
            float ekmdDown = sqrt(pow(k_ordered[iplane][ibin][ib].kappa[2],2) +
                                  pow(kmd_ordered[iplane][ibin][ib].kappa[2],2));
            veyh_kmd[indb].push_back(ekmdUp);            
            veyl_kmd[indb].push_back(ekmdDown);

            //// Data/pseudodata kappas
            vx_mm[indb].push_back(xval+0.1);
            vexl_mm[indb].push_back(0);
            vexh_mm[indb].push_back(0);
            vy_mm[indb].push_back(k_ordered_mm[iplane][ibin][ib].kappa[0]);
            if(mm_scen=="data" || mm_scen=="mc_as_data") {
              veyh_mm[indb].push_back(k_ordered_mm[iplane][ibin][ib].kappa[1]);
              veyl_mm[indb].push_back(k_ordered_mm[iplane][ibin][ib].kappa[2]);         
            } else {     
              veyh_mm[indb].push_back(0);
              veyl_mm[indb].push_back(0);
            }

            //// Printing difference between kappa and kappa_mm
            float kap = k_ordered[iplane][ibin][ib].kappa[0], kap_mm = k_ordered_mm[iplane][ibin][ib].kappa[0];
            TString text = "#Delta_{#kappa} = "+RoundNumber((kap_mm-kap)*100,0,kap)+"%";
            if (mm_scen!="syst_mcstat") syst_values.back().push_back((kap_mm-kap)/kap);
            if(mm_scen=="syst_mcstat") text = "#Delta_{#kappa} = "+RoundNumber((kap-1)*100,0,1)+"%";
            else if(mm_scen=="data") text = "#Delta_{#kappa} = "+RoundNumber((kap_mm-1)*100,0,1)+"%";
            if((abcd.method.Contains("signal")&&iplane>=2) || (abcd.method.Contains("njets1lmet200")&&iplane>=1)
               || (abcd.method.Contains("nb1l")&&iplane>=1) ) 
              klab.SetTextColor(cSignal);
            else klab.SetTextColor(1);
            klab.SetTextSize(0.045);
            if (k_ordered.size()*k_ordered[iplane].size()>8) klab.SetTextSize(0.03);
            if(mm_scen!="mc_as_data") klab.DrawLatex(xval, 0.952*maxy, text);
            //// Printing stat uncertainty of kappa_mm/kappa
            float kapUp = k_ordered[iplane][ibin][ib].kappa[1], kapDown = k_ordered[iplane][ibin][ib].kappa[2];
            float kap_mmUp = k_ordered_mm[iplane][ibin][ib].kappa[1];
            float kap_mmDown = k_ordered_mm[iplane][ibin][ib].kappa[2];
            if(mm_scen=="data" || mm_scen=="mc_as_data") {              
              text = "#sigma_{stat} = ^{+"+RoundNumber(kap_mmUp*100,0, 1)+"%}_{-"+RoundNumber(kap_mmDown*100,0, 1)+"%}";
            } else {
              text = "#sigma_{stat} = ^{+"+RoundNumber(kapUp*100,0, 1)+"%}_{-"+RoundNumber(kapDown*100,0, 1)+"%}";
              if (mm_scen=="syst_mcstat") {
                syst_values[syst_values.size()-2].push_back(kapUp/1);
                syst_values[syst_values.size()-1].push_back(-1*kapDown/1);
              }
            }
            klab.SetTextSize(0.05);
            if (k_ordered.size()*k_ordered[iplane].size()>8) klab.SetTextSize(0.035);
            klab.DrawLatex(xval, 0.888*maxy, text);
           xval += binw;
          }
        } // Loop over nb cuts in ordered TGraphs
      } // Loop over nb cuts in kappa plot
    } // Loop over bin cuts

    // Drawing line separating MET planes
    if (iplane==k_ordered.size()-2 && do_onemet) {line.SetLineStyle(1); line.SetLineWidth(2); line.SetLineColor(kOrange+3);}
    else {line.SetLineStyle(2); line.SetLineWidth(2); line.SetLineColor(kBlack);}
    if (iplane<k_ordered.size()-1) line.DrawLine(bin+0.5, miny, bin+0.5, maxy);
    // Drawing MET labels
    if(label_up) label.DrawLatex((2*bin-k_ordered[iplane].size()+1.)/2., maxy-0.1, CodeToRootTex(abcd.planecuts[iplane].Data()).c_str());
    else label.DrawLatex((2*bin-k_ordered[iplane].size()+1.)/2., -0.10*maxy, CodeToRootTex(abcd.planecuts[iplane].Data()).c_str());
  } // Loop over plane cuts

  //// Drawing legend and TGraphs
  if (debug) cout<<"Building up TGraphs"<<endl;
  int digits_lumi = 1;
  if(lumi < 1) digits_lumi = 3;
  if(lumi-floor(lumi)==0) digits_lumi = 0;
  TString lumi_s = RoundNumber(lumi, digits_lumi);
  double legX(opts.LeftMargin()+0.005), legY(1-0.03), legSingle = 0.05;
  if(label_up) legY = 0.8;
  double legW = 0.35, legH = legSingle*(ind_bcuts.size()+1)/2;
  if(ind_bcuts.size()>3) legH = legSingle*((ind_bcuts.size()+1)/2);
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(opts.LegendEntryHeight()*1.15); leg.SetFillColor(0);
  leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(42);
  leg.SetNColumns(2);
  TGraphAsymmErrors graph[20]; // There's problems with vectors of TGraphs, so using an array
  TGraphAsymmErrors graph_kmd[20]; // There's problems with vectors of TGraphs, so using an array
  TGraphAsymmErrors graph_mm[20]; // There's problems with vectors of TGraphs, so using an array
  for(size_t indb=0; indb<ind_bcuts.size(); indb++){
    graph_kmd[indb] = TGraphAsymmErrors(vx_kmd[indb].size(), &(vx_kmd[indb][0]), &(vy_kmd[indb][0]),
                                        &(vexl_kmd[indb][0]), &(vexh_kmd[indb][0]), &(veyl_kmd[indb][0]), 
                                        &(veyh_kmd[indb][0]));
    graph_kmd[indb].SetMarkerStyle(ind_bcuts[indb].style); graph_kmd[indb].SetMarkerSize(markerSize);
    graph_kmd[indb].SetMarkerColor(ind_bcuts[indb].color);
    graph_kmd[indb].SetLineColor(1); graph_kmd[indb].SetLineWidth(2);
    if(mm_scen=="mc_as_data") graph_kmd[indb].Draw("p0 same");

    graph[indb] = TGraphAsymmErrors(vx[indb].size(), &(vx[indb][0]), &(vy[indb][0]),
                                    &(vexl[indb][0]), &(vexh[indb][0]), &(veyl[indb][0]), &(veyh[indb][0]));
    graph[indb].SetMarkerStyle(ind_bcuts[indb].style); graph[indb].SetMarkerSize(markerSize);
    graph[indb].SetMarkerColor(ind_bcuts[indb].color);
    graph[indb].SetLineColor(ind_bcuts[indb].color); graph[indb].SetLineWidth(2);
    graph[indb].Draw("p0 same");

    graph_mm[indb] = TGraphAsymmErrors(vx_mm[indb].size(), &(vx_mm[indb][0]), &(vy_mm[indb][0]),
                                       &(vexl_mm[indb][0]), &(vexh_mm[indb][0]), &(veyl_mm[indb][0]), 
                                       &(veyh_mm[indb][0]));
    graph_mm[indb].SetMarkerStyle(20); graph_mm[indb].SetMarkerSize(markerSize*1.2);
    graph_mm[indb].SetMarkerColor(1);
    graph_mm[indb].SetLineColor(1); graph_mm[indb].SetLineWidth(2);
    if(mm_scen!="mc_as_data" && mm_scen!="syst_mcstat") graph_mm[indb].Draw("p0 same");

    leg.AddEntry(&graph[indb], "MC", "p");
    TString data_s = (mm_scen=="data"||mm_scen=="off"||mm_scen=="no_mismeasurement"?"Data":"Pseudodata");
    if(mm_scen!="mc_as_data" && mm_scen!="syst_mcstat") leg.AddEntry(&graph_mm[indb], data_s+" "+lumi_s+" fb^{-1}", "p");
    //leg.AddEntry(&graph[indb], CodeToRootTex(ind_bcuts[indb].cut.Data()).c_str(), "p");

  } // Loop over TGraphs
  leg.Draw();
  //if(ind_bcuts.size()>1) leg.Draw();

  //// Drawing CMS labels and line at 1
  TLatex cmslabel;
  cmslabel.SetTextSize(0.06);
  cmslabel.SetNDC(kTRUE);
  cmslabel.SetTextAlign(11);
  //cmslabel.DrawLatex(opts.LeftMargin()+0.005, 1-opts.TopMargin()+0.015,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}");
  cmslabel.SetTextAlign(31);
  //cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.015,"#font[42]{13 TeV}");
  cmslabel.SetTextSize(0.053);
  TString title = "#font[42]{"+abcd.title+"}";
  TString newSignal = "#color["; newSignal += cSignal; newSignal += "]{Signal}";
  title.ReplaceAll("Signal", newSignal);
  cmslabel.DrawLatex(1-opts.RightMargin()-0.005, 1-opts.TopMargin()+0.03, title);

  line.SetLineStyle(3); line.SetLineWidth(1);
  line.DrawLine(minx, 1, maxx, 1);

  TString fname="plots/kappa_"+skim+"_tight_" +abcd.method;
  if (do_loose) fname="plots/kappa_"+skim+"_loose_"+abcd.method;
  if (do_highnb) fname +="_highnb";
  if (do_midnb) fname +="_midnb";
  if(do_ht) fname  += "_ht500";
  lumi_s.ReplaceAll(".","p");
  fname += "_lumi"+lumi_s;
  fname += ".pdf";
  can.SaveAs(fname);
  cout<<endl<<" open "<<fname<<endl; 

}

//// Calculating kappa and Total bkg prediction
// allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
vector<vector<float> > findPreds(abcd_method &abcd, vector<vector<GammaParams> > &allyields,
               vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &kappas_mm, 
               vector<vector<vector<float> > > &kmcdat, vector<vector<vector<float> > > &datapreds,
               vector<vector<vector<float> > > &preds){
  // Powers for kappa:   ({R1, R2, D3, R4})
  vector<float> pow_kappa({ 1, -1, -1,  1});
  vector<float> pow_kk({ 1, -1, -1,  1});
  // Powers for TotBkg pred:({R1, R2, D3,  R1, R2, D3, D4})
  vector<float> pow_totpred( {-1,  1,  1,   1, -1, -1,  1});
  vector<float> pow_datapred( {-1,  1,  1});

  float val(1.), valup(1.), valdown(1.);
  vector<vector<float> > yieldsPlane;
  for(size_t iplane=0; iplane < abcd.planecuts.size(); iplane++) {
    //// Counting yields in plane without double-counting R1/R3 yields when integrated
    GammaParams NdataPlane, NmcPlane;
    for(size_t ibin=0; ibin < abcd.bincuts[iplane].size(); ibin++){
      for(size_t iabcd=0; iabcd < 4; iabcd++){
        if( ! (abcd.int_nbnj && ibin>0 && (iabcd==0||iabcd==2)) ){
          size_t index = abcd.indexBin(iplane, ibin, iabcd);
          NdataPlane += allyields[0][index];
          NmcPlane += allyields[1][index];
        }
      } // Loop over ABCD cuts
    } // Loop over bin cuts
    //cout<<"Plane "<<iplane<<": MC is "<<NmcPlane<<", data is "<<NdataPlane<<endl;
    float Nobs = NdataPlane.Yield(), Nmc = NmcPlane.Yield();
    float dataMC = Nobs/Nmc;
    float edataMC = sqrt(pow(sqrt(Nobs)/Nmc,2) + pow(Nobs*NmcPlane.Uncertainty()/Nmc/Nmc,2));
    yieldsPlane.push_back({dataMC, edataMC});
    

    kappas.push_back(vector<vector<float> >());
    kmcdat.push_back(vector<vector<float> >());
    kappas_mm.push_back(vector<vector<float> >());
    preds.push_back(vector<vector<float> >());
    datapreds.push_back(vector<vector<float> >());
    for(size_t ibin=0; ibin < abcd.bincuts[iplane].size(); ibin++){
      vector<vector<float> > entries;
      vector<vector<float> > weights;
      //// Pushing data yields for predictions
      for(size_t iabcd=0; iabcd < 3; iabcd++){
        size_t index = abcd.indexBin(iplane, ibin, iabcd);
        entries.push_back(vector<float>());
        weights.push_back(vector<float>());
        entries.back().push_back(allyields[0][index].Yield());
        weights.back().push_back(1.);
      } // Loop over ABCD cuts

      // Throwing toys to find predictions with no kappa -> used for data stat. unc. in systematics table
      val = calcKappa(entries, weights, pow_datapred, valdown, valup);
      if(valdown<0) valdown = 0;
      datapreds[iplane].push_back(vector<float>({val, valup, valdown}));

      vector<vector<float> > kentries;
      vector<vector<float> > kweights;
      vector<vector<float> > kkentries;
      vector<vector<float> > kkweights;
      vector<vector<float> > kentries_mm;
      vector<vector<float> > kweights_mm;
      //// Pushing MC yields for predictions and kappas
      for(size_t iabcd=0; iabcd < 4; iabcd++){
        size_t index = abcd.indexBin(iplane, ibin, iabcd);
        // Renormalizing MC to data
        allyields[1][index] *= dataMC;

        // Yields for predictions
        entries.push_back(vector<float>());
        weights.push_back(vector<float>());
        entries.back().push_back(allyields[1][index].NEffective());
        weights.back().push_back(allyields[1][index].Weight());
        // Yields for kappas
        kentries.push_back(vector<float>());
        kweights.push_back(vector<float>());
        kentries.back().push_back(allyields[1][index].NEffective());
        kweights.back().push_back(allyields[1][index].Weight());
        // Yields for kappas on pseudodata
        kentries_mm.push_back(vector<float>());
        kweights_mm.push_back(vector<float>());
        kentries_mm.back().push_back(allyields[0][index].Yield());
        kweights_mm.back().push_back(1.);
        // Yields for kappas_mc normalized to data
        kkentries.push_back(vector<float>());
        kkweights.push_back(vector<float>());
        kkentries.back().push_back(allyields[1][index].Yield());
        kkweights.back().push_back(1.);

      } // Loop over ABCD cuts

      // Throwing toys to find predictions and uncertainties
      val = calcKappa(entries, weights, pow_totpred, valdown, valup);
      if(valdown<0) valdown = 0;
      preds[iplane].push_back(vector<float>({val, valup, valdown}));
      // Throwing toys to find kappas and uncertainties
      val = calcKappa(kentries, kweights, pow_kappa, valdown, valup);
      if(valdown<0) valdown = 0;
      kappas[iplane].push_back(vector<float>({val, valup, valdown}));
      // Throwing toys to find kappas and uncertainties
      val = calcKappa(kentries_mm, kweights_mm, pow_kappa, valdown, valup);
      if(valdown<0) valdown = 0;
      kappas_mm[iplane].push_back(vector<float>({val, valup, valdown}));
      // Throwing toys to find kappas and uncertainties
      val = calcKappa(kkentries, kkweights, pow_kk, valdown, valup);
      if(valdown<0) valdown = 0;
      kmcdat[iplane].push_back(vector<float>({val, valup, valdown}));
    } // Loop over bin cuts
  } // Loop over plane cuts

  return yieldsPlane;
} // findPreds

// allyields: [0] data, [1] bkg, [2] T1tttt(NC), [3] T1tttt(C)
void printDebug(abcd_method &abcd, vector<vector<GammaParams> > &allyields, TString baseline,
                vector<vector<vector<float> > > &kappas, vector<vector<vector<float> > > &kappas_mm, 
                vector<vector<vector<float> > > &preds){

  int digits = 3;
  cout<<endl<<endl<<"=================== Printing cuts for method "<<abcd.method<<" ==================="<<endl;
  cout<<"-- Baseline cuts: "<<baseline<<endl;
  for(size_t iplane=0; iplane < abcd.planecuts.size(); iplane++) {
    cout<<endl<<" **** Plane "<<abcd.planecuts[iplane]<<" ***"<<endl;
    for(size_t ibin=0; ibin < abcd.bincuts[iplane].size(); ibin++){
      for(size_t iabcd=0; iabcd < abcd.abcdcuts.size(); iabcd++){
        size_t index = abcd.indexBin(iplane, ibin, iabcd);
        cout<<"MC: "<<setw(8)<<RoundNumber(allyields[1][index].Yield(),digits)
            <<"  Data: "<<setw(4)<<RoundNumber(allyields[0][index].Yield(), 0)
            <<"  - "<< abcd.allcuts[index]<<endl;
      } // Loop over ABCD cuts
      cout<<"Kappa MC = "<<RoundNumber(kappas[iplane][ibin][0],digits)<<"+"<<RoundNumber(kappas[iplane][ibin][1],digits)
          <<"-"<<RoundNumber(kappas[iplane][ibin][2],digits)
          <<", Kappa Data = "<<RoundNumber(kappas_mm[iplane][ibin][0],digits)
          <<"+"<<RoundNumber(kappas_mm[iplane][ibin][1],digits)
          <<"-"<<RoundNumber(kappas_mm[iplane][ibin][2],digits)<<", Prediction = "
          <<RoundNumber(preds[iplane][ibin][0],digits)<<"+"<<RoundNumber(preds[iplane][ibin][1],digits)
          <<"-"<<RoundNumber(preds[iplane][ibin][2],digits)<<endl;
      cout<<endl;
    } // Loop over bin cuts
  } // Loop over plane cuts

} // printDebug

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"method", required_argument, 0, 'm'},  // Method to run on (if you just want one)
      {"correct", no_argument, 0, 'c'},       // Apply correction
      {"lumi", required_argument, 0, 'l'},    // Luminosity to normalize MC with (no data)
      {"skim", required_argument, 0, 's'},    // Which skim to use: standard, 2015 data
      {"json", required_argument, 0, 'j'},    // Which JSON to use: 0p815, 2p6, 4p0, 7p65, 12p9
      {"split_bkg", no_argument, 0, 'b'},     // Prints Other, tt1l, tt2l contributions
      {"no_signal", no_argument, 0, 'n'},     // Does not print signal columns
      {"do_leptons", no_argument, 0, 'p'},    // Does tables for e/mu/emu as well
      {"do_loose", no_argument, 0, 't'},      // Use tight selection
      {"unblind", no_argument, 0, 'u'},       // Unblinds R4/D4
      {"only_mc", no_argument, 0, 'o'},       // Uses MC as data for the predictions
      {"only_kappa", no_argument, 0, 'k'},    // Only plots kappa (no table)
      {"debug", no_argument, 0, 'd'},         // Debug: prints yields and cuts used
      {"only_dilepton", no_argument, 0, '2'}, // Makes tables only for dilepton tests
      {"ht", no_argument, 0, 0},              // Cuts on ht>500 instead of st>500
      {"mm", required_argument, 0, 0},            // Mismeasurment scenario, 0 for data
      {"quick", no_argument, 0, 0},           // Used inclusive ttbar for quick testing
      {"zbi", no_argument, 0, 0},             // Use Zbi instead of toys
      {"highnb", no_argument, 0, 0},             // Do 3b and 4b for QCD CR
      {"midnb", no_argument, 0, 0},             // Do 3b and 4b for QCD CR
      {"onemet", no_argument, 0, 0},             
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "m:cs:j:udbntl:p2ok", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'm':
      only_method = optarg;
      break;
    case 'c':
      do_correction = true;
      break;
    case 'l':
      mc_lumi = optarg;
      only_mc = true;
      break;
    case 'k':
      only_kappa = true;
      only_mc = true;
      break;
    case 's':
      skim = optarg;
      break;
    case 'j':
      json = optarg;
      break;
    case 'b':
      split_bkg = true;
      break;
    case 'o':
      only_mc = true;
      break;
    case '2':
      only_dilepton = true;
      break;
    case 'p':
      do_leptons = true;
      break;
    case 'n':
      do_signal = false;
      break;
    case 't': //tight
      do_loose = false;
      break;
    case 'u':
      unblind = true;
      break;
    case 'd':
      debug = true;
      break;
    case 0:
      optname = long_options[option_index].name;
      if(optname == "ht"){
        do_ht = true;
      } else if(optname == "mm"){
        mm_scen = optarg;
      }else if(optname == "quick"){
        quick_test = true;
      }else if(optname == "zbi"){
        actualZbi = true;
      }else if(optname == "highnb"){
        do_highnb = true;
      }else if(optname == "midnb"){
        do_midnb = true;
      }else if(optname == "onemet"){
        do_onemet = true;
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
