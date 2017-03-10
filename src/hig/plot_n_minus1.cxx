// study differences in mjj as a function of b-cat for single vs dilepton ttbar

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/plot_opt.hpp"
#include "core/functions.hpp"
#include "hig/hig_functions.hpp"

using namespace std;
using namespace PlotOptTypes;

void GetOptions(int argc, char *argv[]);

namespace{
  string sample = "tt";
  bool paper = true;
}

const NamedFunc fakeb_fromw("fakeb_fromw",[](const Baby &b) -> NamedFunc::ScalarType{
  int nfakesw(0);
  vector<unsigned> mc_match;
  for (unsigned ijet(0); ijet<b.jets_pt()->size(); ijet++){
    // Finding jets used in the di-Higgs reco that are not true b's
    if (!b.jets_h1()->at(ijet) && !b.jets_h2()->at(ijet)) continue;
    if (b.jets_hflavor()->at(ijet)!=5) {
      float mindR=999.;
      unsigned imc_min;
      for (unsigned imc(0); imc<b.mc_id()->size(); imc++){
	bool already_matched = false;
	for(unsigned ind=0; ind<mc_match.size(); ind++)
	  if(imc == mc_match[ind]) already_matched = true;
	if(already_matched) continue;

	if(abs(b.mc_mom()->at(imc)) != 24 || abs(b.mc_id()->at(imc))>5) continue;

	float dR= deltaR(b.jets_eta()->at(ijet), b.jets_phi()->at(ijet), b.mc_eta()->at(imc), b.mc_phi()->at(imc));
	if(dR < mindR) {
	  mindR = dR;
	  imc_min = imc;
	}
      } // Loop over MC
      if(mindR < 0.4) {
	nfakesw++;
	mc_match.push_back(imc_min);
      }
    } // If jet is not matched to b-quark
  } // Loop over jets

  if(nfakesw<2 && false){
    cout<<endl<<endl<<"Nfakes = "<<nfakesw<<endl;
    for (unsigned ijet(0); ijet<b.jets_pt()->size(); ijet++){
      if (!b.jets_h1()->at(ijet) && !b.jets_h2()->at(ijet)) continue;
      if (b.jets_hflavor()->at(ijet)!=5) {
	cout<<"Jet "<<setw(4)<<ijet<<": ("<<setw(8)<<b.jets_pt()->at(ijet)<<", "<<setw(8)<<b.jets_eta()->at(ijet)<<", "<<
	  setw(8)<<b.jets_phi()->at(ijet)<<")"<<endl;
      }
    }
    for (unsigned imc(0); imc<b.mc_id()->size(); imc++){
      if(abs(b.mc_mom()->at(imc)) != 24 || abs(b.mc_id()->at(imc))>5) continue;
	cout<<"MC  "<<setw(4)<<imc<<": ("<<setw(8)<<b.mc_pt()->at(imc)<<", "<<setw(8)<<b.mc_eta()->at(imc)<<", "<<
	  setw(8)<<b.mc_phi()->at(imc)<<")"<<endl;
    }
  }    
  return nfakesw;
});

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  Palette colors("txt/colors.txt", "default");


  /////////////////// PLOT STYLES //////////////////////////////////////
  PlotOpt lin_norm_info("txt/plot_styles.txt", "CMSPaper");
  lin_norm_info.Title(TitleType::info)   
    .Bottom(BottomType::off)
    .YAxis(YAxisType::linear)
    .Stack(StackType::data_norm);
  PlotOpt log_norm_info = lin_norm_info().YAxis(YAxisType::log);
  vector<PlotOpt> plt_norm_info = {lin_norm_info, log_norm_info};

  PlotOpt log_norm = lin_norm_info().YAxis(YAxisType::log).Title(TitleType::info).LogMinimum(.7).Bottom(BottomType::ratio);
  PlotOpt lin_norm = lin_norm_info().YAxis(YAxisType::linear).Title(TitleType::info).Bottom(BottomType::ratio);
  if (paper) {
    lin_norm.Title(TitleType::preliminary);
    log_norm.Title(TitleType::preliminary);
  }
  vector<PlotOpt> plt_norm = {lin_norm, log_norm};

  PlotOpt lin_shapes = lin_norm().Stack(StackType::shapes).Bottom(BottomType::ratio);
  vector<PlotOpt> plt_shapes = {lin_shapes};

  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::off);
  vector<PlotOpt> plt_shapes_info = {lin_shapes_info};



  //////////////////////////////////// PROCESSES /////////////////////////////
  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////// Higgsino //////////////////////////////////////////////
  //string folderhigmc = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higlep1/";
  string folderhigmc = bfolder+"/cms2r0/babymaker/babies/2017_01_27/mc/merged_higmc_higloose/";
  string folderhigdata = bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/merged_higdata_higloose/";
  string foldersig(bfolder+"/cms2r0/babymaker/babies/2017_02_26/TChiHH/merged_higmc_higloose/");

  map<string, set<string>> mctags; 
  mctags["tt"]     = set<string>({"*_TTJets*Lept*.root"});
  mctags["ttx"]     = set<string>({"*_TTJets*Lept*.root", "*_TTZ*.root", "*_TTW*.root",
                                     "*_TTGJets*.root", "*_ttHJetTobb*.root","*_TTTT*.root"});
  mctags["vjets"]   = set<string>({"*_ZJet*.root", "*_WJetsToLNu*.root", "*DYJetsToLL*.root"});
  mctags["singlet"] = set<string>({"*_ST_*.root"});
  mctags["qcd"]     = set<string>({"*QCD_HT*0_Tune*.root", "*QCD_HT*Inf_Tune*.root"});
  mctags["other"]   = set<string>({"*_WH_HToBB*.root", "*_ZH_HToBB*.root",
                                     "*_WWTo*.root", "*_WZ*.root", "*_ZZ_*.root"});
  set<string> allmctags;
  for (auto &iset: mctags) {
      allmctags.insert(iset.second.begin(), iset.second.end());
  }

  NamedFunc wgt = Higfuncs::weight_higd * Higfuncs::eff_higtrig;
  NamedFunc base_func = "pass && pass_ra2_badmu && met/met_calo<5 && nbt>=2 && nvleps==0&&ntks==0&&!low_dphi&&met>150";

  vector<shared_ptr<Process> > procs_hig;
  procs_hig.push_back(Process::MakeShared<Baby_full>("QCD",        
    Process::Type::background, colors("other"),    attach_folder(folderhigmc,mctags["qcd"]),     base_func&&"stitch")); 
  procs_hig.push_back(Process::MakeShared<Baby_full>("t#bar{t}+X", 
    Process::Type::background, colors("tt_1l"),    attach_folder(folderhigmc,mctags["ttx"]),     base_func&&"stitch"));
  procs_hig.push_back(Process::MakeShared<Baby_full>("V+jets",     
    Process::Type::background, kOrange+1,          attach_folder(folderhigmc,mctags["vjets"]),   base_func&&"stitch"));
  procs_hig.push_back(Process::MakeShared<Baby_full>("Single t",   
    Process::Type::background, colors("single_t"), attach_folder(folderhigmc,mctags["singlet"]), base_func&&"stitch"));
  procs_hig.push_back(Process::MakeShared<Baby_full>("Other",      
    Process::Type::background, kGreen+1,           attach_folder(folderhigmc,mctags["other"]),   base_func&&"stitch"));      

  vector<string> sig2m = {"225","400","700"}; 
  vector<int> sig2_colors = {kGreen, kRed, kBlue}; // need sigm.size() >= sig_colors.size()
  for (unsigned isig(0); isig<sig2m.size(); isig++)
    procs_hig.push_back(Process::MakeShared<Baby_full>("TChiHH("+sig2m[isig]+",1)", Process::Type::signal, 
			sig2_colors[isig], {foldersig+"*TChiHH_mGluino-"+sig2m[isig]+"*.root"}, base_func));


  /////////////// MC+DATA PROCESSES
  vector<shared_ptr<Process> > procs_data = procs_hig;
  procs_data.push_back(Process::MakeShared<Baby_full>("Data", Process::Type::data, kBlack,
    {folderhigdata+"*root"},  Higfuncs::trig_hig>0. && base_func)); 


  /////////////// SIGNAL PROCESSES
  vector<shared_ptr<Process> > procs_sig;
  vector<string> sigm = {"225","300", "400","700","1000"}; 
  vector<int> sig_colors = {kGreen+1, kTeal-4, kRed, kBlue, kOrange}; // need sigm.size() >= sig_colors.size()
  for (unsigned isig(0); isig<sigm.size(); isig++)
    procs_sig.push_back(Process::MakeShared<Baby_full>("TChiHH("+sigm[isig]+",1)", Process::Type::signal, 
			sig_colors[isig], {foldersig+"*TChiHH_mGluino-"+sigm[isig]+"*.root"}, base_func));

  /////////////// TTBAR ONLY PROCESSES
  vector<shared_ptr<Process> > procs_tt;
  procs_tt.push_back(Process::MakeShared<Baby_full>("t#bar{t}(1l)", 
    Process::Type::background, colors("tt_1l"), 
    attach_folder(folderhigmc,mctags["tt"]), base_func&&"stitch&&ntruleps==1"));
  ///////////////////////////////////////////////////////////////////////////////////////////////////


  string cuts = "nbm>=2";
  PlotMaker pm;

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////// N-1 //////////////////////////////////////////////
  cuts = "nbdt>=2&&met>150&&higd_drmax<2.2&&higd_dm<40&&njets>=4&&njets<=5";
  pm.Push<Hist1D>(Axis(9, 150., 600., "met", "E^{miss}_{T} [GeV]", {150., 200., 300., 450.}),
    cuts, procs_data, plt_norm).Weight(wgt).Tag("n1");

  // cuts = "nbdt>=2&&nbdm>=3&&met>150&&higd_drmax<2.2&&higd_dm<40&&njets>=4&&njets<=5";
  // pm.Push<Hist1D>(Axis(9, 150., 600., "met", "E^{miss}_{T} [GeV]", {150., 200., 300., 450.}),
  //   cuts, procs_data, plt_norm).Weight(wgt).Tag("n1");

  cuts = "nbdt>=2&&nbdm>=3&&nbdl>=4&&met>150&&higd_drmax<2.2&&higd_dm<40&&njets>=4&&njets<=5";
  pm.Push<Hist1D>(Axis(20, 0., 200., "higd_am", "#LTm#GT [GeV]", {100., 140.}),
    cuts, procs_data, plt_norm).Weight(wgt).Tag("n1");

  cuts = "nbdt>=2&&nbdm>=3&&nbdl>=4&&met>150&&higd_drmax<2.2&&njets>=4&&njets<=5";
  pm.Push<Hist1D>(Axis(15, 0., 120., "higd_dm", "#Deltam [GeV]", {40.}),
    cuts, procs_data, plt_norm).Weight(wgt).Tag("n1");

  cuts = "nbdt>=2&&nbdm>=3&&nbdl>=4&&met>150&&higd_dm<40&&njets>=4&&njets<=5";
  pm.Push<Hist1D>(Axis(20,0,4,"higd_drmax", "#DeltaR_{max} [GeV]", {2.2}),
    cuts, procs_data, plt_norm).Weight(wgt).Tag("n1");

  ///////////////////////////////////////////////////////////////////////////////////////////////////



  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////// MIA //////////////////////////////////////////////
  // cuts = "nbdt>=2&&nbdm>=3&&nbdl>=4&&met>150&& njets>=4&&njets<=5&&higd_dm<20";
  // pm.Push<Hist1D>(Axis(50, 0., 250., "higd_am", "#LTm#GT [GeV]", {100., 140.}),cuts, procs_sig, plt_shapes_info).Weight(wgt).Tag("mia");
  // cuts = "nbdt>=2&&nbdm>=3&&nbdl>=4&&met>150&& njets>=4&&njets<=5&&higd_dm>20&&higd_dm<40";
  // pm.Push<Hist1D>(Axis(50, 0., 250., "higd_am", "#LTm#GT [GeV]", {100., 140.}),cuts, procs_sig, plt_shapes_info).Weight(wgt).Tag("mia");
  // cuts = "nbdt>=2&&nbdm>=3&&nbdl>=4&&met>150&& njets>=4&&njets<=5&&higd_dm>40";
  // pm.Push<Hist1D>(Axis(50, 0., 250., "higd_am", "#LTm#GT [GeV]", {100., 140.}),cuts, procs_sig, plt_shapes_info).Weight(wgt).Tag("mia");

  // cuts = "nbdt>=2&&nbdm>=3&&nbdl>=4&&met>150&&higd_drmax<2.2&&njets>=4&&njets<=5";
  // pm.Push<Hist1D>(Axis(30, 0., 150., "higd_dm", "#Deltam [GeV]", {40.}),cuts, procs_hig, plt_norm).Weight(wgt).Tag("mia");

  // cuts = "nbdt>=2&&nbdm>=3&&nbdl>=4&&met>150&&higd_dm<40&&njets>=4&&njets<=5";
  // pm.Push<Hist1D>(Axis(20,0,4,"higd_drmax", "#DeltaR_{max} [GeV]", {2.2}),cuts, procs_hig, plt_norm).Weight(wgt).Tag("mia");

  // ///////////////////////////////////////////////////////////////////////////////////////////////////



  // ///////////////////////////////////////////////////////////////////////////////////////////////////
  // ////////////////////////////////////// FILIP //////////////////////////////////////////////
  // cuts = "nbdt>=2&&nbdm==3&&nbdl==3&&met>150&&higd_drmax<2.2&&higd_dm<40&& njets>=4&&njets<=5";
  // pm.Push<Hist1D>(Axis(5, -0.5, 4.5, Higfuncs::nb_exci, "b quarks from flavor excitation"),cuts, procs_tt, plt_shapes_info).Weight(wgt).Tag("filip");
  // pm.Push<Hist1D>(Axis(5, -0.5, 4.5, Higfuncs::nb_gs, "b quarks from gluon splitting"),cuts, procs_tt, plt_shapes_info).Weight(wgt).Tag("filip");
  // pm.Push<Hist1D>(Axis(5, -0.5, 4.5, Higfuncs::nb_gs+Higfuncs::nb_exci, "b quarks not from top"),cuts, procs_tt, plt_shapes_info).Weight(wgt).Tag("filip");

  // cuts = "nbdt>=2&&nbdm>=3&&nbdl>=4&&met>150&&higd_drmax<2.2&&higd_dm<40&& njets>=4&&njets<=5";
  // pm.Push<Hist1D>(Axis(5, -0.5, 4.5, Higfuncs::nb_exci, "b quarks from flavor excitation"),cuts, procs_tt, plt_shapes_info).Weight(wgt).Tag("filip");
  // pm.Push<Hist1D>(Axis(5, -0.5, 4.5, Higfuncs::nb_gs, "b quarks from gluon splitting"),cuts, procs_tt, plt_shapes_info).Weight(wgt).Tag("filip");
  // pm.Push<Hist1D>(Axis(5, -0.5, 4.5, Higfuncs::nb_gs+Higfuncs::nb_exci, "b quarks not from top"),cuts, procs_tt, plt_shapes_info).Weight(wgt).Tag("filip");

  // cuts = "nbdt>=2&&nbdm>=3&&nbdl>=4&&met>150&&higd_drmax<2.2&&higd_dm<40&& njets==4";
  // pm.Push<Hist1D>(Axis(5, -0.5, 4.5, fakeb_fromw, "Fake b-tags from W"),cuts, procs_tt, plt_shapes_info).Weight(wgt).Tag("filip");
  // pm.Push<Hist1D>(Axis(5, -0.5, 4.5, Higfuncs::ntrub, "Higgs jets truthmatched to b quarks"),cuts, procs_tt, plt_shapes_info).Weight(wgt).Tag("filip");
  // pm.Push<Hist1D>(Axis(5, -0.5, 4.5, 4.-Higfuncs::ntrub- fakeb_fromw, "Fake b-tags not from W"),cuts, procs_tt, plt_shapes_info).Weight(wgt).Tag("filip");

  // cuts = "nbdt>=2&&nbdm>=3&&nbdl>=4&&met>150&&higd_drmax<2.2&&higd_dm<40&& njets==5";
  // pm.Push<Hist1D>(Axis(5, -0.5, 4.5, fakeb_fromw, "Fake b-tags from W"),cuts, procs_tt, plt_shapes_info).Weight(wgt).Tag("filip");
  // pm.Push<Hist1D>(Axis(5, -0.5, 4.5, Higfuncs::ntrub, "Higgs jets truthmatched to b quarks"),cuts, procs_tt, plt_shapes_info).Weight(wgt).Tag("filip");
  // pm.Push<Hist1D>(Axis(5, -0.5, 4.5, 4.-Higfuncs::ntrub-fakeb_fromw, "Fake b-tags not from W"),cuts, procs_tt, plt_shapes_info).Weight(wgt).Tag("filip");
  // ///////////////////////////////////////////////////////////////////////////////////////////////////



  pm.min_print_ = true;
  pm.MakePlots(35.9);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
} // main

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"sample", required_argument, 0, 's'},    
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 's':
      sample = optarg;
      break;
    case 0:
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
