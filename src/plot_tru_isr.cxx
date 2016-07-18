#include "test.hpp"

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
#include "utilities.hpp"


using namespace std;
using namespace PlotOptTypes;

namespace{
  bool single_thread = false;
  bool verbose = false;
  TString method = "";
}


NamedFunc::ScalarType nisrMatch(const Baby &b);
bool isGoodJet(const Baby &b, size_t ijet);

template<typename T>
shared_ptr<Process> Proc(const string process_name, Process::Type type,
                         int color, const set<string> &files, const string &cut = "1"){
  return make_shared<Process>(process_name, type, color,
                              unique_ptr<Baby>(new T(files)),
                              cut);
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel = 6000;
  GetOptions(argc, argv);

  double lumi = 2.6;

  string bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  string foldermc(bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/unskimmed/");
  string folderb(bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/newalg/");
  string baseline = "1";
  Palette colors("txt/colors.txt", "default");

  auto proc_tt = Proc<Baby_full>("t#bar{t} 2l", Process::Type::signal, colors("tt_2l"),
    {folderb+"*TT*DiLep*.root"}, baseline);//+" && ntruleps==0 && nleps==0");
  auto proc_t1c = Proc<Baby_full>("T1tttt(1200,800)", Process::Type::signal, kRed,
    {folderb+"*SMS-T1tttt_mGluino-1200_mLSP-800*.root"}, baseline);//+" && ntruleps==0 && nleps==0");
  auto proc_t1nc = Proc<Baby_full>("T1tttt(1500,100)", Process::Type::signal, kRed,
    {folderb+"*SMS-T1tttt_mGluino-1500_mLSP-100*.root"}, baseline);
  // auto proc_t1nc_3l = Proc<Baby_full>("T1tttt(1500,100) 3l", Process::Type::signal, kGreen+3,
  //   {folderb+"*SMS-T1tttt_mGluino-1500_mLSP-100*PUSpring16_80*.root"}, baseline+" && ntruleps==4 && nleps==4");
  proc_t1c->SetLineStyle(2);
 
  auto proc_t14qc = Proc<Baby_full>("T1qqqq(1000,800)", Process::Type::signal, kAzure+2 ,
    {folderb+"*SMS-T1qqqq_mGluino-1000_mLSP-800*.root"}, baseline);
  auto proc_t14qnc = Proc<Baby_full>("T1qqqq(1400,100)", Process::Type::signal, kAzure+2 ,
    {folderb+"*SMS-T1qqqq_mGluino-1400_mLSP-100*.root"}, baseline);
  proc_t14qc->SetLineStyle(2);

  auto proc_t2ttvc = Proc<Baby_full>("T2tt(425,325)", Process::Type::signal, kGreen+2,
    {folderb+"*SMS-T2tt_mStop-425_mLSP-325*.root"}, baseline);
  proc_t2ttvc->SetLineStyle(3);
  auto proc_t2ttc = Proc<Baby_full>("T2tt(500,325)", Process::Type::signal, kGreen+2,
    {folderb+"*SMS-T2tt_mStop-500_mLSP-325*.root"}, baseline);
  proc_t2ttc->SetLineStyle(2);
  auto proc_t2ttnc = Proc<Baby_full>("T2tt(850,100)", Process::Type::signal, kGreen+2,
    {folderb+"*SMS-T2tt_mStop-850_mLSP-100*.root"}, baseline);

  vector<shared_ptr<Process> > tt_procs = {proc_tt};
  vector<shared_ptr<Process> > proc_isrpt = {proc_t1c, proc_t1nc, proc_t14qnc, proc_t14qc, proc_t2ttnc, proc_t2ttc, proc_t2ttvc,proc_tt}; 

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::preliminary)
    .Bottom(BottomType::off)
    .YAxis(YAxisType::log)
    .Stack(StackType::shapes);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  vector<PlotOpt> plot_types = {log_lumi_info, lin_lumi_info};

  NamedFunc nisr_match("nisr_match",nisrMatch);

  PlotMaker pm;

  float minx(0), maxx(460);
  int nbins(static_cast<int>((maxx-minx)/10));

  minx = 0; maxx = 800; nbins = static_cast<int>((maxx-minx)/25);
  pm.Push<HistoStack>(HistoDef(nbins, minx, maxx, "isr_tru_pt", "True ISR p_{T} [GeV]",
			       "(ntrutaush+ntrutausl)==0", "weight"), proc_isrpt, plot_types);
 
  minx = -0.5; maxx = 7.5; nbins = static_cast<int>((maxx-minx)/1);
  pm.Push<HistoStack>(HistoDef(nbins, minx, maxx, "nisr", "Number of ISR jets",
			       "1==1", "1"), proc_isrpt, plot_types);
 
  minx = -2.5; maxx =2.5; nbins = static_cast<int>((maxx-minx)/1);
  pm.Push<HistoStack>(HistoDef(nbins, minx, maxx, "nisr-(njets-2)", "true - reco ISR multiplicity",
             "nbm==2&&nleps==2", "1"), tt_procs, plot_types);
  minx = -2; maxx = 2; nbins = static_cast<int>((maxx-minx)/0.1);
  pm.Push<HistoStack>(HistoDef(nbins, minx, maxx, "(jetsys_nob_pt-isr_tru_pt)/isr_tru_pt", "(Reco-True)/True ISR p_{T} [GeV]",
			       "nbm>=2", "weight"), tt_procs, plot_types);

  minx = -300; maxx = 200; nbins = static_cast<int>((maxx-minx)/10);
  pm.Push<HistoStack>(HistoDef(nbins, minx, maxx, "jetsys_nob_pt-isr_tru_pt", "Reco-True ISR p_{T} [GeV]",
			       "nbm>=2", "weight"), tt_procs, plot_types);

  if(single_thread || verbose) pm.multithreaded_ = false;
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


NamedFunc::ScalarType nisrMatch(const Baby &b){
  int Nisr=0;
  for (size_t ijet(0); ijet<b.jets_pt()->size(); ijet++){
    if(!isGoodJet(b, ijet) || b.jets_pt()->at(ijet)<30) continue;
    bool matched=false;
    for (size_t imc(0); imc<b.mc_pt()->size(); imc++){
      if(b.mc_status()->at(imc)!=23 || abs(b.mc_id()->at(imc))>5) continue;
      if(!(abs(b.mc_mom()->at(imc))==6 || abs(b.mc_mom()->at(imc))==23 || 
	   abs(b.mc_mom()->at(imc))==24 || abs(b.mc_mom()->at(imc))==15 || abs(b.mc_mom()->at(imc))>1e6)) continue; // In our ntuples where all taus come from W
      float dR = deltaR(b.jets_eta()->at(ijet), b.jets_phi()->at(ijet), b.mc_eta()->at(imc), b.mc_phi()->at(imc));
      if(dR<0.3){
        if (verbose) cout<<"Jet: ("<<b.jets_pt()->at(ijet)<<", "<<b.jets_eta()->at(ijet)<<", "<<b.jets_phi()->at(ijet)<<"), MC: ("
	      <<b.mc_pt()->at(imc)<<", "<<b.mc_eta()->at(imc)<<", "<<b.mc_phi()->at(imc)<<"), ID "<<b.mc_id()->at(imc)<<". dR "<<dR <<endl;
	matched = true;
	break;
      }
    } // Loop over MC particles
    if(!matched) {
      Nisr++;
      if (verbose) cout<<"Jet: ("<<b.jets_pt()->at(ijet)<<", "<<b.jets_eta()->at(ijet)<<", "<<b.jets_phi()->at(ijet)<<" not matched"<<endl;
    }
  } // Loop over jets
  if (verbose){
    for (size_t imc(0); imc<b.mc_pt()->size(); imc++){
      if(b.mc_status()->at(imc)!=23 || abs(b.mc_id()->at(imc))>5) continue;
      if(!(abs(b.mc_mom()->at(imc))==6 || abs(b.mc_mom()->at(imc))==23 || 
    	 abs(b.mc_mom()->at(imc))==24 || abs(b.mc_mom()->at(imc))==15 || abs(b.mc_mom()->at(imc))>1e6)) continue; // In our ntuples where all taus come from W
      cout<<" MC: ("
    	<<b.mc_pt()->at(imc)<<", "<<b.mc_eta()->at(imc)<<", "<<b.mc_phi()->at(imc)<<"), ID "<<b.mc_id()->at(imc)<<endl;
    }
    cout<<" ======== New event: njets "<<b.njets()<<", Nisr "<<Nisr<<endl<<endl;
  }
  return Nisr;
}

bool isGoodJet(const Baby &b, size_t ijet){
  return ijet<b.jets_pt()->size()
      && fabs(b.jets_eta()->at(ijet))<2.4
      && !b.jets_islep()->at(ijet);
}
