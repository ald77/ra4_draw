#include <set>
#include <string>
#include <iostream>
#include <memory>

#include "TStyle.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TError.h"
#include "TGraphAsymmErrors.h"

#include "timer.hpp"
#include "baby_full.hpp"

using namespace std;

void Fill(bool pass, TH1D &h_pass, TH1D &h_total, double x, double w){
  if(pass) h_pass.Fill(x, w);
  h_total.Fill(x, w);
}

void Format(TH1D &h, bool no_stats = false){
  h.SetLineColor(kBlack);
  h.SetLineStyle(1);
  h.SetLineWidth(5);
  if(no_stats) h.SetStats(false);
}

void Format(TH2D &h, bool no_stats = false){
  //h.SetMarkerStyle(20);
  //h.SetMarkerSize(0.5);
  h.SetMarkerColor(kBlack);
  h.SetMinimum(0.);
  if(no_stats) h.SetStats(false);
}

void Print(TH1D &h, bool no_stats = false){
  Format(h, no_stats);
  TCanvas c;
  h.Draw("e1p");
  c.Print((string("plots/")+h.GetName()+".pdf").c_str());
}

void Print(TH2D &h, bool no_stats = false){
  Format(h, no_stats);
  h.SetTitle(("#rho="+to_string(h.GetCorrelationFactor())).c_str());
  TCanvas c;
  h.Draw("scat");
  c.Print((string("plots/")+h.GetName()+".pdf").c_str());
}

void Print(TH1D &pass, TH1D &total){
  TCanvas c;
  TGraphAsymmErrors g(&pass, &total);
  g.SetTitle((string(";")+pass.GetXaxis()->GetTitle()+";Fraction Mismeasured").c_str());
  g.SetMinimum(0.);
  g.SetMaximum(1.);
  g.Draw("ap0");
  c.Print((string("plots/")+pass.GetName()+".pdf").c_str());
}

int main(){
  gErrorIgnoreLevel = 6000;
  gStyle->SetOptStat("iouMRen");

  double lumi = 2.6;
  int igood = 0;
  int ibad = 2;

  string folder_mc = "/net/cms2/cms2r0/babymaker/babies/mismeasured/2016_06_14/mc/merged_mm_std_nj5mj250/";
  //folder_mc = "/net/cms2/cms2r0/babymaker/babies/mismeasured/2016_06_14/mc/skim_1lh1500met200/";
  folder_mc="/net/cms26/cms26r0/babymaker/babies/mismeasured_v2/2016_06_14/mc/mm_std_nj5mj250/";
  auto baby_nontt = make_shared<Baby_full>(set<string>{
      folder_mc+"*_DYJetsToLL*.root",
	folder_mc+"*_QCD_HT*.root",
	folder_mc+"*_ST*channel*.root",
	folder_mc+"*_TTGJets*.root",
	folder_mc+"*_TTTT*.root",
	folder_mc+"*_TTWJetsTo*.root",
	folder_mc+"*_TTZTo*.root",
	folder_mc+"*_WH_HToBB*.root",
	folder_mc+"*_WJetsToLNu_HT-*.root",
	folder_mc+"*_WWTo*.root",
	folder_mc+"*_WZTo*.root",
	folder_mc+"*_ZH*.root",
	folder_mc+"*_ZZ*.root",
	folder_mc+"*_ttHJetTobb*.root"
	});
  auto baby_tt = make_shared<Baby_full>(set<string>{
      folder_mc+"*_TTJets_*.root",
	});

  TH1D h_1l_mm_lep_pt("h_1l_mm_lep_pt", ";Mismeasured Lepton p_{T} [GeV];Entries", 20, 0., 2000.);
  TH2D h_1l_mt_mj_lep("h_1l_mt_mj_lep", ";M_{J} (with lep) [GeV];m_{T} [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_1l_mt_mj_nolep("h_1l_mt_mj_nolep", ";M_{J} (no lep) [GeV];m_{T} [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_1l_mm_mt_mm_mj_lep("h_1l_mm_mt_mm_mj_lep", ";M_{J} (with lep) [GeV];m_{T} [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_1l_mm_mt_mm_mj_nolep("h_1l_mm_mt_mm_mj_nolep", ";M_{J} (no lep) [GeV];m_{T} [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_1l_mj_nolep_mj_lep("h_1l_mj_nolep_mj_lep", ";M_{J} (with lep) [GeV];M_{J} (no lep) [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_1l_mm_mj_nolep_mm_mj_lep("h_1l_mm_mj_nolep_mm_mj_lep", ";M_{J} (with lep) [GeV];M_{J} (no lep) [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_1l_mm_mt_mt("h_1l_mm_mt_mt", ";Correct m_{T} [GeV];Mismeasured m_{T} [GeV]", 100, -1000., 1000., 100, 0., 2000.);
  TH2D h_1l_mm_mj_nolep_mj_nolep("h_1l_mm_mj_nolep_mj_nolep", ";Correct M_{J} (no lep) [GeV];Mismeasured M_{J} (no lep) [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_1l_mm_mj_lep_mj_lep("h_1l_mm_mj_lep_mj_lep", ";Correct M_{J} (with lep) [GeV];Mismeasured M_{J} (with lep) [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH1D h_2l_mm_lep_pt("h_2l_mm_lep_pt", ";Mismeasured Lepton p_{T} [GeV];Entries", 20, 0., 2000.);
  TH2D h_2l_mt_mj_lep("h_2l_mt_mj_lep", ";M_{J} (with lep) [GeV];m_{T} [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_2l_mt_mj_nolep("h_2l_mt_mj_nolep", ";M_{J} (no lep) [GeV];m_{T} [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_2l_mm_mt_mm_mj_lep("h_2l_mm_mt_mm_mj_lep", ";M_{J} (with lep) [GeV];m_{T} [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_2l_mm_mt_mm_mj_nolep("h_2l_mm_mt_mm_mj_nolep", ";M_{J} (no lep) [GeV];m_{T} [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_2l_mj_nolep_mj_lep("h_2l_mj_nolep_mj_lep", ";M_{J} (with lep) [GeV];M_{J} (no lep) [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_2l_mm_mj_nolep_mm_mj_lep("h_2l_mm_mj_nolep_mm_mj_lep", ";M_{J} (with lep) [GeV];M_{J} (no lep) [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_2l_mm_mt_mt("h_2l_mm_mt_mt", ";Correct m_{T} [GeV];Mismeasured m_{T} [GeV]", 100, -1000., 1000., 100, 0., 2000.);
  TH2D h_2l_mm_mj_nolep_mj_nolep("h_2l_mm_mj_nolep_mj_nolep", ";Correct M_{J} (no lep) [GeV];Mismeasured M_{J} (no lep) [GeV]", 100, 0., 1000., 100, 0., 1000.);
  TH2D h_2l_mm_mj_lep_mj_lep("h_2l_mm_mj_lep_mj_lep", ";Correct M_{J} (with lep) [GeV];Mismeasured M_{J} (with lep) [GeV]", 100, 0., 1000., 100, 0., 1000.);

  TH1D h_num_1l_mt("h_num_1l_mt", ";m_{T} [GeV];Fraction Mismeasured", 20, 0., 1000.);
  TH1D h_den_1l_mt("h_den_1l_mt", ";m_{T} [GeV];Fraction Mismeasured", 20, 0., 1000.);
  TH1D h_num_1l_mj_lep("h_num_1l_mj_lep", ";M_{J} (with lep) [GeV];Fraction Mismeasured", 20, 0., 1000.);
  TH1D h_den_1l_mj_lep("h_den_1l_mj_lep", ";M_{J} (with lep) [GeV];Fraction Mismeasured", 20, 0., 1000.);
  TH1D h_num_1l_mj_nolep("h_num_1l_mj_nolep", ";M_{J} (no lep) [GeV];Fraction Mismeasured", 20, 0., 1000.);
  TH1D h_den_1l_mj_nolep("h_den_1l_mj_nolep", ";M_{J} (no lep) [GeV];Fraction Mismeasured", 20, 0., 1000.);
  TH1D h_num_1l_met("h_num_1l_met", ";MET [GeV];Fraction Mismeasured", 20, 0., 1000.);
  TH1D h_den_1l_met("h_den_1l_met", ";MET [GeV];Fraction Mismeasured", 20, 0., 1000.);
  TH1D h_num_1l_njets("h_num_1l_njets", ";N_{jets};Fraction Mismeasured", 7, 4.5, 11.5);
  TH1D h_den_1l_njets("h_den_1l_njets", ";N_{jets};Fraction Mismeasured", 7, 4.5, 11.5);
  TH1D h_num_1l_nbm("h_num_1l_nbm", ";N_{b};Fraction Mismeasured", 5, -0.5, 4.5);
  TH1D h_den_1l_nbm("h_den_1l_nbm", ";N_{b};Fraction Mismeasured", 5, -0.5, 4.5);
  TH1D h_num_2l_mt("h_num_2l_mt", ";m_{T} [GeV];Fraction Mismeasured", 20, 0., 1000.);
  TH1D h_den_2l_mt("h_den_2l_mt", ";m_{T} [GeV];Fraction Mismeasured", 20, 0., 1000.);
  TH1D h_num_2l_mj_lep("h_num_2l_mj_lep", ";M_{J} (with lep) [GeV];Fraction Mismeasured", 20, 0., 1000.);
  TH1D h_den_2l_mj_lep("h_den_2l_mj_lep", ";M_{J} (with lep) [GeV];Fraction Mismeasured", 20, 0., 1000.);
  TH1D h_num_2l_mj_nolep("h_num_2l_mj_nolep", ";M_{J} (no lep) [GeV];Fraction Mismeasured", 20, 0., 1000.);
  TH1D h_den_2l_mj_nolep("h_den_2l_mj_nolep", ";M_{J} (no lep) [GeV];Fraction Mismeasured", 20, 0., 1000.);
  TH1D h_num_2l_met("h_num_2l_met", ";MET [GeV];Fraction Mismeasured", 20, 0., 1000.);
  TH1D h_den_2l_met("h_den_2l_met", ";MET [GeV];Fraction Mismeasured", 20, 0., 1000.);
  TH1D h_num_2l_njets("h_num_2l_njets", ";N_{jets};Fraction Mismeasured", 7, 4.5, 11.5);
  TH1D h_den_2l_njets("h_den_2l_njets", ";N_{jets};Fraction Mismeasured", 7, 4.5, 11.5);
  TH1D h_num_2l_nbm("h_num_2l_nbm", ";N_{b};Fraction Mismeasured", 5, -0.5, 4.5);
  TH1D h_den_2l_nbm("h_den_2l_nbm", ";N_{b};Fraction Mismeasured", 5, -0.5, 4.5);

  size_t ibaby = 0;
  for(auto &b: {baby_nontt, baby_tt}){
    cout << "Getting entries..." << endl;
    long num_entries = b->GetEntries();
    cout << num_entries << " entries found." << endl;
    Timer timer(num_entries, 1.);
    for(long event = 0; event < num_entries; ++event){
      timer.Iterate();
      b->GetEntry(event);

      if(ibaby == 1 && !b->stitch()) continue;
      double w = 10.*lumi*b->weight();

      bool mm = b->mm()->at(ibad);
      if(mm
	 && (b->mm_nleps()->at(igood)>0 || b->mm_nleps()->at(ibad)>0)
	 && (b->mm_ht()->at(igood)>500 || b->mm_ht()->at(ibad)>500)
	 && (b->mm_met()->at(igood)>200 || b->mm_met()->at(ibad)>200)
	 && (b->mm_njets()->at(igood)>=5 || b->mm_njets()->at(ibad)>=5)){
	if(b->mm_nleps()->at(ibad)==1){
	  h_1l_mm_lep_pt.Fill(b->mm_lep_pt()->at(ibad), w);
	  h_1l_mt_mj_lep.Fill(b->mm_mj14_lep()->at(igood), b->mm_mt()->at(igood), w);
	  h_1l_mt_mj_nolep.Fill(b->mm_mj14_nolep()->at(igood), b->mm_mt()->at(igood), w);
	  h_1l_mm_mt_mm_mj_lep.Fill(b->mm_mj14_lep()->at(ibad), b->mm_mt()->at(ibad), w);
	  h_1l_mm_mt_mm_mj_nolep.Fill(b->mm_mj14_nolep()->at(ibad), b->mm_mt()->at(ibad), w);
	  h_1l_mj_nolep_mj_lep.Fill(b->mm_mj14_lep()->at(igood), b->mm_mj14_nolep()->at(igood), w);
	  h_1l_mm_mj_nolep_mm_mj_lep.Fill(b->mm_mj14_lep()->at(ibad), b->mm_mj14_nolep()->at(ibad), w);
	  h_1l_mm_mt_mt.Fill(b->mm_mt()->at(igood), b->mm_mt()->at(ibad), w);
	  h_1l_mm_mj_nolep_mj_nolep.Fill(b->mm_mj14_nolep()->at(igood), b->mm_mj14_nolep()->at(ibad), w);
	  h_1l_mm_mj_lep_mj_lep.Fill(b->mm_mj14_lep()->at(igood), b->mm_mj14_lep()->at(ibad), w);
	}else if(b->mm_nleps()->at(ibad)==2){
	  h_2l_mm_lep_pt.Fill(b->mm_lep_pt()->at(ibad), w);
	  h_2l_mt_mj_lep.Fill(b->mm_mj14_lep()->at(igood), b->mm_mt()->at(igood), w);
	  h_2l_mt_mj_nolep.Fill(b->mm_mj14_nolep()->at(igood), b->mm_mt()->at(igood), w);
	  h_2l_mm_mt_mm_mj_lep.Fill(b->mm_mj14_lep()->at(ibad), b->mm_mt()->at(ibad), w);
	  h_2l_mm_mt_mm_mj_nolep.Fill(b->mm_mj14_nolep()->at(ibad), b->mm_mt()->at(ibad), w);
	  h_2l_mj_nolep_mj_lep.Fill(b->mm_mj14_lep()->at(igood), b->mm_mj14_nolep()->at(igood), w);
	  h_2l_mm_mj_nolep_mm_mj_lep.Fill(b->mm_mj14_lep()->at(ibad), b->mm_mj14_nolep()->at(ibad), w);
	  h_2l_mm_mt_mt.Fill(b->mm_mt()->at(igood), b->mm_mt()->at(ibad), w);
	  h_2l_mm_mj_nolep_mj_nolep.Fill(b->mm_mj14_nolep()->at(igood), b->mm_mj14_nolep()->at(ibad), w);
	  h_2l_mm_mj_lep_mj_lep.Fill(b->mm_mj14_lep()->at(igood), b->mm_mj14_lep()->at(ibad), w);
	}
      }

      if(b->mm_ht()->at(ibad)>500. && b->mm_met()->at(ibad)>200.
	 && b->mm_njets()->at(ibad)>=5 && b->mm_nbm()->at(ibad)>=1){
	if(b->mm_nleps()->at(ibad)==1){
	  Fill(mm, h_num_1l_mt, h_den_1l_mt, b->mm_mt()->at(ibad), w);
	  Fill(mm, h_num_1l_mj_lep, h_den_1l_mj_lep, b->mm_mj14_lep()->at(ibad), w);
	  Fill(mm, h_num_1l_mj_nolep, h_den_1l_mj_nolep, b->mm_mj14_nolep()->at(ibad), w);
	  Fill(mm, h_num_1l_met, h_den_1l_met, b->mm_met()->at(ibad), w);
	  Fill(mm, h_num_1l_njets, h_den_1l_njets, b->mm_njets()->at(ibad), w);
	  Fill(mm, h_num_1l_nbm, h_den_1l_nbm, b->mm_nbm()->at(ibad), w);
	}else if(b->mm_nleps()->at(ibad)>=2){
	  Fill(mm, h_num_2l_mt, h_den_2l_mt, b->mm_mt()->at(ibad), w);
	  Fill(mm, h_num_2l_mj_lep, h_den_2l_mj_lep, b->mm_mj14_lep()->at(ibad), w);
	  Fill(mm, h_num_2l_mj_nolep, h_den_2l_mj_nolep, b->mm_mj14_nolep()->at(ibad), w);
	  Fill(mm, h_num_2l_met, h_den_2l_met, b->mm_met()->at(ibad), w);
	  Fill(mm, h_num_2l_njets, h_den_2l_njets, b->mm_njets()->at(ibad), w);
	  Fill(mm, h_num_2l_nbm, h_den_2l_nbm, b->mm_nbm()->at(ibad), w);
	}
      }
    }
    ++ibaby;
  }

  Print(h_num_1l_mt, h_den_1l_mt);
  Print(h_num_1l_mj_lep, h_den_1l_mj_lep);
  Print(h_num_1l_mj_nolep, h_den_1l_mj_nolep);
  Print(h_num_1l_met, h_den_1l_met);
  Print(h_num_1l_njets, h_den_1l_njets);
  Print(h_num_1l_nbm, h_den_1l_nbm);
  Print(h_num_2l_mt, h_den_2l_mt);
  Print(h_num_2l_mj_lep, h_den_2l_mj_lep);
  Print(h_num_2l_mj_nolep, h_den_2l_mj_nolep);
  Print(h_num_2l_met, h_den_2l_met);
  Print(h_num_2l_njets, h_den_2l_njets);
  Print(h_num_2l_nbm, h_den_2l_nbm);

  Print(h_1l_mm_lep_pt);
  Print(h_1l_mt_mj_lep);
  Print(h_1l_mt_mj_nolep);
  Print(h_1l_mm_mt_mm_mj_lep);
  Print(h_1l_mm_mt_mm_mj_nolep);
  Print(h_1l_mj_nolep_mj_lep);
  Print(h_1l_mm_mj_nolep_mm_mj_lep);
  Print(h_1l_mm_mt_mt);
  Print(h_1l_mm_mj_nolep_mj_nolep);
  Print(h_1l_mm_mj_lep_mj_lep);
  Print(h_2l_mm_lep_pt);
  Print(h_2l_mt_mj_lep);
  Print(h_2l_mt_mj_nolep);
  Print(h_2l_mm_mt_mm_mj_lep);
  Print(h_2l_mm_mt_mm_mj_nolep);
  Print(h_2l_mj_nolep_mj_lep);
  Print(h_2l_mm_mj_nolep_mm_mj_lep);
  Print(h_2l_mm_mt_mt);
  Print(h_2l_mm_mj_nolep_mj_nolep);
  Print(h_2l_mm_mj_lep_mj_lep);
}
