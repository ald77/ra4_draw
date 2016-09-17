#include "ra4/plot_dilep_angles.hpp"

#include <cmath>

#include "TError.h"
#include "TVector2.h"

#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/hist1d.hpp"
#include "core/event_scan.hpp"

using namespace std;
using namespace PlotOptTypes;

int main(){
  gErrorIgnoreLevel = 6000;

  double lumi = 2.6;

  Palette colors("txt/colors.txt", "default");

  string folder_mc = "/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/merged_standard/";
  auto tt1l = Process::MakeShared<Baby_full>("t#bar{t} (1l)", Process::Type::background, colors("tt_1l"),
    {folder_mc+"*_TTJets*Lept*.root", folder_mc+"*_TTJets_HT*.root"}, "ntruleps<=1&&stitch");
  auto tt2l = Process::MakeShared<Baby_full>("t#bar{t} (2l)", Process::Type::background, colors("tt_2l"),
    {folder_mc+"*_TTJets*Lept*.root", folder_mc+"*_TTJets_HT*.root"}, "ntruleps>=2&&stitch");
  auto wjets = Process::MakeShared<Baby_full>("W+jets", Process::Type::background, colors("wjets"),
    {folder_mc+"*_WJetsToLNu*.root"});
  auto single_t = Process::MakeShared<Baby_full>("Single t", Process::Type::background, colors("single_t"),
    {folder_mc+"*_ST_*.root"});
  auto ttv = Process::MakeShared<Baby_full>("t#bar{t}V", Process::Type::background, colors("ttv"),
    {folder_mc+"*_TTWJets*.root", folder_mc+"*_TTZTo*.root"});
  auto other = Process::MakeShared<Baby_full>("Other", Process::Type::background, colors("other"),
    {folder_mc+"*DYJetsToLL*.root", folder_mc+"*_QCD_HT*.root",
        folder_mc+"*_ZJet*.root", folder_mc+"*_WWTo*.root",
        folder_mc+"*ggZH_HToBB*.root", folder_mc+"*ttHJetTobb*.root",
        folder_mc+"*_TTGJets*.root", folder_mc+"*_TTTT_*.root",
        folder_mc+"*_WH_HToBB*.root", folder_mc+"*_WZTo*.root",
        folder_mc+"*_ZH_HToBB*.root", folder_mc+"_ZZ_*.root"});
  auto data_2016 = Process::MakeShared<Baby_full>("2016 Data", Process::Type::data, kBlack,
    {"/net/cms2/cms2r0/babymaker/babies/2016_06_21/data/skim_standard/*.root"},
    "pass&&(trig[4]||trig[8]||trig[13]||trig[33])");
  auto data_2015 = Process::MakeShared<Baby_full>("2015 Data", Process::Type::data, kRed,
    {"/net/cms2/cms2r0/babymaker/babies/2016_04_29/data/merged_1lht500met200/*.root"},
    "pass&&(trig[4]||trig[8]||trig[13]||trig[28])");

  vector<shared_ptr<Process> > procs = {data_2016, data_2015, tt1l, tt2l, wjets, single_t, ttv, other};

  PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::info)
    .Bottom(BottomType::off)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  vector<PlotOpt> plot_types = {log_lumi, lin_lumi};

  PlotMaker pm;
  double pi = acos(-1.);

  NamedFunc baseline = "met<=500&&nbm<=2&&ht>500&&met>200";
  NamedFunc dilep_all = baseline && "(nleps==2&&njets>=4) || (nleps==1&&nveto==1&&njets>=5&&nbm>=1&&mt>140)";
  NamedFunc dilep_ee = baseline && "nels==2&&nmus==0&&njets>=4";
  NamedFunc dilep_em = baseline && "nels==1&&nmus==1&&njets>=4";
  NamedFunc dilep_mm = baseline && "nels==0&&nmus==2&&njets>=4";
  NamedFunc dilep_ev = baseline && "nels==1&&nmus==0&&nveto==1&&njets>=5&&nbm>=1&&mt>140";
  NamedFunc dilep_mv = baseline && "nels==0&&nmus==1&&nveto==1&&njets>=5&&nbm>=1&&mt>140";

  NamedFunc dphi_ll("dphi_ll", [](const Baby &b) -> NamedFunc::ScalarType{
      double phi1, eta1, phi2, eta2;
      GetAngles(b, phi1, eta1, phi2, eta2);
      if(phi1 == -999 || phi2 == -999) return -1.;
      return fabs(TVector2::Phi_mpi_pi(phi2-phi1));
    });
  NamedFunc min_dphi_lmet("min_dphi_lmet", MinDeltaPhiLepMet);
  NamedFunc max_dphi_lmet("max_dphi_lmet", MaxDeltaPhiLepMet);
  NamedFunc min_dphi_lj("min_dphi_lj", MinDeltaPhiLepJet);
  NamedFunc max_dphi_lj("max_dphi_lj", MaxDeltaPhiLepJet);
  NamedFunc min_dphi_metj("min_dphi_metj", MinDeltaPhiMetJet);
  NamedFunc max_dphi_metj("max_dphi_metj", MaxDeltaPhiMetJet);

  NamedFunc dr_ll("dr_ll", [](const Baby &b) -> NamedFunc::ScalarType{
      double phi1, eta1, phi2, eta2;
      GetAngles(b, phi1, eta1, phi2, eta2);
      if(phi1 == -999 || phi2 == -999) return -1.;
      return hypot(TVector2::Phi_mpi_pi(phi2-phi1), eta2-eta1);
    });
  NamedFunc min_dr_lj("min_dr_lj", MinDeltaRLepJet);
  NamedFunc max_dr_lj("max_dr_lj", MaxDeltaRLepJet);

  vector<pair<string, NamedFunc> > categories = {{"all", dilep_all}, {"ee", dilep_ee}, {"em", dilep_em},
                                                 {"mm", dilep_mm}, {"ev", dilep_ev}, {"mv", dilep_mv}};

  for(const auto &cat: categories){
    std::string tag = "angles_"+cat.first;
    const NamedFunc &cut = cat.second;

    pm.Push<Hist1D>(Axis(20, 0., pi, dphi_ll, "#Delta#phi(l,l)"), cut, procs, plot_types).Tag(tag);
    pm.Push<Hist1D>(Axis(20, 0., pi, min_dphi_lmet, "Min. #Delta#phi(l,MET)"), cut, procs, plot_types).Tag(tag);
    pm.Push<Hist1D>(Axis(20, 0., pi, max_dphi_lmet, "Max. #Delta#phi(l,MET)"), cut, procs, plot_types).Tag(tag);
    pm.Push<Hist1D>(Axis(20, 0., pi, min_dphi_lj, "Min. #Delta#phi(l,jet)"), cut, procs, plot_types).Tag(tag);
    pm.Push<Hist1D>(Axis(20, 0., pi, max_dphi_lj, "Max. #Delta#phi(l,jet)"), cut, procs, plot_types).Tag(tag);
    pm.Push<Hist1D>(Axis(20, 0., pi, min_dphi_metj, "Min. #Delta#phi(MET,jet)"), cut, procs, plot_types).Tag(tag);
    pm.Push<Hist1D>(Axis(20, 0., pi, max_dphi_metj, "Max. #Delta#phi(MET,jet)"), cut, procs, plot_types).Tag(tag);

    pm.Push<Hist1D>(Axis(20, 0., 6, dr_ll, "#Delta R(l,l)"), cut, procs, plot_types).Tag(tag);
    pm.Push<Hist1D>(Axis(20, 0., 6, min_dr_lj, "Min. #Delta R(l,jet)"), cut, procs, plot_types).Tag(tag);
    pm.Push<Hist1D>(Axis(20, 0., 6, max_dr_lj, "Max. #Delta R(l,jet)"), cut, procs, plot_types).Tag(tag);
  }

  pm.Push<EventScan>("dilep_angle", dilep_mm || dilep_em, vector<NamedFunc>{
      "run", "lumiblock", "event", "nmus", "nels", "nveto",
        dphi_ll, dr_ll, min_dphi_lmet, max_dphi_lmet, min_dphi_lj,
        max_dphi_lj, min_dr_lj, max_dr_lj, min_dphi_metj, max_dphi_metj,
        "mj14", "mt", "njets", "nbm", "met", "met_phi", "leps_pt",
        "leps_phi", "leps_eta", "leps_id", "jets_pt", "jets_phi", "jets_eta"
        }, procs, 10);

  //pm.multithreaded_ = false;
  pm.MakePlots(lumi);
}

bool IsGoodMuon(const Baby &b, size_t imu){
  return imu<b.mus_pt()->size()
    && b.mus_pt()->at(imu)>20.
    && fabs(b.mus_eta()->at(imu))<2.4
                                  && b.mus_sigid()->at(imu)
    && b.mus_miniso()->at(imu) >= 0.
    && b.mus_miniso()->at(imu) < 0.2;
}

bool IsGoodElectron(const Baby &b, size_t iel){
  return iel<b.els_pt()->size()
    && b.els_pt()->at(iel)>20.
    && fabs(b.els_sceta()->at(iel))<2.5
                                    && b.els_sigid()->at(iel)
    && b.els_miniso()->at(iel) >= 0.
    && b.els_miniso()->at(iel) < 0.1;
}

bool IsGoodTrack(const Baby &b, size_t itk){
  if(itk >= b.tks_pt()->size()) return false;

  if(fabs(b.tks_pdg()->at(itk))==211  && b.tks_pt()->at(itk)>15. && b.tks_miniso()->at(itk)<0.1 && b.tks_mt2()->at(itk)<60 && b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05 ){
    return true;
  }else if(fabs(b.tks_pdg()->at(itk))==13 && b.tks_pt()->at(itk)>10. && b.tks_miniso()->at(itk)<0.2 && b.tks_mt2()->at(itk)<80 && b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05){
    return true;
  }else if(fabs(b.tks_pdg()->at(itk))==11 && b.tks_pt()->at(itk)>10. && b.tks_miniso()->at(itk)<0.2 && b.tks_mt2()->at(itk)<80 && b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05){
    return true;
  }else{
    return false;
  }
}

bool IsGoodJet(const Baby &b, size_t ijet){
  if(ijet >= b.jets_pt()->size()) return false;

  return b.jets_pt()->at(ijet) > 30.
    && b.jets_eta()->at(ijet) < 2.
                                && !b.jets_islep()->at(ijet);
}

void GetAngles(const Baby &b,
               double &phi1, double &eta1,
               double &phi2, double &eta2){
  phi1 = -999.; eta1 = -999.;
  phi2 = -999.; eta2 = -999.;
  bool h1=false, h2=false;
  if(b.nels()==2 && b.nmus()==0){
    for(size_t iel = 0; iel < b.els_pt()->size() && !(h1&&h2); ++iel){
      if(IsGoodElectron(b, iel)){
        if(!h1){
          phi1 = b.els_phi()->at(iel);
          eta1 = b.els_sceta()->at(iel);
          h1 = true;
        }else if(!h2){
          phi2 = b.els_phi()->at(iel);
          eta2 = b.els_sceta()->at(iel);
          h2 = true;
        }
      }
    }
  }else if(b.nels()==1 && b.nmus()==1){
    for(size_t iel = 0; iel < b.els_pt()->size() && !h1; ++iel){
      if(IsGoodElectron(b, iel)){
        if(!h1){
          phi1 = b.els_phi()->at(iel);
          eta1 = b.els_sceta()->at(iel);
          h1 = true;
        }else if(!h2){
          phi2 = b.els_phi()->at(iel);
          eta2 = b.els_sceta()->at(iel);
          h2 = true;
        }
      }
    }
    for(size_t imu = 0; imu < b.mus_pt()->size() && h1 && !h2; ++imu){
      if(IsGoodMuon(b, imu)){
        if(!h1){
          phi1 = b.mus_phi()->at(imu);
          eta1 = b.mus_eta()->at(imu);
          h1 = true;
        }else if(!h2){
          phi2 = b.mus_phi()->at(imu);
          eta2 = b.mus_eta()->at(imu);
          h2 = true;
        }
      }
    }
  }else if(b.nels()==0 && b.nmus()==2){
    for(size_t imu = 0; imu < b.mus_pt()->size() && !(h1&&h2); ++imu){
      if(IsGoodMuon(b, imu)){
        if(!h1){
          phi1 = b.mus_phi()->at(imu);
          eta1 = b.mus_eta()->at(imu);
          h1 = true;
        }else if(!h2){
          phi2 = b.mus_phi()->at(imu);
          eta2 = b.mus_eta()->at(imu);
          h2 = true;
        }
      }
    }
  }else if(b.nels()==1 && b.nmus()==0 && b.nveto()==1){
    for(size_t iel = 0; iel < b.els_pt()->size() && !h1; ++iel){
      if(IsGoodElectron(b, iel)){
        if(!h1){
          phi1 = b.els_phi()->at(iel);
          eta1 = b.els_sceta()->at(iel);
          h1 = true;
        }else if(!h2){
          phi2 = b.els_phi()->at(iel);
          eta2 = b.els_sceta()->at(iel);
          h2 = true;
        }
      }
    }
    for(size_t itk = 0; itk < b.tks_pt()->size() && h1 && !h2; ++itk){
      if(IsGoodTrack(b, itk)){
        if(!h1){
          phi1 = b.tks_phi()->at(itk);
          eta1 = b.tks_eta()->at(itk);
          h1 = true;
        }else if(!h2){
          phi2 = b.tks_phi()->at(itk);
          eta2 = b.tks_eta()->at(itk);
          h2 = true;
        }
      }
    }
  }else if(b.nels()==0 && b.nmus()==1 && b.nveto()==1){
    for(size_t imu = 0; imu < b.mus_pt()->size() && !h1; ++imu){
      if(IsGoodMuon(b, imu)){
        if(!h1){
          phi1 = b.mus_phi()->at(imu);
          eta1 = b.mus_eta()->at(imu);
          h1 = true;
        }else if(!h2){
          phi2 = b.mus_phi()->at(imu);
          eta2 = b.mus_eta()->at(imu);
          h2 = true;
        }
      }
    }
    for(size_t itk = 0; itk < b.tks_pt()->size() && h1 && !h2; ++itk){
      if(IsGoodTrack(b, itk)){
        if(!h1){
          phi1 = b.tks_phi()->at(itk);
          eta1 = b.tks_eta()->at(itk);
          h1 = true;
        }else if(!h2){
          phi2 = b.tks_phi()->at(itk);
          eta2 = b.tks_eta()->at(itk);
          h2 = true;
        }
      }
    }
  }
}

NamedFunc::ScalarType MinDeltaPhiLepMet(const Baby &b){
  double phi1, eta1, phi2, eta2;
  GetAngles(b, phi1, eta1, phi2, eta2);
  double dphi1 = fabs(TVector2::Phi_mpi_pi(phi1-b.met_phi()));
  double dphi2 = fabs(TVector2::Phi_mpi_pi(phi2-b.met_phi()));
  if(phi1 != -999 && phi2 != -999){
    return std::min(dphi1,dphi2);
  }else if(phi1 != -999 && phi2 == -999){
    return dphi1;
  }else if(phi1 == -999 && phi2 != -999){
    return dphi2;
  }
  return -1;
}

NamedFunc::ScalarType MaxDeltaPhiLepMet(const Baby &b){
  double phi1, eta1, phi2, eta2;
  GetAngles(b, phi1, eta1, phi2, eta2);
  double dphi1 = fabs(TVector2::Phi_mpi_pi(phi1-b.met_phi()));
  double dphi2 = fabs(TVector2::Phi_mpi_pi(phi2-b.met_phi()));
  if(phi1 != -999 && phi2 != -999){
    return std::max(dphi1,dphi2);
  }else if(phi1 != -999 && phi2 == -999){
    return dphi1;
  }else if(phi1 == -999 && phi2 != -999){
    return dphi2;
  }
  return -1;
}

NamedFunc::ScalarType MinDeltaPhiLepJet(const Baby &b){
  double phi1, eta1, phi2, eta2;
  GetAngles(b, phi1, eta1, phi2, eta2);
  double minphi = -1.;
  for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
    if(!IsGoodJet(b,ijet)) continue;
    double dphi1 = fabs(TVector2::Phi_mpi_pi(phi1-b.jets_phi()->at(ijet)));
    double dphi2 = fabs(TVector2::Phi_mpi_pi(phi2-b.jets_phi()->at(ijet)));
    double thisdphi = -1;
    if(phi1 != -999 && phi2 != -999){
      thisdphi = std::min(dphi1, dphi2);
    }else if(phi1 != -999 && phi2 == -999){
      thisdphi = dphi1;
    }else if(phi1 == -999 && phi2 != -999){
      thisdphi = dphi2;
    }
    if(minphi < 0. || thisdphi < minphi){
      minphi = thisdphi;
    }
  }
  return minphi;
}

NamedFunc::ScalarType MaxDeltaPhiLepJet(const Baby &b){
  double phi1, eta1, phi2, eta2;
  GetAngles(b, phi1, eta1, phi2, eta2);
  double maxphi = -1.;
  for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
    if(!IsGoodJet(b,ijet)) continue;
    double dphi1 = fabs(TVector2::Phi_mpi_pi(phi1-b.jets_phi()->at(ijet)));
    double dphi2 = fabs(TVector2::Phi_mpi_pi(phi2-b.jets_phi()->at(ijet)));
    double thisdphi = -1;
    if(phi1 != -999 && phi2 != -999){
      thisdphi = std::max(dphi1, dphi2);
    }else if(phi1 != -999 && phi2 == -999){
      thisdphi = dphi1;
    }else if(phi1 == -999 && phi2 != -999){
      thisdphi = dphi2;
    }
    if(maxphi < 0. || thisdphi > maxphi){
      maxphi = thisdphi;
    }
  }
  return maxphi;
}

NamedFunc::ScalarType MinDeltaPhiMetJet(const Baby &b){
  double minphi = -1.;
  for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
    if(!IsGoodJet(b,ijet)) continue;
    double thisdphi = fabs(TVector2::Phi_mpi_pi(b.met_phi()-b.jets_phi()->at(ijet)));
    if(minphi < 0. || thisdphi < minphi){
      minphi = thisdphi;
    }
  }
  return minphi;
}

NamedFunc::ScalarType MaxDeltaPhiMetJet(const Baby &b){
  double maxphi = -1.;
  for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
    if(!IsGoodJet(b,ijet)) continue;
    double thisdphi = fabs(TVector2::Phi_mpi_pi(b.met_phi()-b.jets_phi()->at(ijet)));
    if(maxphi < 0. || thisdphi > maxphi){
      maxphi = thisdphi;
    }
  }
  return maxphi;
}

NamedFunc::ScalarType MinDeltaRLepJet(const Baby &b){
  double phi1, eta1, phi2, eta2;
  GetAngles(b, phi1, eta1, phi2, eta2);
  double minr = -1.;
  for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
    if(!IsGoodJet(b,ijet)) continue;
    double dr1 = hypot(TVector2::Phi_mpi_pi(phi1-b.jets_phi()->at(ijet)), eta2-eta1);
    double dr2 = hypot(TVector2::Phi_mpi_pi(phi2-b.jets_phi()->at(ijet)), eta2-eta1);
    double thisdr = -1;
    if(phi1 != -999 && phi2 != -999){
      thisdr = std::min(dr1, dr2);
    }else if(phi1 != -999 && phi2 == -999){
      thisdr = dr1;
    }else if(phi1 == -999 && phi2 != -999){
      thisdr = dr2;
    }
    if(minr < 0. || thisdr < minr){
      minr = thisdr;
    }
  }
  return minr;
}

NamedFunc::ScalarType MaxDeltaRLepJet(const Baby &b){
  double phi1, eta1, phi2, eta2;
  GetAngles(b, phi1, eta1, phi2, eta2);
  double maxr = -1.;
  for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
    if(!IsGoodJet(b,ijet)) continue;
    double dr1 = hypot(TVector2::Phi_mpi_pi(phi1-b.jets_phi()->at(ijet)), eta2-eta1);
    double dr2 = hypot(TVector2::Phi_mpi_pi(phi2-b.jets_phi()->at(ijet)), eta2-eta1);
    double thisdr = -1;
    if(phi1 != -999 && phi2 != -999){
      thisdr = std::max(dr1, dr2);
    }else if(phi1 != -999 && phi2 == -999){
      thisdr = dr1;
    }else if(phi1 == -999 && phi2 != -999){
      thisdr = dr2;
    }
    if(maxr < 0. || thisdr > maxr){
      maxr = thisdr;
    }
  }
  return maxr;
}
