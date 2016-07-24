#include "functions.hpp"

#include "TVector2.h"

#include "utilities.hpp"

namespace Functions{
  const NamedFunc n_isr_match("n_isr_match",[](const Baby &b) -> NamedFunc::ScalarType{
      int Nisr=0;
      for (size_t ijet(0); ijet<b.jets_pt()->size(); ijet++){
        if(!IsGoodJet(b, ijet)) continue;
        bool matched=false;
        for (size_t imc(0); imc<b.mc_pt()->size(); imc++){
          if(b.mc_status()->at(imc)!=23 || abs(b.mc_id()->at(imc))>5) continue;
          if(!(abs(b.mc_mom()->at(imc))==6 || abs(b.mc_mom()->at(imc))==23 ||
               abs(b.mc_mom()->at(imc))==24 || abs(b.mc_mom()->at(imc))==15)) continue; // In our ntuples where all taus come from W
          float dR = deltaR(b.jets_eta()->at(ijet), b.jets_phi()->at(ijet), b.mc_eta()->at(imc), b.mc_phi()->at(imc));
          if(dR<0.4){
            matched = true;
            break;
          }
        } // Loop over MC particles
        if(!matched) Nisr++;
      } // Loop over jets

      return Nisr;
    });

  const NamedFunc njets_weights_ttisr("njets_weights_ttisr", [](const Baby &b){
      return NJetsWeights_ttISR(b, false);
    });

  const NamedFunc njets_weights_visr("njets_weights_visr", [](const Baby &b) -> NamedFunc::ScalarType{
      if (b.ntrupv()<0) return 1.; // Do not reweight Data

      float wgt = b.weight()/b.eff_trig()/b.w_toppt();
      if(b.SampleType()<30 && b.SampleType()>=60) return wgt;

      int nisrjets(b.njets());
      // weights derived in DY+jets
      if      (nisrjets==0) return 0.981*wgt; //  +- 0.001
      else if (nisrjets==1) return 1.071*wgt; //  +- 0.001
      else if (nisrjets==2) return 1.169*wgt; //  +- 0.003
      else if (nisrjets==3) return 1.157*wgt; //  +- 0.007
      else if (nisrjets==4) return 1.014*wgt; //  +- 0.013
      else if (nisrjets==5) return 0.920*wgt; //  +- 0.025
      else if (nisrjets==6) return 0.867*wgt; //  +- 0.048
      else if (nisrjets>=7) return 0.935*wgt; //  +- 0.088
      else return wgt;
    });

  const NamedFunc min_dphi_lep_met("min_dphi_lep_met", [](const Baby &b) -> NamedFunc::ScalarType{
      double phi1, eta1, phi2, eta2;
      DileptonAngles(b, eta1, phi1, eta2, phi2);
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
    });

  const NamedFunc max_dphi_lep_met("max_dphi_lep_met", [](const Baby &b) -> NamedFunc::ScalarType{
      double phi1, eta1, phi2, eta2;
      DileptonAngles(b, eta1, phi1, eta2, phi2);
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
    });

  const NamedFunc min_dphi_lep_jet("min_dphi_lep_jet", [](const Baby &b) ->NamedFunc::ScalarType{
      double phi1, eta1, phi2, eta2;
      DileptonAngles(b, eta1, phi1, eta2, phi2);
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
    });

  const NamedFunc max_dphi_lep_jet("max_dphi_lep_jet", [](const Baby &b) ->NamedFunc::ScalarType{
      double phi1, eta1, phi2, eta2;
      DileptonAngles(b, eta1, phi1, eta2, phi2);
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
    });

  const NamedFunc min_dphi_met_jet("min_dphi_met_jet", [](const Baby &b) ->NamedFunc::ScalarType{
      double minphi = -1.;
      for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
        if(!IsGoodJet(b,ijet)) continue;
        double thisdphi = fabs(TVector2::Phi_mpi_pi(b.met_phi()-b.jets_phi()->at(ijet)));
        if(minphi < 0. || thisdphi < minphi){
          minphi = thisdphi;
        }
      }
      return minphi;
    });

  const NamedFunc max_dphi_met_jet("max_dphi_met_jet", [](const Baby &b) ->NamedFunc::ScalarType{
      double maxphi = -1.;
      for(size_t ijet = 0; ijet < b.jets_pt()->size(); ++ijet){
        if(!IsGoodJet(b,ijet)) continue;
        double thisdphi = fabs(TVector2::Phi_mpi_pi(b.met_phi()-b.jets_phi()->at(ijet)));
        if(maxphi < 0. || thisdphi > maxphi){
          maxphi = thisdphi;
        }
      }
      return maxphi;
    });

  const NamedFunc min_dr_lep_jet("min_dr_lep_jet", [](const Baby &b) ->NamedFunc::ScalarType{
      double phi1, eta1, phi2, eta2;
      DileptonAngles(b, eta1, phi1, eta2, phi2);
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
    });

  const NamedFunc max_dr_lep_jet("max_dr_lep_jet", [](const Baby &b) ->NamedFunc::ScalarType{
      double phi1, eta1, phi2, eta2;
      DileptonAngles(b, eta1, phi1, eta2, phi2);
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
    });

  bool IsGoodJet(const Baby &b, size_t ijet){
    return ijet<b.jets_pt()->size()
      && b.jets_pt()->at(ijet) > 30.
      && fabs(b.jets_eta()->at(ijet))<2.4
      && !b.jets_islep()->at(ijet);
  }

  bool IsGoodElectron(const Baby &b, size_t iel){
    return iel<b.els_pt()->size()
      && b.els_pt()->at(iel)>20.
      && fabs(b.els_sceta()->at(iel))<2.5
      && b.els_sigid()->at(iel)
      && b.els_miniso()->at(iel) >= 0.
      && b.els_miniso()->at(iel) < 0.1;
  }

  bool IsGoodMuon(const Baby &b, size_t imu){
    return imu<b.mus_pt()->size()
      && b.mus_pt()->at(imu)>20.
      && fabs(b.mus_eta()->at(imu))<2.4
      && b.mus_sigid()->at(imu)
      && b.mus_miniso()->at(imu) >= 0.
      && b.mus_miniso()->at(imu) < 0.2;
  }

  bool IsGoodTrack(const Baby &b, size_t itk){
    if(itk >= b.tks_pt()->size()) return false;

    if(abs(b.tks_pdg()->at(itk))==211  && b.tks_pt()->at(itk)>15. && b.tks_miniso()->at(itk)<0.1 && b.tks_mt2()->at(itk)<60 && b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05 ){
      return true;
    }else if(abs(b.tks_pdg()->at(itk))==13 && b.tks_pt()->at(itk)>10. && b.tks_miniso()->at(itk)<0.2 && b.tks_mt2()->at(itk)<80 && b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05){
      return true;
    }else if(abs(b.tks_pdg()->at(itk))==11 && b.tks_pt()->at(itk)>10. && b.tks_miniso()->at(itk)<0.2 && b.tks_mt2()->at(itk)<80 && b.tks_dz()->at(itk)<0.07 && b.tks_d0()->at(itk)<0.05){
      return true;
    }else{
      return false;
    }
  }

  NamedFunc::ScalarType NJetsWeights_ttISR(const Baby &b, bool use_baby_nisr){
    if (b.ntrupv()<0) return 1.; // Do not reweight Data

    float wgt = b.weight()/b.eff_trig()/b.w_toppt();

    int nisrjets = b.njets();
    if (b.SampleType()==20) {
      if (use_baby_nisr) nisrjets = b.nisr();
      else nisrjets = b.njets() - 2;
    }
    // weights derived in TTJets and applied using the nisr calculation algorithm
    if      (nisrjets==0) return 1.099*wgt; //  +- 0.012
    else if (nisrjets==1) return 0.969*wgt; //  +- 0.014
    else if (nisrjets==2) return 0.870*wgt; //  +- 0.020
    else if (nisrjets==3) return 0.772*wgt; //  +- 0.031
    else if (nisrjets==4) return 0.712*wgt; //  +- 0.051
    else if (nisrjets==5) return 0.661*wgt; //  +- 0.088
    else if (nisrjets>=6) return 0.566*wgt; //  +- 0.133
    else return wgt;
  }

  void DileptonAngles(const Baby &b,
                      NamedFunc::ScalarType &eta1, NamedFunc::ScalarType &phi1,
                      NamedFunc::ScalarType &eta2, NamedFunc::ScalarType &phi2){
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
}
