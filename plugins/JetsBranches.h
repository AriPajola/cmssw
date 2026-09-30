#ifndef MuonAnalysis_MuonAnalyzer_JetBranches
#define MuonAnalysis_MuonAnalyzer_JetBranches

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <type_traits>
#include "TLorentzVector.h"

#include "NtupleContent.h"
#include "StandAloneNtupleContent.h"
#include "helper.h"

template <typename JET, typename MU>
inline bool CrossClean(const JET &jet, std::vector<MU> &muForJetCleaning) {
  bool muFakeJet = false;
  for (const auto &mu : muForJetCleaning) {
    if (deltaR(jet, mu) < 0.4)
      muFakeJet = true;
  }
  return muFakeJet;
}


// ---------------------------------------------------------
// Generic implementation used by both ntuple formats
// ---------------------------------------------------------

template <typename JET, typename NT>
inline void FillJetBranchesImpl(const JET &jet,
                                const JET &corrJet,
                                NT &nt,
                                std::string era) {
  bool jetTightID = true;
  bool jetTightLepVeto = true;

  if (era.find("2016") != std::string::npos) {
    if (abs(jet.eta()) <= 2.7) {
      if (jet.neutralHadronEnergyFraction() >= 0.9 ||
          jet.neutralEmEnergyFraction() >= 0.9 ||
          jet.numberOfDaughters() <= 1) {
        jetTightID = false;
        jetTightLepVeto = false;
      }

      if (jet.muonEnergyFraction() >= 0.8)
        jetTightLepVeto = false;

      if (abs(jet.eta()) < 2.4) {
        if (jet.chargedHadronEnergyFraction() <= 0 ||
            jet.chargedMultiplicity() <= 0) {
          jetTightID = false;
          jetTightLepVeto = false;
        }

        if (jet.chargedEmEnergyFraction() >= 0.9)
          jetTightLepVeto = false;

        if (jet.chargedEmEnergyFraction() >= 0.99)
          jetTightID = false;
      }
    }
  }

  else if (era.find("2017") != std::string::npos) {
    if (abs(jet.eta()) <= 2.7) {
      if (jet.neutralHadronEnergyFraction() >= 0.9 ||
          jet.neutralEmEnergyFraction() >= 0.9 ||
          jet.numberOfDaughters() <= 1) {
        jetTightID = false;
        jetTightLepVeto = false;
      }

      if (jet.muonEnergyFraction() >= 0.8)
        jetTightLepVeto = false;

      if (abs(jet.eta()) < 2.4) {
        if (jet.chargedHadronEnergyFraction() <= 0 ||
            jet.chargedMultiplicity() <= 0) {
          jetTightID = false;
          jetTightLepVeto = false;
        }

        if (jet.chargedEmEnergyFraction() >= 0.8)
          jetTightLepVeto = false;
      }
    }
  }

  else if (era.find("2018") != std::string::npos) {
    if (abs(jet.eta()) <= 2.6) {
      if (jet.neutralHadronEnergyFraction() >= 0.9 ||
          jet.neutralEmEnergyFraction() >= 0.9 ||
          jet.numberOfDaughters() <= 1 ||
          jet.chargedHadronEnergyFraction() <= 0 ||
          jet.chargedMultiplicity() <= 0) {
        jetTightID = false;
        jetTightLepVeto = false;
      }

      if (jet.muonEnergyFraction() >= 0.8 ||
          jet.chargedEmEnergyFraction() >= 0.8) {
        jetTightLepVeto = false;
      }
    }

    else if (abs(jet.eta()) <= 2.7) {
      if (jet.neutralHadronEnergyFraction() >= 0.9 ||
          jet.neutralEmEnergyFraction() >= 0.99 ||
          jet.chargedMultiplicity() <= 0) {
        jetTightID = false;
        jetTightLepVeto = false;
      }

      if (jet.muonEnergyFraction() >= 0.8 ||
          jet.chargedEmEnergyFraction() >= 0.8) {
        jetTightLepVeto = false;
      }
    }
  }

  if (abs(jet.eta()) > 2.7 && abs(jet.eta()) <= 3.0) {
    if (jet.neutralEmEnergyFraction() <= 0.02 ||
        jet.neutralEmEnergyFraction() >= 0.99 ||
        jet.neutralMultiplicity() <= 2) {
      jetTightID = false;
      jetTightLepVeto = false;
    }
  }
  else if (abs(jet.eta()) > 3.0) {
    if (jet.neutralHadronEnergyFraction() <= 0.2 ||
        jet.neutralEmEnergyFraction() >= 0.9 ||
        jet.neutralMultiplicity() <= 10) {
      jetTightID = false;
      jetTightLepVeto = false;
    }
  }

  if (jetTightID)
    nt.nTightJets++;

  if (jetTightLepVeto)
    nt.nTightLepVetoJets++;

  nt.jets_isTight.push_back(jetTightID);
  nt.jets_isTightLepVeto.push_back(jetTightLepVeto);
  nt.jets_pt.push_back(corrJet.pt());
  nt.jets_eta.push_back(corrJet.eta());
  nt.jets_phi.push_back(corrJet.phi());
  nt.jets_mass.push_back(corrJet.mass());
}


// NtupleContent overload
template <typename JET>
inline void FillJetBranches(const JET &jet,
                            const JET &corrJet,
                            NtupleContent &nt,
                            std::string era) {
  FillJetBranchesImpl(jet, corrJet, nt, era);
}


// StandAloneNtupleContent overload
template <typename JET>
inline void FillJetBranches(const JET &jet,
                            const JET &corrJet,
                            StandAloneNtupleContent &nt,
                            std::string era) {
  FillJetBranchesImpl(jet, corrJet, nt, era);
}


// ---------------------------------------------------------
// Generic FindJetProbePair implementation
// ---------------------------------------------------------

template <typename JET, typename MUON, typename NT>
inline void FindJetProbePairImpl(const std::vector<JET> &jets,
                                 const MUON &mu,
                                 NT &nt) {
  TLorentzVector muonP4;
  muonP4.SetPtEtaPhiM(mu.pt(), mu.eta(), mu.phi(), mu.mass());

  float minDR = 999.f;
  TLorentzVector closestJetP4;

  for (const auto &jet : jets) {
    TLorentzVector jetP4;
    jetP4.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.mass());

    if (jetP4.DeltaR(muonP4) < minDR) {
      minDR = jetP4.DeltaR(muonP4);
      closestJetP4 = jetP4;
    }
  }

  nt.probe_minDR = minDR;
  nt.probe_ptRel_minDR =
      muonP4.P() *
      (closestJetP4.Vect().Cross(muonP4.Vect()).Mag() /
       closestJetP4.P() /
       muonP4.P());
}


template <typename JET, typename MUON>
inline void FindJetProbePair(const std::vector<JET> &jets,
                             const MUON &mu,
                             NtupleContent &nt) {
  FindJetProbePairImpl(jets, mu, nt);
}

template <typename JET, typename MUON>
inline void FindJetProbePair(const std::vector<JET> &jets,
                             const MUON &mu,
                             StandAloneNtupleContent &nt) {
  FindJetProbePairImpl(jets, mu, nt);
}

#endif