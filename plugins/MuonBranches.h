//
// Original Author:
//                george karathanasis
//         Created:  Thu, 23 Mar 2019 17:40:23 GMT
//
// filling functions for aod and miniaod tag/probe

#ifndef MuonAnalysis_MuonAnalyzer_MuonBranches
#define MuonAnalysis_MuonAnalyzer_MuonBranches

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonSimInfo.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuonSetup.h"

#include <type_traits>
#include "NtupleContent.h"
#include "helper.h"

template <typename MUON, typename TRK>
inline void FillTagBranches(const MUON &muon,
                            const std::vector<TRK> &tracks,
                            NtupleContent &nt,
                            const reco::Vertex &vertex) {
  nt.tag_pt = muon.pt();
  nt.tag_eta = muon.eta();
  nt.tag_phi = muon.phi();
  nt.tag_charge = muon.charge();
  nt.tag_pterr = muon.innerTrack()->ptError() / muon.innerTrack()->pt();
  nt.tag_dxy = muon.innerTrack()->dxy(reco::TrackBase::Point(nt.pv_x, nt.pv_y, nt.pv_z));
  nt.tag_dz = muon.innerTrack()->dz(reco::TrackBase::Point(nt.pv_x, nt.pv_y, nt.pv_z));
  nt.tag_isPF = muon.isPFMuon();
  nt.tag_isSA = muon.isStandAloneMuon();
  nt.tag_isTracker = muon.isTrackerMuon();
  nt.tag_isGlobal = muon.isGlobalMuon();
  // Use selectors instead of 'muon.passed' method which is only introduced in CMSSW_9_4_X
  nt.tag_isLoose = muon::isLooseMuon(muon);
  nt.tag_isMedium = muon::isMediumMuon(muon);
  nt.tag_isTight = muon::isTightMuon(muon, vertex);
  nt.tag_isSoft = muon::isSoftMuon(muon, vertex, false);
  nt.tag_isHighPt = muon::isHighPtMuon(muon, vertex);
  float iso04 = (TrackerEnergy04<TRK>(muon.eta(), muon.phi(), tracks) - muon.pt()) / muon.pt();
  nt.tag_relIso04 = (iso04 > 0) ? iso04 : 0;
  /*nt.tag_iso03_sumPt = muon.isolationR03().sumPt;
  nt.tag_pfIso03_charged = muon.pfIsolationR03().sumChargedHadronPt;
  nt.tag_pfIso03_neutral = muon.pfIsolationR03().sumNeutralHadronEt;
  nt.tag_pfIso03_photon = muon.pfIsolationR03().sumPhotonEt;
  nt.tag_pfIso03_sumPU = muon.pfIsolationR03().sumPUPt;*/
  nt.tag_pfIso04_charged = muon.pfIsolationR04().sumChargedHadronPt;
  nt.tag_pfIso04_neutral = muon.pfIsolationR04().sumNeutralHadronEt;
  nt.tag_pfIso04_photon = muon.pfIsolationR04().sumPhotonEt;
  nt.tag_pfIso04_sumPU = muon.pfIsolationR04().sumPUPt;
  if (muon.tunePMuonBestTrack().isNonnull()) {
    nt.tag_tuneP_pt = muon.tunePMuonBestTrack()->pt();
    nt.tag_tuneP_pterr = muon.tunePMuonBestTrack()->ptError();
  } else {
    nt.tag_tuneP_pt = -99.;
    nt.tag_tuneP_pterr = -99.;
  }
  int nsegments = 0;
  for (auto &chamber : muon.matches()) {
    if (chamber.id.det() != DetId::Muon)
      continue;
    if (chamber.id.subdetId() != MuonSubdetId::DT && chamber.id.subdetId() != MuonSubdetId::CSC)
      continue;
    nsegments += chamber.segmentMatches.size();
  }
  nt.tag_nsegments = nsegments;
}

template <typename MUON, typename TRK>
inline void FillProbeBranches(
    const MUON &mu, const std::vector<TRK> &tracks, NtupleContent &nt, bool success, const reco::Vertex &vertex) {
  nt.probe_pt = mu.pt();
  nt.probe_eta = mu.eta();
  nt.probe_phi = mu.phi();
  nt.probe_charge = mu.charge();
  float iso04 = (TrackerEnergy04<TRK>(mu.eta(), mu.phi(), tracks) - mu.pt()) / mu.pt();
  nt.probe_relIso04 = (iso04 > 0) ? iso04 : 0;
  // success --> muon obj and track match in dR
  if (success) {
    // Use selectors instead of 'mu.passed' method which is only introduced in CMSSW_9_4_X
    nt.probe_isLoose = muon::isLooseMuon(mu);
    nt.probe_isMedium = muon::isMediumMuon(mu);
    nt.probe_isTight = muon::isTightMuon(mu, vertex);
    nt.probe_isSoft = muon::isSoftMuon(mu, vertex, false);
    nt.probe_isHighPt = muon::isHighPtMuon(mu, vertex);
    nt.probe_isArbitratedTracker = muon::isGoodMuon(mu, muon::TrackerMuonArbitrated);
    nt.probe_isPF = mu.isPFMuon();
    nt.probe_isSA = mu.isStandAloneMuon();
    nt.probe_isTracker = mu.isTrackerMuon();
    nt.probe_isGlobal = mu.isGlobalMuon();
    /*nt.probe_iso03_sumPt = mu.isolationR03().sumPt;
    nt.probe_pfIso03_charged = mu.pfIsolationR03().sumChargedHadronPt;
    nt.probe_pfIso03_neutral = mu.pfIsolationR03().sumNeutralHadronEt;
    nt.probe_pfIso03_photon = mu.pfIsolationR03().sumPhotonEt;
    nt.probe_pfIso03_sumPU = mu.pfIsolationR03().sumPUPt;
    nt.probe_pfIso04_charged = mu.pfIsolationR04().sumChargedHadronPt;*/
    nt.probe_pfIso04_neutral = mu.pfIsolationR04().sumNeutralHadronEt;
    nt.probe_pfIso04_photon = mu.pfIsolationR04().sumPhotonEt;
    nt.probe_pfIso04_sumPU = mu.pfIsolationR04().sumPUPt;
    nt.probe_matchedStations = mu.numberOfMatchedStations();
    nt.probe_expectedMatchedStations = mu.expectedNnumberOfMatchedStations();
    nt.probe_RPCLayers = mu.numberOfMatchedRPCLayers();
    nt.probe_stationMask = mu.stationMask();
    nt.probe_nShowers = mu.numberOfShowers();
    if (mu.globalTrack().isNonnull()) {
      nt.probe_muonHits = mu.globalTrack()->hitPattern().numberOfValidMuonHits();
      nt.probe_trkChi2 = mu.globalTrack()->normalizedChi2();
    } else if (mu.innerTrack().isNonnull() && mu.innerTrack().isAvailable()) {
      nt.probe_trkChi2 = mu.innerTrack()->normalizedChi2();
      nt.probe_muonHits = mu.innerTrack()->hitPattern().numberOfValidMuonHits();
    } else {
      nt.probe_trkChi2 = -99;
      nt.probe_muonHits = -99;
    }
    if (mu.innerTrack().isNonnull() && mu.innerTrack().isAvailable()) {
      nt.probe_validFraction = mu.innerTrack()->validFraction();
      nt.probe_trackerLayers = mu.innerTrack()->hitPattern().trackerLayersWithMeasurement();
      nt.probe_pixelLayers = mu.innerTrack()->hitPattern().pixelLayersWithMeasurement();
      nt.probe_pterr = mu.innerTrack()->ptError() / mu.innerTrack()->pt();
      nt.probe_dxy = mu.innerTrack()->dxy(reco::TrackBase::Point(nt.pv_x, nt.pv_y, nt.pv_z));
      nt.probe_dz = mu.innerTrack()->dz(reco::TrackBase::Point(nt.pv_x, nt.pv_y, nt.pv_z));
      nt.probe_pixelHits = mu.innerTrack()->hitPattern().numberOfValidPixelHits();
    } else {
      nt.probe_validFraction = -99;
      nt.probe_trackerLayers = -99;
      nt.probe_pixelLayers = -99;
      nt.probe_pterr = -99;
      nt.probe_dxy = -99;
      nt.probe_dz = -99;
      nt.probe_pixelHits = -99;
    }
    if (mu.outerTrack().isNonnull() && mu.outerTrack().isAvailable()) {
      nt.probe_muonStations = mu.outerTrack()->hitPattern().muonStationsWithValidHits();
      nt.probe_DTHits = mu.outerTrack()->hitPattern().numberOfValidMuonDTHits();
      nt.probe_CSCHits = mu.outerTrack()->hitPattern().numberOfValidMuonCSCHits();
    } else {
      nt.probe_muonStations = -99;
      nt.probe_muonHits = -99;
      nt.probe_DTHits = -99;
      nt.probe_CSCHits = -99;
    }
    if (mu.tunePMuonBestTrack().isNonnull()) {
      nt.probe_tuneP_pt = mu.tunePMuonBestTrack()->pt();
      nt.probe_tuneP_pterr = mu.tunePMuonBestTrack()->ptError();
      nt.probe_tuneP_muonHits = mu.tunePMuonBestTrack()->hitPattern().numberOfValidMuonHits();
    } else {
      nt.probe_tuneP_pt = -99;
      nt.probe_tuneP_pterr = -99;
      nt.probe_tuneP_muonHits = -99;
    }
    nt.probe_positionChi2 = mu.combinedQuality().chi2LocalPosition;
    nt.probe_trkKink = mu.combinedQuality().trkKink;
    nt.probe_segmentCompatibility = muon::segmentCompatibility(mu);
    nt.probe_isMuMatched = true;

    int nsegments = 0;
    for (auto &chamber : mu.matches()) {
      if (chamber.id.det() != DetId::Muon)
        continue;
      if (chamber.id.subdetId() != MuonSubdetId::DT && chamber.id.subdetId() != MuonSubdetId::CSC)
        continue;
      nsegments += chamber.segmentMatches.size();
    }
    nt.probe_nsegments = nsegments;
  }
  // no successs (no match)
  else {
    nt.probe_isLoose = false;
    nt.probe_isMedium = false;
    nt.probe_isTight = false;
    nt.probe_isSoft = false;
    nt.probe_isHighPt = false;
    nt.probe_isMuMatched = false;
    nt.probe_isPF = false;
    nt.probe_isSA = false;
    nt.probe_isTracker = false;
    nt.probe_isGlobal = false;
    nt.probe_validFraction = -99;
    nt.probe_trkChi2 = -99;
    nt.probe_positionChi2 = -99;
    nt.probe_trkKink = -99;
    nt.probe_trackerLayers = -99;
    nt.probe_pixelLayers = -99;
    nt.probe_dxy = -99;
    nt.probe_dz = -99;
    nt.probe_muonStations = -99;
    nt.probe_muonHits = -99;
    nt.probe_DTHits = -99;
    nt.probe_CSCHits = -99;
    nt.probe_pterr = -99;
    /*nt.probe_iso03_sumPt = -99;
    nt.probe_pfIso03_charged = -99;
    nt.probe_pfIso03_neutral = -99;
    nt.probe_pfIso03_photon = -99;
    nt.probe_pfIso03_sumPU = -99;*/
    nt.probe_pfIso04_charged = -99;
    nt.probe_pfIso04_neutral = -99;
    nt.probe_pfIso04_photon = -99;
    nt.probe_pfIso04_sumPU = -99;
    nt.probe_pixelHits = -99;
    nt.probe_matchedStations = -99;
    nt.probe_expectedMatchedStations = -99;
    nt.probe_RPCLayers = -99;
    nt.probe_stationMask = 0;
    nt.probe_nShowers = -99;
    nt.probe_tuneP_pt = -99;
    nt.probe_tuneP_pterr = -99;
    nt.probe_tuneP_muonHits = -99;
    nt.probe_nsegments = -99;
    nt.l1pt = -99;
    nt.l1q = -99;
    nt.l1dr = -99;
    nt.l1ptByQ = -99;
    nt.l1qByQ = -99;
    nt.l1drByQ = -99;
  }
}

template <typename MUON>
inline void FillProbeBranchesSelector(const MUON &mu,
                                      NtupleContent &nt,
                                      const std::vector<unsigned> selectorBits,
                                      bool success) {
  for (unsigned int ibit = 0; ibit < selectorBits.size(); ++ibit) {
    if (success && mu.selectors() != 0) {
      unsigned bit = selectorBits.at(ibit);
      nt.probe_selectors[ibit] = mu.passed(1UL << bit);
    } else {
      nt.probe_selectors[ibit] = false;
    }
  }
}

template <typename TRK>
inline void FillProbeBranchesdSA(const TRK &trk, NtupleContent &nt, bool passdSA) {
  nt.probe_isdSA = passdSA;

  nt.probe_dsa_pt = trk.pt();
  nt.probe_dsa_eta = trk.eta();
  nt.probe_dsa_phi = trk.phi();
  nt.probe_dsa_charge = trk.charge();

  if (passdSA) {
    nt.probe_dsa_outerEta = trk.outerEta();
    nt.probe_dsa_outerPhi = trk.outerPhi();
    nt.probe_dsa_dxy = trk.dxy(reco::TrackBase::Point(nt.pv_x, nt.pv_y, nt.pv_z));
    nt.probe_dsa_dz = trk.dz(reco::TrackBase::Point(nt.pv_x, nt.pv_y, nt.pv_z));
    nt.probe_dsa_muonStations = trk.hitPattern().muonStationsWithValidHits();
    nt.probe_dsa_muonHits = trk.hitPattern().numberOfValidMuonHits();
    nt.probe_dsa_DTHits = trk.hitPattern().numberOfValidMuonDTHits();
    nt.probe_dsa_CSCHits = trk.hitPattern().numberOfValidMuonCSCHits();
    nt.probe_dsa_pterr = trk.ptError() / trk.pt();
    nt.probe_dsa_trkChi2 = trk.normalizedChi2();

    // [Adapted from displaced dimuon analysis]
    // Number of DT+CSC segments
    unsigned int nsegments = 0;
    for (auto &hit : trk.recHits()) {
      if (!hit->isValid())
        continue;
      DetId id = hit->geographicalId();
      if (id.det() != DetId::Muon)
        continue;
      if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC)
        nsegments++;
    }
    nt.probe_dsa_nsegments = nsegments;
  } else {
    nt.probe_dsa_outerEta = -99;
    nt.probe_dsa_outerPhi = -99;
    nt.probe_dsa_dxy = -99;
    nt.probe_dsa_dz = -99;
    nt.probe_dsa_muonStations = -99;
    nt.probe_dsa_muonHits = -99;
    nt.probe_dsa_DTHits = -99;
    nt.probe_dsa_CSCHits = -99;
    nt.probe_dsa_pterr = -99;
    nt.probe_dsa_trkChi2 = -99;
    nt.probe_dsa_nsegments = -99;
  }
}

template <typename TRK>
inline void FillTagBranchesdSA(const TRK &trk, NtupleContent &nt, bool passdSA) {
  nt.tag_isdSA = passdSA;

  nt.tag_dsa_pt = trk.pt();
  nt.tag_dsa_eta = trk.eta();
  nt.tag_dsa_phi = trk.phi();
  nt.tag_dsa_charge = trk.charge();

  if (passdSA) {
    nt.tag_dsa_outerEta = trk.outerEta();
    nt.tag_dsa_outerPhi = trk.outerPhi();
    nt.tag_dsa_dxy = trk.dxy(reco::TrackBase::Point(nt.pv_x, nt.pv_y, nt.pv_z));
    nt.tag_dsa_dz = trk.dz(reco::TrackBase::Point(nt.pv_x, nt.pv_y, nt.pv_z));
    nt.tag_dsa_muonStations = trk.hitPattern().muonStationsWithValidHits();
    nt.tag_dsa_muonHits = trk.hitPattern().numberOfValidMuonHits();
    nt.tag_dsa_DTHits = trk.hitPattern().numberOfValidMuonDTHits();
    nt.tag_dsa_CSCHits = trk.hitPattern().numberOfValidMuonCSCHits();
    nt.tag_dsa_pterr = trk.ptError() / trk.pt();
    nt.tag_dsa_trkChi2 = trk.normalizedChi2();
    // [Adapted from displaced dimuon analysis]
    // Number of DT+CSC segments
    unsigned int nsegments = 0;
    for (auto &hit : trk.recHits()) {
      if (!hit->isValid())
        continue;
      DetId id = hit->geographicalId();
      if (id.det() != DetId::Muon)
        continue;
      if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC)
        nsegments++;
    }
    nt.tag_dsa_nsegments = nsegments;
  } else {
    nt.tag_dsa_outerEta = -99;
    nt.tag_dsa_outerPhi = -99;
    nt.tag_dsa_dxy = -99;
    nt.tag_dsa_dz = -99;
    nt.tag_dsa_muonStations = -99;
    nt.tag_dsa_muonHits = -99;
    nt.tag_dsa_DTHits = -99;
    nt.tag_dsa_CSCHits = -99;
    nt.tag_dsa_pterr = -99;
    nt.tag_dsa_trkChi2 = -99;
    nt.tag_dsa_nsegments = -99;
  }
}

template <typename TRK>
inline void FillProbeBranchesdgl(const TRK &trk, NtupleContent &nt, bool passdgl) {
  nt.probe_isdGlobal = passdgl;

  nt.probe_dgl_pt = trk.pt();
  nt.probe_dgl_eta = trk.eta();
  nt.probe_dgl_phi = trk.phi();
  nt.probe_dgl_charge = trk.charge();

  if (passdgl) {
    nt.probe_dgl_dxy = trk.dxy(reco::TrackBase::Point(nt.pv_x, nt.pv_y, nt.pv_z));
    nt.probe_dgl_dz = trk.dz(reco::TrackBase::Point(nt.pv_x, nt.pv_y, nt.pv_z));
    nt.probe_dgl_muonStations = trk.hitPattern().muonStationsWithValidHits();
    nt.probe_dgl_muonHits = trk.hitPattern().numberOfValidMuonHits();
    nt.probe_dgl_outerTrackerHits = trk.hitPattern().numberOfValidStripHits();
    nt.probe_dgl_trackerHits = trk.hitPattern().numberOfValidTrackerHits();
    nt.probe_dgl_totalHits = trk.hitPattern().numberOfValidHits();
    nt.probe_dgl_DTHits = trk.hitPattern().numberOfValidMuonDTHits();
    nt.probe_dgl_CSCHits = trk.hitPattern().numberOfValidMuonCSCHits();
    nt.probe_dgl_pterr = trk.ptError() / trk.pt();
    nt.probe_dgl_trkChi2 = trk.normalizedChi2();
    // [Adapted from displaced dimuon analysis]
    // Number of DT+CSC segments
    unsigned int nsegments = 0;
    for (auto & hit : trk.recHits()) {
      if (!hit->isValid()) continue;
      DetId id = hit->geographicalId();
      if (id.det() != DetId::Muon) continue;
      if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC)
        nsegments++;
    }
    nt.probe_dgl_nsegments = nsegments;
  } else {
    nt.probe_dgl_dxy = -99;
    nt.probe_dgl_dz = -99;
    nt.probe_dgl_muonStations = -99;
    nt.probe_dgl_muonHits = -99;
    nt.probe_dgl_DTHits = -99;
    nt.probe_dgl_CSCHits = -99;
    nt.probe_dgl_pterr = -99;
    nt.probe_dgl_trkChi2 = -99;
    nt.probe_dgl_nsegments = -99;
  }
}

template <typename TRK>
inline void FillProbeBranchesCosmic(const TRK &trk, NtupleContent &nt, bool passcosmic) {
  nt.probe_isCosmic = passcosmic;
}
template <typename MUOTT, typename TRKTT>
inline void FillPairBranches(const MUOTT &muon, const TRKTT &trk, NtupleContent &nt, PropagateToMuon &prop1_) {
  math::PtEtaPhiMLorentzVector mu1(muon.first.pt(), muon.first.eta(), muon.first.phi(), MU_MASS);
  math::PtEtaPhiMLorentzVector mu2(trk.first.pt(), trk.first.eta(), trk.first.phi(), MU_MASS);
  nt.pair_pt = (mu1 + mu2).pt();
  nt.pair_mass = (mu1 + mu2).mass();
  nt.pair_eta = (mu1 + mu2).eta();
  nt.pair_phi = (mu1 + mu2).phi();
  nt.pair_dz = muon.first.vz() - trk.first.vz();
  nt.pair_dR = deltaR(muon.first.eta(), muon.first.phi(), trk.first.eta(), trk.first.phi());

  FreeTrajectoryState trajectory_state_muon = muon.second.impactPointTSCP().theState();
  FreeTrajectoryState trajectory_state_trk = trk.second.impactPointTSCP().theState();

  TrajectoryStateOnSurface prop1_M1 = prop1_.extrapolate(trajectory_state_muon);
  TrajectoryStateOnSurface prop2_M1 = prop1_.extrapolate(trajectory_state_trk);

  if (prop1_M1.isValid() && prop2_M1.isValid()) {
    float dphiM1 = deltaPhi<float>(prop1_M1.globalPosition().phi(), prop2_M1.globalPosition().phi());
    nt.pair_drM1 = hypot(dphiM1, std::abs(prop1_M1.globalPosition().eta() - prop2_M1.globalPosition().eta()));
  }
}

template <typename MUO, typename TRK>
inline void FillTunePPairBranches(const MUO &muon, const TRK &trk, NtupleContent &nt) {
  math::PtEtaPhiMLorentzVector mu1(muon.pt(), muon.eta(), muon.phi(), MU_MASS);
  math::PtEtaPhiMLorentzVector mu2(trk.pt(), trk.eta(), trk.phi(), MU_MASS);
  nt.pair_tuneP_pt = (mu1 + mu2).pt();
  nt.pair_tuneP_mass = (mu1 + mu2).mass();
  nt.pair_tuneP_eta = (mu1 + mu2).eta();
  nt.pair_tuneP_phi = (mu1 + mu2).phi();
  nt.pair_tuneP_dz = muon.vz() - trk.vz();
  nt.pair_tuneP_dR = deltaR(muon.eta(), muon.phi(), trk.eta(), trk.phi());
}

inline void FillTunePPairBranchesDummy(NtupleContent &nt) {
  nt.pair_tuneP_pt = -99;
  nt.pair_tuneP_mass = -99;
  nt.pair_tuneP_eta = -99;
  nt.pair_tuneP_phi = -99;
  nt.pair_tuneP_dz = -99;
  nt.pair_tuneP_dR = -99;
  nt.pair_tuneP_fit_mass = -99;
  nt.pair_tuneP_svprob = -99;
  nt.pair_tuneP_normalchi2 = -99;
}

inline void FillSimMatchingBranches(const pat::Muon &mu, NtupleContent &nt, bool isTag) {
  if (isTag) {
    nt.tag_simType = mu.simType();
    nt.tag_simExtType = mu.simExtType();
    nt.tag_simFlavour = mu.simFlavour();
    nt.tag_simHeaviestMotherFlavour = mu.simHeaviestMotherFlavour();
    nt.tag_simPdgId = mu.simPdgId();
    nt.tag_simMotherPdgId = mu.simMotherPdgId();
    nt.tag_simBX = mu.simBX();
    nt.tag_simProdRho = mu.simProdRho();
    nt.tag_simProdZ = mu.simProdZ();
    nt.tag_simPt = mu.simPt();
    nt.tag_simEta = mu.simEta();
    nt.tag_simPhi = mu.simPhi();
  } else {
    nt.probe_simType = mu.simType();
    nt.probe_simExtType = mu.simExtType();
    nt.probe_simFlavour = mu.simFlavour();
    nt.probe_simHeaviestMotherFlavour = mu.simHeaviestMotherFlavour();
    nt.probe_simPdgId = mu.simPdgId();
    nt.probe_simMotherPdgId = mu.simMotherPdgId();
    nt.probe_simBX = mu.simBX();
    nt.probe_simProdRho = mu.simProdRho();
    nt.probe_simProdZ = mu.simProdZ();
    nt.probe_simPt = mu.simPt();
    nt.probe_simEta = mu.simEta();
    nt.probe_simPhi = mu.simPhi();
  }
}

inline void FillSimMatchingBranchesAOD(const reco::MuonSimInfo &msi, NtupleContent &nt, bool isTag) {
  if (isTag) {
    nt.tag_simType = msi.primaryClass;
    nt.tag_simExtType = msi.extendedClass;
    nt.tag_simFlavour = msi.flavour;
    nt.tag_simHeaviestMotherFlavour = msi.heaviestMotherFlavour;
    nt.tag_simPdgId = msi.pdgId;
    nt.tag_simMotherPdgId = msi.motherPdgId;
    nt.tag_simBX = msi.tpBX;
    nt.tag_simProdRho = msi.vertex.Rho();
    nt.tag_simProdZ = msi.vertex.Z();
    nt.tag_simPt = msi.p4.pt();
    nt.tag_simEta = msi.p4.eta();
    nt.tag_simPhi = msi.p4.phi();
  } else {
    nt.probe_simType = msi.primaryClass;
    nt.probe_simExtType = msi.extendedClass;
    nt.probe_simFlavour = msi.flavour;
    nt.probe_simHeaviestMotherFlavour = msi.heaviestMotherFlavour;
    nt.probe_simPdgId = msi.pdgId;
    nt.probe_simMotherPdgId = msi.motherPdgId;
    nt.probe_simBX = msi.tpBX;
    nt.probe_simProdRho = msi.vertex.Rho();
    nt.probe_simProdZ = msi.vertex.Z();
    nt.probe_simPt = msi.p4.pt();
    nt.probe_simEta = msi.p4.eta();
    nt.probe_simPhi = msi.p4.phi();
  }
}

inline void FillSimMatchingBranchesDummy(NtupleContent &nt, bool isTag) {
  if (isTag) {
    nt.tag_simType = -99;
    nt.tag_simExtType = -99;
    nt.tag_simFlavour = -99;
    nt.tag_simHeaviestMotherFlavour = -99;
    nt.tag_simPdgId = -99;
    nt.tag_simMotherPdgId = -99;
    nt.tag_simBX = -99;
    nt.tag_simProdRho = -99;
    nt.tag_simProdZ = -99;
    nt.tag_simPt = -99;
    nt.tag_simEta = -99;
    nt.tag_simPhi = -99;
  } else {
    nt.probe_simType = -99;
    nt.probe_simExtType = -99;
    nt.probe_simFlavour = -99;
    nt.probe_simHeaviestMotherFlavour = -99;
    nt.probe_simPdgId = -99;
    nt.probe_simMotherPdgId = -99;
    nt.probe_simBX = -99;
    nt.probe_simProdRho = -99;
    nt.probe_simProdZ = -99;
    nt.probe_simPt = -99;
    nt.probe_simEta = -99;
    nt.probe_simPhi = -99;
  }
}

#endif
