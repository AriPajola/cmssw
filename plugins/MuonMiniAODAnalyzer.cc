// Package:    HLTAnalysis/TriggerAnalyzer
// Class:      TriggerAnalyzer
//
/**\class TriggerAnalyzer TriggerAnalyzer.cc
 HLTAnalysis/TriggerAnalyzer/plugins/TriggerAnalyzer.cc

 Description: Ntuplizer for miniAOD files
*/
//
// Original Author:
//                george karathanasis
//         Created:  Thu, 23 Mar 2017 17:40:23 GMT
//
// Modified:
//                Minseok Oh (Feb. 2021)
//
//

// system include files
#include <iostream>
#include <memory>
#include <random>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/DataRecord/interface/JetResolutionRcd.h"
#include "CondFormats/DataRecord/interface/JetResolutionScaleFactorRcd.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "HLTrigger/HLTcore/interface/defaultModuleLabel.h"

#include <iostream>
#include <string>
#include <vector>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "KlFitter.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MuonBranches.h"
#include "MuonGenAnalyzer.h"
#include "NtupleContent.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "helper.h"
#include "MuonMiniIsolation.h"
#include "JetsBranches.h"

using namespace std;
// using namespace edm;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
class MuonMiniAODAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  typedef std::vector<std::pair<pat::Muon, reco::TransientTrack>> PatMuonAndTransientTrkCollection;
  typedef std::pair<pat::Muon, reco::TransientTrack> PatMuonAndTransientTrk;
  typedef std::pair<pat::PackedCandidate, reco::TransientTrack> PatPackedCandAndTransientTrk;
  explicit MuonMiniAODAnalyzer(const edm::ParameterSet&);
  ~MuonMiniAODAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  bool HLTaccept(const edm::Event&, NtupleContent&, std::vector<std::string>&);
  void embedTriggerMatching(const edm::Event&,
                            edm::Handle<edm::TriggerResults>&,
                            const pat::Muon&,
                            NtupleContent&,
                            std::vector<std::string>&,
                            bool);
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupSummaryToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
  edm::EDGetToken muonsToken_;
  edm::EDGetTokenT<edm::View<reco::Muon>> muonsViewToken_;
  edm::EDGetToken PFCands_;
  edm::EDGetToken LostTracks_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneMatch> l1MatchesToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> l1MatchesQualityToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> l1MatchesDeltaRToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneMatch> l1MatchesByQToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> l1MatchesByQQualityToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> l1MatchesByQDeltaRToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> genToken_;
  edm::EDGetTokenT<double> rhoJetsNC_;
  edm::EDGetToken jetsToken_;
  edm::EDGetToken genJetsToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magfieldToken_;
  JME::JetResolution::Token jetResolutionToken_;
  JME::JetResolutionScaleFactor::Token jetResolutionScaleFactorToken_;
  std::vector<std::string> HLTPaths_;
  std::vector<std::string> tagFilters_;
  std::vector<std::string> probeFilters_;
  std::vector<std::string> probeSelectorNames_;
  std::vector<unsigned> probeSelectorBits_;

  const unsigned int tagQual_;
  const StringCutObjectSelector<pat::Muon> tagSelection_;  // kinematic cuts for tag
  const bool HighPurity_;
  const StringCutObjectSelector<pat::PackedCandidate> probeSelection_;  // kinematic cuts for probe
  const bool muonOnly_;
  const StringCutObjectSelector<pat::Muon> probeMuonSelection_;
  const double pairMassMin_;
  const double pairMassMax_;
  const double pairDz_;
  const bool RequireVtxCreation_;  // if true skip pairs that do not create
                                   // that do not have a vertex
  const double minSVtxProb_;       // min probability of a vertex to be kept. If < 0 inactive
  const double maxdz_trk_mu_;
  const double maxpt_relative_dif_trk_mu_;
  const double maxdr_trk_mu_;
  const unsigned momPdgId_;
  const double genRecoDrMatch_;
  PropagateToMuon prop1_;
  PropagateToMuonSetup propSetup1_;

  edm::Service<TFileService> fs;
  TTree* t1;
  NtupleContent nt;

  std::mt19937 m_random_generator = std::mt19937(37428479);
  const bool isMC_, includeJets_;
  const std::string era_;

  // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

MuonMiniAODAnalyzer::MuonMiniAODAnalyzer(const edm::ParameterSet& iConfig)
    : genEventInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))),
      rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("Rho"))),
      pileupSummaryToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupInfo"))),
      beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
      vtxToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
      muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      muonsViewToken_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      PFCands_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCands"))),
      LostTracks_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("lostTracks"))),
      trgresultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
      l1MatchesToken_(consumes<pat::TriggerObjectStandAloneMatch>(iConfig.getParameter<edm::InputTag>("l1Matches"))),
      l1MatchesQualityToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("l1MatchesQuality"))),
      l1MatchesDeltaRToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("l1MatchesDeltaR"))),
      l1MatchesByQToken_(
          consumes<pat::TriggerObjectStandAloneMatch>(iConfig.getParameter<edm::InputTag>("l1MatchesByQ"))),
      l1MatchesByQQualityToken_(
          consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("l1MatchesByQQuality"))),
      l1MatchesByQDeltaRToken_(
          consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("l1MatchesByQDeltaR"))),
      genToken_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("gen"))),
      rhoJetsNC_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoJetsNC"))),
      jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      genJetsToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
      magfieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
      jetResolutionToken_(esConsumes(edm::ESInputTag("", iConfig.getParameter<std::string>("jetResType")))),
      jetResolutionScaleFactorToken_(esConsumes(edm::ESInputTag("", iConfig.getParameter<std::string>("jetResSFType")))),
      HLTPaths_(iConfig.getParameter<std::vector<std::string>>("triggerPaths")),
      tagFilters_(iConfig.getParameter<std::vector<std::string>>("tagFilters")),
      probeFilters_(iConfig.getParameter<std::vector<std::string>>("probeFilters")),
      probeSelectorNames_(iConfig.getParameter<std::vector<std::string>>("probeSelectorNames")),
      probeSelectorBits_(iConfig.getParameter<std::vector<unsigned>>("probeSelectorBits")),
      tagQual_(iConfig.getParameter<unsigned>("tagQuality")),
      tagSelection_(iConfig.getParameter<std::string>("tagSelection")),
      HighPurity_(iConfig.getParameter<bool>("ProbeHPurity")),
      probeSelection_(iConfig.getParameter<std::string>("probeSelection")),
      muonOnly_(iConfig.getParameter<bool>("muonOnly")),
      probeMuonSelection_(iConfig.getParameter<std::string>("probeMuonSelection")),
      pairMassMin_(iConfig.getParameter<double>("pairMassMin")),
      pairMassMax_(iConfig.getParameter<double>("pairMassMax")),
      pairDz_(iConfig.getParameter<double>("pairDz")),
      RequireVtxCreation_(iConfig.getParameter<bool>("RequireVtxCreation")),
      minSVtxProb_(iConfig.getParameter<double>("minSVtxProb")),
      maxdz_trk_mu_(iConfig.getParameter<double>("maxDzProbeTrkMuon")),
      maxpt_relative_dif_trk_mu_(iConfig.getParameter<double>("maxRelPtProbeTrkMuon")),
      maxdr_trk_mu_(iConfig.getParameter<double>("maxDRProbeTrkMuon")),
      momPdgId_(iConfig.getParameter<unsigned>("momPdgId")),
      genRecoDrMatch_(iConfig.getParameter<double>("genRecoDrMatch")),
      //magfieldToken_(iConfig.getParameter<edm::ParameterSet>("propM1")),
      //prop1_(iConfig.getParameter<edm::ParameterSet>("propM1")),
      prop1_(),
      propSetup1_(iConfig, consumesCollector()),
      isMC_(iConfig.getParameter<bool>("isMC")),
      includeJets_(iConfig.getParameter<bool>("includeJets")),
      era_(iConfig.getParameter<std::string>("era")) {
  //  edm::ParameterSet
  //  runParameters=iConfig.getParameter<edm::ParameterSet>("RunParameters");

  if (probeSelectorNames_.size() != probeSelectorBits_.size()) {
    throw cms::Exception("ParameterError")
        << "length of probeSelectorNames and probeSelectorBits should be identical\n";
  }
}

MuonMiniAODAnalyzer::~MuonMiniAODAnalyzer() {
  // cout << "total " << trg_counter << " fires " << fire_counter << " l3"
  // << l3_counter << endl; do anything here that needs to be done at
  // desctruction time
}

//
// member functions
bool MuonMiniAODAnalyzer::HLTaccept(const edm::Event& iEvent, NtupleContent& nt, std::vector<std::string>& HLTPaths) {
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsToken_, trigResults);
  edm::TriggerNames trigName;
  trigName = iEvent.triggerNames(*trigResults);
  bool EvtFire = false;
  unsigned int ipath = 0;
  for (auto path : HLTPaths) {
    bool TrgFire = false;
    for (unsigned int itrg = 0; itrg < trigResults->size(); ++itrg) {
      TString TrigPath = trigName.triggerName(itrg);
      if (!trigResults->accept(itrg))
        continue;
      if (!TrigPath.Contains(path))
        continue;
      EvtFire = true;
      TrgFire = true;
    }
    nt.trigger[ipath] = TrgFire;
    ipath++;
  }
  return EvtFire;
}

void MuonMiniAODAnalyzer::embedTriggerMatching(const edm::Event& iEvent,
                                               edm::Handle<edm::TriggerResults>& trigResults,
                                               const pat::Muon& mu,
                                               NtupleContent& nt,
                                               std::vector<std::string>& Triggers,
                                               bool isTag) {
  for (const auto& trg : Triggers) {
    TString trg_tstr = TString(trg);
    bool matched = false;
    float matched_pt = -99;
    float matched_eta = -99;
    float matched_phi = -99;
    float matched_dr = 99;
    for (auto trigobj : mu.triggerObjectMatches()) {
      trigobj.unpackNamesAndLabels(iEvent, *trigResults);
      float dR_tmp = deltaR(mu.eta(), mu.phi(), trigobj.eta(), trigobj.phi());

      // check path names
      if (trg_tstr.Contains("HLT_")) {
        for (auto path : trigobj.pathNames(true, true)) {
          TString path_tstr = TString(path);
          if (path_tstr.Contains(trg_tstr) && dR_tmp < matched_dr) {
            matched = true;
            matched_pt = trigobj.pt();
            matched_eta = trigobj.eta();
            matched_phi = trigobj.phi();
            matched_dr = dR_tmp;
          }
        }
      }
      // check filters
      else {
        for (auto filter : trigobj.filterLabels()) {
          TString filter_tstr = TString(filter);
          if (filter_tstr.Contains(trg_tstr) && dR_tmp < matched_dr) {
            matched = true;
            matched_pt = trigobj.pt();
            matched_eta = trigobj.eta();
            matched_phi = trigobj.phi();
            matched_dr = dR_tmp;
          }
        }
      }
    }

    if (isTag) {
      nt.tag_trg[&trg - &Triggers[0]] = matched;
      nt.tag_trg_pt[&trg - &Triggers[0]] = matched_pt;
      nt.tag_trg_eta[&trg - &Triggers[0]] = matched_eta;
      nt.tag_trg_phi[&trg - &Triggers[0]] = matched_phi;
      nt.tag_trg_dr[&trg - &Triggers[0]] = matched_dr;
    } else {
      nt.probe_trg[&trg - &Triggers[0]] = matched;
      nt.probe_trg_pt[&trg - &Triggers[0]] = matched_pt;
      nt.probe_trg_eta[&trg - &Triggers[0]] = matched_eta;
      nt.probe_trg_phi[&trg - &Triggers[0]] = matched_phi;
      nt.probe_trg_dr[&trg - &Triggers[0]] = matched_dr;
    }
  }

  return;
}
// ------------ method called for each event  ------------

void MuonMiniAODAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //prop1_.init(iSetup);
  prop1_ = propSetup1_.init(iSetup);

  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_, theBeamSpot);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);

  // Skip evts if there are no vertices
  if (vertices->empty())
    return;
  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonsToken_, muons);
  edm::Handle<edm::View<reco::Muon>> muonsView;
  iEvent.getByToken(muonsViewToken_, muonsView);
  edm::Handle<std::vector<pat::PackedCandidate>> pfcands;
  iEvent.getByToken(PFCands_, pfcands);
  edm::Handle<std::vector<pat::PackedCandidate>> lostTracks;
  iEvent.getByToken(LostTracks_, lostTracks);
  edm::ESHandle<MagneticField> bField;
  bField = iSetup.getHandle(magfieldToken_);
  //iSetup.get<IdealMagneticFieldRecord>().get(bField);
  edm::Handle<double> rhoJetsNC;
  iEvent.getByToken(rhoJetsNC_, rhoJetsNC);
  edm::Handle<std::vector<pat::Jet>> jets;
  iEvent.getByToken(jetsToken_, jets);
  //iSetup.get<IdealMagneticFieldRecord>().get(bField);
  edm::Handle<std::vector<reco::GenJet>> genJets;
  iEvent.getByToken(genJetsToken_, genJets);
  JME::JetResolution resolution;
  //resolution = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
  resolution = iSetup.getData(jetResolutionToken_);
  JME::JetResolutionScaleFactor resolution_sf;
  //resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, "AK4PFchs");
  resolution_sf = iSetup.getData(jetResolutionScaleFactorToken_);

  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsToken_, trigResults);
  edm::Handle<pat::TriggerObjectStandAloneMatch> l1Matches;
  iEvent.getByToken(l1MatchesToken_, l1Matches);
  edm::Handle<edm::ValueMap<int>> l1Qualities;
  iEvent.getByToken(l1MatchesQualityToken_, l1Qualities);
  edm::Handle<edm::ValueMap<float>> l1Drs;
  iEvent.getByToken(l1MatchesDeltaRToken_, l1Drs);
  edm::Handle<pat::TriggerObjectStandAloneMatch> l1MatchesByQ;
  iEvent.getByToken(l1MatchesByQToken_, l1MatchesByQ);
  edm::Handle<edm::ValueMap<int>> l1QualitiesByQ;
  iEvent.getByToken(l1MatchesByQQualityToken_, l1QualitiesByQ);
  edm::Handle<edm::ValueMap<float>> l1DrsByQ;
  iEvent.getByToken(l1MatchesByQDeltaRToken_, l1DrsByQ);

  // Information about run
  nt.ClearBranches();
  nt.run = iEvent.id().run();
  nt.ls = iEvent.luminosityBlock();
  nt.event = iEvent.id().event();
  nt.fromFullAOD = false;
  nt.BSpot_x = theBeamSpot->x0();
  nt.BSpot_y = theBeamSpot->y0();
  nt.BSpot_z = theBeamSpot->z0();
  nt.nvertices = vertices->size();

  // Gen weights
  if (!iEvent.isRealData()) {
    edm::Handle<GenEventInfoProduct> genEventInfoHandle;
    iEvent.getByToken(genEventInfoToken_, genEventInfoHandle);
    nt.genWeight = genEventInfoHandle->weight();
  } else {  // data
    nt.genWeight = 1.;
  }

  // Pileup information
  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_, rhoHandle);
  nt.Rho = *rhoHandle;

  float trueNumInteractions = -1;
  int puNumInteractions = -1;
  if (isMC_) {
    edm::Handle<std::vector<PileupSummaryInfo>> PupInfo;
    iEvent.getByToken(pileupSummaryToken_, PupInfo);

    for (auto PVI : *PupInfo) {
      int BX = PVI.getBunchCrossing();
      if (BX == 0) {
        trueNumInteractions = PVI.getTrueNumInteractions();
        puNumInteractions = PVI.getPU_NumInteractions();
        continue;
      }
    }
  }

  nt.trueNumInteractions = trueNumInteractions;
  nt.puNumInteractions = puNumInteractions;

  bool goodVtx = false;
  reco::Vertex const* pv;
  for (const reco::Vertex& vtx : *vertices) {
    if (vtx.isFake() || !vtx.isValid())
      continue;
    nt.pv_x = vtx.x();
    nt.pv_y = vtx.y();
    nt.pv_z = vtx.z();
    goodVtx = true;
    pv = &vtx;
    break;
  }
  if (!goodVtx)
    return;

  // check if path fired, if so save hlt muons
  if (!HLTaccept(iEvent, nt, HLTPaths_))
    return;
  //  HLTaccept(iEvent, nt.doublemu_trg, DoubleMuPaths_);

  // Gen information
  MuonGenAnalyzer genmu;
  std::vector<unsigned> matched_muon_idx;
  if (!iEvent.isRealData()) {
    genmu.SetInputs(iEvent, genToken_, momPdgId_);
    genmu.FillNtuple(nt);
    auto reco_match_genmu1 =
        MatchReco<pat::Muon>(*muons, nt.genmu1_eta, nt.genmu1_phi, nt.genmu1_charge, genRecoDrMatch_);
    auto reco_match_genmu2 =
        MatchReco<pat::Muon>(*muons, nt.genmu2_eta, nt.genmu2_phi, nt.genmu2_charge, genRecoDrMatch_);
    if (reco_match_genmu1.first)
      matched_muon_idx.push_back(reco_match_genmu1.second);
    if (reco_match_genmu2.first)
      matched_muon_idx.push_back(reco_match_genmu2.second);
  }

  // Find triggering muon
  std::vector<unsigned> tag_muon_map;  // idx of tag muon in muons
  PatMuonAndTransientTrkCollection tag_muon_ttrack;
  std::vector<bool> genmatched_tag;
  for (const pat::Muon& mu : *muons) {
    if (mu.selectors() != 0) {  // Only 9_4_X and later have selector bits
      if (!mu.passed(pow(2, tagQual_)))
        continue;
    } else {  // For 2016, assume loose ID on the tag (can be tightened at spark level)
      if (!muon::isLooseMuon(mu))
        continue;
    }
    bool fired = false;
    for (const auto& path : HLTPaths_) {
      const std::string triggerPattern = path + "*";
      if (!mu.triggered(triggerPattern.c_str()))
        continue;
      fired = true;
      break;
    }
    if (!fired)
      continue;
    if (!tagSelection_(mu))
      continue;
    tag_muon_ttrack.emplace_back(std::make_pair(mu, reco::TransientTrack(*mu.bestTrack(), &(*bField))));
    tag_muon_map.push_back(&mu - &muons->at(0));
    if (std::find(matched_muon_idx.begin(), matched_muon_idx.end(), &mu - &muons->at(0)) != matched_muon_idx.end())
      genmatched_tag.push_back(true);
    else
      genmatched_tag.push_back(false);
  }
  if (tag_muon_ttrack.empty())
    return;
  nt.nmuons = muons->size();
  nt.ntag = tag_muon_ttrack.size();

  // Add Lost Tracks to Packed cands
  std::vector<pat::PackedCandidate> tracks;
  for (const auto& container : {pfcands, lostTracks}) {
    for (const pat::PackedCandidate& trk : *container) {
      if (!probeSelection_(trk))
        continue;
      if (!trk.hasTrackDetails())
        continue;
      if (HighPurity_ && !trk.trackHighPurity())
        continue;
      tracks.emplace_back(trk);
    }
  }
  std::vector<unsigned> matched_track_idx;
  if (!iEvent.isRealData()) {
    auto reco_match_genmu1 =
        MatchReco<pat::PackedCandidate>(tracks, nt.genmu1_eta, nt.genmu1_phi, nt.genmu1_charge, genRecoDrMatch_);
    auto reco_match_genmu2 =
        MatchReco<pat::PackedCandidate>(tracks, nt.genmu2_eta, nt.genmu2_phi, nt.genmu2_charge, genRecoDrMatch_);
    if (reco_match_genmu1.first)
      matched_track_idx.push_back(reco_match_genmu1.second);
    if (reco_match_genmu2.first)
      matched_track_idx.push_back(reco_match_genmu2.second);
  }
  std::pair<std::vector<unsigned>, std::vector<unsigned>> trk_muon_map;
  for (const auto& mu : *muons) {
    if (muonOnly_ && !probeMuonSelection_(mu))
      continue;
    float minDR = 1000;
    unsigned int idx_trk;
    for (const auto& trk : tracks) {
      if (mu.charge() != trk.charge())
        continue;
      if (fabs(mu.vz() - trk.vz()) > maxdz_trk_mu_)
        continue;
      if (fabs(mu.pt() - trk.pt()) / mu.pt() > maxpt_relative_dif_trk_mu_)
        continue;
      float DR = deltaR(mu.eta(), mu.phi(), trk.eta(), trk.phi());
      if (minDR < DR)
        continue;
      minDR = DR;
      idx_trk = &trk - &tracks[0];
    }
    if (minDR > maxdr_trk_mu_)
      continue;
    trk_muon_map.first.push_back(idx_trk);
    trk_muon_map.second.push_back(&mu - &muons->at(0));
  }

  // Muon collection for jet cleaning
  std::vector<reco::Muon> muForJetCleaning;
  for (const auto& mu : *muons) {
    if (!muon::isLooseMuon(mu))
      continue;
    muForJetCleaning.push_back(mu);
  }

  std::vector<pat::Jet> corrJets;
  if (includeJets_) {
    for (const auto& jet : *jets) {
      if (CrossClean(jet, muForJetCleaning))
        continue;
      std::unique_ptr<pat::Jet> corrJet(jet.clone());
      // slimmed jets have corrections applied (L1FastJet, L2, L3) with pT cut at 10 GeV
      double jec = 1.0;
      corrJet->scaleEnergy(jec);
      // JER
      double smearFactor = 1.0;
      if (isMC_) {
        // Gen Jet Info
        for (const auto& genJet : *genJets) {
          nt.genJets_pt.push_back(genJet.pt());
          nt.genJets_eta.push_back(genJet.eta());
          nt.genJets_phi.push_back(genJet.phi());
          nt.genJets_mass.push_back(genJet.mass());
        }

        double jet_resolution = resolution.getResolution({{JME::Binning::JetPt, corrJet->pt()},
                                                          {JME::Binning::JetEta, corrJet->eta()},
                                                          {JME::Binning::Rho, *rhoHandle}});
        double jer_sf = resolution_sf.getScaleFactor({{JME::Binning::JetPt, corrJet->pt()},
                                                      {JME::Binning::JetEta, corrJet->eta()},
                                                      {JME::Binning::Rho, *rhoHandle}},
                                                     Variation::NOMINAL);
        // gen matching
        double min_dR = std::numeric_limits<double>::infinity();
        const reco::GenJet* matched_genJet = nullptr;
        for (const auto& genJet : *genJets) {
          double dR = deltaR(genJet, *corrJet);
          if (dR > min_dR)
            continue;
          if (dR >= 0.2)
            continue;
          min_dR = dR;
          matched_genJet = &genJet;
        }
        if (matched_genJet) {
          double dPt = corrJet->pt() - matched_genJet->pt();
          smearFactor = 1 + (jer_sf - 1.) * dPt / corrJet->pt();
        } else if (jer_sf > 1) {
          double sigma = jet_resolution * std::sqrt(jer_sf * jer_sf - 1);
          std::normal_distribution<> d(0, sigma);
          smearFactor = 1. + d(m_random_generator);
        }
        if (corrJet->pt() * smearFactor < 0.01) {
          smearFactor = 0.01 / corrJet->energy();
        }
      }
      corrJet->scaleEnergy(smearFactor);
      FillJetBranches(jet, *corrJet, nt, era_);
      float deepCSVprobb = -9999., deepCSVprobbb = -9999.;
      float deepFlavprobb = -9999., deepFlavprobbb = -9999.;
      float deepFlavproblepb = -9999.;
      for (const auto& pair : jet.getPairDiscri()) {
        if (pair.first == "pfDeepCSVJetTags:probb") {
          deepCSVprobb = pair.second;
        }
        if (pair.first == "pfDeepCSVJetTags:probbb") {
          deepCSVprobbb = pair.second;
        }
        if (pair.first == "pfDeepFlavourJetTags:probb") {
          deepFlavprobb = pair.second;
        }
        if (pair.first == "pfDeepFlavourJetTags:probbb") {
          deepFlavprobbb = pair.second;
        }
        if (pair.first == "pfDeepFlavourJetTags:problepb") {
          deepFlavproblepb = pair.second;
        }
      }
      if (deepCSVprobb != -9999. && deepCSVprobbb != -9999.) {
        nt.jets_bTag_deepCSV.push_back(deepCSVprobb + deepCSVprobbb);
      } else
        nt.jets_bTag_deepCSV.push_back(-9999.);
      if (deepFlavprobb != -9999. && deepFlavprobbb != -9999. && deepFlavproblepb != -9999.) {
        nt.jets_bTag_deepFlav.push_back(deepFlavprobb + deepFlavprobbb + deepFlavproblepb);
      } else
        nt.jets_bTag_deepFlav.push_back(-9999.);
    }
  }

  // run over tracks and probes once prior to filling tree to determine ordering of pairs
  // this is necessary to use tag-probe pair with highest "quality" later on in spark_tnp
  using t_pair_prob = std::pair<std::pair<int, int>, float>;
  std::vector<t_pair_prob> pair_vtx_probs;
  // loop over tags
  for (const auto& tag : tag_muon_ttrack) {
    auto tag_idx = &tag - &tag_muon_ttrack[0];
    // loop over probes
    for (const auto& probe : tracks) {
      auto probe_idx = &probe - &tracks[0];

      // apply cuts on pairs
      if (tag.first.charge() == probe.charge())
        continue;
      if (fabs(tag.first.vz() - probe.vz()) > pairDz_ && pairDz_ > 0)
        continue;

      float mass = DimuonMass(tag.first.pt(), tag.first.eta(), tag.first.phi(), probe.pt(), probe.eta(), probe.phi());
      if (mass < pairMassMin_ || mass > pairMassMax_)
        continue;

      std::vector<reco::TransientTrack> trk_pair = {tag.second, reco::TransientTrack(probe.pseudoTrack(), &(*bField))};
      KlFitter vtx(trk_pair);
      if (RequireVtxCreation_ && !vtx.status())
        continue;
      if (minSVtxProb_ > 0 && vtx.prob() < minSVtxProb_)
        continue;

      // save vtx prob to sort later
      pair_vtx_probs.emplace_back(std::make_pair(std::make_pair(tag_idx, probe_idx), vtx.prob()));
    }
  }
  nt.npairs = pair_vtx_probs.size();

  // reverse sort vertices by probability
  auto compare_vtx = [=](t_pair_prob& a, t_pair_prob& b) { return a.second > b.second; };
  std::sort(pair_vtx_probs.begin(), pair_vtx_probs.end(), compare_vtx);
  // assign sorted vtx indices to ranking
  map<std::pair<int, int>, int> pair_rank_vtx_prob;
  for (size_t i = 0; i < pair_vtx_probs.size(); i++)
    pair_rank_vtx_prob[pair_vtx_probs[i].first] = i;

  // Final pair selection
  // loop over tags
  for (const auto& tag : tag_muon_ttrack) {
    // loop over probe tracks
    for (const auto& probe : tracks) {
      // apply cuts on pairs
      if (tag.first.charge() == probe.charge())
        continue;
      if (fabs(tag.first.vz() - probe.vz()) > pairDz_ && pairDz_ > 0)
        continue;

      float mass = DimuonMass(tag.first.pt(), tag.first.eta(), tag.first.phi(), probe.pt(), probe.eta(), probe.phi());
      if (mass < pairMassMin_ || mass > pairMassMax_)
        continue;

      std::vector<reco::TransientTrack> trk_pair = {tag.second, reco::TransientTrack(probe.pseudoTrack(), &(*bField))};
      KlFitter vtx(trk_pair);
      if (RequireVtxCreation_ && !vtx.status())
        continue;
      if (minSVtxProb_ > 0 && vtx.prob() < minSVtxProb_)
        continue;

      auto it = std::find(trk_muon_map.first.begin(), trk_muon_map.first.end(), &probe - &tracks[0]);
      if (muonOnly_ && it == trk_muon_map.first.end())
        continue;

      FillTagBranches<pat::Muon, pat::PackedCandidate>(tag.first, tracks, nt, *pv);
      nt.tag_isMatchedGen = genmatched_tag[&tag - &tag_muon_ttrack[0]];
      FillMiniIso<pat::Muon, pat::PackedCandidate>(*pfcands, tag.first, *rhoJetsNC, nt, true);

      // Tag-trigger matching
      auto tagRef = muonsView->refAt(tag_muon_map[&tag - &tag_muon_ttrack[0]]);
      pat::TriggerObjectStandAloneRef tagl1Match = (*l1Matches)[tagRef];
      if (tagl1Match.isNonnull()) {
        nt.tag_l1pt = tagl1Match->pt();
        nt.tag_l1q = (*l1Qualities)[tagRef];
        nt.tag_l1dr = (*l1Drs)[tagRef];
      } else {
        nt.tag_l1pt = -99.;
        nt.tag_l1q = -99;
        nt.tag_l1dr = 99.;
      }

      pat::TriggerObjectStandAloneRef tagl1MatchByQ = (*l1MatchesByQ)[tagRef];
      if (tagl1MatchByQ.isNonnull()) {
        nt.tag_l1ptByQ = tagl1MatchByQ->pt();
        nt.tag_l1qByQ = (*l1QualitiesByQ)[tagRef];
        nt.tag_l1drByQ = (*l1DrsByQ)[tagRef];
      } else {
        nt.tag_l1ptByQ = -99.;
        nt.tag_l1qByQ = -99;
        nt.tag_l1drByQ = 99.;
      }

      embedTriggerMatching(iEvent, trigResults, tag.first, nt, tagFilters_, true);

      if (iEvent.isRealData())
        FillSimMatchingBranchesDummy(nt, true);
      else
        FillSimMatchingBranches(tag.first, nt, true);

      math::PtEtaPhiMLorentzVector mu1(tag.first.pt(), tag.first.eta(), tag.first.phi(), MU_MASS);
      math::PtEtaPhiMLorentzVector mu2(probe.pt(), probe.eta(), probe.phi(), MU_MASS);

      if (it != trk_muon_map.first.end()) {
        unsigned idx = std::distance(trk_muon_map.first.begin(), it);
        FillProbeBranches<pat::Muon, pat::PackedCandidate>(muons->at(trk_muon_map.second[idx]), tracks, nt, true, *pv);
        FillProbeBranchesSelector<pat::Muon>(muons->at(trk_muon_map.second[idx]), nt, probeSelectorBits_, true);
        FillMiniIso<pat::Muon, pat::PackedCandidate>(
            *pfcands, muons->at(trk_muon_map.second[idx]), *rhoJetsNC, nt, false);
        if (includeJets_)
          FindJetProbePair<pat::Jet, pat::Muon>(*jets, muons->at(trk_muon_map.second[idx]), nt);

        // Probe-trigger matching
        auto muRef = muonsView->refAt(trk_muon_map.second[idx]);
        pat::TriggerObjectStandAloneRef l1Match = (*l1Matches)[muRef];
        if (l1Match.isNonnull()) {
          nt.l1pt = l1Match->pt();
          nt.l1q = (*l1Qualities)[muRef];
          nt.l1dr = (*l1Drs)[muRef];
        } else {
          nt.l1pt = -99.;
          nt.l1q = -99;
          nt.l1dr = 99.;
        }

        pat::TriggerObjectStandAloneRef l1MatchByQ = (*l1MatchesByQ)[muRef];
        if (l1MatchByQ.isNonnull()) {
          nt.l1ptByQ = l1MatchByQ->pt();
          nt.l1qByQ = (*l1QualitiesByQ)[muRef];
          nt.l1drByQ = (*l1DrsByQ)[muRef];
        } else {
          nt.l1ptByQ = -99.;
          nt.l1qByQ = -99;
          nt.l1drByQ = 99.;
        }

        // Probe-trigger matching
        embedTriggerMatching(iEvent, trigResults, muons->at(trk_muon_map.second[idx]), nt, probeFilters_, false);

        if (iEvent.isRealData())
          FillSimMatchingBranchesDummy(nt, false);
        else
          FillSimMatchingBranches(muons->at(trk_muon_map.second[idx]), nt, false);

        // TuneP pair branches
        if (tag.first.tunePMuonBestTrack().isNonnull() &&
            muons->at(trk_muon_map.second[idx]).tunePMuonBestTrack().isNonnull()) {
          const reco::TrackRef tag_tuneP = tag.first.tunePMuonBestTrack();
          const reco::TrackRef probe_tuneP = muons->at(trk_muon_map.second[idx]).tunePMuonBestTrack();
          FillTunePPairBranches<reco::Track, reco::Track>(*tag_tuneP, *probe_tuneP, nt);

          std::vector<reco::TransientTrack> ttrk_pair_tuneP = {reco::TransientTrack(*tag_tuneP, &(*bField)),
                                                               reco::TransientTrack(*probe_tuneP, &(*bField))};
          KlFitter vtx_tuneP(ttrk_pair_tuneP);
          vtx_tuneP.fillNtuple(nt, true);
        } else {
          FillTunePPairBranchesDummy(nt);
        }
      } else {
        reco::Muon fakeMuon;
        fakeMuon.setP4(mu2);
        fakeMuon.setCharge(probe.charge());
        FillProbeBranches<reco::Muon, pat::PackedCandidate>(fakeMuon, tracks, nt, false, *pv);
        FillProbeBranchesSelector<reco::Muon>(fakeMuon, nt, probeSelectorBits_, false);
        FillMiniIso<pat::Muon, pat::PackedCandidate>(*pfcands, fakeMuon, *rhoJetsNC, nt, false);
        if (includeJets_)
          FindJetProbePair<pat::Jet, pat::Muon>(*jets, fakeMuon, nt);

        // store dummy trigger variables if offline muon is not found
        for (const auto& path : probeFilters_) {
          nt.probe_trg[&path - &probeFilters_[0]] = false;
          nt.probe_trg_pt[&path - &probeFilters_[0]] = -99;
          nt.probe_trg_eta[&path - &probeFilters_[0]] = -99;
          nt.probe_trg_phi[&path - &probeFilters_[0]] = -99;
          nt.probe_trg_dr[&path - &probeFilters_[0]] = 99;
        }
        nt.l1pt = -99.;
        nt.l1q = -99;
        nt.l1dr = 99.;
        nt.l1ptByQ = -99.;
        nt.l1qByQ = -99;
        nt.l1drByQ = 99.;

        FillSimMatchingBranchesDummy(nt, false);

        FillTunePPairBranchesDummy(nt);
      }

      PatPackedCandAndTransientTrk probe_pair = std::make_pair(probe, reco::TransientTrack(probe.pseudoTrack(), &(*bField)));
      
      FillPairBranches<PatMuonAndTransientTrk, PatPackedCandAndTransientTrk>(tag, probe_pair, nt, prop1_);

      vtx.fillNtuple(nt);

      auto it_match = std::find(matched_track_idx.begin(), matched_track_idx.end(), &probe - &tracks[0]);
      nt.probe_isMatchedGen = (it_match != matched_track_idx.end());

      nt.iprobe++;
      nt.pair_rank_vtx_prob = pair_rank_vtx_prob[{&tag - &tag_muon_ttrack[0], &probe - &tracks[0]}];
      nt.probe_isHighPurity = probe.trackHighPurity();

      t1->Fill();
    }
  }
}

// ------------ method called once each job just before starting event loop
// ------------
void MuonMiniAODAnalyzer::beginJob() {
  t1 = fs->make<TTree>("Events", "Events");
  nt.SetTree(t1);
  nt.CreateBranches(HLTPaths_, probeSelectorNames_);
  if (!tagFilters_.empty())
    nt.CreateExtraTrgBranches(tagFilters_, true);
  if (!probeFilters_.empty())
    nt.CreateExtraTrgBranches(probeFilters_, false);
}

// ------------ method called once each job just after ending the event loop
// ------------
void MuonMiniAODAnalyzer::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void MuonMiniAODAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no
  // validation
  // Please change this to state exactly what you do use, even if it is no
  // parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

///////////////////////

// define this as a plug-in
DEFINE_FWK_MODULE(MuonMiniAODAnalyzer);
