#ifndef _HCALCLOSURETEST_ANALYZERS_CALCRESPCORRPHOJETS_H_
#define _HCALCLOSURETEST_ANALYZERS_CALCRESPCORRPHOJETS_H_

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
//#include "RecoEgamma/EgammaTools/interface/ggPFClusters.h"
//#include "RecoEgamma/EgammaTools/interface/ggPFESClusters.h"
//#include "RecoEgamma/EgammaTools/interface/ggPFTracks.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// forward declarations
class TH1D;
class TH2D;
class TFile;
class TTree;

//
// class declarations
//

class PhotonPair : protected std::pair<const reco::Photon*, double> {
 public:
  PhotonPair() {
    first=0;
    second=0.0;
    fIdx=-1;
  }
  PhotonPair(const reco::Photon* ph, double pt, int setIdx=-1) {
    first=ph;
    second=pt;
    fIdx=setIdx;
  }
  ~PhotonPair() {}

  inline const reco::Photon* photon(void) const { return first; }
  inline void photon(const reco::Photon* ph) { first=ph; return; }
  inline double pt(void) const { return second; }
  inline void pt(double d) { second=d; return; }
  void idx(int set_idx) { fIdx=set_idx; };
  int idx() const { return fIdx; }
  bool isValid() const { return (first!=NULL) ? true:false; }

 private:
  int fIdx; // index in the photon collection
};


class CaloJetCorretPair : protected std::pair<const reco::CaloJet*, double> {
 public:
  CaloJetCorretPair() {
    first=0;
    second=1.0;
  }
  CaloJetCorretPair(const reco::CaloJet* j, double s) {
    first=j;
    second=s;
  }
  ~CaloJetCorretPair() {}

  inline const reco::CaloJet* jet(void) const { return first; }
  inline void jet(const reco::CaloJet* j) { first=j; return; }
  inline double scale(void) const { return second; }
  inline void scale(double d) { second=d; return; }
  double scaledEt() const { return first->et() * second; }
  bool isValid() const { return (first!=NULL) ? true:false; }

 private:
  
};

class PFJetCorretPair : protected std::pair<const reco::PFJet*, double> {
 public:
  PFJetCorretPair() {
    first=0;
    second=1.0;
  }
  PFJetCorretPair(const reco::PFJet* j, double s) {
    first=j;
    second=s;
  }
  ~PFJetCorretPair() {}

  inline const reco::PFJet* jet(void) const { return first; }
  inline void jet(const reco::PFJet* j) { first=j; return; }
  inline double scale(void) const { return second; }
  inline void scale(double d) { second=d; return; }
  double scaledEt() const { return first->et() * second; }
  bool isValid() const { return (first!=NULL) ? true:false; }

 private:
  
};

// --------------------------------------------
// Main class
// --------------------------------------------

class CalcRespCorrPhotonPlusJet : public edm::EDAnalyzer {
 public:
  explicit CalcRespCorrPhotonPlusJet(const edm::ParameterSet&);
  ~CalcRespCorrPhotonPlusJet();

  //float  pfEcalIso(const reco::Photon*, edm::Handle<reco::PFCandidateCollection>, float, float, float, float, float, float, float,       reco::PFCandidate::ParticleType);
  float pfEcalIso(const reco::Photon* localPho1, edm::Handle<reco::PFCandidateCollection> pfHandle, float dRmax, float dRVetoBarrel, float dRVetoEndcap, float etaStripBarrel, float etaStripEndcap, float energyBarrel, float energyEndcap, reco::PFCandidate::ParticleType pfToUse);

  //float  pfHcalIso(const reco::Photon*, edm::Handle<reco::PFCandidateCollection>, float, float, reco::PFCandidate::ParticleType);
  float pfHcalIso(const reco::Photon* localPho,edm::Handle<reco::PFCandidateCollection> pfHandle,float dRmax, float dRveto,reco::PFCandidate::ParticleType pfToUse);

  //std::vector<float> pfTkIsoWithVertex(const reco::Photon*, edm::Handle<reco::PFCandidateCollection>, edm::Handle<reco::VertexCollection>, float, float, float, float, float, float, reco::PFCandidate::ParticleType);
  std::vector<float> pfTkIsoWithVertex(const reco::Photon* localPho1, edm::Handle<reco::PFCandidateCollection> pfHandle, edm::Handle<reco::VertexCollection> vtxHandle, float dRmax, float dRvetoBarrel, float dRvetoEndcap, float ptMin, float dzMax, float dxyMax, reco::PFCandidate::ParticleType pfToUse);

 private:
  virtual void beginJob();//(const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  void beginRun(const edm::Run&, const edm::EventSetup&);

  
  // parameters
  int debug_;                      // print debug statements
  edm::InputTag rhoCollection_;
  std::string photonCollName_;      // label for the photon collection
  std::string caloJetCollName_;     // label for the calo jet collection
  std::string caloJetCorrName_;     // label for the calo jet correction service
  std::string pfJetCollName_;       // label for the PF jet collection
  std::string pfJetCorrName_;       // label for the PF jet correction service
  std::string genJetCollName_;      // label for the genjet collection
  std::string genParticleCollName_; // label for the genparticle collection
  std::string genEventInfoName_;    // label for the generator event info collection
  std::string hbheRecHitName_;      // label for HBHERecHits collection
  std::string hfRecHitName_;        // label for HFRecHit collection
  std::string hoRecHitName_;        // label for HORecHit collection
  std::string rootHistFilename_;    // name of the histogram file

  bool allowNoPhoton_; // whether module is used for dijet analysis
  double photonPtMin_;   // lowest value of the leading photon pT
  double photonJetDPhiMin_; // phi angle between the leading photon and the leading jet
  double jetEtMin_;      // lowest value of the leading jet ET
  double jet2EtMax_;     // largest value of the subleading jet ET
  double jet3EtMax_;     // largest value of the third jet ET
  std::vector<std::string> photonTrigNamesV_; // photon trigger names
  std::vector<std::string> jetTrigNamesV_; // jet trigger names
  bool writeTriggerPrescale_; // whether attempt to record the prescale

  //bool doPhotons_;                  // use Photons
  bool doCaloJets_;                 // use CaloJets
  bool doPFJets_;                   // use PFJets
  bool doGenJets_;                  // use GenJets
  bool ignoreHLT_;

  // root file/histograms
  TFile* rootfile_;

  TH1D* h_types_;
  TH1D* h_ntypes_;
  TH1D* h_ietaHCAL_;
  TH1D* h_etaHFHAD_;
  TH1D* h_etaHFEM_;
  TH1D* h_ietaHO_;
  TH1D* h_HFHAD_n_;
  TH1D* h_HFEM_n_;
  TH1D* h_HFHAD_type_;
  TH1D* h_HFEM_type_;
  TH1D* h_HBHE_n_;
  TH1D* h_HF_n_;
  TH1D* h_HO_n_;
  TH1D* h_twrietas_;
  TH2D* h_rechitspos_;
  TH1D* h_hbherecoieta_;

  TH1D* hcount;
  TH1D* h_nPho;
  TH1D* h_nPFJets;
  TH1D* h_leadJetCorrEt;
  TH1D* h_2ndleadJetCorrEt;
  TH1D* h_3rdleadJetCorrEt;
  TH1D* h_dRphojet;
  TH1D* h_dphi;
  TH1D*   h_pfrecoOgen_et;

  TTree* misc_tree_; // misc.information. Will be filled only once
  TTree* calo_tree_;
  TTree* pf_tree_;

  // trigger info
  HLTConfigProvider hltConfig_; // variable for the access
  std::vector<int> photonTrigFired_;
  std::vector<int> photonTrigPrescale_;
  std::vector<int> jetTrigFired_;
  std::vector<int> jetTrigPrescale_;

  // Event info
  int runNumber_, lumiBlock_, eventNumber_;
  float eventWeight_;
  int nPhotons_, nGenJets_;
  int nCaloJets_, nPFJets_;

  // photon info
  float rho2012_;
  float tagPho_pt_, pho_2nd_pt_, tagPho_energy_, tagPho_eta_, tagPho_phi_, tagPho_sieie_;
  float tagPho_HoE_, tagPho_r9_, tagPho_EcalIsoDR04_, tagPho_HcalIsoDR04_, tagPho_HcalIsoDR0412_, tagPho_TrkIsoHollowDR04_, tagPho_pfiso_myphoton03_;
  float tagPho_pfiso_myneutral03_;
  std::vector<std::vector<float> >  tagPho_pfiso_mycharged03 ;
  int tagPho_pixelSeed_;
  int tagPho_ConvSafeEleVeto_;
  int tagPho_idTight_, tagPho_idLoose_;
  float tagPho_genPt_, tagPho_genEnergy_, tagPho_genEta_, tagPho_genPhi_;
  float tagPho_genDeltaR_;

  // Calo jets
  float tcalojet_et_;
  float tcalojet_pt_, tcalojet_p_, tcalojet_eta_, tcalojet_phi_, tcalojet_emf_, tcalojet_scale_;
  float tcalojet_gendr_, tcalojet_genpt_, tcalojet_genp_;
  float tcalojet_EBE_, tcalojet_EEE_, tcalojet_HBE_, tcalojet_HEE_, tcalojet_HFEtot_;
  float tcalojet_HFEcalE_, tcalojet_HFHadE_;
  int tcalojet_ntwrs_;
  int tcalojet_twr_ieta_[100];
  float tcalojet_twr_totE_[100], tcalojet_twr_eme_[100], tcalojet_twr_hade_[100], tcalojet_twr_hadoe_[100];
  //float tcalojet_twr_numBadEcalCells[100], tcalojet_twr_numBadHcalCells[100];
  //float tcalojet_twr_numProblematicEcalCells[100], tcalojet_twr_numProblematicHcalCells[100];
  //float tcalojet_twr_numRecoveredEcalCells[100], tcalojet_twr_numRecoveredHcalCells[100];

  float pcalojet_et_;
  float pcalojet_pt_, pcalojet_p_, pcalojet_eta_, pcalojet_phi_, pcalojet_emf_, pcalojet_scale_;
  float pcalojet_gendr_, pcalojet_genpt_, pcalojet_genp_;
  float pcalojet_EBE_, pcalojet_EEE_, pcalojet_HBE_, pcalojet_HEE_, pcalojet_HFEtot_;
  int pcalojet_ntwrs_;
  int pcalojet_twr_ieta_[100];
  float pcalojet_twr_eme_[100], pcalojet_twr_hade_[100];

  float calo_thirdjet_et_, calo_thirdjet_pt_, calo_thirdjet_p_, calo_thirdjet_px_, calo_thirdjet_py_;
  float calo_thirdjet_E_, calo_thirdjet_eta_, calo_thirdjet_phi_, calo_thirdjet_scale_;

  // Particle-flow jets
  // leading Et jet info
  float ppfjet_pt_, ppfjet_p_, ppfjet_E_, ppfjet_eta_, ppfjet_phi_, ppfjet_scale_;
  float ppfjet_NeutralHadronFrac_, ppfjet_NeutralEMFrac_;
  int ppfjet_nConstituents_;
  float ppfjet_ChargedHadronFrac_, ppfjet_ChargedMultiplicity_, ppfjet_ChargedEMFrac_;
  float ppfjet_gendr_, ppfjet_genpt_, ppfjet_genp_, ppfjet_genE_;
  //float ppfjet_EBE_, ppfjet_EEE_, ppfjet_HBE_, ppfjet_HEE_, ppfjet_HFE_;
  float ppfjet_unkown_E_, ppfjet_unkown_px_, ppfjet_unkown_py_, ppfjet_unkown_pz_, ppfjet_unkown_EcalE_;
  float ppfjet_electron_E_, ppfjet_electron_px_, ppfjet_electron_py_, ppfjet_electron_pz_, ppfjet_electron_EcalE_;
  float ppfjet_muon_E_, ppfjet_muon_px_, ppfjet_muon_py_, ppfjet_muon_pz_, ppfjet_muon_EcalE_;
  float ppfjet_photon_E_, ppfjet_photon_px_, ppfjet_photon_py_, ppfjet_photon_pz_, ppfjet_photon_EcalE_;
  int ppfjet_unkown_n_, ppfjet_electron_n_, ppfjet_muon_n_, ppfjet_photon_n_;
  int ppfjet_had_n_;
  std::vector<float> ppfjet_had_E_, ppfjet_had_px_, ppfjet_had_py_, ppfjet_had_pz_, ppfjet_had_EcalE_, ppfjet_had_rawHcalE_, ppfjet_had_emf_, ppfjet_had_E_mctruth_;
  std::vector<int> ppfjet_had_id_, ppfjet_had_candtrackind_, ppfjet_had_mcpdgId_, ppfjet_had_ntwrs_;
  int ppfjet_ntwrs_;
  std::vector<int> ppfjet_twr_ieta_, ppfjet_twr_iphi_, ppfjet_twr_depth_, ppfjet_twr_subdet_, ppfjet_twr_candtrackind_, ppfjet_twr_hadind_, ppfjet_twr_elmttype_, ppfjet_twr_clusterind_;
  std::vector<float> ppfjet_twr_hade_, ppfjet_twr_frac_, ppfjet_twr_dR_;
  int ppfjet_cluster_n_;
  std::vector<float> ppfjet_cluster_eta_, ppfjet_cluster_phi_, ppfjet_cluster_dR_;
  int ppfjet_ncandtracks_;
  std::vector<float> ppfjet_candtrack_px_, ppfjet_candtrack_py_, ppfjet_candtrack_pz_, ppfjet_candtrack_EcalE_;

  // subleading Et jet info
  float pfjet2_pt_, pfjet2_p_, pfjet2_E_, pfjet2_eta_, pfjet2_phi_, pfjet2_scale_;
  float pfjet2_NeutralHadronFrac_, pfjet2_NeutralEMFrac_;
  int pfjet2_nConstituents_;
  float pfjet2_ChargedHadronFrac_, pfjet2_ChargedMultiplicity_, pfjet2_ChargedEMFrac_;
  float pfjet2_gendr_, pfjet2_genpt_, pfjet2_genp_, pfjet2_genE_;
  //float pfjet2_EBE_, pfjet2_EEE_, pfjet2_HBE_, pfjet2_HEE_, pfjet2_HFE_;
  float pfjet2_unkown_E_, pfjet2_unkown_px_, pfjet2_unkown_py_, pfjet2_unkown_pz_, pfjet2_unkown_EcalE_;
  float pfjet2_electron_E_, pfjet2_electron_px_, pfjet2_electron_py_, pfjet2_electron_pz_, pfjet2_electron_EcalE_;
  float pfjet2_muon_E_, pfjet2_muon_px_, pfjet2_muon_py_, pfjet2_muon_pz_, pfjet2_muon_EcalE_;
  float pfjet2_photon_E_, pfjet2_photon_px_, pfjet2_photon_py_, pfjet2_photon_pz_, pfjet2_photon_EcalE_;
  int pfjet2_unkown_n_, pfjet2_electron_n_, pfjet2_muon_n_, pfjet2_photon_n_;
  int pfjet2_had_n_;
  std::vector<float> pfjet2_had_E_, pfjet2_had_px_, pfjet2_had_py_, pfjet2_had_pz_, pfjet2_had_EcalE_, pfjet2_had_rawHcalE_, pfjet2_had_emf_, pfjet2_had_E_mctruth_;
  std::vector<int> pfjet2_had_id_, pfjet2_had_candtrackind_, pfjet2_had_mcpdgId_, pfjet2_had_ntwrs_;
  int pfjet2_ntwrs_;
  std::vector<int> pfjet2_twr_ieta_, pfjet2_twr_iphi_, pfjet2_twr_depth_, pfjet2_twr_subdet_, pfjet2_twr_candtrackind_, pfjet2_twr_hadind_, pfjet2_twr_elmttype_, pfjet2_twr_clusterind_;
  std::vector<float> pfjet2_twr_hade_, pfjet2_twr_frac_, pfjet2_twr_dR_;
  int pfjet2_cluster_n_;
  std::vector<float> pfjet2_cluster_eta_, pfjet2_cluster_phi_, pfjet2_cluster_dR_;
  int pfjet2_ncandtracks_;
  std::vector<float> pfjet2_candtrack_px_, pfjet2_candtrack_py_, pfjet2_candtrack_pz_, pfjet2_candtrack_EcalE_;

  float pf_thirdjet_et_;
  float pf_thirdjet_pt_, pf_thirdjet_p_, pf_thirdjet_px_, pf_thirdjet_py_;
  float pf_thirdjet_E_, pf_thirdjet_eta_, pf_thirdjet_phi_, pf_thirdjet_scale_;

  // ------------------------------
  // helper functions

  template<class JetPair_type>
  float calc_dPhi(const PhotonPair &pho, const JetPair_type &jet) {
    if (!pho.isValid() || !jet.isValid()) return 9999.;
    float phi1=pho.photon()->phi();
    float phi2=jet.jet()->phi();
    float dphi=fabs(phi1-phi2);
    const float cPi= 4*atan(1);
    while (dphi>cPi) dphi = fabs(2*cPi - dphi);
    return dphi;
  }

  double deltaR(const reco::Jet* j1, const reco::Jet* j2);
  double deltaR(const double eta1, const double phi1, const double eta2, const double phi2);
  int getEtaPhi(const DetId id);
  int getEtaPhi(const HcalDetId id);

  void clear_leadingPfJetVars();
  void copy_leadingPfJetVars_to_pfJet2();

  template<class Jet_type>
  double deltaR(const PhotonPair &photon, const Jet_type *jet) {
    if (!photon.isValid()) return 9999.;
    return deltaR(photon.photon()->eta(),photon.photon()->phi(),
		  jet->eta(), jet->phi());
  }


  struct CaloJetCorretPairComp {
    inline bool operator() ( const CaloJetCorretPair& a, const CaloJetCorretPair& b) {
      return (a.jet()->pt()*a.scale()) > (b.jet()->pt()*b.scale());
    }
  };

  struct PFJetCorretPairComp {
    inline bool operator() ( const PFJetCorretPair& a, const PFJetCorretPair& b) {
      return (a.jet()->pt()*a.scale()) > (b.jet()->pt()*b.scale());
    }
  };

  struct PhotonPairComp {
    inline bool operator() ( const PhotonPair& a, const PhotonPair& b) {
      return ( (a.photon()->pt()) > (b.photon()->pt()) );
    }
  };


};


#endif
