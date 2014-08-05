#ifndef _HCALCLOSURETEST_ANALYZERS_CALCRESPCORRDIJETS_H_
#define _HCALCLOSURETEST_ANALYZERS_CALCRESPCORRDIJETS_H_

//
// CalcRespCorrDiJets.h
//
//    description: Makes plots to calculate the response correction using dijets.
//
//    author: J.P. Chou, Brown
//
//

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

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


// forward declarations
class TH1D;
class TH2D;
class TFile;
class TTree;

//
// class declarations
//

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

 private:
  
};

class CalcRespCorrDiJets : public edm::EDAnalyzer {
 public:
  explicit CalcRespCorrDiJets(const edm::ParameterSet&);
  ~CalcRespCorrDiJets();
  
  
 private:
  virtual void beginJob();//(const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  
  // parameters
  bool debug_;                   // print debug statements
  std::string caloJetCollName_;  // label for the calo jet collection
  std::string caloJetCorrName_;  // label for the calo jet correction service
  std::string pfJetCollName_;    // label for the PF jet collection
  std::string pfJetCorrName_;    // label for the PF jet correction service
  std::string genJetCollName_;   // label for the genjet collection
  std::string RecHitLabelName_;  // label for the rechits
  std::string hbheRecHitInstance_; // instance for HBHERecHits
  std::string hfRecHitInstance_; // instance for HFRecHit
  std::string hoRecHitInstance_; // instance for HORecHit
  std::string rootHistFilename_; // name of the histogram file
  double maxDeltaEta_;           // maximum delta-|Eta| between Jets
  double minTagJetEta_;          // minimum |eta| of the tag jet
  double maxTagJetEta_;          // maximum |eta| of the tag jet
  double minSumJetEt_;           // minimum Sum of the tag and probe jet Et
  double minJetEt_;              // minimum Jet Et
  double maxThirdJetEt_;         // maximum 3rd jet Et
  double maxJetEMF_;             // maximum EMF of the tag and probe jets
  bool doCaloJets_;              // use CaloJets
  bool doPFJets_;                // use PFJets

  // root file/histograms
  TFile* rootfile_;

  TH1D* hPassSelCalo_;
  TH1D* hPassSelPF_;
  TH1D* h_types_;
  TH1D* h_ntypes_;
  TTree* calo_tree_;
  TTree* pf_tree_;
  float tcalojet_pt_, tcalojet_p_, tcalojet_eta_, tcalojet_phi_, tcalojet_emf_, tcalojet_scale_;
  float tcalojet_gendr_, tcalojet_genpt_, tcalojet_genp_;
  float tcalojet_EBE_, tcalojet_EEE_, tcalojet_HBE_, tcalojet_HEE_, tcalojet_HFE_;
  int tcalojet_ntwrs_;
  int tcalojet_twr_ieta_[100];
  float tcalojet_twr_eme_[100], tcalojet_twr_hade_[100];
  float pcalojet_pt_, pcalojet_p_, pcalojet_eta_, pcalojet_phi_, pcalojet_emf_, pcalojet_scale_;
  float pcalojet_gendr_, pcalojet_genpt_, pcalojet_genp_;
  float pcalojet_EBE_, pcalojet_EEE_, pcalojet_HBE_, pcalojet_HEE_, pcalojet_HFE_;
  int pcalojet_ntwrs_;
  int pcalojet_twr_ieta_[100];
  float pcalojet_twr_eme_[100], pcalojet_twr_hade_[100];
  float calo_dijet_deta_, calo_dijet_dphi_, calo_dijet_balance_;
  float calo_thirdjet_px_, calo_thirdjet_py_;
  float tpfjet_pt_, tpfjet_p_, tpfjet_eta_, tpfjet_phi_, tpfjet_scale_;
  float tpfjet_gendr_, tpfjet_genpt_, tpfjet_genp_;
  float tpfjet_EBE_, tpfjet_EEE_, tpfjet_HBE_, tpfjet_HEE_, tpfjet_HFE_;
  float tpfjet_unkown_E_, tpfjet_unkown_px_, tpfjet_unkown_py_, tpfjet_unkown_pz_;
  float tpfjet_chHad_E_, tpfjet_chHad_px_, tpfjet_chHad_py_, tpfjet_chHad_pz_;
  float tpfjet_electron_E_, tpfjet_electron_px_, tpfjet_electron_py_, tpfjet_electron_pz_;
  float tpfjet_muon_E_, tpfjet_muon_px_, tpfjet_muon_py_, tpfjet_muon_pz_;
  float tpfjet_photon_E_, tpfjet_photon_px_, tpfjet_photon_py_, tpfjet_photon_pz_;
  float tpfjet_Had0_E_, tpfjet_Had0_px_, tpfjet_Had0_py_, tpfjet_Had0_pz_;
  float tpfjet_HFHad_E_, tpfjet_HFHad_px_, tpfjet_HFHad_py_, tpfjet_HFHad_pz_;
  float tpfjet_HFEM_E_, tpfjet_HFEM_px_, tpfjet_HFEM_py_, tpfjet_HFEM_pz_;
  int tpfjet_unkown_n_, tpfjet_chHad_n_, tpfjet_electron_n_, tpfjet_muon_n_;
  int tpfjet_photon_n_, tpfjet_Had0_n_, tpfjet_HFHad_n_, tpfjet_HFEM_n_;
  float tpfjet_chHad_EcalE_, tpfjet_HFHad_EcalE_;
  int tpfjet_ntwrs_;
  int tpfjet_twr_ieta_[1000];
  float tpfjet_twr_eme_[1000], tpfjet_twr_hade_[1000], tpfjet_twr_frac_[1000];
  float ppfjet_pt_, ppfjet_p_, ppfjet_eta_, ppfjet_phi_, ppfjet_scale_;
  float ppfjet_gendr_, ppfjet_genpt_, ppfjet_genp_;
  float ppfjet_EBE_, ppfjet_EEE_, ppfjet_HBE_, ppfjet_HEE_, ppfjet_HFE_;
  float ppfjet_unkown_E_, ppfjet_unkown_px_, ppfjet_unkown_py_, ppfjet_unkown_pz_;
  float ppfjet_chHad_E_, ppfjet_chHad_px_, ppfjet_chHad_py_, ppfjet_chHad_pz_;
  float ppfjet_electron_E_, ppfjet_electron_px_, ppfjet_electron_py_, ppfjet_electron_pz_;
  float ppfjet_muon_E_, ppfjet_muon_px_, ppfjet_muon_py_, ppfjet_muon_pz_;
  float ppfjet_photon_E_, ppfjet_photon_px_, ppfjet_photon_py_, ppfjet_photon_pz_;
  float ppfjet_Had0_E_, ppfjet_Had0_px_, ppfjet_Had0_py_, ppfjet_Had0_pz_;
  float ppfjet_HFHad_E_, ppfjet_HFHad_px_, ppfjet_HFHad_py_, ppfjet_HFHad_pz_;
  float ppfjet_HFEM_E_, ppfjet_HFEM_px_, ppfjet_HFEM_py_, ppfjet_HFEM_pz_;
  int ppfjet_unkown_n_, ppfjet_chHad_n_, ppfjet_electron_n_, ppfjet_muon_n_;
  int ppfjet_photon_n_, ppfjet_Had0_n_, ppfjet_HFHad_n_, ppfjet_HFEM_n_;
  float ppfjet_chHad_EcalE_, ppfjet_HFHad_EcalE_;
  int ppfjet_ntwrs_;
  int ppfjet_twr_ieta_[1000];
  float ppfjet_twr_eme_[1000], ppfjet_twr_hade_[1000], ppfjet_twr_frac_[1000];
  float pf_dijet_deta_, pf_dijet_dphi_, pf_dijet_balance_;
  float pf_thirdjet_px_, pf_thirdjet_py_;

  float hfieta[13];
  //float hfieta[13] = {2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191};
  int maxhfiphi;

  // helper functions
  double deltaR(const reco::Jet* j1, const reco::Jet* j2);

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

};


#endif
