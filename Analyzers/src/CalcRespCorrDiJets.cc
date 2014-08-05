//
// CalcRespCorr.cc
//
//   description: Calculation of dijet response corrections
//
//   author: J.P. Chou, Brown
//
//

#include "HcalClosureTest/Analyzers/interface/CalcRespCorrDiJets.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TClonesArray.h"

#include <vector>
#include <set>
//using namespace reco;

CalcRespCorrDiJets::CalcRespCorrDiJets(const edm::ParameterSet& iConfig)
{
  // set parameters
  caloJetCollName_   = iConfig.getParameter<std::string>("caloJetCollName");
  caloJetCorrName_   = iConfig.getParameter<std::string>("caloJetCorrName");
  pfJetCollName_     = iConfig.getParameter<std::string>("pfJetCollName");
  pfJetCorrName_     = iConfig.getParameter<std::string>("pfJetCorrName");
  genJetCollName_    = iConfig.getParameter<std::string>("genJetCollName");
  RecHitLabelName_   = iConfig.getParameter<std::string>("RecHitLabelName");
  hbheRecHitInstance_  = iConfig.getParameter<std::string>("hbheRecHitInstance");
  hfRecHitInstance_  = iConfig.getParameter<std::string>("hfRecHitInstance");
  hoRecHitInstance_  = iConfig.getParameter<std::string>("hoRecHitInstance");
  rootHistFilename_  = iConfig.getParameter<std::string>("rootHistFilename");
  maxDeltaEta_       = iConfig.getParameter<double>("maxDeltaEta");
  minTagJetEta_      = iConfig.getParameter<double>("minTagJetEta");
  maxTagJetEta_      = iConfig.getParameter<double>("maxTagJetEta");
  minSumJetEt_       = iConfig.getParameter<double>("minSumJetEt");
  minJetEt_          = iConfig.getParameter<double>("minJetEt");
  maxThirdJetEt_     = iConfig.getParameter<double>("maxThirdJetEt");
  maxJetEMF_         = iConfig.getParameter<double>("maxJetEMF");
  doCaloJets_        = iConfig.getParameter<bool>("doCaloJets");
  doPFJets_          = iConfig.getParameter<bool>("doPFJets");
  debug_             = iConfig.getUntrackedParameter<bool>("debug", false);
}

CalcRespCorrDiJets::~CalcRespCorrDiJets()
{
}
  
//
// member functions
//
  
// ------------ method called to for each event  ------------
void
CalcRespCorrDiJets::analyze(const edm::Event& iEvent, const edm::EventSetup& evSetup)
{ 
  /*~~~edm::Handle<reco::GenJetCollection> genjets;
  iEvent.getByLabel(genJetCollName_,genjets);
  if(!genjets.isValid()) {
    throw edm::Exception(edm::errors::ProductNotFound)
      << " could not find GenJetCollection named " << genJetCollName_ << ".\n";
    return;
    }*/
  edm::Handle<std::vector<reco::GenJet>> genjets;
  iEvent.getByLabel(genJetCollName_,genjets);
  if(!genjets.isValid()) {
    throw edm::Exception(edm::errors::ProductNotFound)
      << " could not find GenJet vector named " << genJetCollName_ << ".\n";
    return;
  }

  unsigned int testEvent = 600000;
  if(iEvent.id().event() == testEvent){
    std::cout << "Event " << iEvent.id().event() << std::endl;
  }
  
  //std::cout << "Event " << iEvent.id().event() << std::endl;

  if(doCaloJets_){
    edm::Handle<reco::CaloJetCollection> calojets;
    iEvent.getByLabel(caloJetCollName_,calojets);
    if(!calojets.isValid()) {
      throw edm::Exception(edm::errors::ProductNotFound)
	<< " could not find CaloJetCollection named " << caloJetCollName_ << ".\n";
      return;
    }

    const JetCorrector* correctorCalo = JetCorrector::getJetCorrector(caloJetCorrName_,evSetup);

    //////////////////////////////
    // Event Selection
    //////////////////////////////
    
    // determine which cut results in failure
    int passSelCalo=0;
    
    // sort jets by corrected et
    std::set<CaloJetCorretPair, CaloJetCorretPairComp> calojetcorretpairset;
    for(reco::CaloJetCollection::const_iterator it=calojets->begin(); it!=calojets->end(); ++it) {
      const reco::CaloJet* jet=&(*it);
      calojetcorretpairset.insert( CaloJetCorretPair(jet, correctorCalo->correction(jet->p4())) );
    }
    
    // find highest two (corrected) et jets
    CaloJetCorretPair calo_tag, calo_probe;
    calo_thirdjet_px_=calo_thirdjet_py_=0.0;
    int cntr=0;
    for(std::set<CaloJetCorretPair, CaloJetCorretPairComp>::const_iterator it=calojetcorretpairset.begin(); it!=calojetcorretpairset.end(); ++it) {
      CaloJetCorretPair jet=(*it);
      ++cntr;
      if(cntr==1) calo_tag=jet;
      else if(cntr==2) calo_probe=jet;
      else {
	calo_thirdjet_px_ += jet.scale()*jet.jet()->px();
	calo_thirdjet_py_ += jet.scale()*jet.jet()->py();
      }
    }
    
    if(calo_tag.jet() && calo_probe.jet()){
      // require that the first two jets are above some minimum,
      // and the rest are below some maximum
      if((calo_tag.jet()->et()+calo_probe.jet()->et())<minSumJetEt_) passSelCalo |= 0x1;
      if(calo_tag.jet()->et()<minJetEt_ || calo_probe.jet()->et()<minJetEt_) passSelCalo |= 0x2;
      if(sqrt(calo_thirdjet_px_*calo_thirdjet_px_ + calo_thirdjet_py_*calo_thirdjet_py_)>maxThirdJetEt_) passSelCalo |= 0x4;
      
      // force the tag jet to have the smaller |eta|
      if(std::fabs(calo_tag.jet()->eta())>std::fabs(calo_probe.jet()->eta())) {
	CaloJetCorretPair temp=calo_tag;
	calo_tag=calo_probe;
	calo_probe=temp;
      }
      
      // eta cuts
      double dAbsEta=std::fabs(std::fabs(calo_tag.jet()->eta())-std::fabs(calo_probe.jet()->eta()));
      if(dAbsEta>maxDeltaEta_) passSelCalo |= 0x8;
      if(fabs(calo_tag.jet()->eta())<minTagJetEta_) passSelCalo |= 0x10;
      if(fabs(calo_tag.jet()->eta())>maxTagJetEta_) passSelCalo |= 0x10;
      
      // emf cuts
      if(calo_tag.jet()->emEnergyFraction()>maxJetEMF_) passSelCalo |= 0x20;
      if(calo_probe.jet()->emEnergyFraction()>maxJetEMF_) passSelCalo |= 0x20;
    }
    else{
      passSelCalo = 0x40;
    }

    // make the cuts
    hPassSelCalo_->Fill(passSelCalo);
    if(!passSelCalo){
      
      // dump
      if(debug_) {
	std::cout << "Run: " << iEvent.id().run() << "; Event: " << iEvent.id().event() << std::endl;
	for(reco::CaloJetCollection::const_iterator it=calojets->begin(); it!=calojets->end(); ++it) {
	  const reco::CaloJet *jet=&(*it);
	  std::cout << "istag=" << (jet==calo_tag.jet()) << "; isprobe=" << (jet==calo_probe.jet()) << "; et=" << jet->et() << "; eta=" << jet->eta() << std::endl;
	}
      }
      
      // fill tag jet variables
      tcalojet_pt_    = calo_tag.jet()->pt();
      tcalojet_p_     = calo_tag.jet()->p();
      tcalojet_eta_   = calo_tag.jet()->eta();
      tcalojet_phi_   = calo_tag.jet()->phi();
      tcalojet_emf_   = calo_tag.jet()->emEnergyFraction();
      tcalojet_scale_ = calo_tag.scale();
      tcalojet_EBE_   = calo_tag.jet()->emEnergyInEB();
      tcalojet_EEE_   = calo_tag.jet()->emEnergyInEE();
      tcalojet_HBE_   = calo_tag.jet()->hadEnergyInHB();
      tcalojet_HEE_   = calo_tag.jet()->hadEnergyInHE();
      tcalojet_HFE_   = calo_tag.jet()->emEnergyInHF() + calo_tag.jet()->hadEnergyInHF();
      tcalojet_ntwrs_=0;
      std::vector<CaloTowerPtr> tagconst=calo_tag.jet()->getCaloConstituents();
      for(std::vector<CaloTowerPtr>::const_iterator it=tagconst.begin(); it!=tagconst.end(); ++it) {
	int ieta=(*it)->id().ieta();
	int ietaAbs=(*it)->id().ietaAbs();
	tcalojet_twr_ieta_[tcalojet_ntwrs_]=ieta;
	if(ietaAbs<=29) {
	  tcalojet_twr_eme_[tcalojet_ntwrs_] = (*it)->emEnergy();
	  tcalojet_twr_hade_[tcalojet_ntwrs_] = (*it)->hadEnergy();
	} else {
	  tcalojet_twr_eme_[tcalojet_ntwrs_] = 0;
	  tcalojet_twr_hade_[tcalojet_ntwrs_] = (*it)->emEnergy()+(*it)->hadEnergy();
	}
	++tcalojet_ntwrs_;
      }
      
      // fill probe jet variables
      pcalojet_pt_    = calo_probe.jet()->pt();
      pcalojet_p_     = calo_probe.jet()->p();
      pcalojet_eta_   = calo_probe.jet()->eta();
      pcalojet_phi_   = calo_probe.jet()->phi();
      pcalojet_emf_   = calo_probe.jet()->emEnergyFraction();
      pcalojet_scale_ = calo_probe.scale();
      pcalojet_EBE_   = calo_probe.jet()->emEnergyInEB();
      pcalojet_EEE_   = calo_probe.jet()->emEnergyInEE();
      pcalojet_HBE_   = calo_probe.jet()->hadEnergyInHB();
      pcalojet_HEE_   = calo_probe.jet()->hadEnergyInHE();
      pcalojet_HFE_   = calo_probe.jet()->emEnergyInHF() + calo_probe.jet()->hadEnergyInHF();
      pcalojet_ntwrs_=0;
      std::vector<CaloTowerPtr> probeconst=calo_probe.jet()->getCaloConstituents();
      for(std::vector<CaloTowerPtr>::const_iterator it=probeconst.begin(); it!=probeconst.end(); ++it) {
	int ieta=(*it)->id().ieta();
	int ietaAbs=(*it)->id().ietaAbs();
	pcalojet_twr_ieta_[pcalojet_ntwrs_]=ieta;
	if(ietaAbs<=29) {
	  pcalojet_twr_eme_[pcalojet_ntwrs_] = (*it)->emEnergy();
	  pcalojet_twr_hade_[pcalojet_ntwrs_] = (*it)->hadEnergy();
	} else {
	  pcalojet_twr_eme_[pcalojet_ntwrs_] = 0;
	  pcalojet_twr_hade_[pcalojet_ntwrs_] = (*it)->emEnergy()+(*it)->hadEnergy();
	}
	++pcalojet_ntwrs_;
      }
      
      // fill genjet tag/probe variables
      tcalojet_gendr_ = 99999.;
      tcalojet_genpt_ = 0;
      tcalojet_genp_  = 0;
      pcalojet_gendr_ = 99999.;
      pcalojet_genpt_ = 0;
      pcalojet_genp_  = 0;
      for(std::vector<reco::GenJet>::const_iterator it=genjets->begin(); it!=genjets->end(); ++it){
	const reco::GenJet* jet=&(*it);
	double dr=deltaR(jet, calo_probe.jet());
	if(dr<pcalojet_gendr_) {
	  pcalojet_gendr_ = dr;
	  pcalojet_genpt_ = jet->pt();
	  pcalojet_genp_ = jet->p();
	}
	dr=deltaR(jet, calo_tag.jet());
	if(dr<tcalojet_gendr_) {
	  tcalojet_gendr_ = dr;
	  tcalojet_genpt_ = jet->pt();
	  tcalojet_genp_ = jet->p();
	}
      }
      /*~~~for(reco::GenJetCollection::const_iterator it=genjets->begin(); it!=genjets->end(); ++it) {
	const reco::GenJet* jet=&(*it);
	double dr=deltaR(jet, calo_probe.jet());
	if(dr<pcalojet_gendr_) {
	  pcalojet_gendr_ = dr;
	  pcalojet_genpt_ = jet->pt();
	  pcalojet_genp_ = jet->p();
	}
	dr=deltaR(jet, calo_tag.jet());
	if(dr<tcalojet_gendr_) {
	  tcalojet_gendr_ = dr;
	  tcalojet_genpt_ = jet->pt();
	  tcalojet_genp_ = jet->p();
	}
	}*/
      
      // fill dijet variables
      calo_dijet_deta_=std::fabs(std::fabs(calo_tag.jet()->eta())-std::fabs(calo_probe.jet()->eta()));
      calo_dijet_dphi_=calo_tag.jet()->phi()-calo_probe.jet()->phi();
      if(calo_dijet_dphi_>3.1415) calo_dijet_dphi_ = 6.2832-calo_dijet_dphi_;
      calo_dijet_balance_ = (tcalojet_pt_-pcalojet_pt_)/(tcalojet_pt_+pcalojet_pt_);
      
      calo_tree_->Fill();
    }
  }
  
  if(doPFJets_){
    edm::Handle<reco::PFJetCollection> pfjets;
    iEvent.getByLabel(pfJetCollName_,pfjets);
    if(!pfjets.isValid()) {
      throw edm::Exception(edm::errors::ProductNotFound)
	<< " could not find PFJetCollection named " << pfJetCollName_ << ".\n";
      return;
    }

    edm::Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>> hbhereco;
    //iEvent.getByLabel(RecHitLabelName_,hbheRecHitInstance_,hbhereco);
    iEvent.getByLabel(hbheRecHitInstance_,hbhereco);
    if(!hbhereco.isValid()) {
      throw edm::Exception(edm::errors::ProductNotFound)
	<< " could not find HBHERecHit named " << RecHitLabelName_ << ":" << hbheRecHitInstance_ << ".\n";
      return;
    }
    
    edm::Handle<edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>> hfreco;
    //iEvent.getByLabel(RecHitLabelName_,hfRecHitInstance_,hfreco);
    iEvent.getByLabel(hfRecHitInstance_,hfreco);
    if(!hfreco.isValid()) {
      throw edm::Exception(edm::errors::ProductNotFound)
	<< " could not find HFRecHit named " << RecHitLabelName_ << ":" << hfRecHitInstance_ << ".\n";
      return;
    }

    edm::Handle<edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>> horeco;
    //iEvent.getByLabel(RecHitLabelName_,hoRecHitInstance_,horeco);
    iEvent.getByLabel(hoRecHitInstance_,horeco);
    if(!horeco.isValid()) {
      throw edm::Exception(edm::errors::ProductNotFound)
	<< " could not find HORecHit named " << RecHitLabelName_ << ":" << hoRecHitInstance_ << ".\n";
      return;
    }

    edm::Handle<edm::SortedCollection<CaloTower,edm::StrictWeakOrdering<CaloTower>>> towerMaker;
    iEvent.getByLabel("towerMaker",towerMaker);
    if(!horeco.isValid()) {
      throw edm::Exception(edm::errors::ProductNotFound)
	<< " could not find CaloTower named " << "towerMaker" << ".\n";
      return;
    }
    
    const JetCorrector* correctorPF = JetCorrector::getJetCorrector(pfJetCorrName_,evSetup);
    
    //////////////////////////////
    // Event Selection
    //////////////////////////////
    
    // determine which cut results in failure
    int passSelPF=0;

    // sort jets by corrected et
    std::set<PFJetCorretPair, PFJetCorretPairComp> pfjetcorretpairset;
    for(reco::PFJetCollection::const_iterator it=pfjets->begin(); it!=pfjets->end(); ++it) {
      const reco::PFJet* jet=&(*it);
      pfjetcorretpairset.insert( PFJetCorretPair(jet, correctorPF->correction(jet->p4())) );
    }

    PFJetCorretPair pf_tag, pf_probe;
    pf_thirdjet_px_=pf_thirdjet_py_=0.0;
    int cntr=0;
    for(std::set<PFJetCorretPair, PFJetCorretPairComp>::const_iterator it=pfjetcorretpairset.begin(); it!=pfjetcorretpairset.end(); ++it) {
      PFJetCorretPair jet=(*it);
      ++cntr;
      if(cntr==1) pf_tag=jet;
      else if(cntr==2) pf_probe=jet;
      else {
	pf_thirdjet_px_ += jet.scale()*jet.jet()->px();
	pf_thirdjet_py_ += jet.scale()*jet.jet()->py();
      }
    }
    
    if(pf_tag.jet() && pf_probe.jet()){
    // require that the first two jets are above some minimum,
    // and the rest are below some maximum
      if((pf_tag.jet()->et()+pf_probe.jet()->et())<minSumJetEt_) passSelPF |= 0x1;
      if(pf_tag.jet()->et()<minJetEt_ || pf_probe.jet()->et()<minJetEt_) passSelPF |= 0x2;
      if(sqrt(pf_thirdjet_px_*pf_thirdjet_px_ + pf_thirdjet_py_*pf_thirdjet_py_)>maxThirdJetEt_) passSelPF |= 0x4;
      
      // force the tag jet to have the smaller |eta|
      if(std::fabs(pf_tag.jet()->eta())>std::fabs(pf_probe.jet()->eta())) {
	PFJetCorretPair temp=pf_tag;
	pf_tag=pf_probe;
	pf_probe=temp;
      }
      
      // eta cuts
      double dAbsEta=std::fabs(std::fabs(pf_tag.jet()->eta())-std::fabs(pf_probe.jet()->eta()));
      if(dAbsEta>maxDeltaEta_) passSelPF |= 0x8;
      if(fabs(pf_tag.jet()->eta())<minTagJetEta_) passSelPF |= 0x10;
      if(fabs(pf_tag.jet()->eta())>maxTagJetEta_) passSelPF |= 0x10;
    }
    else{
      passSelPF = 0x40;
    }
    
    hPassSelPF_->Fill(passSelPF);
    if(!passSelPF){
      // dump
      if(debug_) {
	std::cout << "Run: " << iEvent.id().run() << "; Event: " << iEvent.id().event() << std::endl;
	for(reco::PFJetCollection::const_iterator it=pfjets->begin(); it!=pfjets->end(); ++it) {
	  const reco::PFJet *jet=&(*it);
	  std::cout << "istag=" << (jet==pf_tag.jet()) << "; isprobe=" << (jet==pf_probe.jet()) << "; et=" << jet->et() << "; eta=" << jet->eta() << std::endl;
	}
      }

      // Reset particle variables
      tpfjet_unkown_E_ = tpfjet_unkown_px_ = tpfjet_unkown_py_ = tpfjet_unkown_pz_ = 0.0;
      tpfjet_chHad_E_ = tpfjet_chHad_px_ = tpfjet_chHad_py_ = tpfjet_chHad_pz_ = 0.0;
      tpfjet_electron_E_ = tpfjet_electron_px_ = tpfjet_electron_py_ = tpfjet_electron_pz_ = 0.0;
      tpfjet_muon_E_ = tpfjet_muon_px_ = tpfjet_muon_py_ = tpfjet_muon_pz_ = 0.0;
      tpfjet_photon_E_ = tpfjet_photon_px_ = tpfjet_photon_py_ = tpfjet_photon_pz_ = 0.0;
      tpfjet_Had0_E_ = tpfjet_Had0_px_ = tpfjet_Had0_py_ = tpfjet_Had0_pz_ = 0.0;
      tpfjet_HFHad_E_ = tpfjet_HFHad_px_ = tpfjet_HFHad_py_ = tpfjet_HFHad_pz_ = 0.0;
      tpfjet_HFEM_E_ = tpfjet_HFEM_px_ = tpfjet_HFEM_py_ = tpfjet_HFEM_pz_ = 0.0;
      tpfjet_unkown_n_ = tpfjet_chHad_n_ = tpfjet_electron_n_ = tpfjet_muon_n_ = 0;
      tpfjet_photon_n_ = tpfjet_Had0_n_ = tpfjet_HFHad_n_ = tpfjet_HFEM_n_ = 0;
      tpfjet_chHad_EcalE_ = tpfjet_HFHad_EcalE_ = 0.0;
      ppfjet_unkown_E_ = ppfjet_unkown_px_ = ppfjet_unkown_py_ = ppfjet_unkown_pz_ = 0.0;
      ppfjet_chHad_E_ = ppfjet_chHad_px_ = ppfjet_chHad_py_ = ppfjet_chHad_pz_ = 0.0;
      ppfjet_electron_E_ = ppfjet_electron_px_ = ppfjet_electron_py_ = ppfjet_electron_pz_ = 0.0;
      ppfjet_muon_E_ = ppfjet_muon_px_ = ppfjet_muon_py_ = ppfjet_muon_pz_ = 0.0;
      ppfjet_photon_E_ = ppfjet_photon_px_ = ppfjet_photon_py_ = ppfjet_photon_pz_ = 0.0;
      ppfjet_Had0_E_ = ppfjet_Had0_px_ = ppfjet_Had0_py_ = ppfjet_Had0_pz_ = 0.0;
      ppfjet_HFHad_E_ = ppfjet_HFHad_px_ = ppfjet_HFHad_py_ = ppfjet_HFHad_pz_ = 0.0;
      ppfjet_HFEM_E_ = ppfjet_HFEM_px_ = ppfjet_HFEM_py_ = ppfjet_HFEM_pz_ = 0.0;
      ppfjet_unkown_n_ = ppfjet_chHad_n_ = ppfjet_electron_n_ = ppfjet_muon_n_ = 0;
      ppfjet_photon_n_ = ppfjet_Had0_n_ = ppfjet_HFHad_n_ = ppfjet_HFEM_n_ = 0;
      ppfjet_chHad_EcalE_ = ppfjet_HFHad_EcalE_ = 0.0;
      
      
      // fill tag jet variables
      tpfjet_pt_    = pf_tag.jet()->pt();
      tpfjet_p_     = pf_tag.jet()->p();
      tpfjet_eta_   = pf_tag.jet()->eta();
      tpfjet_phi_   = pf_tag.jet()->phi();
      //tpfjet_photonE_   = pf_tag.jet()->photonEnergy();
      //tpfjet_photonN_   = pf_tag.jet()->photonMultiplicity();
      //tpfjet_electronE_ = pf_tag.jet()->electronEnergy();
      //tpfjet_electronN_ = pf_tag.jet()->electronMultiplicity();
      //tpfjet_muonE_     = pf_tag.jet()->muonEnergy();
      //tpfjet_muonN_     = pf_tag.jet()->muonMultiplicity();
      //tpfjet_HFEME_     = pf_tag.jet()->HFEMEnergy();
      //tpfjet_HFEMN_     = pf_tag.jet()->HFEMMultiplicity();
      tpfjet_scale_ = pf_tag.scale();
      tpfjet_ntwrs_=0;

      //std::cout << pf_tag.jet()->print() << std::endl;
      //bool isMatchedt = false;
      int types = 0;
      int ntypes = 0;
      
      // Get PF constituents and fill HCAL towers
      std::vector<reco::PFCandidatePtr> tagconst=pf_tag.jet()->getPFConstituents();
      if(iEvent.id().event() == testEvent){
	std::cout << "do PFCandidatePtr loop" << std::endl;
      }
      for(std::vector<reco::PFCandidatePtr>::const_iterator it=tagconst.begin(); it!=tagconst.end(); ++it){
	//std::cout << (*it) << std::endl;
	reco::PFCandidate::ParticleType candidateType = (*it)->particleId();
	if(candidateType == 0){
	  tpfjet_unkown_E_ += (*it)->energy();
	  tpfjet_unkown_px_ += (*it)->px();
	  tpfjet_unkown_py_ += (*it)->py();
	  tpfjet_unkown_pz_ += (*it)->pz();
	  tpfjet_unkown_n_++;
	  continue;
	}
	if(candidateType == 1){
	  tpfjet_chHad_E_ += (*it)->energy();
	  tpfjet_chHad_px_ += (*it)->px();
	  tpfjet_chHad_py_ += (*it)->py();
	  tpfjet_chHad_pz_ += (*it)->pz();
	  tpfjet_chHad_EcalE_ += (*it)->ecalEnergy();
	  tpfjet_chHad_n_++;

	  /*reco::TrackRef trackRef = (*it)->trackRef();
	  if(trackRef.isNonnull()){
	    //reco::TrackCollection tracks = *trackRef;
	    reco::Track track = *trackRef;
	    float track_p = track.p();
	    float track_E = sqrt(track_p*track_p + 0.13957*0.13957);
	    float track_E_m = sqrt(track_p*track_p + (*it)->mass()*(*it)->mass());
	    std::cout << "Candidate E: " << (*it)->energy() << " Track E: " << track_E << " Candidate mass: " << (*it)->mass() << " Track E with candidate mass: " << track_E_m << std::endl;
	  }
	  else{
	    std::cout << "NULL" << std::endl;
	  }*/
	}
	if(candidateType == 2){
	  tpfjet_electron_E_ += (*it)->energy();
	  tpfjet_electron_px_ += (*it)->px();
	  tpfjet_electron_py_ += (*it)->py();
	  tpfjet_electron_pz_ += (*it)->pz();
	  tpfjet_electron_n_++;
	  continue;
	}
	if(candidateType == 3){
	  tpfjet_muon_E_ += (*it)->energy();
	  tpfjet_muon_px_ += (*it)->px();
	  tpfjet_muon_py_ += (*it)->py();
	  tpfjet_muon_pz_ += (*it)->pz();
	  tpfjet_muon_n_++;
	  continue;
	}
	if(candidateType == 4){
	  tpfjet_photon_E_ += (*it)->energy();
	  tpfjet_photon_px_ += (*it)->px();
	  tpfjet_photon_py_ += (*it)->py();
	  tpfjet_photon_pz_ += (*it)->pz();
	  tpfjet_photon_n_++;
	  continue;
	}
	if(candidateType == 5){
	  tpfjet_Had0_E_ += (*it)->energy();
	  tpfjet_Had0_px_ += (*it)->px();
	  tpfjet_Had0_py_ += (*it)->py();
	  tpfjet_Had0_pz_ += (*it)->pz();
	  tpfjet_Had0_n_++;
	}
	if(candidateType == 6){
	  tpfjet_HFHad_E_ += (*it)->energy();
	  tpfjet_HFHad_px_ += (*it)->px();
	  tpfjet_HFHad_py_ += (*it)->py();
	  tpfjet_HFHad_pz_ += (*it)->pz();
	  tpfjet_HFHad_EcalE_ += (*it)->ecalEnergy();
	  tpfjet_HFHad_n_++;

	  reco::TrackRef trackRef = (*it)->trackRef();
	  if(trackRef.isNonnull()){
	    //reco::TrackCollection tracks = *trackRef;
	    reco::Track track = *trackRef;
	    float track_p = track.p();
	    float track_E = sqrt(track_p*track_p + 0.13957*0.13957);
	    float track_E_m = sqrt(track_p*track_p + (*it)->mass()*(*it)->mass());
	    std::cout << "Candidate E: " << (*it)->energy() << " Track E: " << track_E << " Candidate mass: " << (*it)->mass() << " Track E with candidate mass: " << track_E_m << std::endl;
	  }
	  else{
	    //std::cout << "NULL" << std::endl;
	  }
	}
	if(candidateType == 7){
	  tpfjet_HFEM_E_ += (*it)->energy();
	  tpfjet_HFEM_px_ += (*it)->px();
	  tpfjet_HFEM_py_ += (*it)->py();
	  tpfjet_HFEM_pz_ += (*it)->pz();
	  tpfjet_HFEM_n_++;
	}
	
	//std::cout << candidateType << ": " << (*it)->energy() << " " << (*it)->px() << " " << (*it)->py() << " " << (*it)->pz() << " m= " << (*it)->mass() << " m2= " << (*it)->massSqr() << " " << sqrt((*it)->energy()*(*it)->energy() - (*it)->px()*(*it)->px() - (*it)->py()*(*it)->py() - (*it)->pz()*(*it)->pz()) << std::endl;
	
	int maxElement=(*it)->elementsInBlocks().size();
	if(iEvent.id().event() == testEvent){
	  std::cout << "do e loop" << std::endl;
	}
	for(int e=0; e<maxElement; ++e){
	  reco::PFBlockRef blockRef = (*it)->elementsInBlocks()[e].first;
	  const edm::OwnVector<reco::PFBlockElement>& elements = blockRef->elements();
	  if(iEvent.id().event() == testEvent){
	    std::cout << "do iEle loop" << std::endl;
	  }
	  for(unsigned iEle=0; iEle<elements.size(); iEle++) {
	    reco::PFClusterRef tmpclusterref = elements[iEle].clusterRef();
	    if(elements[iEle].index() == (*it)->elementsInBlocks()[e].second){
	      if(iEvent.id().event() == testEvent){
		std::cout << "passed event == testEvent" << std::endl;
		std::cout << "=====Type " << elements[iEle].type() << std::endl;
	      }

	      ////////////////////////////////////////
	      //
	      // This is the section with the matching issues.
	      // I'm not sure exactly what these PFBlockElement types are
	      // because type HO will give clusters with ieta > 15.
	      // And there's the question of what rechits to use for
	      // PFBlockElement types HFHAD and HFEM.
	      //
	      ////////////////////////////////////////

	      if(elements[iEle].type() == reco::PFBlockElement::HCAL){
		types |= 0x1;
		ntypes++;
		reco::PFClusterRef clusterref = elements[iEle].clusterRef();
		reco::PFCluster cluster = *clusterref;

		//	std::vector<reco::PFRecHitFraction> pfRecHits = cluster.recHitFractions();

		//int NpfRecHit = pfRecHits.size();
		//for(int jkl; jkl<NpfRecHit; jkl++){
		/*
		for(std::vector<reco::PFRecHitFraction>::const_iterator itRecHits=pfRecHits.begin(); itRecHits!=pfRecHits.end(); ++itRecHits){
		  //reco::PFRecHitRef pfRecHitRef = (*itRecHits).recHitRef();
		  //reco::PFRecHit pfRecHit = *pfRecHitRef;
		  std::cerr << itRecHits->recHitRef()->layer() << std::endl;
		}
		*/
		std::vector<std::pair<DetId,float>> hitsAndFracs = cluster.hitsAndFractions();
		if(iEvent.id().event() == testEvent){
		  std::cout << "Tag hits" << std::endl;
		}
		int nHits = hitsAndFracs.size();
		for(int iHit=0; iHit<nHits; iHit++){
		  HcalDetId hDet(hitsAndFracs[iHit].first.rawId());
		  int etaPhiPF = hitsAndFracs[iHit].first.rawId() & 0x3FFF;

		  /*int tmpzside = ((hitsAndFracs[iHit].first.rawId() >> 13) & 0x1) ? 1 : -1;
		  int tmpieta = ((hitsAndFracs[iHit].first.rawId() >> 7) & 0x3F);
		  if(tmpieta == 29){
		    std::cout << iEvent.id().event() << ": " << tmpzside*tmpieta << std::endl;
		  }
		  else if(tmpieta > 29){
		    std::cout << iEvent.id().event() << ": " << tmpzside*tmpieta << " ==================" << std::endl;
		  }*/

		  if(iEvent.id().event() == testEvent){
		    std::cout << hitsAndFracs[iHit].first.rawId() <<  " ";
		    for(int ijk=31; ijk>=0; ijk--){
		      int x = (hitsAndFracs[iHit].first.rawId() >> ijk) & 0x1;
		      std::cout << x;
		    }
		    std::cout << " " << (hitsAndFracs[iHit].first.rawId() >> 28) << " " << ((hitsAndFracs[iHit].first.rawId() >> 25) & 0x7) << " " << hDet << " " << etaPhiPF <<std::endl;
		  }

		  for(edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>::const_iterator ith=hbhereco->begin(); ith!=hbhereco->end(); ++ith){
		    int etaPhiRecHit = (*ith).id().rawId() & 0x3FFF;
		    /*if(etaPhiPF == 3587){
		      std::cout << "   +++ " << etaPhiPF << " " << etaPhiRecHit << std::endl;
		      if(etaPhiRecHit == 3587){
			std::cout << "     -----" << etaPhiPF << " " << etaPhiRecHit << " ";
			if(etaPhiPF == etaPhiRecHit){
			  std::cout << "1";
			}else{
			  std::cout << "0";
			}
			std::cout << std::endl;
		      }
		      }*/
		    if(etaPhiPF == etaPhiRecHit){
		      //std::cout << "    +-----" << etaPhiPF << " " << etaPhiRecHit << " " << std::endl;
		      //isMatchedt = true;
		      //std::cout<<"Matches Rechit SD "<<(*ith).subdetId()<<std::endl;
		      //if((*ith).energy() > 0.0 && hitsAndFracs[iHit].second > 0.0){
		      //std::cout << "HBHE phi=" << (*it)->phi() << " iphi=" << (hitsAndFracs[iHit].first.rawId() & 0x3F) << std::endl;
		      if(true){
			tpfjet_twr_ieta_[tpfjet_ntwrs_] = (*ith).id().ieta();
			//int tmpieta = (*ith).id().ieta();
			//std::cout << tmpieta << std::endl;
			tpfjet_twr_hade_[tpfjet_ntwrs_] = (*ith).energy();
			tpfjet_twr_frac_[tpfjet_ntwrs_] = hitsAndFracs[iHit].second;
			++tpfjet_ntwrs_;
		      }
		    }
		  }

		  for(edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>::const_iterator ith=hfreco->begin(); ith!=hfreco->end(); ++ith){
		    int etaPhiRecHit = (*ith).id().rawId() & 0x3FFF;
		    if(etaPhiPF == etaPhiRecHit){
		      //isMatchedt = true;
		      //std::cout << "HF" << std::endl;
		      //if((*ith).energy() > 0.0 && hitsAndFracs[iHit].second > 0.0){
		      if(true){
			tpfjet_twr_ieta_[tpfjet_ntwrs_] = (*ith).id().ieta();
			tpfjet_twr_hade_[tpfjet_ntwrs_] = (*ith).energy();
			tpfjet_twr_frac_[tpfjet_ntwrs_] = hitsAndFracs[iHit].second;
			++tpfjet_ntwrs_;
		      }
		    }
		  }

		  for(edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>::const_iterator ith=horeco->begin(); ith!=horeco->end(); ++ith){
		    int etaPhiRecHit = (*ith).id().rawId() & 0x3FFF;
		    if(etaPhiPF == etaPhiRecHit){
		      //isMatchedt = true;
		      //std::cout << "HO" << std::endl;
		      //if((*ith).energy() > 0.0 && hitsAndFracs[iHit].second > 0.0){
		      if(true){
			tpfjet_twr_ieta_[tpfjet_ntwrs_] = (*ith).id().ieta();
			tpfjet_twr_hade_[tpfjet_ntwrs_] = (*ith).energy();
			tpfjet_twr_frac_[tpfjet_ntwrs_] = hitsAndFracs[iHit].second;
			++tpfjet_ntwrs_;
		      }
		    }
		    
		    /*int tmpzside = ((hitsAndFracs[iHit].first.rawId() >> 13) & 0x1) ? 1 : -1;
		    int tmpieta = ((hitsAndFracs[iHit].first.rawId() >> 7) & 0x3F);
		    if(tmpieta > 15){
		      std::cout << iEvent.id().event() << ": " << tmpzside*tmpieta << " " << hitsAndFracs[iHit].first.rawId() << std::endl;
		      }*/
		  }
		} // Loop over hits
	      } // Test if element is from HCAL
	      else if(elements[iEle].type() == reco::PFBlockElement::HFHAD){
		types |= 0x2;
		ntypes++;
		//reco::PFClusterRef clusterref = elements[iEle].clusterRef();
		//reco::PFCluster cluster = *clusterref;

		/*int candieta = 63;
		if(fabs((*it)->eta()) > hfieta[0]){
		  for(int itieta = 1; itieta<13; itieta++){
		    if(fabs((*it)->eta()) < hfieta[itieta]){
		      candieta = itieta + 29;
		      if((*it)->eta() > 0.0){
			candieta |= 0x40;
		      }
		    }
		  }
		}
		int candiphi;
		if((*it)->phi() < 0.0){
		  candiphi = floor((6.28319 - (*it)->phi())/0.087);
		  if(candiphi % 2 == 0){
		    candiphi++;
		  }
		}
		else{
		  candiphi = floor((*it)->phi()/0.087);
		  if(candiphi % 2 ==0){
		    candiphi++;
		  }
		}

		int etaPhiPF = (candieta << 7) | candiphi;

		for(edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>::const_iterator ith=hfreco->begin(); ith!=hfreco->end(); ++ith){
		  int etaPhiRecHit = (*ith).id().rawId() & 0x3FFF;
		  if(etaPhiPF == etaPhiRecHit){
		    if(true){
		      tpfjet_twr_ieta_[tpfjet_ntwrs_] = (*ith).id().ieta();
		      tpfjet_twr_hade_[tpfjet_ntwrs_] = (*ith).energy();
		      //tpfjet_twr_frac_[tpfjet_ntwrs_] = hitsAndFracs[iHit].second;
		      ++tpfjet_ntwrs_;
		      std::cout << "candidate energy: " << (*it)->energy() << " RecHit energy: " << (*ith).energy() << std::endl;
		    }
		  }
		  }*/

		//std::cout << "HF phi=" << (*it)->phi() << " iphi=" << 
	      }
	      else if(elements[iEle].type() == reco::PFBlockElement::HFEM){
		types |= 0x4;
		ntypes++;
		reco::PFClusterRef clusterref = elements[iEle].clusterRef();
		//reco::PFCluster cluster = *clusterref;
	      }
	      else if(elements[iEle].type() == reco::PFBlockElement::HO){
		types |= 0x8;
		ntypes++;
		reco::PFClusterRef clusterref = elements[iEle].clusterRef();
		reco::PFCluster cluster = *clusterref;

		std::vector<std::pair<DetId,float>> hitsAndFracs = cluster.hitsAndFractions();
		if(iEvent.id().event() == testEvent){
		  std::cout << "HO Tag hits" << std::endl;
		}
		int nHits = hitsAndFracs.size();
		for(int iHit=0; iHit<nHits; iHit++){
		  HcalDetId hDet(hitsAndFracs[iHit].first.rawId());
		  int etaPhiPF = hitsAndFracs[iHit].first.rawId() & 0x3FFF;

		  if(iEvent.id().event() == testEvent){
		    std::cout << hitsAndFracs[iHit].first.rawId() <<  " ";
		    for(int ijk=31; ijk>=0; ijk--){
		      int x = (hitsAndFracs[iHit].first.rawId() >> ijk) & 0x1;
		      std::cout << x;
		    }
		    std::cout << " " << (hitsAndFracs[iHit].first.rawId() >> 28) << " " << ((hitsAndFracs[iHit].first.rawId() >> 25) & 0x7) << " " << hDet << " " << etaPhiPF <<std::endl;
		  }

		  for(edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>::const_iterator ith=hbhereco->begin(); ith!=hbhereco->end(); ++ith){
		    int etaPhiRecHit = (*ith).id().rawId() & 0x3FFF;
		    if(etaPhiPF == etaPhiRecHit){
		      //isMatchedt = true;
		      if(true){
			tpfjet_twr_ieta_[tpfjet_ntwrs_] = (*ith).id().ieta();
			tpfjet_twr_hade_[tpfjet_ntwrs_] = (*ith).energy();
			tpfjet_twr_frac_[tpfjet_ntwrs_] = hitsAndFracs[iHit].second;
			++tpfjet_ntwrs_;
		      }
		    }
		  }

		  for(edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>::const_iterator ith=hfreco->begin(); ith!=hfreco->end(); ++ith){
		    int etaPhiRecHit = (*ith).id().rawId() & 0x3FFF;
		    if(etaPhiPF == etaPhiRecHit){
		      //isMatchedt = true;
		      if(true){
			tpfjet_twr_ieta_[tpfjet_ntwrs_] = (*ith).id().ieta();
			tpfjet_twr_hade_[tpfjet_ntwrs_] = (*ith).energy();
			tpfjet_twr_frac_[tpfjet_ntwrs_] = hitsAndFracs[iHit].second;
			++tpfjet_ntwrs_;
		      }
		    }
		  }

		  for(edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>::const_iterator ith=horeco->begin(); ith!=horeco->end(); ++ith){
		    int etaPhiRecHit = (*ith).id().rawId() & 0x3FFF;
		    if(etaPhiPF == etaPhiRecHit){
		      //isMatchedt = true;
		      if(true){
			tpfjet_twr_ieta_[tpfjet_ntwrs_] = (*ith).id().ieta();
			tpfjet_twr_hade_[tpfjet_ntwrs_] = (*ith).energy();
			tpfjet_twr_frac_[tpfjet_ntwrs_] = hitsAndFracs[iHit].second;
			++tpfjet_ntwrs_;
		      }
		    }
		  }
		} // Loop over hits
	      }
	    } // Test for right element index
	  } // Loop over elements
	} // Loop over elements in blocks
      } // Loop over PF constitutents

      h_types_->Fill(types);
      h_ntypes_->Fill(ntypes);
      
      // fill probe jet variables
      ppfjet_pt_    = pf_probe.jet()->pt();
      ppfjet_p_     = pf_probe.jet()->p();
      ppfjet_eta_   = pf_probe.jet()->eta();
      ppfjet_phi_   = pf_probe.jet()->phi();
      //ppfjet_photonE_   = pf_probe.jet()->photonEnergy();
      //ppfjet_photonN_   = pf_probe.jet()->photonMultiplicity();
      //ppfjet_electronE_ = pf_probe.jet()->electronEnergy();
      //ppfjet_electronN_ = pf_probe.jet()->electronMultiplicity();
      //ppfjet_muonE_     = pf_probe.jet()->muonEnergy();
      //ppfjet_muonN_     = pf_probe.jet()->muonMultiplicity();
      //ppfjet_HFEME_     = pf_probe.jet()->HFEMEnergy();
      //ppfjet_HFEMN_     = pf_probe.jet()->HFEMMultiplicity();
      ppfjet_scale_ = pf_probe.scale();
      ppfjet_ntwrs_=0;

      //bool isMatchedp = false;

      // Get PF constituents and fill HCAL towers
      std::vector<reco::PFCandidatePtr> probeconst=pf_probe.jet()->getPFConstituents();
      for(std::vector<reco::PFCandidatePtr>::const_iterator it=probeconst.begin(); it!=probeconst.end(); ++it){
	reco::PFCandidate::ParticleType candidateType = (*it)->particleId();
	if(candidateType == 0){
	  ppfjet_unkown_E_ += (*it)->energy();
	  ppfjet_unkown_px_ += (*it)->px();
	  ppfjet_unkown_py_ += (*it)->py();
	  ppfjet_unkown_pz_ += (*it)->pz();
	  ppfjet_unkown_n_++;
	  continue;
	}
	if(candidateType == 1){
	  ppfjet_chHad_E_ += (*it)->energy();
	  ppfjet_chHad_px_ += (*it)->px();
	  ppfjet_chHad_py_ += (*it)->py();
	  ppfjet_chHad_pz_ += (*it)->pz();
	  ppfjet_chHad_EcalE_ += (*it)->ecalEnergy();
	  ppfjet_chHad_n_++;
	}
	if(candidateType == 2){
	  ppfjet_electron_E_ += (*it)->energy();
	  ppfjet_electron_px_ += (*it)->px();
	  ppfjet_electron_py_ += (*it)->py();
	  ppfjet_electron_pz_ += (*it)->pz();
	  ppfjet_electron_n_++;
	  continue;
	}
	if(candidateType == 3){
	  ppfjet_muon_E_ += (*it)->energy();
	  ppfjet_muon_px_ += (*it)->px();
	  ppfjet_muon_py_ += (*it)->py();
	  ppfjet_muon_pz_ += (*it)->pz();
	  ppfjet_muon_n_++;
	  continue;
	}
	if(candidateType == 4){
	  ppfjet_photon_E_ += (*it)->energy();
	  ppfjet_photon_px_ += (*it)->px();
	  ppfjet_photon_py_ += (*it)->py();
	  ppfjet_photon_pz_ += (*it)->pz();
	  ppfjet_photon_n_++;
	  continue;
	}
	if(candidateType == 5){
	  ppfjet_Had0_E_ += (*it)->energy();
	  ppfjet_Had0_px_ += (*it)->px();
	  ppfjet_Had0_py_ += (*it)->py();
	  ppfjet_Had0_pz_ += (*it)->pz();
	  ppfjet_Had0_n_++;
	}
	if(candidateType == 6){
	  ppfjet_HFHad_E_ += (*it)->energy();
	  ppfjet_HFHad_px_ += (*it)->px();
	  ppfjet_HFHad_py_ += (*it)->py();
	  ppfjet_HFHad_pz_ += (*it)->pz();
	  ppfjet_HFHad_EcalE_ += (*it)->ecalEnergy();
	  ppfjet_HFHad_n_++;
	}
	if(candidateType == 7){
	  ppfjet_HFEM_E_ += (*it)->energy();
	  ppfjet_HFEM_px_ += (*it)->px();
	  ppfjet_HFEM_py_ += (*it)->py();
	  ppfjet_HFEM_pz_ += (*it)->pz();
	  ppfjet_HFEM_n_++;
	}
	
	int maxElement=(*it)->elementsInBlocks().size();
	for(int e=0; e<maxElement; ++e){
	  reco::PFBlockRef blockRef = (*it)->elementsInBlocks()[e].first;
	  const edm::OwnVector<reco::PFBlockElement>& elements = blockRef->elements();
	  for(unsigned iEle=0; iEle<elements.size(); iEle++) {
	    if(elements[iEle].index() == (*it)->elementsInBlocks()[e].second){
	      if(elements[iEle].type() == reco::PFBlockElement::HCAL){
		reco::PFClusterRef clusterref = elements[iEle].clusterRef();
		reco::PFCluster cluster = *clusterref;
		
		std::vector<std::pair<DetId,float>> hitsAndFracs = cluster.hitsAndFractions();
		if(iEvent.id().event() == testEvent){
		  std::cout << "Probe hits" << std::endl;
		}
		int nHits = hitsAndFracs.size();
		for(int iHit=0; iHit<nHits; iHit++){
		  HcalDetId hDet(hitsAndFracs[iHit].first.rawId());
		  int etaPhiPF = hitsAndFracs[iHit].first.rawId() & 0x3FFF;
		  
		  if(iEvent.id().event() == testEvent){
		    std::cout << hitsAndFracs[iHit].first.rawId() <<  " ";
		    for(int ijk=31; ijk>=0; ijk--){
		      int x = (hitsAndFracs[iHit].first.rawId() >> ijk) & 0x1;
		      std::cout << x;
		    }
		    std::cout << " " << (hitsAndFracs[iHit].first.rawId() >> 28) << " " << ((hitsAndFracs[iHit].first.rawId() >> 25) & 0x7) << " " << hDet << " " << etaPhiPF << std::endl;
		  }

		  for(edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>::const_iterator ith=hbhereco->begin(); ith!=hbhereco->end(); ++ith){
		    int etaPhiRecHit = (*ith).id().rawId() & 0x3FFF;
		    if(etaPhiPF == etaPhiRecHit){
		      //isMatchedp = true;
		      //if((*ith).energy() > 0.0 && hitsAndFracs[iHit].second > 0.0){
		      if(true){
			ppfjet_twr_ieta_[ppfjet_ntwrs_] = (*ith).id().ieta();
			//std::cout << ppfjet_twr_ieta_[ppfjet_ntwrs_] << std::endl;
			ppfjet_twr_hade_[ppfjet_ntwrs_] = (*ith).energy();
			ppfjet_twr_frac_[ppfjet_ntwrs_] = hitsAndFracs[iHit].second;
			++ppfjet_ntwrs_;
		      }
		    }
		  }

		  for(edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>::const_iterator ith=hfreco->begin(); ith!=hfreco->end(); ++ith){
		    int etaPhiRecHit = (*ith).id().rawId() & 0x3FFF;
		    if(etaPhiPF == etaPhiRecHit){
		      //isMatchedp = true;
		      //if((*ith).energy() > 0.0 && hitsAndFracs[iHit].second > 0.0){
		      if(true){
			ppfjet_twr_ieta_[ppfjet_ntwrs_] = (*ith).id().ieta();
			//int tmpieta = (*ith).id().ieta();
			//std::cout << ppfjet_twr_ieta_[ppfjet_ntwrs_] << std::endl;
			ppfjet_twr_hade_[ppfjet_ntwrs_] = (*ith).energy();
			ppfjet_twr_frac_[ppfjet_ntwrs_] = hitsAndFracs[iHit].second;
			++ppfjet_ntwrs_;
		      }
		    }
		  }

		  for(edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>::const_iterator ith=horeco->begin(); ith!=horeco->end(); ++ith){
		    int etaPhiRecHit = (*ith).id().rawId() & 0x3FFF;
		    if(etaPhiPF == etaPhiRecHit){
		      //isMatchedp = true;
		      //if((*ith).energy() > 0.0 && hitsAndFracs[iHit].second > 0.0){
		      if(true){
			ppfjet_twr_ieta_[ppfjet_ntwrs_] = (*ith).id().ieta();
			//int tmpieta = (*ith).id().ieta();
			//std::cout << ppfjet_twr_ieta_[ppfjet_ntwrs_] << std::endl;
			ppfjet_twr_hade_[ppfjet_ntwrs_] = (*ith).energy();
			ppfjet_twr_frac_[ppfjet_ntwrs_] = hitsAndFracs[iHit].second;
			++ppfjet_ntwrs_;
		      }
		    }
		  }
		} // Loop over hits
	      } // Test if element is from HCAL
	    } // Test for right element index
	  } // Loop over elements
	} // Loop over elements in blocks
      } // Loop over PF constitutents

      //std::cout << "Matched " << isMatchedt << " " << isMatchedp << std::endl;
      
      // fill genjet tag/probe variables
      tpfjet_gendr_ = 99999.;
      tpfjet_genpt_ = 0;
      tpfjet_genp_  = 0;
      ppfjet_gendr_ = 99999.;
      ppfjet_genpt_ = 0;
      ppfjet_genp_  = 0;
      for(std::vector<reco::GenJet>::const_iterator it=genjets->begin(); it!=genjets->end(); ++it){
	const reco::GenJet* jet=&(*it);
	double dr=deltaR(jet, pf_probe.jet());
	if(dr<ppfjet_gendr_) {
	  ppfjet_gendr_ = dr;
	  ppfjet_genpt_ = jet->pt();
	  ppfjet_genp_ = jet->p();
	}
	dr=deltaR(jet, pf_tag.jet());
	if(dr<tpfjet_gendr_) {
	  tpfjet_gendr_ = dr;
	  tpfjet_genpt_ = jet->pt();
	  tpfjet_genp_ = jet->p();
	}
      }
      /*~~~for(reco::GenJetCollection::const_iterator it=genjets->begin(); it!=genjets->end(); ++it) {
	const reco::GenJet* jet=&(*it);
	double dr=deltaR(jet, pf_probe.jet());
	if(dr<ppfjet_gendr_) {
	  ppfjet_gendr_ = dr;
	  ppfjet_genpt_ = jet->pt();
	  ppfjet_genp_ = jet->p();
	}
	dr=deltaR(jet, pf_tag.jet());
	if(dr<tpfjet_gendr_) {
	  tpfjet_gendr_ = dr;
	  tpfjet_genpt_ = jet->pt();
	  tpfjet_genp_ = jet->p();
	}
	}*/
      
      // fill dijet variables
      pf_dijet_deta_=std::fabs(std::fabs(pf_tag.jet()->eta())-std::fabs(pf_probe.jet()->eta()));
      pf_dijet_dphi_=pf_tag.jet()->phi()-pf_probe.jet()->phi();
      if(pf_dijet_dphi_>3.1415) pf_dijet_dphi_ = 6.2832-pf_dijet_dphi_;
      pf_dijet_balance_ = (tpfjet_pt_-ppfjet_pt_)/(tpfjet_pt_+ppfjet_pt_);

      if(iEvent.id().event() == testEvent){
	std::cout << "RecHits" << std::endl;
	std::cout << "HBHE" << std::endl;
	for(edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>::const_iterator it=hbhereco->begin(); it!=hbhereco->end(); ++it){
	  std::cout << (*it).id().rawId() <<  " ";
	  for(int ijk=31; ijk>=0; ijk--){
	    int x = ((*it).id().rawId() >> ijk) & 0x1;
	    std::cout << x;
	  }
	  std::cout << " " << ((*it).id().rawId() >> 28) << " " << (((*it).id().rawId() >> 25) & 0x7) << " " << (*it).id() << " " << ((*it).id().rawId() & 0x3FFF) << std::endl;
	}
	std::cout << "HF" << std::endl;
	for(edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>::const_iterator it=hfreco->begin(); it!=hfreco->end(); ++it){
	  std::cout << (*it).id().rawId() <<  " ";
	  for(int ijk=31; ijk>=0; ijk--){
	    int x = ((*it).id().rawId() >> ijk) & 0x1;
	    std::cout << x;
	  }
	  std::cout << " " << ((*it).id().rawId() >> 28) << " " << (((*it).id().rawId() >> 25) & 0x7) << " " << (*it).id() << " " << ((*it).id().rawId() & 0x3FFF) << std::endl;
	}
	std::cout << "HO" << std::endl;
	for(edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>::const_iterator it=horeco->begin(); it!=horeco->end(); ++it){
	  std::cout << (*it).id().rawId() <<  " ";
	  for(int ijk=31; ijk>=0; ijk--){
	    int x = ((*it).id().rawId() >> ijk) & 0x1;
	    std::cout << x;
	  }
	  std::cout << " " << ((*it).id().rawId() >> 28) << " " << (((*it).id().rawId() >> 25) & 0x7) << " " << (*it).id() << " " << ((*it).id().rawId() & 0x3FFF) << std::endl;
	}
      }


      /*std::cout << "towerMaker" << std::endl;
      for(edm::SortedCollection<CaloTower,edm::StrictWeakOrdering<CaloTower>>::const_iterator it=towerMaker->begin(); it!=towerMaker->end(); ++it){
	std::cout << (*it).id() << " " << (*it).constituentsSize() << std::endl;
	std::vector<DetId> caloTowerConstituents = (*it).constituents();
	for(unsigned int i=0; i<(*it).constituentsSize(); i++){
	  if((caloTowerConstituents[i].rawId() >> 28) == 4){
	    std::cout << caloTowerConstituents[i].rawId() << " ";
	    for(int ijk=31; ijk>=0; ijk--){
	      int x = (caloTowerConstituents[i].rawId() >> ijk) & 0x1;
	      std::cout << x;
	    }
	    int tmpzside = ((caloTowerConstituents[i].rawId() >> 13) & 0x1) ? 1 : -1;
	    int tmpieta = ((caloTowerConstituents[i].rawId() >> 7) & 0x3F);
	    std::cout << " " << (caloTowerConstituents[i].rawId() >> 28) << " " << ((caloTowerConstituents[i].rawId() >> 25) & 0x7) << " " << (caloTowerConstituents[i].rawId() & 0x3FFF) << " (" << tmpzside*tmpieta << "," << (caloTowerConstituents[i].rawId() & 0x7F) << ")" << std::endl; 
	  }
	}
	}*/


      for(edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>::const_iterator it=hfreco->begin(); it!=hfreco->end(); ++it){
	int tmpiphi = ((*it).id().rawId() & 0x7F);
	if(tmpiphi % 2 == 0){
	  std::cout << "HF iphi " << tmpiphi << std::endl;
	}
	if(tmpiphi > maxhfiphi){
	  maxhfiphi = tmpiphi;
	}
      }
      
      pf_tree_->Fill();
    }
  }
  
  return;
}

// ------------ method called once each job just before starting event loop  ------------
void CalcRespCorrDiJets::beginJob()
{
  // book histograms
  rootfile_ = new TFile(rootHistFilename_.c_str(), "RECREATE");

  if(doCaloJets_){
    hPassSelCalo_ = new TH1D("hPassSelectionCalo", "Selection Pass Failures CaloJets",200,-0.5,199.5);
    
    calo_tree_ = new TTree("calo_dijettree", "tree for dijet balancing using CaloJets");
    
    calo_tree_->Branch("tcalojet_pt",&tcalojet_pt_, "tcalojet_pt/F");
    calo_tree_->Branch("tcalojet_p",&tcalojet_p_, "tcalojet_p/F");
    calo_tree_->Branch("tcalojet_eta",&tcalojet_eta_, "tcalojet_eta/F");
    calo_tree_->Branch("tcalojet_phi",&tcalojet_phi_, "tcalojet_phi/F");
    calo_tree_->Branch("tcalojet_emf",&tcalojet_emf_, "tcalojet_emf/F");
    calo_tree_->Branch("tcalojet_scale",&tcalojet_scale_, "tcalojet_scale/F");
    calo_tree_->Branch("tcalojet_genpt",&tcalojet_genpt_, "tcalojet_genpt/F");
    calo_tree_->Branch("tcalojet_genp",&tcalojet_genp_, "tcalojet_genp/F");
    calo_tree_->Branch("tcalojet_gendr",&tcalojet_gendr_, "tcalojet_gendr/F");
    calo_tree_->Branch("tcalojet_EBE",&tcalojet_EBE_, "tcalojet_EBE/F");
    calo_tree_->Branch("tcalojet_EEE",&tcalojet_EEE_, "tcalojet_EEE/F");
    calo_tree_->Branch("tcalojet_HBE",&tcalojet_HBE_, "tcalojet_HBE/F");
    calo_tree_->Branch("tcalojet_HEE",&tcalojet_HEE_, "tcalojet_HEE/F");
    calo_tree_->Branch("tcalojet_HFE",&tcalojet_HFE_, "tcalojet_HFE/F");
    calo_tree_->Branch("tcalojet_ntwrs",&tcalojet_ntwrs_, "tcalojet_ntwrs/I");
    calo_tree_->Branch("tcalojet_twr_ieta",tcalojet_twr_ieta_, "tcalojet_twr_ieta[tcalojet_ntwrs]/I");
    calo_tree_->Branch("tcalojet_twr_eme",tcalojet_twr_eme_, "tcalojet_twr_eme[tcalojet_ntwrs]/F");
    calo_tree_->Branch("tcalojet_twr_hade",tcalojet_twr_hade_, "tcalojet_twr_hade[tcalojet_ntwrs]/F");
    calo_tree_->Branch("pcalojet_pt",&pcalojet_pt_, "pcalojet_pt/F");
    calo_tree_->Branch("pcalojet_p",&pcalojet_p_, "pcalojet_p/F");
    calo_tree_->Branch("pcalojet_eta",&pcalojet_eta_, "pcalojet_eta/F");
    calo_tree_->Branch("pcalojet_phi",&pcalojet_phi_, "pcalojet_phi/F");
    calo_tree_->Branch("pcalojet_emf",&pcalojet_emf_, "pcalojet_emf/F");
    calo_tree_->Branch("pcalojet_scale",&pcalojet_scale_, "pcalojet_scale/F");
    calo_tree_->Branch("pcalojet_genpt",&pcalojet_genpt_, "pcalojet_genpt/F");
    calo_tree_->Branch("pcalojet_genp",&pcalojet_genp_, "pcalojet_genp/F");
    calo_tree_->Branch("pcalojet_gendr",&pcalojet_gendr_, "pcalojet_gendr/F");
    calo_tree_->Branch("pcalojet_EBE",&pcalojet_EBE_, "pcalojet_EBE/F");
    calo_tree_->Branch("pcalojet_EEE",&pcalojet_EEE_, "pcalojet_EEE/F");
    calo_tree_->Branch("pcalojet_HBE",&pcalojet_HBE_, "pcalojet_HBE/F");
    calo_tree_->Branch("pcalojet_HEE",&pcalojet_HEE_, "pcalojet_HEE/F");
    calo_tree_->Branch("pcalojet_HFE",&pcalojet_HFE_, "pcalojet_HFE/F");
    calo_tree_->Branch("pcalojet_ntwrs",&pcalojet_ntwrs_, "pcalojet_ntwrs/I");
    calo_tree_->Branch("pcalojet_twr_ieta",pcalojet_twr_ieta_, "pcalojet_twr_ieta[pcalojet_ntwrs]/I");
    calo_tree_->Branch("pcalojet_twr_eme",pcalojet_twr_eme_, "pcalojet_twr_eme[pcalojet_ntwrs]/F");
    calo_tree_->Branch("pcalojet_twr_hade",pcalojet_twr_hade_, "pcalojet_twr_hade[pcalojet_ntwrs]/F");
    calo_tree_->Branch("calo_dijet_deta",&calo_dijet_deta_, "calo_dijet_deta/F");
    calo_tree_->Branch("calo_dijet_dphi",&calo_dijet_dphi_, "calo_dijet_dphi/F");
    calo_tree_->Branch("calo_dijet_balance",&calo_dijet_balance_, "calo_dijet_balance/F");
    calo_tree_->Branch("calo_thirdjet_px",&calo_thirdjet_px_, "calo_thirdjet_px/F");
    calo_tree_->Branch("calo_thirdjet_py",&calo_thirdjet_py_, "calo_thirdjet_py/F");
  }

  if(doPFJets_){
    h_types_ = new TH1D("h_types","h_types",16,0,16);
    h_ntypes_ = new TH1D("h_ntypes","h_ntypes",50,0,50);
    hPassSelPF_ = new TH1D("hPassSelectionPF", "Selection Pass Failures PFJets",200,-0.5,199.5);

    pf_tree_ = new TTree("pf_dijettree", "tree for dijet balancing using PFJets");

    pf_tree_->Branch("tpfjet_pt",&tpfjet_pt_, "tpfjet_pt/F");
    pf_tree_->Branch("tpfjet_p",&tpfjet_p_, "tpfjet_p/F");
    pf_tree_->Branch("tpfjet_eta",&tpfjet_eta_, "tpfjet_eta/F");
    pf_tree_->Branch("tpfjet_phi",&tpfjet_phi_, "tpfjet_phi/F");
    pf_tree_->Branch("tpfjet_scale",&tpfjet_scale_, "tpfjet_scale/F");
    pf_tree_->Branch("tpfjet_genpt",&tpfjet_genpt_, "tpfjet_genpt/F");
    pf_tree_->Branch("tpfjet_genp",&tpfjet_genp_, "tpfjet_genp/F");
    pf_tree_->Branch("tpfjet_gendr",&tpfjet_gendr_, "tpfjet_gendr/F");
    pf_tree_->Branch("tpfjet_unkown_E",&tpfjet_unkown_E_, "tpfjet_unkown_E/F");
    pf_tree_->Branch("tpfjet_chHad_E",&tpfjet_chHad_E_, "tpfjet_chHad_E/F");
    pf_tree_->Branch("tpfjet_electron_E",&tpfjet_electron_E_, "tpfjet_electron_E/F");
    pf_tree_->Branch("tpfjet_muon_E",&tpfjet_muon_E_, "tpfjet_muon_E/F");
    pf_tree_->Branch("tpfjet_photon_E",&tpfjet_photon_E_, "tpfjet_photon_E/F");
    pf_tree_->Branch("tpfjet_Had0_E",&tpfjet_Had0_E_, "tpfjet_Had0_E/F");
    pf_tree_->Branch("tpfjet_HFHad_E",&tpfjet_HFHad_E_, "tpfjet_HFHad_E/F");
    pf_tree_->Branch("tpfjet_HFEM_E",&tpfjet_HFEM_E_, "tpfjet_HFEM_E/F");
    pf_tree_->Branch("tpfjet_unkown_px",&tpfjet_unkown_px_, "tpfjet_unkown_px/F");
    pf_tree_->Branch("tpfjet_chHad_px",&tpfjet_chHad_px_, "tpfjet_chHad_px/F");
    pf_tree_->Branch("tpfjet_electron_px",&tpfjet_electron_px_, "tpfjet_electron_px/F");
    pf_tree_->Branch("tpfjet_muon_px",&tpfjet_muon_px_, "tpfjet_muon_px/F");
    pf_tree_->Branch("tpfjet_photon_px",&tpfjet_photon_px_, "tpfjet_photon_px/F");
    pf_tree_->Branch("tpfjet_Had0_px",&tpfjet_Had0_px_, "tpfjet_Had0_px/F");
    pf_tree_->Branch("tpfjet_HFHad_px",&tpfjet_HFHad_px_, "tpfjet_HFHad_px/F");
    pf_tree_->Branch("tpfjet_HFEM_px",&tpfjet_HFEM_px_, "tpfjet_HFEM_px/F");
    pf_tree_->Branch("tpfjet_unkown_py",&tpfjet_unkown_py_, "tpfjet_unkown_py/F");
    pf_tree_->Branch("tpfjet_chHad_py",&tpfjet_chHad_py_, "tpfjet_chHad_py/F");
    pf_tree_->Branch("tpfjet_electron_py",&tpfjet_electron_py_, "tpfjet_electron_py/F");
    pf_tree_->Branch("tpfjet_muon_py",&tpfjet_muon_py_, "tpfjet_muon_py/F");
    pf_tree_->Branch("tpfjet_photon_py",&tpfjet_photon_py_, "tpfjet_photon_py/F");
    pf_tree_->Branch("tpfjet_Had0_py",&tpfjet_Had0_py_, "tpfjet_Had0_py/F");
    pf_tree_->Branch("tpfjet_HFHad_py",&tpfjet_HFHad_py_, "tpfjet_HFHad_py/F");
    pf_tree_->Branch("tpfjet_HFEM_py",&tpfjet_HFEM_py_, "tpfjet_HFEM_py/F");
    pf_tree_->Branch("tpfjet_unkown_pz",&tpfjet_unkown_pz_, "tpfjet_unkown_pz/F");
    pf_tree_->Branch("tpfjet_chHad_pz",&tpfjet_chHad_pz_, "tpfjet_chHad_pz/F");
    pf_tree_->Branch("tpfjet_electron_pz",&tpfjet_electron_pz_, "tpfjet_electron_pz/F");
    pf_tree_->Branch("tpfjet_muon_pz",&tpfjet_muon_pz_, "tpfjet_muon_pz/F");
    pf_tree_->Branch("tpfjet_photon_pz",&tpfjet_photon_pz_, "tpfjet_photon_pz/F");
    pf_tree_->Branch("tpfjet_Had0_pz",&tpfjet_Had0_pz_, "tpfjet_Had0_pz/F");
    pf_tree_->Branch("tpfjet_HFHad_pz",&tpfjet_HFHad_pz_, "tpfjet_HFHad_pz/F");
    pf_tree_->Branch("tpfjet_HFEM_pz",&tpfjet_HFEM_pz_, "tpfjet_HFEM_pz/F");
    pf_tree_->Branch("tpfjet_unkown_n",&tpfjet_unkown_n_, "tpfjet_unkown_n/I");
    pf_tree_->Branch("tpfjet_chHad_n",&tpfjet_chHad_n_, "tpfjet_chHad_n/I");
    pf_tree_->Branch("tpfjet_electron_n",&tpfjet_electron_n_, "tpfjet_electron_n/I");
    pf_tree_->Branch("tpfjet_muon_n",&tpfjet_muon_n_, "tpfjet_muon_n/I");
    pf_tree_->Branch("tpfjet_photon_n",&tpfjet_photon_n_, "tpfjet_photon_n/I");
    pf_tree_->Branch("tpfjet_Had0_n",&tpfjet_Had0_n_, "tpfjet_Had0_n/I");
    pf_tree_->Branch("tpfjet_HFHad_n",&tpfjet_HFHad_n_, "tpfjet_HFHad_n/I");
    pf_tree_->Branch("tpfjet_HFEM_n",&tpfjet_HFEM_n_, "tpfjet_HFEM_n/I");
    pf_tree_->Branch("tpfjet_ntwrs",&tpfjet_ntwrs_, "tpfjet_ntwrs/I");
    pf_tree_->Branch("tpfjet_twr_ieta",tpfjet_twr_ieta_, "tpfjet_twr_ieta[tpfjet_ntwrs]/I");
    pf_tree_->Branch("tpfjet_twr_hade",tpfjet_twr_hade_, "tpfjet_twr_hade[tpfjet_ntwrs]/F");
    pf_tree_->Branch("tpfjet_twr_frac",tpfjet_twr_frac_, "tpfjet_twr_frac[tpfjet_ntwrs]/F");
    pf_tree_->Branch("ppfjet_pt",&ppfjet_pt_, "ppfjet_pt/F");
    pf_tree_->Branch("ppfjet_p",&ppfjet_p_, "ppfjet_p/F");
    pf_tree_->Branch("ppfjet_eta",&ppfjet_eta_, "ppfjet_eta/F");
    pf_tree_->Branch("ppfjet_phi",&ppfjet_phi_, "ppfjet_phi/F");
    pf_tree_->Branch("ppfjet_scale",&ppfjet_scale_, "ppfjet_scale/F");
    pf_tree_->Branch("ppfjet_genpt",&ppfjet_genpt_, "ppfjet_genpt/F");
    pf_tree_->Branch("ppfjet_genp",&ppfjet_genp_, "ppfjet_genp/F");
    pf_tree_->Branch("ppfjet_gendr",&ppfjet_gendr_, "ppfjet_gendr/F");
    pf_tree_->Branch("ppfjet_unkown_E",&ppfjet_unkown_E_, "ppfjet_unkown_E/F");
    pf_tree_->Branch("ppfjet_chHad_E",&ppfjet_chHad_E_, "ppfjet_chHad_E/F");
    pf_tree_->Branch("ppfjet_electron_E",&ppfjet_electron_E_, "ppfjet_electron_E/F");
    pf_tree_->Branch("ppfjet_muon_E",&ppfjet_muon_E_, "ppfjet_muon_E/F");
    pf_tree_->Branch("ppfjet_photon_E",&ppfjet_photon_E_, "ppfjet_photon_E/F");
    pf_tree_->Branch("ppfjet_Had0_E",&ppfjet_Had0_E_, "ppfjet_Had0_E/F");
    pf_tree_->Branch("ppfjet_HFHad_E",&ppfjet_HFHad_E_, "ppfjet_HFHad_E/F");
    pf_tree_->Branch("ppfjet_HFEM_E",&ppfjet_HFEM_E_, "ppfjet_HFEM_E/F");
    pf_tree_->Branch("ppfjet_unkown_px",&ppfjet_unkown_px_, "ppfjet_unkown_px/F");
    pf_tree_->Branch("ppfjet_chHad_px",&ppfjet_chHad_px_, "ppfjet_chHad_px/F");
    pf_tree_->Branch("ppfjet_electron_px",&ppfjet_electron_px_, "ppfjet_electron_px/F");
    pf_tree_->Branch("ppfjet_muon_px",&ppfjet_muon_px_, "ppfjet_muon_px/F");
    pf_tree_->Branch("ppfjet_photon_px",&ppfjet_photon_px_, "ppfjet_photon_px/F");
    pf_tree_->Branch("ppfjet_Had0_px",&ppfjet_Had0_px_, "ppfjet_Had0_px/F");
    pf_tree_->Branch("ppfjet_HFHad_px",&ppfjet_HFHad_px_, "ppfjet_HFHad_px/F");
    pf_tree_->Branch("ppfjet_HFEM_px",&ppfjet_HFEM_px_, "ppfjet_HFEM_px/F");
    pf_tree_->Branch("ppfjet_unkown_py",&ppfjet_unkown_py_, "ppfjet_unkown_py/F");
    pf_tree_->Branch("ppfjet_chHad_py",&ppfjet_chHad_py_, "ppfjet_chHad_py/F");
    pf_tree_->Branch("ppfjet_electron_py",&ppfjet_electron_py_, "ppfjet_electron_py/F");
    pf_tree_->Branch("ppfjet_muon_py",&ppfjet_muon_py_, "ppfjet_muon_py/F");
    pf_tree_->Branch("ppfjet_photon_py",&ppfjet_photon_py_, "ppfjet_photon_py/F");
    pf_tree_->Branch("ppfjet_Had0_py",&ppfjet_Had0_py_, "ppfjet_Had0_py/F");
    pf_tree_->Branch("ppfjet_HFHad_py",&ppfjet_HFHad_py_, "ppfjet_HFHad_py/F");
    pf_tree_->Branch("ppfjet_HFEM_py",&ppfjet_HFEM_py_, "ppfjet_HFEM_py/F");
    pf_tree_->Branch("ppfjet_unkown_pz",&ppfjet_unkown_pz_, "ppfjet_unkown_pz/F");
    pf_tree_->Branch("ppfjet_chHad_pz",&ppfjet_chHad_pz_, "ppfjet_chHad_pz/F");
    pf_tree_->Branch("ppfjet_electron_pz",&ppfjet_electron_pz_, "ppfjet_electron_pz/F");
    pf_tree_->Branch("ppfjet_muon_pz",&ppfjet_muon_pz_, "ppfjet_muon_pz/F");
    pf_tree_->Branch("ppfjet_photon_pz",&ppfjet_photon_pz_, "ppfjet_photon_pz/F");
    pf_tree_->Branch("ppfjet_Had0_pz",&ppfjet_Had0_pz_, "ppfjet_Had0_pz/F");
    pf_tree_->Branch("ppfjet_HFHad_pz",&ppfjet_HFHad_pz_, "ppfjet_HFHad_pz/F");
    pf_tree_->Branch("ppfjet_HFEM_pz",&ppfjet_HFEM_pz_, "ppfjet_HFEM_pz/F");
    pf_tree_->Branch("ppfjet_unkown_n",&ppfjet_unkown_n_, "ppfjet_unkown_n/I");
    pf_tree_->Branch("ppfjet_chHad_n",&ppfjet_chHad_n_, "ppfjet_chHad_n/I");
    pf_tree_->Branch("ppfjet_electron_n",&ppfjet_electron_n_, "ppfjet_electron_n/I");
    pf_tree_->Branch("ppfjet_muon_n",&ppfjet_muon_n_, "ppfjet_muon_n/I");
    pf_tree_->Branch("ppfjet_photon_n",&ppfjet_photon_n_, "ppfjet_photon_n/I");
    pf_tree_->Branch("ppfjet_Had0_n",&ppfjet_Had0_n_, "ppfjet_Had0_n/I");
    pf_tree_->Branch("ppfjet_HFHad_n",&ppfjet_HFHad_n_, "ppfjet_HFHad_n/I");
    pf_tree_->Branch("ppfjet_HFEM_n",&ppfjet_HFEM_n_, "ppfjet_HFEM_n/I");
    pf_tree_->Branch("ppfjet_ntwrs",&ppfjet_ntwrs_, "ppfjet_ntwrs/I");
    pf_tree_->Branch("ppfjet_twr_ieta",ppfjet_twr_ieta_, "ppfjet_twr_ieta[ppfjet_ntwrs]/I");
    pf_tree_->Branch("ppfjet_twr_hade",ppfjet_twr_hade_, "ppfjet_twr_hade[ppfjet_ntwrs]/F");
    pf_tree_->Branch("ppfjet_twr_frac",ppfjet_twr_frac_, "ppfjet_twr_frac[ppfjet_ntwrs]/F");
    pf_tree_->Branch("pf_dijet_deta",&pf_dijet_deta_, "pf_dijet_deta/F");
    pf_tree_->Branch("pf_dijet_dphi",&pf_dijet_dphi_, "pf_dijet_dphi/F");
    pf_tree_->Branch("pf_dijet_balance",&pf_dijet_balance_, "pf_dijet_balance/F");
    pf_tree_->Branch("pf_thirdjet_px",&pf_thirdjet_px_, "pf_thirdjet_px/F");
    pf_tree_->Branch("pf_thirdjet_py",&pf_thirdjet_py_, "pf_thirdjet_py/F");
  }

  hfieta[0]  = 2.964;
  hfieta[1]  = 3.139;
  hfieta[2]  = 3.314;
  hfieta[3]  = 3.489;
  hfieta[4]  = 3.664;
  hfieta[5]  = 3.839;
  hfieta[6]  = 4.013;
  hfieta[7]  = 4.191;
  hfieta[8]  = 4.363;
  hfieta[9]  = 4.538;
  hfieta[10] = 4.716;
  hfieta[11] = 4.889;
  hfieta[12] = 5.191;

  maxhfiphi = 0;

  return;
}  

// ------------ method called once each job just after ending the event loop  ------------
void 
CalcRespCorrDiJets::endJob() {
  // write histograms
  rootfile_->cd();
  if(doCaloJets_){
    hPassSelCalo_->Write();
    calo_tree_->Write();
  }
  if(doPFJets_){
    hPassSelPF_->Write();
    pf_tree_->Write();
  }
  h_types_->Write();
  h_ntypes_->Write();
  rootfile_->Close();

  std::cout << "maxhfiphi " << maxhfiphi << std::endl;
}

// helper function

double CalcRespCorrDiJets::deltaR(const reco::Jet* j1, const reco::Jet* j2)
{
  double deta = j1->eta()-j2->eta();
  double dphi = std::fabs(j1->phi()-j2->phi());
  if(dphi>3.1415927) dphi = 2*3.1415927 - dphi;
  return std::sqrt(deta*deta + dphi*dphi);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CalcRespCorrDiJets);
