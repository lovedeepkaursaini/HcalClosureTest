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
//
// CalcrespCorrDiJets
//

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
  edm::Handle<reco::GenJetCollection> genjets;
  iEvent.getByLabel(genJetCollName_,genjets);
  if(!genjets.isValid()) {
    throw edm::Exception(edm::errors::ProductNotFound)
      << " could not find GenJetCollection named " << genJetCollName_ << ".\n";
    return;
  }

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
      for(reco::GenJetCollection::const_iterator it=genjets->begin(); it!=genjets->end(); ++it) {
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
    //iEvent.get("HBHERecHitsSorted_reducedHcalRecHits_hbhereco_RECO",hbreco);
    //iEvent.getManyByType(hbreco);
    //iEvent.getByLabel("reducedHcalRecHits","hbhereco",hbreco);
    iEvent.getByLabel(RecHitLabelName_,hbheRecHitInstance_,hbhereco);
    //iEvent.getByLabel(hbRecHitName_,hbreco);
    if(!hbhereco.isValid()) {
      throw edm::Exception(edm::errors::ProductNotFound)
	<< " could not find HBHERecHit named " << RecHitLabelName_ << ":" << hbheRecHitInstance_ << ".\n";
      return;
    }

    edm::Handle<edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>> hfreco;
    iEvent.getByLabel(RecHitLabelName_,hfRecHitInstance_,hfreco);
    if(!hfreco.isValid()) {
      throw edm::Exception(edm::errors::ProductNotFound)
	<< " could not find HFRecHit named " << RecHitLabelName_ << ":" << hfRecHitInstance_ << ".\n";
      return;
    }

    edm::Handle<edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>> horeco;
    iEvent.getByLabel(RecHitLabelName_,hoRecHitInstance_,horeco);
    if(!horeco.isValid()) {
      throw edm::Exception(edm::errors::ProductNotFound)
	<< " could not find HORecHit named " << RecHitLabelName_ << ":" << hoRecHitInstance_ << ".\n";
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
      
      // emf cuts
      /*double sumEMtag = pf_tag.jet()->chargedEmEnergy() + pf_tag.jet()->neutralEmEnergy() + pf_tag.jet()->HFEMEnergy();
	double sumHadtag = pf_tag.jet()->chargedHadronEnergy() + pf_tag.jet()->neutralHadronEnergy() + pf_tag.jet()->HFHadronEnergy();
	if(sumEMtag/(sumEMtag + sumHadtag) > maxJetEMF_) passSelPF |= 0x20;
	double sumEMprobe = pf_probe.jet()->chargedEmEnergy() + pf_probe.jet()->neutralEmEnergy() + pf_probe.jet()->HFEMEnergy();
	double sumHadprobe = pf_probe.jet()->chargedHadronEnergy() + pf_probe.jet()->neutralHadronEnergy() + pf_probe.jet()->HFHadronEnergy();
	if(sumEMprobe/(sumEMprobe + sumHadprobe) > maxJetEMF_) passSelPF |= 0x20;*/
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
      
      // fill tag jet variables
      tpfjet_pt_    = pf_tag.jet()->pt();
      tpfjet_p_     = pf_tag.jet()->p();
      tpfjet_eta_   = pf_tag.jet()->eta();
      tpfjet_phi_   = pf_tag.jet()->phi();
      /*double sumEM = pf_tag.jet()->chargedEmEnergy() + pf_tag.jet()->neutralEmEnergy() + pf_tag.jet()->HFEMEnergy();
	double sumHad = pf_tag.jet()->chargedHadronEnergy() + pf_tag.jet()->neutralHadronEnergy() + pf_tag.jet()->HFHadronEnergy();
	tpfjet_emf_   = sumEM/(sumEM + sumHad);*/
      tpfjet_scale_ = pf_tag.scale();
      //tpfjet_EBE_   = pf_tag.jet()->emEnergyInEB();
      //tpfjet_EEE_   = pf_tag.jet()->emEnergyInEE();
      //tpfjet_HBE_   = pf_tag.jet()->hadEnergyInHB();
      //tpfjet_HEE_   = pf_tag.jet()->hadEnergyInHE();
      //tpfjet_HFE_   = pf_tag.jet()->emEnergyInHF() + tag.jet()->hadEnergyInHF();
      tpfjet_ntwrs_=0;
      
      //std::cout << pf_tag.jet()->print() << std::endl;
      
      std::vector<reco::PFCandidatePtr> tagconst=pf_tag.jet()->getPFConstituents();
      for(std::vector<reco::PFCandidatePtr>::const_iterator it=tagconst.begin(); it!=tagconst.end(); ++it){
	//std::cout << (*it)->particleId() << " " << (*it)->pdgId() << " " << (*it)->eta() << " " << (*it)->phi() << std::endl;
	
	//std::cout<<(*it)->elementsInBlocks().size()<<std::endl;
	int maxElement=(*it)->elementsInBlocks().size();
	for(int e=0; e<maxElement; ++e){
	  reco::PFBlockRef blockRef = (*it)->elementsInBlocks()[e].first;
	  //if(!blockRef.isNull())std::cout<<"NOT NULL "<<std::endl;
	  const edm::OwnVector<reco::PFBlockElement>& elements = blockRef->elements();
	  //std::cout << "size: " << elements.size() << std::endl;
	  for(unsigned iEle=0; iEle<elements.size(); iEle++) {
	    if(elements[iEle].index() == (*it)->elementsInBlocks()[e].second){//matched to your candidate //elements[iEle].trackRef() or elements[iEle].clusterRef()
	      //std::cout << "matched element " << iEle << "   " << elements[iEle].type() << std::endl;
	      if(elements[iEle].type() == reco::PFBlockElement::HCAL){
		std::cout << "HCal Cluster" << std::endl;
		reco::PFClusterRef clusterref = elements[iEle].clusterRef();
		/*for(reco::PFClusterRef::const_iterator itt=clusterref.begin(); itt!=clusterref.end(); ++itt){
		  std::cout << (*itt) << std::endl;
		  }*/
		reco::PFCluster cluster = *clusterref;
		std::vector<std::pair<DetId,float>> hitsAndFracs = cluster.hitsAndFractions();
		int nHits = hitsAndFracs.size();
		for(int iHit=0; iHit<nHits; iHit++){
		  //std::cout << hitsAndFracs[iHit].first << ": " << hitsAndFracs[iHit].second << std::endl;
		  HcalDetId hDet(hitsAndFracs[iHit].first.rawId());
		  std::cout << "detector: " << hitsAndFracs[iHit].first.det() << " " << hitsAndFracs[iHit].first.subdetId() << " " << hitsAndFracs[iHit].first.rawId() << " " << hDet << " " << hitsAndFracs[iHit].second << std::endl;

		  for(edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>::const_iterator ithb=hbhereco->begin(); ithb!=hbhereco->end(); ++ithb){
		    if((*ithb).id() == hitsAndFracs[iHit].first.rawId()){
		      std::cout << "hit matches " << (*ithb).id() << std::endl;
		    }
		    //std::cout << (*ithb).id() << ": " << (*it).energy() << std::endl;
		    //std::cout << (*it).energy() << std::endl;
		  }
		}
	      }
	    }
	    /*else{
	    std::cout<<"elements "<<iEle<<std::endl;
	    }*/
	  }
	  //	std::cout<<"Looking for Element "<<elements[(*it)->elementsInBlocks()[e].second]<<std::endl;
	  //std::cout<<"\t"<<(*it)->elementsInBlocks()[e].second<<std::endl;
	  //std::cout<<"type "<<(*it)->elementsInBlocks()[e].type()<<std::endl;
	  
	}
	//reco::SuperClusterRef cluster = (*it)->superClusterRef();
	
	//std::cout << cluster << std::endl;
	//const auto& hitsAndFracs = cluster.hitsAndFractions();
      }
      
      /*for(unsigned i=0; i<tagconst.size(); ++i){
	if(tagconst[i].get()){
	std::cout << i << ": " << *(tagconst[i].id()) << std::endl;
	}
	else{
	std::cout << i << ": none" << std::endl;
	}
	}*/
      
      /*std::vector<CaloTowerPtr> tagconst=tag.jet()->getCaloConstituents();
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
	}*/
      
      // fill probe jet variables
      ppfjet_pt_    = pf_probe.jet()->pt();
      ppfjet_p_     = pf_probe.jet()->p();
      ppfjet_eta_   = pf_probe.jet()->eta();
      ppfjet_phi_   = pf_probe.jet()->phi();
      //ppfjet_emf_   = pf_probe.jet()->emEnergyFraction();
      ppfjet_scale_ = pf_probe.scale();
      //ppfjet_EBE_   = pf_probe.jet()->emEnergyInEB();
      //ppfjet_EEE_   = pf_probe.jet()->emEnergyInEE();
      //ppfjet_HBE_   = pf_probe.jet()->hadEnergyInHB();
      //ppfjet_HEE_   = pf_probe.jet()->hadEnergyInHE();
      //ppfjet_HFE_   = pf_probe.jet()->emEnergyInHF() + probe.jet()->hadEnergyInHF();
      ppfjet_ntwrs_=0;
      /*std::vector<CaloTowerPtr> probeconst=probe.jet()->getCaloConstituents();
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
	}*/
      
      // fill genjet tag/probe variables
      tpfjet_gendr_ = 99999.;
      tpfjet_genpt_ = 0;
      tpfjet_genp_  = 0;
      ppfjet_gendr_ = 99999.;
      ppfjet_genpt_ = 0;
      ppfjet_genp_  = 0;
      for(reco::GenJetCollection::const_iterator it=genjets->begin(); it!=genjets->end(); ++it) {
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
      
      // fill dijet variables
      pf_dijet_deta_=std::fabs(std::fabs(pf_tag.jet()->eta())-std::fabs(pf_probe.jet()->eta()));
      pf_dijet_dphi_=pf_tag.jet()->phi()-pf_probe.jet()->phi();
      if(pf_dijet_dphi_>3.1415) pf_dijet_dphi_ = 6.2832-pf_dijet_dphi_;
      pf_dijet_balance_ = (tpfjet_pt_-ppfjet_pt_)/(tpfjet_pt_+ppfjet_pt_);
      
      pf_tree_->Fill();
    }
    
    std::cout << "RecHits " << hbhereco->size() << " " << hfreco->size() << " " << horeco->size() << std::endl;
    //edm::Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>> hbreco;
    for(edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>::const_iterator it=hbhereco->begin(); it!=hbhereco->end(); ++it){
      std::cout << (*it).id().rawId() << " " << (*it).id() << ": " << (*it).energy() << std::endl;
      //std::cout << (*it).energy() << std::endl;
    }
    for(edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>::const_iterator it=hfreco->begin(); it!=hfreco->end(); ++it){
      std::cout << (*it).id().rawId() << " " << (*it).id() << ": " << (*it).energy() << std::endl;
    }
    for(edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>::const_iterator it=horeco->begin(); it!=horeco->end(); ++it){
      std::cout << (*it).id().rawId() << " " << (*it).id() << ": " << (*it).energy() << std::endl;
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
    hPassSelPF_ = new TH1D("hPassSelectionPF", "Selection Pass Failures PFJets",200,-0.5,199.5);

    pf_tree_ = new TTree("pf_dijettree", "tree for dijet balancing using PFJets");
  }

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
  rootfile_->Close();
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
