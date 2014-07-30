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
  unsigned int testEvent = 100000;
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
      
      // fill tag jet variables
      tpfjet_pt_    = pf_tag.jet()->pt();
      tpfjet_p_     = pf_tag.jet()->p();
      tpfjet_eta_   = pf_tag.jet()->eta();
      tpfjet_phi_   = pf_tag.jet()->phi();
      tpfjet_photonE_   = pf_tag.jet()->photonEnergy();
      tpfjet_photonN_   = pf_tag.jet()->photonMultiplicity();
      tpfjet_electronE_ = pf_tag.jet()->electronEnergy();
      tpfjet_electronN_ = pf_tag.jet()->electronMultiplicity();
      tpfjet_muonE_     = pf_tag.jet()->muonEnergy();
      tpfjet_muonN_     = pf_tag.jet()->muonMultiplicity();
      tpfjet_HFEME_     = pf_tag.jet()->HFEMEnergy();
      tpfjet_HFEMN_     = pf_tag.jet()->HFEMMultiplicity();
      tpfjet_scale_ = pf_tag.scale();
      tpfjet_ntwrs_=0;

      //std::cout << pf_tag.jet()->print() << std::endl;
      
      // Get PF constituents and fill HCAL towers
      std::vector<reco::PFCandidatePtr> tagconst=pf_tag.jet()->getPFConstituents();
      if(iEvent.id().event() == testEvent){
	std::cout << "do PFCandidatePtr loop" << std::endl;
      }
      for(std::vector<reco::PFCandidatePtr>::const_iterator it=tagconst.begin(); it!=tagconst.end(); ++it){
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
	    //std::cout << "++++++ " << tmpclusterref.id() << std::endl;
	    if(elements[iEle].index() == (*it)->elementsInBlocks()[e].second){
	      if(iEvent.id().event() == testEvent){
		std::cout << "passed event == testEvent" << std::endl;
	      }
	      if(elements[iEle].type() == reco::PFBlockElement::HCAL// ||
		 //		 elements[iEle].type() == reco::PFBlockElement::HFHAD// ||
		 //		 elements[iEle].type() == reco::PFBlockElement::HO
){
		reco::PFClusterRef clusterref = elements[iEle].clusterRef();
		//std::cout << "------ " << clusterref.id() << std::endl;
		reco::PFCluster cluster = *clusterref;
		
		std::vector<std::pair<DetId,float>> hitsAndFracs = cluster.hitsAndFractions();
		if(iEvent.id().event() == testEvent){
		  std::cout << "Tag hits" << std::endl;
		}
		int nHits = hitsAndFracs.size();
		for(int iHit=0; iHit<nHits; iHit++){
		  HcalDetId hDet(hitsAndFracs[iHit].first.rawId());
		  int etaPhiPF = hitsAndFracs[iHit].first.rawId() & 0x1FFF;

		  if(iEvent.id().event() == testEvent){
		    std::cout << hitsAndFracs[iHit].first.rawId() <<  " ";
		    for(int ijk=31; ijk>=0; ijk--){
		      int x = (hitsAndFracs[iHit].first.rawId() >> ijk) & 0x1;
		      std::cout << x;
		    }
		    std::cout << " " << (hitsAndFracs[iHit].first.rawId() >> 28) << " " << ((hitsAndFracs[iHit].first.rawId() >> 25) & 0x7) << std::endl;
		  }

		  for(edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>::const_iterator ith=hbhereco->begin(); ith!=hbhereco->end(); ++ith){
		    int etaPhiRecHit = (*ith).id().rawId() & 0x1FFF;
		    if(etaPhiPF == etaPhiRecHit){
		      //std::cout << "hit matches " << (*ith).id() << " " << (*ith).id().ieta() << " " << (*ith).energy()*hitsAndFracs[iHit].second << std::endl;

		      if((*ith).energy() > 0.0 && hitsAndFracs[iHit].second > 0.0){
			tpfjet_twr_ieta_[tpfjet_ntwrs_] = (*ith).id().ieta();
			int tmpieta = (*ith).id().ieta();
			if(tmpieta > 41 || tmpieta < -41 || tmpieta == 0){
			  std::cout << tmpieta << std::endl;
			}
			tpfjet_twr_hade_[tpfjet_ntwrs_] = (*ith).energy();
			tpfjet_twr_frac_[tpfjet_ntwrs_] = hitsAndFracs[iHit].second;
			//std::cout << tpfjet_twr_frac_[tpfjet_ntwrs_] << std::endl;
		      

		      //if((*ith).energy() < 0.0 || hitsAndFracs[iHit].second < 0.0){
		      //std::cout << iEvent.id().event() << " tag: " << (*ith).energy() << " " << hitsAndFracs[iHit].second << std::endl;
		      //}

			++tpfjet_ntwrs_;
		      }
		      
		      /*std::cout << hitsAndFracs[iHit].first.rawId() <<  " ";
		      for(int ijk=31; ijk>=0; ijk--){
			int x = (hitsAndFracs[iHit].first.rawId() >> ijk) & 0x1;
			std::cout << x;
		      }
		      std::cout << " " <<  hDet << " " << (hitsAndFracs[iHit].first.rawId() >> 28) << " " << ((hitsAndFracs[iHit].first.rawId() >> 25) & 0x7) << std::endl;
		      std::cout << (*ith).id().rawId() << " ";
		      for(int ijk=31; ijk>=0; ijk--){
			int x = ((*ith).id().rawId() >> ijk) & 0x1;
			std::cout << x;
		      }
		      std::cout << " " << (*ith).id() << " " << ((*ith).id().rawId() >> 28) << " " << (((*ith).id().rawId() >> 25) & 0x7) << std::endl;

		      int  tmpxor = hitsAndFracs[iHit].first.rawId() ^ (*ith).id().rawId();
		      std::cout << "           ";
		      for(int ijk=31; ijk>=0; ijk--){
			int x = (tmpxor >> ijk) & 0x1;
			std::cout << x;
		      }
		      std::cout << std::endl;*/
		    }
		  }

		  for(edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>::const_iterator ith=hfreco->begin(); ith!=hfreco->end(); ++ith){
		    int etaPhiRecHit = (*ith).id().rawId() & 0x1FFF;
		    if(etaPhiPF == etaPhiRecHit){
		      //std::cout << "hit matches " << (*ith).id() << " " << (*ith).id().ieta() << " " << (*ith).energy()*hitsAndFracs[iHit].second << std::endl;

		      if((*ith).energy() > 0.0 && hitsAndFracs[iHit].second > 0.0){
			tpfjet_twr_ieta_[tpfjet_ntwrs_] = (*ith).id().ieta();
			tpfjet_twr_hade_[tpfjet_ntwrs_] = (*ith).energy();
			tpfjet_twr_frac_[tpfjet_ntwrs_] = hitsAndFracs[iHit].second;
			//std::cout << tpfjet_twr_frac_[tpfjet_ntwrs_] << std::endl;

		      //if((*ith).energy() < 0.0 || hitsAndFracs[iHit].second < 0.0){
		      //std::cout << iEvent.id().event() << " tag: " << (*ith).energy() << " " << hitsAndFracs[iHit].second << std::endl;
		      //}

			++tpfjet_ntwrs_;
		      }
		      
		      /*std::cout << hitsAndFracs[iHit].first.rawId() <<  " ";
		      for(int ijk=31; ijk>=0; ijk--){
			int x = (hitsAndFracs[iHit].first.rawId() >> ijk) & 0x1;
			std::cout << x;
		      }
		      std::cout << " " <<  hDet << std::endl;
		      std::cout << (*ith).id().rawId() << " ";
		      for(int ijk=31; ijk>=0; ijk--){
			int x = ((*ith).id().rawId() >> ijk) & 0x1;
			std::cout << x;
		      }
		      std::cout << " " << (*ith).id() << std::endl;

		      int  tmpxor = hitsAndFracs[iHit].first.rawId() ^ (*ith).id().rawId();
		      std::cout << "           ";
		      for(int ijk=31; ijk>=0; ijk--){
			int x = (tmpxor >> ijk) & 0x1;
			std::cout << x;
		      }
		      std::cout << std::endl;*/
		    }
		  }

		  for(edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>::const_iterator ith=horeco->begin(); ith!=horeco->end(); ++ith){
		    int etaPhiRecHit = (*ith).id().rawId() & 0x1FFF;
		    if(etaPhiPF == etaPhiRecHit){
		      //std::cout << "hit matches " << (*ith).id() << " " << (*ith).id().ieta() << " " << (*ith).energy()*hitsAndFracs[iHit].second << std::endl;

		      if((*ith).energy() > 0.0 && hitsAndFracs[iHit].second > 0.0){
			tpfjet_twr_ieta_[tpfjet_ntwrs_] = (*ith).id().ieta();
			tpfjet_twr_hade_[tpfjet_ntwrs_] = (*ith).energy();
			tpfjet_twr_frac_[tpfjet_ntwrs_] = hitsAndFracs[iHit].second;
			//std::cout << tpfjet_twr_frac_[tpfjet_ntwrs_] << std::endl;
		      
		      //if((*ith).energy() < 0.0 || hitsAndFracs[iHit].second < 0.0){
		      //std::cout << iEvent.id().event() << " tag: " << (*ith).energy() << " " << hitsAndFracs[iHit].second << std::endl;
		      //}

			++tpfjet_ntwrs_;
		      }
		      
		      /*std::cout << hitsAndFracs[iHit].first.rawId() <<  " ";
		      for(int ijk=31; ijk>=0; ijk--){
			int x = (hitsAndFracs[iHit].first.rawId() >> ijk) & 0x1;
			std::cout << x;
		      }
		      std::cout << " " <<  hDet << std::endl;
		      std::cout << (*ith).id().rawId() << " ";
		      for(int ijk=31; ijk>=0; ijk--){
			int x = ((*ith).id().rawId() >> ijk) & 0x1;
			std::cout << x;
		      }
		      std::cout << " " << (*ith).id() << std::endl;

		      int  tmpxor = hitsAndFracs[iHit].first.rawId() ^ (*ith).id().rawId();
		      std::cout << "           ";
		      for(int ijk=31; ijk>=0; ijk--){
			int x = (tmpxor >> ijk) & 0x1;
			std::cout << x;
		      }
		      std::cout << std::endl;*/
		    }
		  }
		} // Loop over hits
	      } // Test if element is from HCAL
	      else if(elements[iEle].type() == reco::PFBlockElement::HFHAD){
		//std::cout << "HFHAD " << iEle << " " << (*it)->elementsInBlocks()[e].second << " " << elements[iEle] << std::endl;
		reco::PFClusterRef clusterref = elements[iEle].clusterRef();
		//std::cout << "id " << clusterref.id() << std::endl;
		//reco::PFCluster cluster = *clusterref;
	      }
	    } // Test for right element index
	  } // Loop over elements
	} // Loop over elements in blocks
      } // Loop over PF constitutents
      
      // fill probe jet variables
      ppfjet_pt_    = pf_probe.jet()->pt();
      ppfjet_p_     = pf_probe.jet()->p();
      ppfjet_eta_   = pf_probe.jet()->eta();
      ppfjet_phi_   = pf_probe.jet()->phi();
      ppfjet_photonE_   = pf_probe.jet()->photonEnergy();
      ppfjet_photonN_   = pf_probe.jet()->photonMultiplicity();
      ppfjet_electronE_ = pf_probe.jet()->electronEnergy();
      ppfjet_electronN_ = pf_probe.jet()->electronMultiplicity();
      ppfjet_muonE_     = pf_probe.jet()->muonEnergy();
      ppfjet_muonN_     = pf_probe.jet()->muonMultiplicity();
      ppfjet_HFEME_     = pf_probe.jet()->HFEMEnergy();
      ppfjet_HFEMN_     = pf_probe.jet()->HFEMMultiplicity();
      ppfjet_scale_ = pf_probe.scale();
      ppfjet_ntwrs_=0;

      // Get PF constituents and fill HCAL towers
      std::vector<reco::PFCandidatePtr> probeconst=pf_probe.jet()->getPFConstituents();
      for(std::vector<reco::PFCandidatePtr>::const_iterator it=probeconst.begin(); it!=probeconst.end(); ++it){
	int maxElement=(*it)->elementsInBlocks().size();
	for(int e=0; e<maxElement; ++e){
	  reco::PFBlockRef blockRef = (*it)->elementsInBlocks()[e].first;
	  const edm::OwnVector<reco::PFBlockElement>& elements = blockRef->elements();
	  for(unsigned iEle=0; iEle<elements.size(); iEle++) {
	    if(elements[iEle].index() == (*it)->elementsInBlocks()[e].second){
	      if(elements[iEle].type() == reco::PFBlockElement::HCAL// ||
		 //		 elements[iEle].type() == reco::PFBlockElement::HFHAD// ||
		 //		 elements[iEle].type() == reco::PFBlockElement::HO
){
		reco::PFClusterRef clusterref = elements[iEle].clusterRef();
		reco::PFCluster cluster = *clusterref;
		
		std::vector<std::pair<DetId,float>> hitsAndFracs = cluster.hitsAndFractions();
		if(iEvent.id().event() == testEvent){
		  std::cout << "Probe hits" << std::endl;
		}
		int nHits = hitsAndFracs.size();
		for(int iHit=0; iHit<nHits; iHit++){
		  HcalDetId hDet(hitsAndFracs[iHit].first.rawId());
		  int etaPhiPF = hitsAndFracs[iHit].first.rawId() & 0x1FFF;
		  
		  if(iEvent.id().event() == testEvent){
		    std::cout << hitsAndFracs[iHit].first.rawId() <<  " ";
		    for(int ijk=31; ijk>=0; ijk--){
		      int x = (hitsAndFracs[iHit].first.rawId() >> ijk) & 0x1;
		      std::cout << x;
		    }
		    std::cout << " " << (hitsAndFracs[iHit].first.rawId() >> 28) << " " << ((hitsAndFracs[iHit].first.rawId() >> 25) & 0x7) << std::endl;
		  }

		  for(edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>::const_iterator ith=hbhereco->begin(); ith!=hbhereco->end(); ++ith){
		    int etaPhiRecHit = (*ith).id().rawId() & 0x1FFF;
		    if(etaPhiPF == etaPhiRecHit){
		      //std::cout << "hit matches " << (*ith).id() << " " << (*ith).id().ieta() << " " << (*ith).energy()*hitsAndFracs[iHit].second << std::endl;

		      if((*ith).energy() > 0.0 && hitsAndFracs[iHit].second > 0.0){
			ppfjet_twr_ieta_[ppfjet_ntwrs_] = (*ith).id().ieta();
			int tmpieta = (*ith).id().ieta();
			if(tmpieta > 41 || tmpieta < -41 || tmpieta == 0){
			  std::cout << tmpieta << std::endl;
			}
			ppfjet_twr_hade_[ppfjet_ntwrs_] = (*ith).energy();
			ppfjet_twr_frac_[ppfjet_ntwrs_] = hitsAndFracs[iHit].second;
		      }
		      //if((*ith).energy() < 0.0 || hitsAndFracs[iHit].second < 0.0){
		      //std::cout << iEvent.id().event() << " probe: " << (*ith).energy() << " " << hitsAndFracs[iHit].second << std::endl;
		      //}

		      ++ppfjet_ntwrs_;
		      
		      /*std::cout << hitsAndFracs[iHit].first.rawId() <<  " ";
		      for(int ijk=31; ijk>=0; ijk--){
			int x = (hitsAndFracs[iHit].first.rawId() >> ijk) & 0x1;
			std::cout << x;
		      }
		      std::cout << " " <<  hDet << " " << (hitsAndFracs[iHit].first.rawId() >> 28) << " " << ((hitsAndFracs[iHit].first.rawId() >> 25) & 0x7) << std::endl;
		      std::cout << (*ith).id().rawId() << " ";
		      for(int ijk=31; ijk>=0; ijk--){
			int x = ((*ith).id().rawId() >> ijk) & 0x1;
			std::cout << x;
		      }
		      std::cout << " " << (*ith).id() << " " << ((*ith).id().rawId() >> 28) << " " << (((*ith).id().rawId() >> 25) & 0x7) << std::endl;

		      int  tmpxor = hitsAndFracs[iHit].first.rawId() ^ (*ith).id().rawId();
		      std::cout << "           ";
		      for(int ijk=31; ijk>=0; ijk--){
			int x = (tmpxor >> ijk) & 0x1;
			std::cout << x;
		      }
		      std::cout << std::endl;*/
		    }
		  }

		  for(edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>::const_iterator ith=hfreco->begin(); ith!=hfreco->end(); ++ith){
		    int etaPhiRecHit = (*ith).id().rawId() & 0x1FFF;
		    if(etaPhiPF == etaPhiRecHit){
		      //std::cout << "hit matches " << (*ith).id() << " " << (*ith).id().ieta() << " " << (*ith).energy()*hitsAndFracs[iHit].second << std::endl;

		      if((*ith).energy() > 0.0 && hitsAndFracs[iHit].second > 0.0){
			ppfjet_twr_ieta_[ppfjet_ntwrs_] = (*ith).id().ieta();
			ppfjet_twr_hade_[ppfjet_ntwrs_] = (*ith).energy();
			ppfjet_twr_frac_[ppfjet_ntwrs_] = hitsAndFracs[iHit].second;
		      }

		      //if((*ith).energy() < 0.0 || hitsAndFracs[iHit].second < 0.0){
		      //std::cout << iEvent.id().event() << " probe: " << (*ith).energy() << " " << hitsAndFracs[iHit].second << std::endl;
		      //}

		      ++ppfjet_ntwrs_;
		      
		      /*std::cout << hitsAndFracs[iHit].first.rawId() <<  " ";
		      for(int ijk=31; ijk>=0; ijk--){
			int x = (hitsAndFracs[iHit].first.rawId() >> ijk) & 0x1;
			std::cout << x;
		      }
		      std::cout << " " <<  hDet << std::endl;
		      std::cout << (*ith).id().rawId() << " ";
		      for(int ijk=31; ijk>=0; ijk--){
			int x = ((*ith).id().rawId() >> ijk) & 0x1;
			std::cout << x;
		      }
		      std::cout << " " << (*ith).id() << std::endl;

		      int  tmpxor = hitsAndFracs[iHit].first.rawId() ^ (*ith).id().rawId();
		      std::cout << "           ";
		      for(int ijk=31; ijk>=0; ijk--){
			int x = (tmpxor >> ijk) & 0x1;
			std::cout << x;
		      }
		      std::cout << std::endl;*/
		    }
		  }

		  for(edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>::const_iterator ith=horeco->begin(); ith!=horeco->end(); ++ith){
		    int etaPhiRecHit = (*ith).id().rawId() & 0x1FFF;
		    if(etaPhiPF == etaPhiRecHit){
		      //std::cout << "hit matches " << (*ith).id() << " " << (*ith).id().ieta() << " " << (*ith).energy()*hitsAndFracs[iHit].second << std::endl;

		      if((*ith).energy() > 0.0 && hitsAndFracs[iHit].second > 0.0){
			ppfjet_twr_ieta_[ppfjet_ntwrs_] = (*ith).id().ieta();
			ppfjet_twr_hade_[ppfjet_ntwrs_] = (*ith).energy();
			ppfjet_twr_frac_[ppfjet_ntwrs_] = hitsAndFracs[iHit].second;

			//if((*ith).energy() < 0.0){
			//std::cout << iEvent.id().event() << ": " << (*ith).energy() << std::endl;
			//}
			//if(hitsAndFracs[iHit].second < 0.0){
			//std::cout << iEvent.id().event() << ": " << hitsAndFracs[iHit].second << std::endl;
			//}

			++ppfjet_ntwrs_;
		      }
		      
		      /*std::cout << hitsAndFracs[iHit].first.rawId() <<  " ";
		      for(int ijk=31; ijk>=0; ijk--){
			int x = (hitsAndFracs[iHit].first.rawId() >> ijk) & 0x1;
			std::cout << x;
		      }
		      std::cout << " " <<  hDet << std::endl;
		      std::cout << (*ith).id().rawId() << " ";
		      for(int ijk=31; ijk>=0; ijk--){
			int x = ((*ith).id().rawId() >> ijk) & 0x1;
			std::cout << x;
		      }
		      std::cout << " " << (*ith).id() << std::endl;

		      int  tmpxor = hitsAndFracs[iHit].first.rawId() ^ (*ith).id().rawId();
		      std::cout << "           ";
		      for(int ijk=31; ijk>=0; ijk--){
			int x = (tmpxor >> ijk) & 0x1;
			std::cout << x;
		      }
		      std::cout << std::endl;*/
		    }
		  }
		} // Loop over hits
	      } // Test if element is from HCAL
	    } // Test for right element index
	  } // Loop over elements
	} // Loop over elements in blocks
      } // Loop over PF constitutents
      
      // fill genjet tag/probe variables
      tpfjet_gendr_ = 99999.;
      tpfjet_genpt_ = 0;
      tpfjet_genp_  = 0;
      ppfjet_gendr_ = 99999.;
      ppfjet_genpt_ = 0;
      ppfjet_genp_  = 0;
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
	for(edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>::const_iterator it=hbhereco->begin(); it!=hbhereco->end(); ++it){
	  std::cout << (*it).id().rawId() <<  " ";
	  for(int ijk=31; ijk>=0; ijk--){
	    int x = ((*it).id().rawId() >> ijk) & 0x1;
	    std::cout << x;
	  }
	  std::cout << " " << ((*it).id().rawId() >> 28) << " " << (((*it).id().rawId() >> 25) & 0x7) << std::endl;
	}
	for(edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit>>::const_iterator it=hfreco->begin(); it!=hfreco->end(); ++it){
	  std::cout << (*it).id().rawId() <<  " ";
	  for(int ijk=31; ijk>=0; ijk--){
	    int x = ((*it).id().rawId() >> ijk) & 0x1;
	    std::cout << x;
	  }
	  std::cout << " " << ((*it).id().rawId() >> 28) << " " << (((*it).id().rawId() >> 25) & 0x7) << std::endl;
	}
	for(edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit>>::const_iterator it=horeco->begin(); it!=horeco->end(); ++it){
	  std::cout << (*it).id().rawId() <<  " ";
	  for(int ijk=31; ijk>=0; ijk--){
	    int x = ((*it).id().rawId() >> ijk) & 0x1;
	    std::cout << x;
	  }
	  std::cout << " " << ((*it).id().rawId() >> 28) << " " << (((*it).id().rawId() >> 25) & 0x7) << std::endl;
	}
      }
      
      /*if(tpfjet_ntwrs_ == 0 || ppfjet_ntwrs_ == 0){
	std::cout << "Event " << iEvent.id().event() << ": " << tpfjet_ntwrs_ << " " << ppfjet_ntwrs_ << std::endl;
      }
      else{
	//std::cout << iEvent.id().event() << ": " << ppfjet_ntwrs_ << std::endl;
	}*/
      
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
    //calo_tree_->Branch("tcalojet_genpt",&tcalojet_genpt_, "tcalojet_genpt/F");
    //calo_tree_->Branch("tcalojet_genp",&tcalojet_genp_, "tcalojet_genp/F");
    //calo_tree_->Branch("tcalojet_gendr",&tcalojet_gendr_, "tcalojet_gendr/F");
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
    //calo_tree_->Branch("pcalojet_genpt",&pcalojet_genpt_, "pcalojet_genpt/F");
    //calo_tree_->Branch("pcalojet_genp",&pcalojet_genp_, "pcalojet_genp/F");
    //calo_tree_->Branch("pcalojet_gendr",&pcalojet_gendr_, "pcalojet_gendr/F");
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

    pf_tree_->Branch("tpfjet_pt",&tpfjet_pt_, "tpfjet_pt/F");
    pf_tree_->Branch("tpfjet_p",&tpfjet_p_, "tpfjet_p/F");
    pf_tree_->Branch("tpfjet_eta",&tpfjet_eta_, "tpfjet_eta/F");
    pf_tree_->Branch("tpfjet_phi",&tpfjet_phi_, "tpfjet_phi/F");
    pf_tree_->Branch("tpfjet_scale",&tpfjet_scale_, "tpfjet_scale/F");
    //pf_tree_->Branch("tpfjet_genpt",&tpfjet_genpt_, "tpfjet_genpt/F");
    //pf_tree_->Branch("tpfjet_genp",&tpfjet_genp_, "tpfjet_genp/F");
    //pf_tree_->Branch("tpfjet_gendr",&tpfjet_gendr_, "tpfjet_gendr/F");
    pf_tree_->Branch("tpfjet_photonE",&tpfjet_photonE_, "tpfjet_photonE/F");
    pf_tree_->Branch("tpfjet_photonN",&tpfjet_photonN_, "tpfjet_photonN/I");
    pf_tree_->Branch("tpfjet_electronE",&tpfjet_electronE_, "tpfjet_electronE/F");
    pf_tree_->Branch("tpfjet_electronN",&tpfjet_electronN_, "tpfjet_electronN/I");
    pf_tree_->Branch("tpfjet_muonE",&tpfjet_muonE_, "tpfjet_muonE/F");
    pf_tree_->Branch("tpfjet_muonN",&tpfjet_muonN_, "tpfjet_muonN/I");
    pf_tree_->Branch("tpfjet_HFEME",&tpfjet_HFEME_, "tpfjet_HFEME/F");
    pf_tree_->Branch("tpfjet_HFEMN",&tpfjet_HFEMN_, "tpfjet_HFEMN/I");
    pf_tree_->Branch("tpfjet_ntwrs",&tpfjet_ntwrs_, "tpfjet_ntwrs/I");
    pf_tree_->Branch("tpfjet_twr_ieta",tpfjet_twr_ieta_, "tpfjet_twr_ieta[tpfjet_ntwrs]/I");
    pf_tree_->Branch("tpfjet_twr_hade",tpfjet_twr_hade_, "tpfjet_twr_hade[tpfjet_ntwrs]/F");
    pf_tree_->Branch("tpfjet_twr_frac",tpfjet_twr_frac_, "tpfjet_twr_frac[tpfjet_ntwrs]/F");
    pf_tree_->Branch("ppfjet_pt",&ppfjet_pt_, "ppfjet_pt/F");
    pf_tree_->Branch("ppfjet_p",&ppfjet_p_, "ppfjet_p/F");
    pf_tree_->Branch("ppfjet_eta",&ppfjet_eta_, "ppfjet_eta/F");
    pf_tree_->Branch("ppfjet_phi",&ppfjet_phi_, "ppfjet_phi/F");
    pf_tree_->Branch("ppfjet_scale",&ppfjet_scale_, "ppfjet_scale/F");
    //pf_tree_->Branch("ppfjet_genpt",&ppfjet_genpt_, "ppfjet_genpt/F");
    //pf_tree_->Branch("ppfjet_genp",&ppfjet_genp_, "ppfjet_genp/F");
    //pf_tree_->Branch("ppfjet_gendr",&ppfjet_gendr_, "ppfjet_gendr/F");
    pf_tree_->Branch("ppfjet_photonE",&ppfjet_photonE_, "ppfjet_photonE/F");
    pf_tree_->Branch("ppfjet_photonN",&ppfjet_photonN_, "ppfjet_photonN/I");
    pf_tree_->Branch("ppfjet_electronE",&ppfjet_electronE_, "ppfjet_electronE/F");
    pf_tree_->Branch("ppfjet_electronN",&ppfjet_electronN_, "ppfjet_electronN/I");
    pf_tree_->Branch("ppfjet_muonE",&ppfjet_muonE_, "ppfjet_muonE/F");
    pf_tree_->Branch("ppfjet_muonN",&ppfjet_muonN_, "ppfjet_muonN/I");
    pf_tree_->Branch("ppfjet_HFEME",&ppfjet_HFEME_, "ppfjet_HFEME/F");
    pf_tree_->Branch("ppfjet_HFEMN",&ppfjet_HFEMN_, "ppfjet_HFEMN/I");
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
