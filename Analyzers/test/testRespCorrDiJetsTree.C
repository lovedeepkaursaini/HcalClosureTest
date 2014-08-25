#include <vector>
#include <map>

void testRespCorrDiJetsTree()
{
  //gROOT->ProcessLine(".L loader.C+");
  
  TChain* tree = new TChain("pf_dijettree");
  TString input = "/uscms_data/d3/dgsheffi/HCal/Trees/QCD_Pt-15to3000_0030487D5E5F.root";
  tree->Add(input);

  TString output = "/uscms_data/d3/dgsheffi/HCal/Trees/validation/QCD_Pt-15to3000_0030487D5E5F.root";

  float tpfjet_E_;
  float tpfjet_eta_;
  float tpfjet_genE_;
  float tpfjet_unkown_E_, tpfjet_electron_E_, tpfjet_muon_E_, tpfjet_photon_E_;
  int tpfjet_had_n_;
  vector<float>* tpfjet_had_E_;
  vector<float>* tpfjet_had_EcalE_;
  vector<int>* tpfjet_had_id_;
  vector<int>* tpfjet_had_candtrackind_;
  int tpfjet_ntwrs_;
  vector<int>* tpfjet_twr_hadind_;
  vector<float>* tpfjet_twr_hade_;
  vector<float>* tpfjet_twr_frac_;
  int tpfjet_ncandtracks_;
  vector<float>* tpfjet_candtrack_px_;
  vector<float>* tpfjet_candtrack_py_;
  int pf_Run_;
  int pf_Lumi_;
  int pf_Event_;

  tree->SetBranchAddress("tpfjet_E",&tpfjet_E_);
  tree->SetBranchAddress("tpfjet_eta",&tpfjet_eta_);
  tree->SetBranchAddress("tpfjet_E",&tpfjet_genE_);
  tree->SetBranchAddress("tpfjet_unkown_E",&tpfjet_unkown_E_);
  tree->SetBranchAddress("tpfjet_electron_E",&tpfjet_electron_E_);
  tree->SetBranchAddress("tpfjet_muon_E",&tpfjet_muon_E_);
  tree->SetBranchAddress("tpfjet_photon_E",&tpfjet_photon_E_);
  tree->SetBranchAddress("tpfjet_had_n",&tpfjet_had_n_);
  tree->SetBranchAddress("tpfjet_had_E",&tpfjet_had_E_);
  tree->SetBranchAddress("tpfjet_had_EcalE",&tpfjet_had_EcalE_);
  tree->SetBranchAddress("tpfjet_had_id",&tpfjet_had_id_);
  tree->SetBranchAddress("tpfjet_had_candtrackind",&tpfjet_had_candtrackind_);
  tree->SetBranchAddress("tpfjet_ntwrs",&tpfjet_ntwrs_);
  tree->SetBranchAddress("tpfjet_twr_hadind",&tpfjet_twr_hadind_);
  tree->SetBranchAddress("tpfjet_twr_hade",&tpfjet_twr_hade_);
  tree->SetBranchAddress("tpfjet_twr_frac",&tpfjet_twr_frac_);
  tree->SetBranchAddress("tpfjet_ncandtracks",&tpfjet_ncandtracks_);
  tree->SetBranchAddress("tpfjet_candtrack_px",&tpfjet_candtrack_px_);
  tree->SetBranchAddress("tpfjet_candtrack_py",&tpfjet_candtrack_py_);
  tree->SetBranchAddress("pf_Run",&pf_Run_);
  tree->SetBranchAddress("pf_Lumi",&pf_Lumi_);
  tree->SetBranchAddress("pf_Event",&pf_Event_);

  TH1D* h_tpfjet_jet_Ediff_ = new TH1D("h_tpfjet_jet_Ediff_","tag (rechits - pfjet)/pfjet",200,-1,8);
  TH1D* h_tpfjet_genjet_Ediff_ = new TH1D("h_tpfjet_genjet_Ediff_","tag (rechits - genjet)/genjet",200,-1,8);
  TH1D* h_tpfjet_Ediff_rechit_ = new TH1D("h_tpfjet_Ediff_rechit","tag (rechits - cand)/cand",200,-1,8);
  TH1D* h_tpfjet_h_Ediff_rechit_ = new TH1D("h_tpfjet_h_Ediff_rechit","tag h (rechits - cand)/cand",200,-1,8);
  TH1D* h_tpfjet_h0_Ediff_rechit_ = new TH1D("h_tpfjet_h0_Ediff_rechit","tag h0 (rechits - cand)/cand",200,-1,8);
  TH1D* h_tpfjet_HFhad_Ediff_rechit_ = new TH1D("h_tpfjet_HFhad_Ediff_rechit","tag HFhad (rechits - cand)/cand",200,-1,8);
  TH1D* h_tpfjet_egammaHF_Ediff_rechit_ = new TH1D("h_tpfjet_egammaHF_Ediff_rechit","tag egammaHF (rechits - cand)/cand",200,-1,8);
  TH1D* h_tpfjet_h_Ediff_allrechit_ = new TH1D("h_tpfjet_h_Ediff_allrechit","tag h (rechits - cand)/cand all rechits",200,-1,8);
  TH1D* h_tpfjet_h_Ediff_rechit_HB_ = new TH1D("h_tpfjet_h_Ediff_rechit_HB","tag h (rechits - cand)/cand HB",200,-1,8);
  TH1D* h_tpfjet_h_Ediff_rechit_HE_ = new TH1D("h_tpfjet_h_Ediff_rechit_HE","tag h (rechits - cand)/cand HB",200,-1,8);
  TH1D* h_tpfjet_h_rechits_neg_ = new TH1D("h_tpfjet_h_rechits_neg","tag h rechits negative hits",100,0,100);
  TH1D* h_tpfjet_h_rechits_pt_ = new TH1D("h_tpfjet_h_rechits_pt","tag h rechits pt",100,0,100);
  TH1D* h_tpfjet_h_rechits_EcalEfrac_ = new TH1D("h_tpfjet_h_rechits_EcalEfrac","tag h rechits EcalEfrac",100,0,1);
  TH1D* h_tpfjet_h_norechits_neg_ = new TH1D("h_tpfjet_h_norechits_neg","tag h norechits negative hits",100,0,100);
  TH1D* h_tpfjet_h_norechits_pt_ = new TH1D("h_tpfjet_h_norechits_pt","tag h norechits pt",100,0,100);
  TH1D* h_tpfjet_h_norechits_EcalEfrac_ = new TH1D("h_tpfjet_h_norechits_EcalEfrac","tag h norechits EcalEfrac",100,0,1);
  TH1D* h_tpfjet_h_norechits_EcalEfrac_lowPt_ = new TH1D("h_tpfjet_h_norechits_EcalEfrac_lowPt","tag h norechits EcalEfrac low Pt",100,0,1);
  TH1D* h_tpfjet_h_norechits_EcalEfrac_midPt_ = new TH1D("h_tpfjet_h_norechits_EcalEfrac_midPt","tag h norechits EcalEfrac mid Pt",100,0,1);
  TH1D* h_tpfjet_h_norechits_EcalEfrac_highPt_ = new TH1D("h_tpfjet_h_norechits_EcalEfrac_highPt","tag h norechits EcalEfrac high Pt",100,0,1);
  TH2D* h_tpfjet_h_EdiffvsE_rechit_ = new TH2D("h_tpfjet_h_EdiffvsE_rechit","tag h (rechits - cand)/cand vs E",200,-1,8,200,0,500);
  TH2D* h_tpfjet_h_EdiffvsE_rechit_lowE_ = new TH2D("h_tpfjet_h_EdiffvsE_rechit_lowE","tag h (rechits - cand)/cand vs E",200,-1,8,200,0,50);
  TH2D* h_tpfjet_h0_EdiffvsE_rechit_ = new TH2D("h_tpfjet_h0_EdiffvsE_rechit","tag h0 (rechits - cand)/cand vs E",200,-1,8,200,0,50);
  TH1D* h_tpfjet_rechitE_h_highEdiff_ = new TH1D("h_tpfjet_rechitE_h_highEdiff","rechit E for high tag h (rechits - cand)/cand",200,0,1000);
  TH1D* h_tpfjet_rechitE_h_lowEdiff_ = new TH1D("h_tpfjet_rechitE_h_lowEdiff","rechit E for low tag h (rechits - cand)/cand",200,0,1000);
  TH1D* h_tpfjet_candtwrs_h_highEdiff_ = new TH1D("h_tpfjet_candtwrs_h_highEdiff","towers for high tag h (rechits - cand)/cand",200,0,400);
  TH1D* h_tpfjet_candtwrs_h_lowEdiff_ = new TH1D("h_tpfjet_candtwrs_h_lowEdiff","towers E for low tag h (rechits - cand)/cand",200,0,400);
  TH1D* h_tpfjet_twrE_h_highEdiff_ = new TH1D("h_tpfjet_twrE_h_highEdiff","tower E for high tag h (rechits - cand)/cand",200,0,100);
  TH1D* h_tpfjet_twrE_h_lowEdiff_ = new TH1D("h_tpfjet_twrE_h_lowEdiff","tower E for low tag h (rechits - cand)/cand",200,0,100);

  //TH1D* h_pf_balance_ = new TH1D("h_pf_balance","dijet balance using PF jets",200,-1,1);
  //TH1D* h_pfrechit_balance_ = new TH1D("h_pfrechit_balance","dijet balance using rechits",200,-1,1);

  int nEvents = tree->GetEntries();
  cout << "Running over " << nEvents << " events" << endl;
  //nEvents = 5;
  for(int iEvent=0; iEvent<nEvents; iEvent++){
    if(iEvent % 100 == 0){
      cout << "Processing event " << iEvent << endl;
    }
    tree->GetEntry(iEvent);

    //////////////////////////
    // Fill tag histograms
    //////////////////////////
    
    float tjet_rechit_E = 0;
    for(int i=0; i<tpfjet_ntwrs_; i++){
      if(tpfjet_twr_hade_->at(i) > 0.0) jet_rechit_E += tpfjet_twr_hade_->at(i)*tpfjet_twr_frac_->at(i);
    }
    float tjet_had_EcalE = 0;
    for(int i=0; i<tpfjet_had_n_; i++){
      tjet_had_EcalE += tpfjet_had_EcalE_->at(i);
    }
    float pfjet_total_E = jet_rechit_E + tpfjet_unkown_E_ + tpfjet_electron_E_ + tpfjet_muon_E_ + tpfjet_photon_E_ + tjet_had_EcalE;

    for(int i=0; i<tpfjet_had_n_; i++){
      int ncandtwrs = 0;
      int nNeg = 0;
      float rechit_E = 0;
      float allrechit_E = 0;
      float totalrechit_E = 0;
      for(int j=0; j<tpfjet_ntwrs_; j++){
	totalrechit_E += tpfjet_twr_hade_->at(j)*tpfjet_twr_frac_->at(j);
	if(tpfjet_twr_hadind_->at(j) == i){
	  allrechit_E += tpfjet_twr_hade_->at(j)*tpfjet_twr_frac_->at(j);
	  if(tpfjet_twr_hade_->at(j) > 0.0){
	    rechit_E += tpfjet_twr_hade_->at(j)*tpfjet_twr_frac_->at(j);
	    //if(tpfjet_had_id_->at(i) == 2) cout << i << ": " << tpfjet_twr_hade_->at(j) << endl;
	    ncandtwrs++;
	  }
	  else{
	    nNeg++;
	  }
	}
      } // Loop over rechits
      h_tpfjet_Ediff_rechit_->Fill((rechit_E - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
      
      if(tpfjet_had_id_->at(i) == 0){
	if(tpfjet_had_candtrackind_->at(i) > -1){
	  int ind = tpfjet_had_candtrackind_->at(i);
	  float track_pt = sqrt(tpfjet_candtrack_px_->at(ind)*tpfjet_candtrack_px_->at(ind) + tpfjet_candtrack_py_->at(ind)*tpfjet_candtrack_py_->at(ind));
	}
	
	h_tpfjet_h_Ediff_rechit_->Fill((rechit_E - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
	h_tpfjet_h_EdiffvsE_rechit_->Fill((rechit_E - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i),tpfjet_had_E_->at(i));
	h_tpfjet_h_EdiffvsE_rechit_lowE_->Fill((rechit_E - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i),tpfjet_had_E_->at(i));
	//if((rechit_E - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i) > 7.0 && tpfjet_had_E_->at(i) > 20.0) std::cout << pf_Run_ << ":" << pf_Lumi_ << ":" << pf_Event_ << std::endl;
	h_tpfjet_h_Ediff_allrechit_->Fill((allrechit_E - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
	if(tpfjet_eta_ < 1.4) h_tpfjet_h_Ediff_rechit_HB_->Fill((rechit_E - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
	else if(tpfjet_eta_ < 3.0) h_tpfjet_h_Ediff_rechit_HE_->Fill((rechit_E - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));

	if((rechit_E - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i) > 3.0){
	  h_tpfjet_rechitE_h_highEdiff_->Fill(totalrechit_E);
	  h_tpfjet_candtwrs_h_highEdiff_->Fill(ncandtwrs);
	  for(int j=0; j<tpfjet_ntwrs_; j++){
	    if(tpfjet_twr_hadind_->at(j) == i && tpfjet_twr_hade_->at(j) > 0.0){
	      h_tpfjet_twrE_h_highEdiff_->Fill(tpfjet_twr_hade_->at(j)*tpfjet_twr_frac_->at(j));
	    }
	  }
	}
	else if(rechit_E > 0){
	  h_tpfjet_rechitE_h_lowEdiff_->Fill(totalrechit_E);
	  h_tpfjet_candtwrs_h_lowEdiff_->Fill(ncandtwrs);
	  for(int j=0; j<tpfjet_ntwrs_; j++){
	    if(tpfjet_twr_hadind_->at(j) == i && tpfjet_twr_hade_->at(j) > 0.0){
	      h_tpfjet_twrE_h_lowEdiff_->Fill(tpfjet_twr_hade_->at(j)*tpfjet_twr_frac_->at(j));
	    }
	  }
	}

	if(rechit_E == 0){
	  h_tpfjet_h_norechits_neg_->Fill(nNeg);
	  h_tpfjet_h_norechits_pt_->Fill(track_pt);
	  h_tpfjet_h_norechits_EcalEfrac_->Fill(tpfjet_had_EcalE_->at(i)/tpfjet_had_E_->at(i));
	  if(track_pt < 1.0) h_tpfjet_h_norechits_EcalEfrac_lowPt_->Fill(tpfjet_had_EcalE_->at(i)/tpfjet_had_E_->at(i));
	  else if(track_pt > 1.0 && track_pt < 10.0) h_tpfjet_h_norechits_EcalEfrac_midPt_->Fill(tpfjet_had_EcalE_->at(i)/tpfjet_had_E_->at(i));
	  else h_tpfjet_h_norechits_EcalEfrac_highPt_->Fill(tpfjet_had_EcalE_->at(i)/tpfjet_had_E_->at(i));
	}
	else{
	  h_tpfjet_h_rechits_neg_->Fill(nNeg);
	  h_tpfjet_h_rechits_pt_->Fill(track_pt);
	  h_tpfjet_h_rechits_EcalEfrac_->Fill(tpfjet_had_EcalE_->at(i)/tpfjet_had_E_->at(i));
	}
      }
      else if(tpfjet_had_id_->at(i) == 1){
	h_tpfjet_h0_Ediff_rechit_->Fill((rechit_E - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
	h_tpfjet_h0_EdiffvsE_rechit_->Fill((rechit_E - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i),tpfjet_had_E_->at(i));
      }
      else if(tpfjet_had_id_->at(i) == 2) h_tpfjet_HFhad_Ediff_rechit_->Fill((rechit_E - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
      else if(tpfjet_had_id_->at(i) == 3) h_tpfjet_egammaHF_Ediff_rechit_->Fill((rechit_E - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
    } // Loop over candidates
  
  }

  TFile* fout = new TFile(output,"RECREATE");
  fout->cd();

  h_tpfjet_jet_Ediff_->Write();
  h_tpfjet_genjet_Ediff_->Write();
  h_tpfjet_Ediff_rechit_->Write();
  h_tpfjet_h_Ediff_rechit_->Write();
  h_tpfjet_h0_Ediff_rechit_->Write();
  h_tpfjet_HFhad_Ediff_rechit_->Write();
  h_tpfjet_egammaHF_Ediff_rechit_->Write();
  h_tpfjet_h_Ediff_allrechit_->Write();
  h_tpfjet_h_Ediff_rechit_HB_->Write();
  h_tpfjet_h_Ediff_rechit_HE_->Write();
  h_tpfjet_h_rechits_neg_->Write();
  h_tpfjet_h_rechits_pt_->Write();
  h_tpfjet_h_rechits_EcalEfrac_->Write();
  h_tpfjet_h_norechits_neg_->Write();
  h_tpfjet_h_norechits_pt_->Write();
  h_tpfjet_h_norechits_EcalEfrac_->Write();
  h_tpfjet_h_norechits_EcalEfrac_lowPt_->Write();
  h_tpfjet_h_norechits_EcalEfrac_midPt_->Write();
  h_tpfjet_h_norechits_EcalEfrac_highPt_->Write();
  h_tpfjet_h_EdiffvsE_rechit_->Write();
  h_tpfjet_h_EdiffvsE_rechit_lowE_->Write();
  h_tpfjet_h0_EdiffvsE_rechit_->Write();
  h_tpfjet_rechitE_h_highEdiff_->Write();
  h_tpfjet_rechitE_h_lowEdiff_->Write();
  h_tpfjet_candtwrs_h_highEdiff_->Write();
  h_tpfjet_candtwrs_h_lowEdiff_->Write();
  h_tpfjet_twrE_h_highEdiff_->Write();
  h_tpfjet_twrE_h_lowEdiff_->Write();

  fout->Close();
  
  return;
}
