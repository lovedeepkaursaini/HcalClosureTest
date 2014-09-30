#include <vector>

void runPFJetCorr()
{
  gSystem->Load("../../../../lib/slc5_amd64_gcc462/libHcalClosureTestDataFormat.so");
  gROOT->ProcessLine(".L loader.C+");
  

  TChain* tree = new TChain("pf_dijettree");
  TString input = "/eos/uscms/store/user/dgsheffi/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/DijetCalibration_dEta-1p5_Et-10_3rdEt-50/e02441adc4b1f61e7a01cc47fa7cba8d/tree_*.root";
  //  TString input = "/eos/uscms/store/user/dgsheffi/QCD_Pt-15to3000_TuneD6T_Flat_8TeV_pythia6/DijetCalibration_dEta-1p5_Et-10_3rdEt-50/e02441adc4b1f61e7a01cc47fa7cba8d/tree_10_1_WgX.root";
  cout << "Opening file: " << input << endl;
  tree->Add(input);
  cout << "File opened." << endl;

  TString output = "/uscms_data/d3/dgsheffi/HCal/corrections/QCD_Pt-15to3000_TuneD6R_Flat_8TeV_pythia6_dEta-0p5_Et-20_3rdEt-15_weight0p035.root";

  DijetRespCorrData data;

  const int MAXIETA = 41;
  const int NUMTOWERS = 83;

  float tjet_eta_, tjet_phi_, tjet_E_;
  float tjet_unkown_E_, tjet_electron_E_, tjet_muon_E_, tjet_photon_E_;
  int tjet_had_n_;
  vector<float>* tjet_had_EcalE_;
  vector<int>* tjet_had_id_;
  vector<int>* tjet_had_ntwrs_;
  vector<int>* tjet_had_candtrackind_;
  int tjet_ntwrs_;
  vector<int>* tjet_twr_ieta_;
  vector<int>* tjet_twr_candtrackind_;
  vector<int>* tjet_twr_clusterind_;
  vector<float>* tjet_twr_hade_;
  vector<float>* tjet_twr_frac_;
  int tjet_cluster_n_;
  vector<float>* tjet_cluster_dR_;
  int tjet_ncandtracks_;
  vector<float>* tjet_candtrack_px_;
  vector<float>* tjet_candtrack_py_;
  vector<float>* tjet_candtrack_pz_;
  vector<float>* tjet_candtrack_EcalE_;
  float pjet_eta_, pjet_phi_, pjet_E_;
  float pjet_unkown_E_, pjet_electron_E_, pjet_muon_E_, pjet_photon_E_;
  int pjet_had_n_;
  vector<float>* pjet_had_EcalE_;
  vector<int>* pjet_had_id_;
  vector<int>* pjet_had_ntwrs_;
  vector<int>* pjet_had_candtrackind_;
  int pjet_ntwrs_;
  vector<int>* pjet_twr_ieta_;
  vector<int>* pjet_twr_candtrackind_;
  vector<int>* pjet_twr_clusterind_;
  vector<float>* pjet_twr_hade_;
  vector<float>* pjet_twr_frac_;
  int pjet_cluster_n_;
  vector<float>* pjet_cluster_dR_;
  int pjet_ncandtracks_;
  vector<float>* pjet_candtrack_px_;
  vector<float>* pjet_candtrack_py_;
  vector<float>* pjet_candtrack_pz_;
  vector<float>* pjet_candtrack_EcalE_;
  float thirdjet_px_, thirdjet_py_;
  float dijet_deta_;
  float weight_;

  tree->SetBranchAddress("tpfjet_eta",&tjet_eta_);
  tree->SetBranchAddress("tpfjet_phi",&tjet_phi_);
  tree->SetBranchAddress("tpfjet_E",&tjet_E_);
  tree->SetBranchAddress("tpfjet_unkown_E",&tjet_unkown_E_);
  tree->SetBranchAddress("tpfjet_electron_E",&tjet_electron_E_);
  tree->SetBranchAddress("tpfjet_muon_E",&tjet_muon_E_);
  tree->SetBranchAddress("tpfjet_photon_E",&tjet_photon_E_);
  tree->SetBranchAddress("tpfjet_had_n",&tjet_had_n_);
  tree->SetBranchAddress("tpfjet_had_EcalE",&tjet_had_EcalE_);
  tree->SetBranchAddress("tpfjet_had_id",&tjet_had_id_);
  tree->SetBranchAddress("tpfjet_had_ntwrs",&tjet_had_ntwrs_);
  tree->SetBranchAddress("tpfjet_had_candtrackind",&tjet_had_candtrackind_);
  tree->SetBranchAddress("tpfjet_ntwrs",&tjet_ntwrs_);
  tree->SetBranchAddress("tpfjet_twr_ieta",&tjet_twr_ieta_);
  tree->SetBranchAddress("tpfjet_twr_hade",&tjet_twr_hade_);
  tree->SetBranchAddress("tpfjet_twr_frac",&tjet_twr_frac_);
  tree->SetBranchAddress("tpfjet_twr_candtrackind",&tjet_twr_candtrackind_);
  tree->SetBranchAddress("tpfjet_twr_clusterind",&tjet_twr_clusterind_);
  tree->SetBranchAddress("tpfjet_cluster_n",&tjet_cluster_n_);
  tree->SetBranchAddress("tpfjet_cluster_dR",&tjet_cluster_dR_);
  tree->SetBranchAddress("tpfjet_ncandtracks",&tjet_ncandtracks_);
  tree->SetBranchAddress("tpfjet_candtrack_px",&tjet_candtrack_px_);
  tree->SetBranchAddress("tpfjet_candtrack_py",&tjet_candtrack_py_);
  tree->SetBranchAddress("tpfjet_candtrack_pz",&tjet_candtrack_pz_);
  tree->SetBranchAddress("tpfjet_candtrack_EcalE",&tjet_candtrack_EcalE_);
  tree->SetBranchAddress("ppfjet_eta",&pjet_eta_);
  tree->SetBranchAddress("ppfjet_phi",&pjet_phi_);
  tree->SetBranchAddress("ppfjet_E",&pjet_E_);
  tree->SetBranchAddress("ppfjet_unkown_E",&pjet_unkown_E_);
  tree->SetBranchAddress("ppfjet_electron_E",&pjet_electron_E_);
  tree->SetBranchAddress("ppfjet_muon_E",&pjet_muon_E_);
  tree->SetBranchAddress("ppfjet_photon_E",&pjet_photon_E_);
  tree->SetBranchAddress("ppfjet_had_n",&pjet_had_n_);
  tree->SetBranchAddress("ppfjet_had_EcalE",&pjet_had_EcalE_);
  tree->SetBranchAddress("ppfjet_had_id",&pjet_had_id_);
  tree->SetBranchAddress("ppfjet_had_ntwrs",&pjet_had_ntwrs_);
  tree->SetBranchAddress("ppfjet_had_candtrackind",&pjet_had_candtrackind_);
  tree->SetBranchAddress("ppfjet_ntwrs",&pjet_ntwrs_);
  tree->SetBranchAddress("ppfjet_twr_ieta",&pjet_twr_ieta_);
  tree->SetBranchAddress("ppfjet_twr_hade",&pjet_twr_hade_);
  tree->SetBranchAddress("ppfjet_twr_frac",&pjet_twr_frac_);
  tree->SetBranchAddress("ppfjet_twr_candtrackind",&pjet_twr_candtrackind_);
  tree->SetBranchAddress("ppfjet_twr_clusterind",&pjet_twr_clusterind_);
  tree->SetBranchAddress("ppfjet_cluster_n",&pjet_cluster_n_);
  tree->SetBranchAddress("ppfjet_cluster_dR",&pjet_cluster_dR_);
  tree->SetBranchAddress("ppfjet_ncandtracks",&pjet_ncandtracks_);
  tree->SetBranchAddress("ppfjet_candtrack_px",&pjet_candtrack_px_);
  tree->SetBranchAddress("ppfjet_candtrack_py",&pjet_candtrack_py_);
  tree->SetBranchAddress("ppfjet_candtrack_pz",&pjet_candtrack_pz_);
  tree->SetBranchAddress("ppfjet_candtrack_EcalE",&pjet_candtrack_EcalE_);
  tree->SetBranchAddress("pf_thirdjet_px",&thirdjet_px_);
  tree->SetBranchAddress("pf_thirdjet_py",&thirdjet_py_);
  tree->SetBranchAddress("pf_dijet_deta",&dijet_deta_);
  tree->SetBranchAddress("pf_weight",&weight_);

  TH1D* h_PassSel_ = new TH1D("h_PassSelection", "Selection Pass Failures",256,-0.5,255.5);
  int fails = 0;

  int nEvents = tree->GetEntries();
  //nEvents = 100;
  cout << "Running over " << nEvents << " events" << endl;
  for(int iEvent=0; iEvent<nEvents; iEvent++){
    if(iEvent % 1000 == 0){
      cout << "Processing event " << iEvent << endl;
    }
    tree->GetEntry(iEvent);

    int passSel = 0;

    if(tjet_ntwrs_ == 0 || pjet_ntwrs_ == 0){
      fails++;
      passSel | 0x80;
      //cout << "Fails: " << iEvent << " " << tjet_ntwrs_ << " " << pjet_ntwrs_ << endl;
      //continue;
    }
    float tjet_Et = tjet_E_/cosh(tjet_eta_);
    float pjet_Et = pjet_E_/cosh(pjet_eta_);
    float minSumJetEt_ = 40.0;//40.0;
    float minJetEt_ = 20.0;//20.0;
    float maxThirdJetEt_ = 15.0;//15.0;
    float maxDeltaEta_ = 0.5;//0.5;
    if(tjet_Et + pjet_Et < minSumJetEt_) passSel |= 0x1;
    if(tjet_Et < minJetEt_ || pjet_Et < minJetEt_) passSel |= 0x2;
    if(sqrt(thirdjet_px_*thirdjet_px_ + thirdjet_py_*thirdjet_py_) > maxThirdJetEt_) passSel |= 0x4;
    if(dijet_deta_ > maxDeltaEta_) passSel |= 0x8;
    
    h_PassSel_->Fill(passSel);
    if(passSel) continue;

    DijetRespCorrDatum datum;
    
    // Fill datum
    datum.SetWeight(weight_);

    float sumt = 0;
    datum.SetTagEta(tjet_eta_);
    datum.SetTagPhi(tjet_phi_);
    for(int i=0; i<tjet_ntwrs_; i++){
      if(tjet_twr_hade_->at(i) > 0.0 && (tjet_twr_clusterind_->at(i) || tjet_cluster_dR_->at(tjet_twr_clusterind_->at(i)))){
	datum.AddTagHcalE(tjet_twr_hade_->at(i)*tjet_twr_frac_->at(i),tjet_twr_ieta_->at(i));
	sumt += tjet_twr_hade_->at(i)*tjet_twr_frac_->at(i);
      }
    }
    /*datum.SetCandTrackN(tjet_ncandtracks_);
    for(int i=0; i<tjet_ncandtracks_; i++){
      datum.AddCandTrackP(sqrt(tjet_candtrack_px_->at(i)*tjet_candtrack_px_->at(i) + tjet_candtrack_py_->at(i)*tjet_candtrack_py_->at(i) + tjet_candtrack_pz_->at(i)*tjet_candtrack_pz_->at(i)));
      datum.AddCandTrackEcalE(tjet_candtrack_EcalE_->at(i));
      map<Int_t, Double_t> clusterEnergies;
      for(int j=0; j<tjet_ntwrs_; j++){
	if(tjet_twr_candtrackind_->at(j) == i){
	  if(tjet_twr_hade_->at(j) > 0.0){
	    assert(tjet_twr_ieta_->at(j)<=41 && tjet_twr_ieta_->at(j)>=-41 && tjet_twr_ieta_->at(j)!=0);
	    //clusterEnergies[tjet_twr_ieta_[j]] = tjet_twr_hade_[j];
	    clusterEnergies[tjet_twr_ieta_->at(j)] = tjet_twr_hade_->at(j)*tjet_twr_frac_->at(j);
	  }
	}
      }
      datum.AddCandTrackHcalE(clusterEnergies);
      }*/
    float tjet_had_EcalE_total = 0;
    float tjet_had_candNoRecHits_E = 0;
    for(int iHad=0; iHad<tjet_had_n_; iHad++){
      if(tjet_had_id_->at(iHad) < 2) tjet_had_EcalE_total += tjet_had_EcalE_->at(iHad);
      if(tjet_had_ntwrs_->at(iHad) == 0 && tjet_had_candtrackind_->at(iHad) > -1){
	tjet_had_candNoRecHits_E += sqrt(tjet_candtrack_px_->at(tjet_had_candtrackind_->at(iHad))*tjet_candtrack_px_->at(tjet_had_candtrackind_->at(iHad)) + tjet_candtrack_py_->at(tjet_had_candtrackind_->at(iHad))*tjet_candtrack_py_->at(tjet_had_candtrackind_->at(iHad)) + tjet_candtrack_pz_->at(tjet_had_candtrackind_->at(iHad))*tjet_candtrack_pz_->at(tjet_had_candtrackind_->at(iHad))) - tjet_had_EcalE_->at(iHad);
      }
    }
    datum.SetTagEcalE(tjet_unkown_E_ + tjet_electron_E_ + tjet_muon_E_ + tjet_photon_E_ + tjet_had_EcalE_total + tjet_had_candNoRecHits_E);

    float sump = 0;
    datum.SetProbeEta(pjet_eta_);
    datum.SetProbePhi(pjet_phi_);
    for(int i=0; i<pjet_ntwrs_; i++){
      //cout << pjet_twr_clusterind_->size() << " " << pjet_twr_clusterind_->at(i) << " " << pjet_cluster_n_ << endl;
      if(pjet_twr_hade_->at(i) > 0.0 && (pjet_twr_clusterind_->at(i) || pjet_cluster_dR_->at(pjet_twr_clusterind_->at(i)))){
	datum.AddProbeHcalE(pjet_twr_hade_->at(i)*pjet_twr_frac_->at(i),pjet_twr_ieta_->at(i));
	sump += pjet_twr_hade_->at(i)*pjet_twr_frac_->at(i);
      }
    }
    /*datum.SetCandTrackN(pjet_ncandtracks_);
    for(int i=0; i<pjet_ncandtracks_; i++){
      datum.AddCandTrackP(sqrt(pjet_candtrack_px_->at(i)*pjet_candtrack_px_->at(i) + pjet_candtrack_py_->at(i)*pjet_candtrack_py_->at(i) + pjet_candtrack_pz_->at(i)*pjet_candtrack_pz_->at(i)));
      datum.AddCandTrackEcalE(pjet_candtrack_EcalE_->at(i));
      map<Int_t, Double_t> clusterEnergies;
      for(int j=0; j<pjet_ntwrs_; j++){
	if(pjet_twr_candtrackind_->at(j) == i){
	  if(pjet_twr_hade_->at(j) > 0.0){
	    assert(pjet_twr_ieta_->at(j)<=41 && pjet_twr_ieta_->at(j)>=-41 && pjet_twr_ieta_->at(j)!=0);
	    //clusterEnergies[tjet_twr_ieta_->at(j)] = tjet_twr_hade_->at(j);
	    clusterEnergies[pjet_twr_ieta_->at(j)] = pjet_twr_hade_->at(j)*pjet_twr_frac_->at(j);
	  }
	}
      }
      datum.AddCandTrackHcalE(clusterEnergies);
      }*/
    float pjet_had_EcalE_total = 0;
    float pjet_had_candNoRecHits_E = 0;
    for(int iHad=0; iHad<pjet_had_n_; iHad++){
      if(pjet_had_id_->at(iHad) < 2) pjet_had_EcalE_total += pjet_had_EcalE_->at(iHad);
      if(pjet_had_ntwrs_->at(iHad) == 0 && pjet_had_candtrackind_->at(iHad) > -1){
	pjet_had_candNoRecHits_E += sqrt(pjet_candtrack_px_->at(pjet_had_candtrackind_->at(iHad))*pjet_candtrack_px_->at(pjet_had_candtrackind_->at(iHad)) + pjet_candtrack_py_->at(pjet_had_candtrackind_->at(iHad))*pjet_candtrack_py_->at(pjet_had_candtrackind_->at(iHad)) + pjet_candtrack_pz_->at(pjet_had_candtrackind_->at(iHad))*pjet_candtrack_pz_->at(pjet_had_candtrackind_->at(iHad))) - pjet_had_EcalE_->at(iHad);
      }
    }
    datum.SetProbeEcalE(pjet_unkown_E_ + pjet_electron_E_ + pjet_muon_E_ + pjet_photon_E_ + pjet_had_EcalE_total);

    datum.SetThirdJetPx(thirdjet_px_);
    datum.SetThirdJetPy(thirdjet_py_);

    if(sumt == 0 || sump == 0){
      fails++;
      continue;
    }
    

    data.push_back(datum);
  }

  cout << data.GetSize() << " data" << endl;
  
  cout << "Passes: " << nEvents - fails << " Fails: " << fails << endl;
  cout << "Do CandTrack? " << data.GetDoCandTrackEnergyDiff() << endl;
  data.SetDoCandTrackEnergyDiff(false);
  cout << "Do CandTrack? " << data.GetDoCandTrackEnergyDiff() << endl;

  //return;
  
  TH1D* hist = data.doFit("h_corr","Response Corrections");
  hist->GetXaxis()->SetTitle("ieta");
  hist->GetYaxis()->SetTitle("response corrections");

  TFile* fout = new TFile(output,"RECREATE");
  fout->cd();
  hist->Write();
  h_PassSel_->Write();
  fout->Close();

  cout << "Passes: " << nEvents - fails << " Fails: " << fails << endl;
  cout << "Events that passed cuts: " << h_PassSel_->GetBinContent(1) << endl;

  return;
}
