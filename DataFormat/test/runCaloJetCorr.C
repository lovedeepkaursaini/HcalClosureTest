void runCaloJetCorr()
{
  gSystem->Load("../../../../lib/slc5_amd64_gcc462/libHcalClosureTestDataFormat.so");

  TChain* tree = new TChain("calo_dijettree");
  TString input = "/uscms_data/d3/dgsheffi/HCal/Trees/Pion_Pt-50_*_fullRecHits.root";
  tree->Add(input);

  TString output = "/uscms_data/d3/dgsheffi/HCal/caloJetCorr_fullRecHits.root";

  DijetRespCorrData data;

  const int MAXIETA = 41;
  const int NUMTOWERS = 83;

  float tjet_eta_, tjet_phi_;
  float tjet_EBE_, tjet_EEE_;
  int tjet_ntwrs_;
  int tjet_twr_ieta_[100];
  float tjet_twr_hade_[100];
  float pjet_eta_, pjet_phi_;
  float pjet_EBE_, pjet_EEE_;
  int pjet_ntwrs_;
  int pjet_twr_ieta_[100];
  float pjet_twr_hade_[100];
  float thirdjet_px_, thirdjet_py_;

  tree->SetBranchAddress("tjet_eta",&tjet_eta_);
  tree->SetBranchAddress("tjet_phi",&tjet_phi_);
  tree->SetBranchAddress("tjet_EBE",&tjet_EBE_);
  tree->SetBranchAddress("tjet_EEE",&tjet_EEE_);
  tree->SetBranchAddress("tjet_ntwrs",&tjet_ntwrs_);
  tree->SetBranchAddress("tjet_twr_ieta",&tjet_twr_ieta_);
  tree->SetBranchAddress("tjet_twr_hade",&tjet_twr_hade_);
  tree->SetBranchAddress("pjet_eta",&pjet_eta_);
  tree->SetBranchAddress("pjet_phi",&pjet_phi_);
  tree->SetBranchAddress("pjet_EBE",&pjet_EBE_);
  tree->SetBranchAddress("pjet_EEE",&pjet_EEE_);
  tree->SetBranchAddress("pjet_ntwrs",&pjet_ntwrs_);
  tree->SetBranchAddress("pjet_twr_ieta",&pjet_twr_ieta_);
  tree->SetBranchAddress("pjet_twr_hade",&pjet_twr_hade_);
  tree->SetBranchAddress("thirdjet_px",&thirdjet_px_);
  tree->SetBranchAddress("thirdjet_py",&thirdjet_py_);

  int nEvents = tree->GetEntries();
  //nEvents = 10;
  cout << "Running over " << nEvents << " events" << endl;
  for(int iEvent=0; iEvent<nEvents; iEvent++){
    if(iEvent % 100 == 0){
      cout << "Processing event " << iEvent << endl;
    }
    tree->GetEntry(iEvent);

    DijetRespCorrDatum datum;
    
    // Fill datum

    datum.SetTagEta(tjet_eta_);
    datum.SetTagPhi(tjet_phi_);
    for(int i=0; i<tjet_ntwrs_; i++){
      datum.AddTagHcalE(tjet_twr_hade_[i],tjet_twr_ieta_[i]);
    }
    datum.SetTagEcalE(tjet_EBE_ + tjet_EEE_);

    datum.SetProbeEta(pjet_eta_);
    datum.SetProbePhi(pjet_phi_);
    for(int i=0; i<pjet_ntwrs_; i++){
      datum.AddProbeHcalE(pjet_twr_hade_[i],pjet_twr_ieta_[i]);
    }
    datum.SetProbeEcalE(pjet_EBE_ + pjet_EEE_);

    datum.SetThirdJetPx(thirdjet_px_);
    datum.SetThirdJetPy(thirdjet_py_);
    
    data.push_back(datum);
  }

  cout << data.GetSize() << " data" << endl;
    
  TH1D* hist = data.doFit("hcorr","Response Corrections");

  TFile* fout = new TFile(output,"RECREATE");
  fout->cd();
  hist->Write();
  fout->Close();

  return;
}
