void runPFJetCorr()
{
  gSystem->Load("../../../../lib/slc5_amd64_gcc462/libHcalClosureTestDataFormat.so");

  TChain* tree = new TChain("pf_dijettree");
  TString input = "/uscms_data/d3/dgsheffi/HCal/Trees/Pion_Pt-50_noHF.root";
  //TString input = "/uscms_data/d3/dgsheffi/HCal/Trees/Pion_Pt-50_*_fullRecHits.root";
  tree->Add(input);

  //TString output = "/uscms_data/d3/dgsheffi/HCal/pfJetCorr.root";
  TString output = "/uscms_data/d3/dgsheffi/HCal/pfJetCorr_noHF_freeHF.root";

  DijetRespCorrData data;

  const int MAXIETA = 41;
  const int NUMTOWERS = 83;

  float tjet_eta_, tjet_phi_;
  float tjet_unkown_E_, tjet_chHad_E_, tjet_electron_E_, tjet_muon_E_;
  float tjet_photon_E_, tjet_Had0_E_, tjet_HFHad_E_, tjet_HFEM_E_;
  int tjet_ntwrs_;
  int tjet_twr_ieta_[1000];
  float tjet_twr_hade_[1000], tjet_twr_frac_[1000];
  float pjet_eta_, pjet_phi_;
  float pjet_unkown_E_, pjet_chHad_E_, pjet_electron_E_, pjet_muon_E_;
  float pjet_photon_E_, pjet_Had0_E_, pjet_HFHad_E_, pjet_HFEM_E_;
  int pjet_ntwrs_;
  int pjet_twr_ieta_[1000];
  float pjet_twr_hade_[1000], pjet_twr_frac_[1000];
  float thirdjet_px_, thirdjet_py_;

  tree->SetBranchAddress("tpfjet_eta",&tjet_eta_);
  tree->SetBranchAddress("tpfjet_phi",&tjet_phi_);
  tree->SetBranchAddress("tpfjet_unkown_E",&tjet_unkown_E_);
  tree->SetBranchAddress("tpfjet_chHad_E",&tjet_chHad_E_);
  tree->SetBranchAddress("tpfjet_electron_E",&tjet_electron_E_);
  tree->SetBranchAddress("tpfjet_muon_E",&tjet_muon_E_);
  tree->SetBranchAddress("tpfjet_photon_E",&tjet_photon_E_);
  tree->SetBranchAddress("tpfjet_Had0_E",&tjet_Had0_E_);
  tree->SetBranchAddress("tpfjet_HFHad_E",&tjet_HFHad_E_);
  tree->SetBranchAddress("tpfjet_HFEM_E",&tjet_HFEM_E_);
  tree->SetBranchAddress("tpfjet_ntwrs",&tjet_ntwrs_);
  tree->SetBranchAddress("tpfjet_twr_ieta",&tjet_twr_ieta_);
  tree->SetBranchAddress("tpfjet_twr_hade",&tjet_twr_hade_);
  tree->SetBranchAddress("tpfjet_twr_frac",&tjet_twr_frac_);
  tree->SetBranchAddress("ppfjet_eta",&pjet_eta_);
  tree->SetBranchAddress("ppfjet_phi",&pjet_phi_);
  tree->SetBranchAddress("ppfjet_unkown_E",&pjet_unkown_E_);
  tree->SetBranchAddress("ppfjet_chHad_E",&pjet_chHad_E_);
  tree->SetBranchAddress("ppfjet_electron_E",&pjet_electron_E_);
  tree->SetBranchAddress("ppfjet_muon_E",&pjet_muon_E_);
  tree->SetBranchAddress("ppfjet_photon_E",&pjet_photon_E_);
  tree->SetBranchAddress("ppfjet_Had0_E",&pjet_Had0_E_);
  tree->SetBranchAddress("ppfjet_HFHad_E",&pjet_HFHad_E_);
  tree->SetBranchAddress("ppfjet_HFEM_E",&pjet_HFEM_E_);
  tree->SetBranchAddress("ppfjet_ntwrs",&pjet_ntwrs_);
  tree->SetBranchAddress("ppfjet_twr_ieta",&pjet_twr_ieta_);
  tree->SetBranchAddress("ppfjet_twr_hade",&pjet_twr_hade_);
  tree->SetBranchAddress("ppfjet_twr_frac",&pjet_twr_frac_);
  tree->SetBranchAddress("pf_thirdjet_px",&thirdjet_px_);
  tree->SetBranchAddress("pf_thirdjet_py",&thirdjet_py_);

  int fails = 0;

  int nEvents = tree->GetEntries();
  //nEvents = 10;
  cout << "Running over " << nEvents << " events" << endl;
  for(int iEvent=0; iEvent<nEvents; iEvent++){
    if(iEvent % 1000 == 0){
      cout << "Processing event " << iEvent << endl;
    }
    tree->GetEntry(iEvent);

    if(tjet_ntwrs_ == 0 || pjet_ntwrs_ == 0){
      fails++;
      //cout << "Fails: " << iEvent << " " << tjet_ntwrs_ << " " << pjet_ntwrs_ << endl;
      continue;
    }

    DijetRespCorrDatum datum;
    
    // Fill datum

    float sumt = 0;
    datum.SetTagEta(tjet_eta_);
    datum.SetTagPhi(tjet_phi_);
    for(int i=0; i<tjet_ntwrs_; i++){
      if(tjet_twr_hade_[i] > 0.0){
	datum.AddTagHcalE(tjet_twr_hade_[i]*tjet_twr_frac_[i],tjet_twr_ieta_[i]);
	sumt += tjet_twr_hade_[i]*tjet_twr_frac_[i];
      }
    }
    datum.SetTagEcalE(tjet_unkown_E_ + tjet_electron_E_ + tjet_muon_E_ + tjet_photon_E_);

    float sump = 0;
    datum.SetProbeEta(pjet_eta_);
    datum.SetProbePhi(pjet_phi_);
    for(int i=0; i<pjet_ntwrs_; i++){
      if(pjet_twr_hade_[i] > 0.0){
	datum.AddProbeHcalE(pjet_twr_hade_[i]*pjet_twr_frac_[i],pjet_twr_ieta_[i]);
	sump += pjet_twr_hade_[i]*pjet_twr_frac_[i];
      }
    }
    datum.SetProbeEcalE(pjet_unkown_E_ + pjet_electron_E_ + pjet_muon_E_ + pjet_photon_E_);

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
  //return;
  
  TH1D* hist = data.doFit("hcorr","Response Corrections");
  hist->GetXaxis()->SetTitle("ieta");
  hist->GetYaxis()->SetTitle("response corrections");

  TFile* fout = new TFile(output,"RECREATE");
  fout->cd();
  hist->Write();
  fout->Close();

  cout << "Passes: " << nEvents - fails << " Fails: " << fails << endl;

  return;
}
