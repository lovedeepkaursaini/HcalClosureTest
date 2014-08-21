#include <vector>
#include <map>

void testRespCorrDiJetsTree()
{
  //gROOT->ProcessLine(".L loader.C+");
  
  TChain* tree = new TChain("pf_dijettree");
  TString input = "/uscms_data/d3/dgsheffi/HCal/Trees/QCD_Pt-15to3000_0030487D5E5F.root";
  //TString input = "/uscms_data/d3/dgsheffi/HCal/Trees/tmpQCD.root";
  //TString input = "/uscms_data/d3/dgsheffi/HCal/test.root";
  //TString input = "/uscms_data/d3/dgsheffi/HCal/Trees/Pion_Pt-50.root";
  tree->Add(input);

  TString output = "/uscms_data/d3/dgsheffi/HCal/Trees/validation/QCD_Pt-15to3000_0030487D5E5F.root";
  //TString output = "/uscms_data/d3/dgsheffi/HCal/Trees/validation/test.root";
  //TString output = "/uscms_data/d3/dgsheffi/HCal/Trees/validation/Pion_Pt-50.root";

  float tpfjet_E_, tpfjet_eta_, tpfjet_phi_;
  float tpfjet_unkown_E_, tpfjet_electron_E_, tpfjet_muon_E_, tpfjet_photon_E_;
  int tpfjet_unkown_n_, tpfjet_electron_n_, tpfjet_muon_n_, tpfjet_photon_n_;
  int tpfjet_had_n_;
  vector<float>* tpfjet_had_E_;
  vector<float>* tpfjet_had_px_;
  vector<float>* tpfjet_had_py_;
  vector<float>* tpfjet_had_pz_;
  vector<float>* tpfjet_had_EcalE_;
  vector<int>* tpfjet_had_id_;
  vector<float>* tpfjet_had_E_mctruth_;
  vector<int>* tpfjet_had_mcpdgId_;
  int tpfjet_ntwrs_;
  vector<int>* tpfjet_twr_ieta_;
  vector<int>* tpfjet_twr_candtrackind_;
  vector<int>* tpfjet_twr_hadind_;
  vector<float>* tpfjet_twr_hade_;
  vector<float>* tpfjet_twr_frac_;
  int tpfjet_ncandtracks_;
  vector<float>* tpfjet_candtrack_px_;
  vector<float>* tpfjet_candtrack_py_;
  vector<float>* tpfjet_candtrack_pz_;
  int pf_Event_;

  float ppfjet_E_, ppfjet_eta_, ppfjet_phi_;
  float pf_dijet_deta_;
  float ppfjet_unkown_E_, ppfjet_electron_E_, ppfjet_muon_E_, ppfjet_photon_E_;
  int ppfjet_unkown_n_, ppfjet_electron_n_, ppfjet_muon_n_, ppfjet_photon_n_;
  int ppfjet_had_n_;
  vector<float>* ppfjet_had_E_;
  vector<float>* ppfjet_had_px_;
  vector<float>* ppfjet_had_py_;
  vector<float>* ppfjet_had_pz_;
  vector<float>* ppfjet_had_EcalE_;
  vector<int>* ppfjet_had_id_;
  vector<float>* ppfjet_had_E_mctruth_;
  vector<int>* ppfjet_had_mcpdgId_;
  int ppfjet_ntwrs_;
  vector<int>* ppfjet_twr_ieta_;
  vector<int>* ppfjet_twr_candtrackind_;
  vector<int>* ppfjet_twr_hadind_;
  vector<float>* ppfjet_twr_hade_;
  vector<float>* ppfjet_twr_frac_;
  int ppfjet_ncandtracks_;
  vector<float>* ppfjet_candtrack_px_;
  vector<float>* ppfjet_candtrack_py_;
  vector<float>* ppfjet_candtrack_pz_;

  tree->SetBranchAddress("tpfjet_E",&tpfjet_E_);
  tree->SetBranchAddress("tpfjet_eta",&tpfjet_eta_);
  tree->SetBranchAddress("tpfjet_phi",&tpfjet_phi_);
  tree->SetBranchAddress("tpfjet_unkown_E",&tpfjet_unkown_E_);
  tree->SetBranchAddress("tpfjet_electron_E",&tpfjet_electron_E_);
  tree->SetBranchAddress("tpfjet_muon_E",&tpfjet_muon_E_);
  tree->SetBranchAddress("tpfjet_photon_E",&tpfjet_photon_E_);
  tree->SetBranchAddress("tpfjet_unkown_n",&tpfjet_unkown_n_);
  tree->SetBranchAddress("tpfjet_electron_n",&tpfjet_electron_n_);
  tree->SetBranchAddress("tpfjet_muon_n",&tpfjet_muon_n_);
  tree->SetBranchAddress("tpfjet_photon_n",&tpfjet_photon_n_);
  tree->SetBranchAddress("tpfjet_had_n",&tpfjet_had_n_);
  tree->SetBranchAddress("tpfjet_had_E",&tpfjet_had_E_);
  tree->SetBranchAddress("tpfjet_had_px",&tpfjet_had_px_);
  tree->SetBranchAddress("tpfjet_had_py",&tpfjet_had_py_);
  tree->SetBranchAddress("tpfjet_had_pz",&tpfjet_had_pz_);
  tree->SetBranchAddress("tpfjet_had_EcalE",&tpfjet_had_EcalE_);
  tree->SetBranchAddress("tpfjet_had_id",&tpfjet_had_id_);
  tree->SetBranchAddress("tpfjet_had_E_mctruth",&tpfjet_had_E_mctruth_);
  tree->SetBranchAddress("tpfjet_had_mcpdgId",&tpfjet_had_mcpdgId_);
  tree->SetBranchAddress("tpfjet_ntwrs",&tpfjet_ntwrs_);
  tree->SetBranchAddress("tpfjet_twr_ieta",&tpfjet_twr_ieta_);
  tree->SetBranchAddress("tpfjet_twr_candtrackind",&tpfjet_twr_candtrackind_);
  tree->SetBranchAddress("tpfjet_twr_hadind",&tpfjet_twr_hadind_);
  tree->SetBranchAddress("tpfjet_twr_hade",&tpfjet_twr_hade_);
  tree->SetBranchAddress("tpfjet_twr_frac",&tpfjet_twr_frac_);
  tree->SetBranchAddress("tpfjet_ncandtracks",&tpfjet_ncandtracks_);
  tree->SetBranchAddress("tpfjet_candtrack_px",&tpfjet_candtrack_px_);
  tree->SetBranchAddress("tpfjet_candtrack_py",&tpfjet_candtrack_py_);
  tree->SetBranchAddress("tpfjet_candtrack_pz",&tpfjet_candtrack_pz_);
  tree->SetBranchAddress("pf_Event",&pf_Event_);
  
  tree->SetBranchAddress("ppfjet_E",&ppfjet_E_);
  tree->SetBranchAddress("ppfjet_eta",&ppfjet_eta_);
  tree->SetBranchAddress("ppfjet_phi",&ppfjet_phi_);
  tree->SetBranchAddress("pf_dijet_deta",&pf_dijet_deta_);
  tree->SetBranchAddress("ppfjet_unkown_E",&ppfjet_unkown_E_);
  tree->SetBranchAddress("ppfjet_electron_E",&ppfjet_electron_E_);
  tree->SetBranchAddress("ppfjet_muon_E",&ppfjet_muon_E_);
  tree->SetBranchAddress("ppfjet_photon_E",&ppfjet_photon_E_);
  tree->SetBranchAddress("ppfjet_unkown_n",&ppfjet_unkown_n_);
  tree->SetBranchAddress("ppfjet_electron_n",&ppfjet_electron_n_);
  tree->SetBranchAddress("ppfjet_muon_n",&ppfjet_muon_n_);
  tree->SetBranchAddress("ppfjet_photon_n",&ppfjet_photon_n_);
  tree->SetBranchAddress("ppfjet_had_n",&ppfjet_had_n_);
  tree->SetBranchAddress("ppfjet_had_E",&ppfjet_had_E_);
  tree->SetBranchAddress("ppfjet_had_px",&ppfjet_had_px_);
  tree->SetBranchAddress("ppfjet_had_py",&ppfjet_had_py_);
  tree->SetBranchAddress("ppfjet_had_pz",&ppfjet_had_pz_);
  tree->SetBranchAddress("ppfjet_had_EcalE",&ppfjet_had_EcalE_);
  tree->SetBranchAddress("ppfjet_had_id",&ppfjet_had_id_);
  tree->SetBranchAddress("ppfjet_had_E_mctruth",&ppfjet_had_E_mctruth_);
  tree->SetBranchAddress("ppfjet_had_mcpdgId",&ppfjet_had_mcpdgId_);
  tree->SetBranchAddress("ppfjet_ntwrs",&ppfjet_ntwrs_);
  tree->SetBranchAddress("ppfjet_twr_ieta",&ppfjet_twr_ieta_);
  tree->SetBranchAddress("ppfjet_twr_candtrackind",&ppfjet_twr_candtrackind_);
  tree->SetBranchAddress("ppfjet_twr_hadind",&ppfjet_twr_hadind_);
  tree->SetBranchAddress("ppfjet_twr_hade",&ppfjet_twr_hade_);
  tree->SetBranchAddress("ppfjet_twr_frac",&ppfjet_twr_frac_);
  tree->SetBranchAddress("ppfjet_ncandtracks",&ppfjet_ncandtracks_);
  tree->SetBranchAddress("ppfjet_candtrack_px",&ppfjet_candtrack_px_);
  tree->SetBranchAddress("ppfjet_candtrack_py",&ppfjet_candtrack_py_);
  tree->SetBranchAddress("ppfjet_candtrack_pz",&ppfjet_candtrack_pz_);

  TH1D* h_tpfjet_eta_ = new TH1D("h_tpfjet_eta","tag PF jet #eta",100,-5.5,5.5);
  TH1D* h_tpfjet_phi_ = new TH1D("h_tpfjet_phi","tag PF jet #phi",100,-3.141593,3.141593);
  TH1D* h_tpfjet_candn_ = new TH1D("h_tpfjet_candn","tag PF jet number of candidates",150,0,150);
  TH1D* h_tpfjet_E_ = new TH1D("h_tpfjet_E","tag PF jet energy from jet",200,0,500);
  TH1D* h_tpfjet_candE_ = new TH1D("h_tpfjet_candE","tag PF jet energy from candidates",200,0,500);
  TH1D* h_tpfjet_candE_diff_ = new TH1D("h_tpfjet_candE_diff","difference between tag PF jet energy from candidates and jet",200,-0.5,0.5);
  TH1D* h_tpfjet_rechitE_ = new TH1D("h_tpfjet_rechitE","tag PF jet energy from rechits",200,0,500);
  TH1D* h_tpfjet_rechitE_diff_ = new TH1D("h_tpfjet_rechitE_diff","difference between tag PF jet energy from rechits and jet",200,-1,4);
  TH1D* h_tpfjet_emf_ = new TH1D("h_tpfjet_emf","tag PF jet emf from candidates",100,0,1);
  TH1D* h_tpfjet_negE_ = new TH1D("h_tpfjet_negE","number of negative rechits per tag jet",200,0,1000);
  TH1D* h_tpfjet_negE_frac_ = new TH1D("h_tpfjet_negE_frac","fraction of negative rechits per tag jet",200,0,1);
  TH1D* h_tpfjet_rechit_had_E_diff_ = new TH1D("h_tpfjet_rechit_had_E_diff","difference between tag jet had and rechit",200,-1,4);
  TH1D* h_tpfjet_rechit_had_total_E_diff_ = new TH1D("h_tpfjet_rechit_had_total_E_diff","difference between tag jet had and rechit",200,-1.5,1.5);
  TH1D* h_tpfjet_rechit_had_totalind_E_diff_ = new TH1D("h_tpfjet_rechit_had_totalind_E_diff","difference between tag jet had and rechit",200,-1.5,1.5);
  TH1D* h_tpfjet_rechit_had_h_E_diff_ = new TH1D("h_tpfjet_rechit_had_h_E_diff","difference between tag jet h and rechit",200,-1,8);
  TH1D* h_tpfjet_rechit_had_h_fraccut_E_diff_ = new TH1D("h_tpfjet_rechit_had_h_fraccut_E_diff","difference between tag jet h and rechit cut on fraction",200,-1,8);
  TH1D* h_tpfjet_rechit_had_h_Ecut_E_diff_ = new TH1D("h_tpfjet_rechit_had_h_Ecut_E_diff","difference between tag jet h and rechit cut on E",200,-1,8);
  TH1D* h_tpfjet_rechit_had_h_Emc_diff_ = new TH1D("h_tpfjet_rechit_had_h_Emc_diff","difference between tag jet h mc and rechit",200,-1,8);
  TH1D* h_tpfjet_rechit_had_h_mcPi_Emc_diff_ = new TH1D("h_tpfjet_rechit_had_h_mcPi_Emc_diff","difference between tag jet h mc and rechit for MC pions",200,-1,8);
  TH1D* h_tpfjet_rechit_had_h_mcNotPi_Emc_diff_ = new TH1D("h_tpfjet_rechit_had_h_mcNotPi_Emc_diff","difference between tag jet h mc and rechit for not MC pions",200,-1,8);
  TH2D* h_tpfjet_rechit_had_h_EdiffvsE_ = new TH2D("h_tpfjet_rechit_had_h_EdiffvsE","difference between tag jet h and rechit vs h E",200,-1,8,200,0,500);
  TH1D* h_tpfjet_rechit_had_h0_E_diff_ = new TH1D("h_tpfjet_rechit_had_h0_E_diff","difference between tag jet h0 and rechit",200,-1,4);
  TH1D* h_tpfjet_rechit_had_HFhad_E_diff_ = new TH1D("h_tpfjet_rechit_had_HFhad_E_diff","difference between tag jet HF_had and rechit",200,-1,4);
  TH1D* h_tpfjet_rechit_had_egammaHF_E_diff_ = new TH1D("h_tpfjet_rechit_had_egammaHF_E_diff","difference between tag jet egammaHF and rechit",200,-1,4);
  TH1D* h_tpfjet_HFmatches_ = new TH1D("h_tpfjet_HFmatches","number of HF matches per candidate for tag jet",10,0,10);
  TH1D* h_tpfjet_h_EoverP_ = new TH1D("h_tpfjet_h_EoverP","rechit E over track P for tag jet",200,0,8);
  TH1D* h_tpfjet_h_lowE_ = new TH1D("h_tpfjet_h_lowE","rechit low E charged hadron for tag jet",200,0,1);
  TH1D* h_tpfjet_h_norechits_EcalEfrac_ = new TH1D("h_tpfjet_h_norechits_EcalEfrac","Ecal E over candidate E for no rechits for tag jet",200,0,1);
  TH1D* h_tpfjet_h_norechits_EcalEfrac_lowpt_ = new TH1D("h_tpfjet_h_norechits_EcalEfrac_lowpt","Ecal E over candidate E for no rechits for tag jet low pt",200,0,1);
  TH1D* h_tpfjet_h_norechits_EcalEfrac_highpt_ = new TH1D("h_tpfjet_h_norechits_EcalEfrac_highpt","Ecal E over candidate E for no rechits for tag jet high pt",200,0,1);
  TH1D* h_tpfjet_h_norechits_pt_ = new TH1D("h_tpfjet_h_norechits_pt","pf for no rechits for tag jet",200,0,200);
  TH1D* h_tpfjet_h_rechits_pt_ = new TH1D("h_tpfjet_h_rechits_pt","pf for rechits for tag jet",200,0,200);
  TH1D* h_tpfjet_h_norechits_candE_ = new TH1D("h_tpfjet_h_rechits_candE","cand E for no rechits for tag jet",200,0,200);
  
  TH1D* h_ppfjet_eta_ = new TH1D("h_ppfjet_eta","probe PF jet #eta",100,-5.5,5.5);
  TH1D* h_ppfjet_phi_ = new TH1D("h_ppfjet_phi","probe PF jet #phi",100,-3.141593,3.141593);
  TH1D* h_ppfjet_candn_ = new TH1D("h_ppfjet_candn","probe PF jet number of candidates",150,0,150);
  TH1D* h_ppfjet_E_ = new TH1D("h_ppfjet_E","probe PF jet energy from jet",200,0,500);
  TH1D* h_ppfjet_candE_ = new TH1D("h_ppfjet_candE","probe PF jet energy from candidates",200,0,500);
  TH1D* h_ppfjet_candE_diff_ = new TH1D("h_ppfjet_candE_diff","difference between probe PF jet energy from candidates and jet",200,-0.5,0.5);
  TH1D* h_ppfjet_rechitE_ = new TH1D("h_ppfjet_rechitE","probe PF jet energy from rechits",200,0,500);
  TH1D* h_ppfjet_rechitE_diff_ = new TH1D("h_ppfjet_rechitE_diff","difference between probe PF jet energy from rechits and jet",200,-1,4);
  TH1D* h_ppfjet_emf_ = new TH1D("h_ppfjet_emf","probe PF jet emf from candidates",100,0,1);
  TH1D* h_ppfjet_negE_ = new TH1D("h_ppfjet_negE","number of negative rechits per probe jet",200,0,1000);
  TH1D* h_ppfjet_negE_frac_ = new TH1D("h_ppfjet_negE_frac","fraction of negative rechits per probe jet",200,0,1);
  TH1D* h_ppfjet_rechit_had_E_diff_ = new TH1D("h_ppfjet_rechit_had_E_diff","difference between probe jet had and rechit",200,-1,4);
  TH1D* h_ppfjet_rechit_had_h_E_diff_ = new TH1D("h_ppfjet_rechit_had_h_E_diff","difference between probe jet h and rechit",200,-1,8);
  TH1D* h_ppfjet_rechit_had_h_fraccut_E_diff_ = new TH1D("h_ppfjet_rechit_had_h_fraccut_E_diff","difference between probe jet h and rechit cut on fraction",200,-1,8);
  TH1D* h_ppfjet_rechit_had_h_Ecut_E_diff_ = new TH1D("h_ppfjet_rechit_had_h_Ecut_E_diff","difference between probe jet h and rechit cut on E",200,-1,8);
  TH1D* h_ppfjet_rechit_had_h_Emc_diff_ = new TH1D("h_ppfjet_rechit_had_h_Emc_diff","difference between probe jet h mc and rechit",200,-1,8);
  TH1D* h_ppfjet_rechit_had_h_mcPi_Emc_diff_ = new TH1D("h_ppfjet_rechit_had_h_mcPi_Emc_diff","difference between probe jet h mc and rechit for MC pions",200,-1,8);
  TH1D* h_ppfjet_rechit_had_h_mcNotPi_Emc_diff_ = new TH1D("h_ppfjet_rechit_had_h_mcNotPi_Emc_diff","difference between probe jet h mc and rechit for not MC pions",200,-1,8);
  TH2D* h_ppfjet_rechit_had_h_EdiffvsE_ = new TH2D("h_ppfjet_rechit_had_h_EdiffvsE","difference between tag jet h and rechit vs h E",200,-1,8,200,0,500);
  TH1D* h_ppfjet_rechit_had_h0_E_diff_ = new TH1D("h_ppfjet_rechit_had_h0_E_diff","difference between probe jet h0 and rechit",200,-1,4);
  TH1D* h_ppfjet_rechit_had_HFhad_E_diff_ = new TH1D("h_ppfjet_rechit_had_HFhad_E_diff","difference between probe jet HF_had and rechit",200,-1,4);
  TH1D* h_ppfjet_rechit_had_egammaHF_E_diff_ = new TH1D("h_ppfjet_rechit_had_egammaHF_E_diff","difference between probe jet egammaHF and rechit",200,-1,4);
  TH1D* h_ppfjet_HFmatches_ = new TH1D("h_ppfjet_HFmatches","number of HF matches per candidate for probe jet",10,0,10);
  TH1D* h_ppfjet_h_EoverP_ = new TH1D("h_ppfjet_h_EoverP","rechit E over track P for probe jet",200,0,8);
  TH1D* h_ppfjet_h_lowE_ = new TH1D("h_ppfjet_h_lowE","rechit low E charged hadron for probe jet",200,0,1);
  TH1D* h_ppfjet_h_norechits_EcalEfrac_ = new TH1D("h_ppfjet_h_norechits_EcalEfrac","Ecal E over candidate E for no rechits for probe jet",200,0,1);
  TH1D* h_ppfjet_h_norechits_EcalEfrac_lowpt_ = new TH1D("h_ppfjet_h_norechits_EcalEfrac_lowpt","Ecal E over candidate E for no rechits for probe jet low pt",200,0,1);
  TH1D* h_ppfjet_h_norechits_EcalEfrac_highpt_ = new TH1D("h_ppfjet_h_norechits_EcalEfrac_highpt","Ecal E over candidate E for no rechits for probe jet high pt",200,0,1);
  TH1D* h_ppfjet_h_norechits_pt_ = new TH1D("h_ppfjet_h_norechits_pt","pf for no rechits for probe jet",200,0,200);
  TH1D* h_ppfjet_h_rechits_pt_ = new TH1D("h_ppfjet_h_rechits_pt","pf for rechits for probe jet",200,0,200);
  TH1D* h_ppfjet_h_norechits_candE_ = new TH1D("h_ppfjet_h_rechits_candE","cand E for no rechits for probe jet",200,0,200);

  TH1D* h_pf_dijet_deta_ = new TH1D("h_pf_dijet_deta","PF dijet #Delta|eta|",100,0,2);

  int tpfjet_norechits = 0;
  int ppfjet_norechits = 0;
  int dijet_norechits  = 0;

  map<Int_t,Int_t> tpfjet_pdgIdMap;
  map<Int_t,Int_t> ppfjet_pdgIdMap;

  int nEvents = tree->GetEntries();
  cout << "Running over " << nEvents << " events" << endl;
  //nEvents = 5;
  for(int iEvent=0; iEvent<nEvents; iEvent++){
    if(iEvent % 100 == 0){
      cout << "Processing event " << iEvent << endl;
    }
    tree->GetEntry(iEvent);

    if(pf_Event_ == 9404912){
      cout << tpfjet_eta_ << " " << tpfjet_phi_ << endl;
      cout << ppfjet_eta_ << " " << ppfjet_phi_ << endl;
    }

    if(tpfjet_had_n_ != tpfjet_had_E_->size()) cout << "tpfjet_had_n_ != tpfjet_had_E_->size() " << tpfjet_had_n_ << " " << tpfjet_had_E_->size() << endl;
    if(tpfjet_had_n_ != tpfjet_had_px_->size()) cout << "tpfjet_had_n_ != tpfjet_had_px_->size() " << tpfjet_had_n_ << " " << tpfjet_had_px_->size() << endl;
    if(tpfjet_had_n_ != tpfjet_had_py_->size()) cout << "tpfjet_had_n_ != tpfjet_had_py_->size() " << tpfjet_had_n_ << " " << tpfjet_had_py_->size() << endl;
    if(tpfjet_had_n_ != tpfjet_had_pz_->size()) cout << "tpfjet_had_n_ != tpfjet_had_pz_->size() " << tpfjet_had_n_ << " " << tpfjet_had_pz_->size() << endl;
    if(tpfjet_had_n_ != tpfjet_had_EcalE_->size()) cout << "tpfjet_had_n_ != tpfjet_had_EcalE_->size() " << tpfjet_had_n_ << " " << tpfjet_had_EcalE_->size() << endl;
    if(tpfjet_had_n_ != tpfjet_had_id_->size()) cout << "tpfjet_had_n_ != tpfjet_had_id_->size() " << tpfjet_had_n_ << " " << tpfjet_had_id_->size() << endl;
    if(tpfjet_had_n_ != tpfjet_had_E_mctruth_->size()) cout << "tpfjet_had_n_ != tpfjet_had_E_mctruth_->size() " << tpfjet_had_n_ << " " << tpfjet_had_E_mctruth_->size() << endl;
    if(tpfjet_had_n_ != tpfjet_had_mcpdgId_->size()) cout << "tpfjet_had_n_ != tpfjet_had_mcpdgId_->size() " << tpfjet_had_n_ << " " << tpfjet_had_mcpdgId_->size() << endl;
    if(tpfjet_ntwrs_ != tpfjet_twr_ieta_->size()) cout << "tpfjet_ntwrs_ != tpfjet_twr_ieta_->size() " << tpfjet_ntwrs_ << " " << tpfjet_twr_ieta_->size() << endl;
    if(tpfjet_ntwrs_ != tpfjet_twr_candtrackind_->size()) cout << "tpfjet_ntwrs_ != tpfjet_twr_candtrackind_->size() " << tpfjet_ntwrs_ << " " << tpfjet_twr_candtrackind_->size() << endl;
    if(tpfjet_ntwrs_ != tpfjet_twr_hadind_->size()) cout << "tpfjet_ntwrs_ != tpfjet_twr_hadind_->size() " << tpfjet_ntwrs_ << " " << tpfjet_twr_hadind_->size() << endl;
    if(tpfjet_ntwrs_ != tpfjet_twr_hade_->size()) cout << "tpfjet_ntwrs_ != tpfjet_twr_hade_->size() " << tpfjet_ntwrs_ << " " << tpfjet_twr_hade_->size() << endl;
    if(tpfjet_ntwrs_ != tpfjet_twr_frac_->size()) cout << "tpfjet_ntwrs_ != tpfjet_twr_frac_->size() " << tpfjet_ntwrs_ << " " << tpfjet_twr_frac_->size() << endl;
    if(tpfjet_ncandtracks_ != tpfjet_candtrack_px_->size()) cout << "tpfjet_ncandtracks_ != tpfjet_candtrack_px_->size() " << tpfjet_ncandtracks_ << " " << tpfjet_candtrack_px_->size() << endl;
    if(tpfjet_ncandtracks_ != tpfjet_candtrack_py_->size()) cout << "tpfjet_ncandtracks_ != tpfjet_candtrack_py_->size() " << tpfjet_ncandtracks_ << " " << tpfjet_candtrack_py_->size() << endl;
    if(tpfjet_ncandtracks_ != tpfjet_candtrack_pz_->size()) cout << "tpfjet_ncandtracks_ != tpfjet_candtrack_pz_->size() " << tpfjet_ncandtracks_ << " " << tpfjet_candtrack_pz_->size() << endl;

    if(ppfjet_had_n_ != ppfjet_had_E_->size()) cout << "ppfjet_had_n_ != ppfjet_had_E_->size() " << ppfjet_had_n_ << " " << ppfjet_had_E_->size() << endl;
    if(ppfjet_had_n_ != ppfjet_had_px_->size()) cout << "ppfjet_had_n_ != ppfjet_had_px_->size() " << ppfjet_had_n_ << " " << ppfjet_had_px_->size() << endl;
    if(ppfjet_had_n_ != ppfjet_had_py_->size()) cout << "ppfjet_had_n_ != ppfjet_had_py_->size() " << ppfjet_had_n_ << " " << ppfjet_had_py_->size() << endl;
    if(ppfjet_had_n_ != ppfjet_had_pz_->size()) cout << "ppfjet_had_n_ != ppfjet_had_pz_->size() " << ppfjet_had_n_ << " " << ppfjet_had_pz_->size() << endl;
    if(ppfjet_had_n_ != ppfjet_had_EcalE_->size()) cout << "ppfjet_had_n_ != ppfjet_had_EcalE_->size() " << ppfjet_had_n_ << " " << ppfjet_had_EcalE_->size() << endl;
    if(ppfjet_had_n_ != ppfjet_had_id_->size()) cout << "ppfjet_had_n_ != ppfjet_had_id_->size() " << ppfjet_had_n_ << " " << ppfjet_had_id_->size() << endl;
    if(ppfjet_had_n_ != ppfjet_had_E_mctruth_->size()) cout << "ppfjet_had_n_ != ppfjet_had_E_mctruth_->size() " << ppfjet_had_n_ << " " << ppfjet_had_E_mctruth_->size() << endl;
    if(ppfjet_had_n_ != ppfjet_had_mcpdgId_->size()) cout << "ppfjet_had_n_ != ppfjet_had_mcpdgId_->size() " << ppfjet_had_n_ << " " << ppfjet_had_mcpdgId_->size() << endl;
    if(ppfjet_ntwrs_ != ppfjet_twr_ieta_->size()) cout << "ppfjet_ntwrs_ != ppfjet_twr_ieta_->size() " << ppfjet_ntwrs_ << " " << ppfjet_twr_ieta_->size() << endl;
    if(ppfjet_ntwrs_ != ppfjet_twr_candtrackind_->size()) cout << "ppfjet_ntwrs_ != ppfjet_twr_candtrackind_->size() " << ppfjet_ntwrs_ << " " << ppfjet_twr_candtrackind_->size() << endl;
    if(ppfjet_ntwrs_ != ppfjet_twr_hadind_->size()) cout << "ppfjet_ntwrs_ != ppfjet_twr_hadind_->size() " << ppfjet_ntwrs_ << " " << ppfjet_twr_hadind_->size() << endl;
    if(ppfjet_ntwrs_ != ppfjet_twr_hade_->size()) cout << "ppfjet_ntwrs_ != ppfjet_twr_hade_->size() " << ppfjet_ntwrs_ << " " << ppfjet_twr_hade_->size() << endl;
    if(ppfjet_ntwrs_ != ppfjet_twr_frac_->size()) cout << "ppfjet_ntwrs_ != ppfjet_twr_frac_->size() " << ppfjet_ntwrs_ << " " << ppfjet_twr_frac_->size() << endl;
    if(ppfjet_ncandtracks_ != ppfjet_candtrack_px_->size()) cout << "ppfjet_ncandtracks_ != ppfjet_candtrack_px_->size() " << ppfjet_ncandtracks_ << " " << ppfjet_candtrack_px_->size() << endl;
    if(ppfjet_ncandtracks_ != ppfjet_candtrack_py_->size()) cout << "ppfjet_ncandtracks_ != ppfjet_candtrack_py_->size() " << ppfjet_ncandtracks_ << " " << ppfjet_candtrack_py_->size() << endl;
    if(ppfjet_ncandtracks_ != ppfjet_candtrack_pz_->size()) cout << "ppfjet_ncandtracks_ != ppfjet_candtrack_pz_->size() " << ppfjet_ncandtracks_ << " " << ppfjet_candtrack_pz_->size() << endl;

    int tpfjet_negE = 0;
    int ppfjet_negE = 0;
    
    //////////////////////////
    // Fill tag histograms
    //////////////////////////
    
    h_tpfjet_eta_->Fill(tpfjet_eta_);
    h_tpfjet_phi_->Fill(tpfjet_phi_);
    h_tpfjet_candn_->Fill(tpfjet_unkown_n_ + tpfjet_electron_n_ + tpfjet_muon_n_ + tpfjet_photon_n_ + tpfjet_had_n_ );
    h_tpfjet_E_->Fill(tpfjet_E_);
    
    float tpfjet_had_E = 0;
    float tpfjet_egammaHF_E = 0;
    float tpfjet_had_EcalE = 0;
    for(int i=0; i<tpfjet_had_n_; i++){
      if(tpfjet_had_id_->at(i) != 3){
	tpfjet_had_E += tpfjet_had_E_->at(i);
      }
      else{
	tpfjet_egammaHF_E += tpfjet_had_E_->at(i);
      }
      tpfjet_had_EcalE += tpfjet_had_EcalE_->at(i);
    }
    float tpfjet_candE = tpfjet_unkown_E_ + tpfjet_electron_E_ + tpfjet_muon_E_ + tpfjet_photon_E_ + tpfjet_had_E + tpfjet_egammaHF_E;
    float tpfjet_emf = (tpfjet_unkown_E_ + tpfjet_electron_E_ + tpfjet_muon_E_ + tpfjet_photon_E_ + tpfjet_egammaHF_E)/tpfjet_candE;
    h_tpfjet_candE_->Fill(tpfjet_candE);
    h_tpfjet_candE_diff_->Fill((tpfjet_candE - tpfjet_E_)/tpfjet_E_);
    h_tpfjet_emf_->Fill(tpfjet_emf);
    
    float tpfjet_rechitonlyE = 0;
    for(int i=0; i<tpfjet_ntwrs_; i++){
      if(tpfjet_twr_hade_->at(i) > 0.0){
	tpfjet_rechitonlyE += tpfjet_twr_hade_->at(i)*tpfjet_twr_frac_->at(i);
      }
      else{
	tpfjet_negE++;
      }
    }
    float tpfjet_rechitE = tpfjet_unkown_E_ + tpfjet_electron_E_ + tpfjet_muon_E_ + tpfjet_photon_E_ + tpfjet_rechitonlyE + tpfjet_had_EcalE;
    h_tpfjet_rechitE_->Fill(tpfjet_rechitE);
    h_tpfjet_rechitE_diff_->Fill((tpfjet_rechitE - tpfjet_E_)/tpfjet_E_);
    if(tpfjet_rechitonlyE == 0) tpfjet_norechits++;
    h_tpfjet_negE_->Fill(tpfjet_negE);
    h_tpfjet_negE_frac_->Fill((float)tpfjet_negE/(float)tpfjet_ntwrs_);

    h_tpfjet_rechit_had_total_E_diff_->Fill((tpfjet_rechitonlyE - tpfjet_had_E)/tpfjet_had_E);
    
    float had_rechit_totalind_E = 0;
    for(int i=0; i<tpfjet_had_n_; i++){
      int HFmatches = 0;
      float had_rechitE = 0;
      float had_rechitE_fraccut = 0;
      float had_rechitE_Ecut = 0;
      for(int j=0; j<tpfjet_ntwrs_; j++){
	if(tpfjet_twr_hadind_->at(j) == i && tpfjet_twr_hade_->at(j) > 0.0){
	  had_rechitE += tpfjet_twr_hade_->at(j)*tpfjet_twr_frac_->at(j);
	  had_rechit_totalind_E += tpfjet_twr_hade_->at(j)*tpfjet_twr_frac_->at(j);
	  if(tpfjet_twr_frac_->at(j) > 0.25) had_rechitE_fraccut += tpfjet_twr_hade_->at(j)*tpfjet_twr_frac_->at(j);
	  if(tpfjet_twr_hade_->at(j) > 0.25) had_rechitE_Ecut += tpfjet_twr_hade_->at(j)*tpfjet_twr_frac_->at(j);
	  if(tpfjet_had_id_->at(i) == 2 || tpfjet_had_id_->at(i) == 3) HFmatches++;
	}
      }
      h_tpfjet_rechit_had_E_diff_->Fill((had_rechitE - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
      if(tpfjet_had_id_->at(i) == 2 || tpfjet_had_id_->at(i) == 3) h_tpfjet_HFmatches_->Fill(HFmatches);
      
      if(tpfjet_had_id_->at(i) == 0){
	if(tpfjet_had_E_->at(i) > 1.0) h_tpfjet_rechit_had_h_E_diff_->Fill((had_rechitE - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
	h_tpfjet_rechit_had_h_EdiffvsE_->Fill((had_rechitE - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i),tpfjet_had_E_->at(i));
	//if((had_rechitE - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i) > 6.0) cout << pf_Event_ << endl;
	h_tpfjet_rechit_had_h_Emc_diff_->Fill((had_rechitE - tpfjet_had_E_mctruth_->at(i))/tpfjet_had_E_mctruth_->at(i));
	
	if(abs(tpfjet_had_mcpdgId_->at(i)) == 211)  h_tpfjet_rechit_had_h_mcPi_Emc_diff_->Fill((had_rechitE - tpfjet_had_E_mctruth_->at(i))/tpfjet_had_E_mctruth_->at(i));
	else h_tpfjet_rechit_had_h_mcNotPi_Emc_diff_->Fill((had_rechitE - tpfjet_had_E_mctruth_->at(i))/tpfjet_had_E_mctruth_->at(i));
	
	if((had_rechitE - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i) < -0.9) h_tpfjet_h_lowE_->Fill(had_rechitE);
	if(had_rechitE == 0){
	  h_tpfjet_h_norechits_EcalEfrac_->Fill(tpfjet_had_EcalE_->at(i)/tpfjet_had_E_->at(i));
	  h_tpfjet_h_norechits_candE_->Fill(tpfjet_had_E_->at(i));
	}
	tpfjet_pdgIdMap[tpfjet_had_mcpdgId_->at(i)]++;

	h_tpfjet_rechit_had_h_fraccut_E_diff_->Fill((had_rechitE_fraccut - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
	h_tpfjet_rechit_had_h_Ecut_E_diff_->Fill((had_rechitE_Ecut - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
      }
      else if(tpfjet_had_id_->at(i) == 1) h_tpfjet_rechit_had_h0_E_diff_->Fill((had_rechitE - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
      else if(tpfjet_had_id_->at(i) == 2) h_tpfjet_rechit_had_HFhad_E_diff_->Fill((had_rechitE - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
      else if(tpfjet_had_id_->at(i) == 3) h_tpfjet_rechit_had_egammaHF_E_diff_->Fill((had_rechitE - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
    }

    h_tpfjet_rechit_had_totalind_E_diff_->Fill((had_rechit_totalind_E - tpfjet_had_E)/tpfjet_had_E);
    
    if(tpfjet_twr_candtrackind_->size() == tpfjet_ntwrs_){
      for(int i=0; i<tpfjet_ncandtracks_; i++){
	float trackp = sqrt(tpfjet_candtrack_px_->at(i)*tpfjet_candtrack_px_->at(i) + tpfjet_candtrack_py_->at(i)*tpfjet_candtrack_py_->at(i) + tpfjet_candtrack_pz_->at(i)*tpfjet_candtrack_pz_->at(i));
	float trackpt = sqrt(tpfjet_candtrack_px_->at(i)*tpfjet_candtrack_px_->at(i) + tpfjet_candtrack_py_->at(i)*tpfjet_candtrack_py_->at(i));
	float had_rechitE = 0;
	for(int j=0; j<tpfjet_ntwrs_; j++){
	  //cout << i << "/" << tpfjet_ncandtracks_ << " " << j << "/" << tpfjet_ntwrs_ << endl;
	  if(tpfjet_twr_candtrackind_->at(j) == i && tpfjet_twr_hade_->at(j) > 0.0){
	    had_rechitE += tpfjet_twr_hade_->at(j)*tpfjet_twr_frac_->at(j);
	  }
	}
	h_tpfjet_h_EoverP_->Fill(had_rechitE/trackp);
	if(had_rechitE == 0) h_tpfjet_h_norechits_pt_->Fill(trackpt);
	else h_tpfjet_h_rechits_pt_->Fill(trackpt);
      }
    }

    //////////////////////////
    // Fill probe histograms
    //////////////////////////
    
    h_ppfjet_eta_->Fill(ppfjet_eta_);
    h_ppfjet_phi_->Fill(ppfjet_phi_);
    h_ppfjet_candn_->Fill(ppfjet_unkown_n_ + ppfjet_electron_n_ + ppfjet_muon_n_ + ppfjet_photon_n_ + ppfjet_had_n_ );
    h_ppfjet_E_->Fill(ppfjet_E_);
    
    float ppfjet_had_E = 0;
    float ppfjet_egammaHF_E = 0;
    float ppfjet_had_EcalE = 0;
    for(int i=0; i<ppfjet_had_n_; i++){
      if(ppfjet_had_id_->at(i) != 3){
	ppfjet_had_E += ppfjet_had_E_->at(i);
      }
      else{
	ppfjet_egammaHF_E += ppfjet_had_E_->at(i);
      }
      ppfjet_had_EcalE += ppfjet_had_EcalE_->at(i);
    }
    float ppfjet_candE = ppfjet_unkown_E_ + ppfjet_electron_E_ + ppfjet_muon_E_ + ppfjet_photon_E_ + ppfjet_had_E + ppfjet_egammaHF_E;
    float ppfjet_emf = (ppfjet_unkown_E_ + ppfjet_electron_E_ + ppfjet_muon_E_ + ppfjet_photon_E_ + ppfjet_egammaHF_E)/ppfjet_candE;
    h_ppfjet_candE_->Fill(ppfjet_candE);
    h_ppfjet_candE_diff_->Fill((ppfjet_candE - ppfjet_E_)/ppfjet_E_);
    h_ppfjet_emf_->Fill(ppfjet_emf);

    float ppfjet_rechitonlyE = 0;
    for(int i=0; i<ppfjet_ntwrs_; i++){
      if(ppfjet_twr_hade_->at(i) > 0.0){
	ppfjet_rechitonlyE += ppfjet_twr_hade_->at(i)*ppfjet_twr_frac_->at(i);
      }
      else{
	ppfjet_negE++;
      }
    }
    float ppfjet_rechitE = ppfjet_unkown_E_ + ppfjet_electron_E_ + ppfjet_muon_E_ + ppfjet_photon_E_ + ppfjet_rechitonlyE + ppfjet_had_EcalE;
    h_ppfjet_rechitE_->Fill(ppfjet_rechitE);
    h_ppfjet_rechitE_diff_->Fill((ppfjet_rechitE - ppfjet_E_)/ppfjet_E_);
    if(ppfjet_rechitonlyE == 0) ppfjet_norechits++;
    h_ppfjet_negE_->Fill(ppfjet_negE);
    h_ppfjet_negE_frac_->Fill((float)ppfjet_negE/(float)ppfjet_ntwrs_);

    for(int i=0; i<ppfjet_had_n_; i++){
      int HFmatches = 0;
      float had_rechitE = 0;
      float had_rechitE_fraccut = 0;
      float had_rechitE_Ecut = 0;
      for(int j=0; j<ppfjet_ntwrs_; j++){
	if(ppfjet_twr_hadind_->at(j) == i && ppfjet_twr_hade_->at(j) > 0.0){
	  had_rechitE += ppfjet_twr_hade_->at(j)*ppfjet_twr_frac_->at(j);
	  if(ppfjet_twr_frac_->at(j) > 0.25) had_rechitE_fraccut += ppfjet_twr_hade_->at(j)*ppfjet_twr_frac_->at(j);
	  if(ppfjet_twr_hade_->at(j) > 0.25) had_rechitE_Ecut += ppfjet_twr_hade_->at(j)*ppfjet_twr_frac_->at(j);
	  if(ppfjet_had_id_->at(i) == 2 || ppfjet_had_id_->at(i) == 3) HFmatches++;
	}
      }
      h_ppfjet_rechit_had_E_diff_->Fill((had_rechitE - ppfjet_had_E_->at(i))/ppfjet_had_E_->at(i));
      if(ppfjet_had_id_->at(i) == 2 || ppfjet_had_id_->at(i) == 3) h_ppfjet_HFmatches_->Fill(HFmatches);

      if(ppfjet_had_id_->at(i) == 0){
	if(ppfjet_had_E_->at(i) > 1) h_ppfjet_rechit_had_h_E_diff_->Fill((had_rechitE - ppfjet_had_E_->at(i))/ppfjet_had_E_->at(i));
	h_ppfjet_rechit_had_h_EdiffvsE_->Fill((had_rechitE - ppfjet_had_E_->at(i))/ppfjet_had_E_->at(i),ppfjet_had_E_->at(i));
	h_ppfjet_rechit_had_h_Emc_diff_->Fill((had_rechitE - ppfjet_had_E_mctruth_->at(i))/ppfjet_had_E_mctruth_->at(i));
	
	if(abs(ppfjet_had_mcpdgId_->at(i)) == 211) h_ppfjet_rechit_had_h_mcPi_Emc_diff_->Fill((had_rechitE - ppfjet_had_E_mctruth_->at(i))/ppfjet_had_E_mctruth_->at(i));
	else h_ppfjet_rechit_had_h_mcNotPi_Emc_diff_->Fill((had_rechitE - ppfjet_had_E_mctruth_->at(i))/ppfjet_had_E_mctruth_->at(i));

	if((had_rechitE - ppfjet_had_E_->at(i))/ppfjet_had_E_->at(i) < -0.9) h_ppfjet_h_lowE_->Fill(had_rechitE);
	if(had_rechitE == 0){
	  h_ppfjet_h_norechits_EcalEfrac_->Fill(ppfjet_had_EcalE_->at(i)/ppfjet_had_E_->at(i));
	  h_ppfjet_h_norechits_candE_->Fill(ppfjet_had_E_->at(i));
	}
	ppfjet_pdgIdMap[ppfjet_had_mcpdgId_->at(i)]++;

	h_ppfjet_rechit_had_h_fraccut_E_diff_->Fill((had_rechitE_fraccut - ppfjet_had_E_->at(i))/ppfjet_had_E_->at(i));
	h_ppfjet_rechit_had_h_Ecut_E_diff_->Fill((had_rechitE_Ecut - ppfjet_had_E_->at(i))/ppfjet_had_E_->at(i));
      }
      else if(ppfjet_had_id_->at(i) == 1) h_ppfjet_rechit_had_h0_E_diff_->Fill((had_rechitE - ppfjet_had_E_->at(i))/ppfjet_had_E_->at(i));
      else if(ppfjet_had_id_->at(i) == 2) h_ppfjet_rechit_had_HFhad_E_diff_->Fill((had_rechitE - ppfjet_had_E_->at(i))/ppfjet_had_E_->at(i));
      else if(ppfjet_had_id_->at(i) == 3) h_ppfjet_rechit_had_egammaHF_E_diff_->Fill((had_rechitE - ppfjet_had_E_->at(i))/ppfjet_had_E_->at(i));
    }

    if(ppfjet_twr_candtrackind_->size() == ppfjet_ntwrs_){
      for(int i=0; i<ppfjet_ncandtracks_; i++){
	float trackp = sqrt(ppfjet_candtrack_px_->at(i)*ppfjet_candtrack_px_->at(i) + ppfjet_candtrack_py_->at(i)*ppfjet_candtrack_py_->at(i) + ppfjet_candtrack_pz_->at(i)*ppfjet_candtrack_pz_->at(i));
	float trackp = sqrt(ppfjet_candtrack_px_->at(i)*ppfjet_candtrack_px_->at(i) + ppfjet_candtrack_py_->at(i)*ppfjet_candtrack_py_->at(i));
	float had_rechitE = 0;
	for(int j=0; j<ppfjet_ntwrs_; j++){
	  if(ppfjet_twr_candtrackind_->at(j) == i && ppfjet_twr_hade_->at(j) > 0.0){
	    had_rechitE += ppfjet_twr_hade_->at(j)*ppfjet_twr_frac_->at(j);
	  }
	}
	h_ppfjet_h_EoverP_->Fill(had_rechitE/trackp);
	if(had_rechitE == 0) h_ppfjet_h_norechits_pt_->Fill(trackpt);
	else h_ppfjet_h_rechits_pt_->Fill(trackpt);
      }
    }

    //////////////////////////
    // Fill dijet histograms
    //////////////////////////
    
    h_pf_dijet_deta_->Fill(pf_dijet_deta_);
    
    if(tpfjet_rechitonlyE == 0 || ppfjet_rechitonlyE == 0) dijet_norechits++;
  }

  TFile* fout = new TFile(output,"RECREATE");
  fout->cd();
  h_tpfjet_eta_->Write();
  h_tpfjet_phi_->Write();
  h_tpfjet_candn_->Write();
  h_tpfjet_E_->Write();
  h_tpfjet_candE_->Write();
  h_tpfjet_candE_diff_->Write();
  h_tpfjet_rechitE_->Write();
  h_tpfjet_rechitE_diff_->Write();
  h_tpfjet_emf_->Write();
  h_tpfjet_negE_->Write();
  h_tpfjet_negE_frac_->Write();
  h_tpfjet_rechit_had_E_diff_->Write();
  h_tpfjet_rechit_had_total_E_diff_->Write();
  h_tpfjet_rechit_had_totalind_E_diff_->Write();
  h_tpfjet_rechit_had_h_E_diff_->Write();
  h_tpfjet_rechit_had_h_fraccut_E_diff_->Write();
  h_tpfjet_rechit_had_h_Ecut_E_diff_->Write();
  h_tpfjet_rechit_had_h_Emc_diff_->Write();
  h_tpfjet_rechit_had_h_mcPi_Emc_diff_->Write();
  h_tpfjet_rechit_had_h_mcNotPi_Emc_diff_->Write();
  h_tpfjet_rechit_had_h_EdiffvsE_->Write();
  h_tpfjet_rechit_had_h0_E_diff_->Write();
  h_tpfjet_rechit_had_HFhad_E_diff_->Write();
  h_tpfjet_rechit_had_egammaHF_E_diff_->Write();
  h_tpfjet_HFmatches_->Write();
  h_tpfjet_h_EoverP_->Write();
  h_tpfjet_h_lowE_->Write();
  h_tpfjet_h_norechits_EcalEfrac_->Write();
  h_tpfjet_h_norechits_EcalEfrac_lowpt_->Write();
  h_tpfjet_h_norechits_EcalEfrac_highpt_->Write();
  h_tpfjet_h_norechits_pt_->Write();
  h_tpfjet_h_rechits_pt_->Write();
  h_tpfjet_h_norechits_candE_->Write();

  h_ppfjet_eta_->Write();
  h_ppfjet_phi_->Write();
  h_ppfjet_candn_->Write();
  h_ppfjet_E_->Write();
  h_ppfjet_candE_->Write();
  h_ppfjet_candE_diff_->Write();
  h_ppfjet_rechitE_->Write();
  h_ppfjet_rechitE_diff_->Write();
  h_ppfjet_emf_->Write();
  h_ppfjet_negE_->Write();
  h_ppfjet_negE_frac_->Write();
  h_ppfjet_rechit_had_E_diff_->Write();
  h_ppfjet_rechit_had_h_E_diff_->Write();
  h_ppfjet_rechit_had_h_fraccut_E_diff_->Write();
  h_ppfjet_rechit_had_h_Ecut_E_diff_->Write();
  h_ppfjet_rechit_had_h_Emc_diff_->Write();
  h_ppfjet_rechit_had_h_mcPi_Emc_diff_->Write();
  h_ppfjet_rechit_had_h_mcNotPi_Emc_diff_->Write();
  h_ppfjet_rechit_had_h_EdiffvsE_->Write();
  h_ppfjet_rechit_had_h0_E_diff_->Write();
  h_ppfjet_rechit_had_HFhad_E_diff_->Write();
  h_ppfjet_rechit_had_egammaHF_E_diff_->Write();
  h_ppfjet_HFmatches_->Write();
  h_ppfjet_h_EoverP_->Write();
  h_ppfjet_h_lowE_->Write();
  h_ppfjet_h_norechits_EcalEfrac_->Write();
  h_ppfjet_h_norechits_EcalEfrac_lowpt_->Write();
  h_ppfjet_h_norechits_EcalEfrac_highpt_->Write();
  h_ppfjet_h_norechits_pt_->Write();
  h_ppfjet_h_rechits_pt_->Write();
  h_ppfjet_h_norechits_candE_->Write();
  
  h_pf_dijet_deta_->Write();

  fout->Close();

  cout << "tag jet MC PDG IDs" << endl;
  for(map<Int_t,Int_t>::const_iterator it=tpfjet_pdgIdMap.begin(); it!=tpfjet_pdgIdMap.end(); it++){
    if(it->first > 0){
      cout << it->first << " " << it->second + tpfjet_pdgIdMap[-it->first] << endl;
    }
  }
  cout << "probe jet MC PDG IDs" << endl;
  for(map<Int_t,Int_t>::const_iterator it=ppfjet_pdgIdMap.begin(); it!=ppfjet_pdgIdMap.end(); it++){
    if(it->first > 0){
      cout << it->first << " " << it->second + ppfjet_pdgIdMap[-it->first] << endl;
    }
  }

  cout << "Events in tree: " << nEvents << endl;
  cout << "No rechits: " << dijet_norechits << "=" << (float)dijet_norechits/(float)nEvents << " tag:" << tpfjet_norechits << "=" << (float)tpfjet_norechits/(float)nEvents << " probe: " << ppfjet_norechits << "=" << (float)ppfjet_norechits/(float)nEvents << endl;
  
  return;
}
