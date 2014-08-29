#include <vector>
#include <map>
/*#include <iostream>
#include "TROOT.h"
#include "TChain.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace std;


double deltaR(const double px1, const double py1, const double pz1, const double px2, const double py2, const double pz2){
  TLorentzVector p1(px1,py1,pz1,0.0);
  TLorentzVector p2(px2,py2,pz2,0.0);
  double deta = p1.Eta() - p2.Eta();
  double dphi = fabs(p1.Phi() - p2.Phi());
  if(dphi>3.1415927) dphi = 2*3.1415927 - dphi;
  return sqrt(deta*deta + dphi*dphi);
  }*/

void testRespCorrDiJetsTree()
{
  //gROOT->ProcessLine(".L loader.C+");
  gROOT->ProcessLine(".L deltaR.C+");

  TChain* tree = new TChain("pf_dijettree");
  TString input = "/uscms_data/d3/dgsheffi/HCal/Trees/QCD_Pt-15to3000_0030487D5E5F.root";
  //TString input = "/uscms_data/d3/dgsheffi/HCal/Trees/Pion_Pt-50.root";
  tree->Add(input);

  TString output = "/uscms_data/d3/dgsheffi/HCal/Trees/validation/QCD_Pt-15to3000_0030487D5E5F.root";
  //TString output = "/uscms_data/d3/dgsheffi/HCal/Trees/validation/Pion_Pt-50.root";

  float tpfjet_pt_, tpfjet_p_, tpfjet_E_, tpfjet_eta_, tpfjet_phi_, tpfjet_scale_;
  float tpfjet_gendr_, tpfjet_genpt_, tpfjet_genp_, tpfjet_genE_;
  float tpfjet_EBE_, tpfjet_EEE_, tpfjet_HBE_, tpfjet_HEE_, tpfjet_HFE_;
  float tpfjet_unkown_E_, tpfjet_unkown_px_, tpfjet_unkown_py_, tpfjet_unkown_pz_, tpfjet_unkown_EcalE_;
  float tpfjet_electron_E_, tpfjet_electron_px_, tpfjet_electron_py_, tpfjet_electron_pz_, tpfjet_electron_EcalE_;
  float tpfjet_muon_E_, tpfjet_muon_px_, tpfjet_muon_py_, tpfjet_muon_pz_, tpfjet_muon_EcalE_;
  float tpfjet_photon_E_, tpfjet_photon_px_, tpfjet_photon_py_, tpfjet_photon_pz_, tpfjet_photon_EcalE_;
  int tpfjet_unkown_n_, tpfjet_electron_n_, tpfjet_muon_n_, tpfjet_photon_n_;
  int tpfjet_had_n_;
  vector<float>* tpfjet_had_E_;
  vector<float>* tpfjet_had_px_;
  vector<float>* tpfjet_had_py_;
  vector<float>* tpfjet_had_pz_;
  vector<float>* tpfjet_had_EcalE_;
  vector<float>* tpfjet_had_rawHcalE_;
  vector<float>* tpfjet_had_emf_;
  vector<float>* tpfjet_had_E_mctruth_;
  vector<int>* tpfjet_had_id_;
  vector<int>* tpfjet_had_candtrackind_;
  vector<int>* tpfjet_had_mcpdgId_;
  int tpfjet_ntwrs_;
  vector<int>* tpfjet_twr_ieta_;
  vector<int>* tpfjet_twr_candtrackind_;
  vector<int>* tpfjet_twr_hadind_;
  vector<int>* tpfjet_twr_elmttype_;
  vector<float>* tpfjet_twr_hade_;
  vector<float>* tpfjet_twr_frac_;
  int tpfjet_ncandtracks_;
  vector<float>* tpfjet_candtrack_px_;
  vector<float>* tpfjet_candtrack_py_;
  vector<float>* tpfjet_candtrack_pz_;
  vector<float>* tpfjet_candtrack_EcalE_;
  float ppfjet_pt_, ppfjet_p_, ppfjet_E_, ppfjet_eta_, ppfjet_phi_, ppfjet_scale_;
  float ppfjet_gendr_, ppfjet_genpt_, ppfjet_genp_, ppfjet_genE_;
  float ppfjet_EBE_, ppfjet_EEE_, ppfjet_HBE_, ppfjet_HEE_, ppfjet_HFE_;
  float ppfjet_unkown_E_, ppfjet_unkown_px_, ppfjet_unkown_py_, ppfjet_unkown_pz_, ppfjet_unkown_EcalE_;
  float ppfjet_electron_E_, ppfjet_electron_px_, ppfjet_electron_py_, ppfjet_electron_pz_, ppfjet_electron_EcalE_;
  float ppfjet_muon_E_, ppfjet_muon_px_, ppfjet_muon_py_, ppfjet_muon_pz_, ppfjet_muon_EcalE_;
  float ppfjet_photon_E_, ppfjet_photon_px_, ppfjet_photon_py_, ppfjet_photon_pz_, ppfjet_photon_EcalE_;
  int ppfjet_unkown_n_, ppfjet_electron_n_, ppfjet_muon_n_, ppfjet_photon_n_;
  int ppfjet_had_n_;
  vector<float>* ppfjet_had_E_;
  vector<float>* ppfjet_had_px_;
  vector<float>* ppfjet_had_py_;
  vector<float>* ppfjet_had_pz_;
  vector<float>* ppfjet_had_EcalE_;
  vector<float>* ppfjet_had_rawHcalE_;
  vector<float>* ppfjet_had_emf_;
  vector<float>* ppfjet_had_E_mctruth_;
  vector<int>* ppfjet_had_id_;
  vector<int>* ppfjet_had_candtrackind_;
  vector<int>* ppfjet_had_mcpdgId_;
  int ppfjet_ntwrs_;
  vector<int>* ppfjet_twr_ieta_;
  vector<float>* ppfjet_twr_candtrackind_;
  vector<float>* ppfjet_twr_hadind_;
  vector<float>* ppfjet_twr_elmttype_;
  vector<float>* ppfjet_twr_hade_;
  vector<float>* ppfjet_twr_frac_;
  int ppfjet_ncandtracks_;
  vector<float>* ppfjet_candtrack_px_;
  vector<float>* ppfjet_candtrack_py_;
  vector<float>* ppfjet_candtrack_pz_;
  vector<float>* ppfjet_candtrack_EcalE_;
  float pf_dijet_deta_, pf_dijet_dphi_, pf_dijet_balance_;
  float pf_thirdjet_px_, pf_thirdjet_py_;
  int pf_Run_, pf_Lumi_, pf_Event_;

  tree->SetBranchAddress("tpfjet_E",&tpfjet_E_);
  tree->SetBranchAddress("tpfjet_pt",&tpfjet_pt_);
  tree->SetBranchAddress("tpfjet_p",&tpfjet_p_);
  tree->SetBranchAddress("tpfjet_E",&tpfjet_E_);
  tree->SetBranchAddress("tpfjet_eta",&tpfjet_eta_);
  tree->SetBranchAddress("tpfjet_phi",&tpfjet_phi_);
  tree->SetBranchAddress("tpfjet_scale",&tpfjet_scale_);
  tree->SetBranchAddress("tpfjet_genpt",&tpfjet_genpt_);
  tree->SetBranchAddress("tpfjet_genp",&tpfjet_genp_);
  tree->SetBranchAddress("tpfjet_genE",&tpfjet_genE_);
  tree->SetBranchAddress("tpfjet_gendr",&tpfjet_gendr_);
  tree->SetBranchAddress("tpfjet_unkown_E",&tpfjet_unkown_E_);
  tree->SetBranchAddress("tpfjet_electron_E",&tpfjet_electron_E_);
  tree->SetBranchAddress("tpfjet_muon_E",&tpfjet_muon_E_);
  tree->SetBranchAddress("tpfjet_photon_E",&tpfjet_photon_E_);
  tree->SetBranchAddress("tpfjet_unkown_px",&tpfjet_unkown_px_);
  tree->SetBranchAddress("tpfjet_electron_px",&tpfjet_electron_px_);
  tree->SetBranchAddress("tpfjet_muon_px",&tpfjet_muon_px_);
  tree->SetBranchAddress("tpfjet_photon_px",&tpfjet_photon_px_);
  tree->SetBranchAddress("tpfjet_unkown_py",&tpfjet_unkown_py_);
  tree->SetBranchAddress("tpfjet_electron_py",&tpfjet_electron_py_);
  tree->SetBranchAddress("tpfjet_muon_py",&tpfjet_muon_py_);
  tree->SetBranchAddress("tpfjet_photon_py",&tpfjet_photon_py_);
  tree->SetBranchAddress("tpfjet_unkown_pz",&tpfjet_unkown_pz_);
  tree->SetBranchAddress("tpfjet_electron_pz",&tpfjet_electron_pz_);
  tree->SetBranchAddress("tpfjet_muon_pz",&tpfjet_muon_pz_);
  tree->SetBranchAddress("tpfjet_photon_pz",&tpfjet_photon_pz_);
  tree->SetBranchAddress("tpfjet_unkown_EcalE",&tpfjet_unkown_EcalE_);
  tree->SetBranchAddress("tpfjet_electron_EcalE",&tpfjet_electron_EcalE_);
  tree->SetBranchAddress("tpfjet_muon_EcalE",&tpfjet_muon_EcalE_);
  tree->SetBranchAddress("tpfjet_photon_EcalE",&tpfjet_photon_EcalE_);
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
  tree->SetBranchAddress("tpfjet_had_rawHcalE",&tpfjet_had_rawHcalE_);
  tree->SetBranchAddress("tpfjet_had_emf",&tpfjet_had_emf_);
  tree->SetBranchAddress("tpfjet_had_id",&tpfjet_had_id_);
  tree->SetBranchAddress("tpfjet_had_candtrackind",&tpfjet_had_candtrackind_);
  tree->SetBranchAddress("tpfjet_had_E_mctruth",&tpfjet_had_E_mctruth_);
  tree->SetBranchAddress("tpfjet_had_mcpdgId",&tpfjet_had_mcpdgId_);
  tree->SetBranchAddress("tpfjet_ntwrs",&tpfjet_ntwrs_);
  tree->SetBranchAddress("tpfjet_twr_ieta",&tpfjet_twr_ieta_);
  tree->SetBranchAddress("tpfjet_twr_hade",&tpfjet_twr_hade_);
  tree->SetBranchAddress("tpfjet_twr_frac",&tpfjet_twr_frac_);
  tree->SetBranchAddress("tpfjet_twr_candtrackind",&tpfjet_twr_candtrackind_);
  tree->SetBranchAddress("tpfjet_twr_hadind",&tpfjet_twr_hadind_);
  tree->SetBranchAddress("tpfjet_twr_elmttype",&tpfjet_twr_elmttype_);
  tree->SetBranchAddress("tpfjet_ncandtracks",&tpfjet_ncandtracks_);
  tree->SetBranchAddress("tpfjet_candtrack_px",&tpfjet_candtrack_px_);
  tree->SetBranchAddress("tpfjet_candtrack_py",&tpfjet_candtrack_py_);
  tree->SetBranchAddress("tpfjet_candtrack_pz",&tpfjet_candtrack_pz_);
  tree->SetBranchAddress("tpfjet_candtrack_EcalE",&tpfjet_candtrack_EcalE_);
  tree->SetBranchAddress("ppfjet_pt",&ppfjet_pt_);
  tree->SetBranchAddress("ppfjet_p",&ppfjet_p_);
  tree->SetBranchAddress("ppfjet_E",&ppfjet_E_);
  tree->SetBranchAddress("ppfjet_eta",&ppfjet_eta_);
  tree->SetBranchAddress("ppfjet_phi",&ppfjet_phi_);
  tree->SetBranchAddress("ppfjet_scale",&ppfjet_scale_);
  tree->SetBranchAddress("ppfjet_genpt",&ppfjet_genpt_);
  tree->SetBranchAddress("ppfjet_genp",&ppfjet_genp_);
  tree->SetBranchAddress("ppfjet_genE",&ppfjet_genE_);
  tree->SetBranchAddress("ppfjet_gendr",&ppfjet_gendr_);
  tree->SetBranchAddress("ppfjet_unkown_E",&ppfjet_unkown_E_);
  tree->SetBranchAddress("ppfjet_electron_E",&ppfjet_electron_E_);
  tree->SetBranchAddress("ppfjet_muon_E",&ppfjet_muon_E_);
  tree->SetBranchAddress("ppfjet_photon_E",&ppfjet_photon_E_);
  tree->SetBranchAddress("ppfjet_unkown_px",&ppfjet_unkown_px_);
  tree->SetBranchAddress("ppfjet_electron_px",&ppfjet_electron_px_);
  tree->SetBranchAddress("ppfjet_muon_px",&ppfjet_muon_px_);
  tree->SetBranchAddress("ppfjet_photon_px",&ppfjet_photon_px_);
  tree->SetBranchAddress("ppfjet_unkown_py",&ppfjet_unkown_py_);
  tree->SetBranchAddress("ppfjet_electron_py",&ppfjet_electron_py_);
  tree->SetBranchAddress("ppfjet_muon_py",&ppfjet_muon_py_);
  tree->SetBranchAddress("ppfjet_photon_py",&ppfjet_photon_py_);
  tree->SetBranchAddress("ppfjet_unkown_pz",&ppfjet_unkown_pz_);
  tree->SetBranchAddress("ppfjet_electron_pz",&ppfjet_electron_pz_);
  tree->SetBranchAddress("ppfjet_muon_pz",&ppfjet_muon_pz_);
  tree->SetBranchAddress("ppfjet_photon_pz",&ppfjet_photon_pz_);
  tree->SetBranchAddress("ppfjet_unkown_EcalE",&ppfjet_unkown_EcalE_);
  tree->SetBranchAddress("ppfjet_electron_EcalE",&ppfjet_electron_EcalE_);
  tree->SetBranchAddress("ppfjet_muon_EcalE",&ppfjet_muon_EcalE_);
  tree->SetBranchAddress("ppfjet_photon_EcalE",&ppfjet_photon_EcalE_);
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
  tree->SetBranchAddress("ppfjet_had_rawHcalE",&ppfjet_had_rawHcalE_);
  tree->SetBranchAddress("ppfjet_had_emf",&ppfjet_had_emf_);
  tree->SetBranchAddress("ppfjet_had_id",&ppfjet_had_id_);
  tree->SetBranchAddress("ppfjet_had_candtrackind",&ppfjet_had_candtrackind_);
  tree->SetBranchAddress("ppfjet_had_E_mctruth",&ppfjet_had_E_mctruth_);
  tree->SetBranchAddress("ppfjet_had_mcpdgId",&ppfjet_had_mcpdgId_);
  tree->SetBranchAddress("ppfjet_ntwrs",&ppfjet_ntwrs_);
  tree->SetBranchAddress("ppfjet_twr_ieta",&ppfjet_twr_ieta_);
  tree->SetBranchAddress("ppfjet_twr_hade",&ppfjet_twr_hade_);
  tree->SetBranchAddress("ppfjet_twr_frac",&ppfjet_twr_frac_);
  tree->SetBranchAddress("ppfjet_twr_candtrackind",&ppfjet_twr_candtrackind_);
  tree->SetBranchAddress("ppfjet_twr_hadind",&ppfjet_twr_hadind_);
  tree->SetBranchAddress("ppfjet_twr_elmttype",&ppfjet_twr_elmttype_);
  tree->SetBranchAddress("ppfjet_ncandtracks",&ppfjet_ncandtracks_);
  tree->SetBranchAddress("ppfjet_candtrack_px",&ppfjet_candtrack_px_);
  tree->SetBranchAddress("ppfjet_candtrack_py",&ppfjet_candtrack_py_);
  tree->SetBranchAddress("ppfjet_candtrack_pz",&ppfjet_candtrack_pz_);
  tree->SetBranchAddress("ppfjet_candtrack_EcalE",&ppfjet_candtrack_EcalE_);
  tree->SetBranchAddress("pf_dijet_deta",&pf_dijet_deta_);
  tree->SetBranchAddress("pf_dijet_dphi",&pf_dijet_dphi_);
  tree->SetBranchAddress("pf_dijet_balance",&pf_dijet_balance_);
  tree->SetBranchAddress("pf_thirdjet_px",&pf_thirdjet_px_);
  tree->SetBranchAddress("pf_thirdjet_py",&pf_thirdjet_py_);
  tree->SetBranchAddress("pf_Run",&pf_Run_);
  tree->SetBranchAddress("pf_Lumi",&pf_Lumi_);
  tree->SetBranchAddress("pf_Event",&pf_Event_);
  
  // Jet
  TH1D* h_tag_jet_Ediff_ = new TH1D("h_tag_jet_Ediff","tag (rechits - pfjet)/pfjet",200,-1,8);
  TH1D* h_tag_jet_genEdiff_ = new TH1D("h_tag_jet_genEdiff","tag (rechits - genjet)/genjet",200,-1,8);
  TH1D* h_tag_jet_Ediff_cut_ = new TH1D("h_tag_jet_Ediff_cut","tag (rechits - pfjet)/pfjet with cuts",200,-1,8);
  TH1D* h_tag_jet_negEraw_ = new TH1D("h_tag_jet_negEraw","tag negative rechit energies",200,-2,0);
  TH1D* h_tag_jet_negEtimesFrac_ = new TH1D("h_tag_jet_negEtimesFrac","tag negative rechit energies times fraction",200,-2,0);
  TH1D* h_tag_jet_negEfrac_ = new TH1D("h_tag_jet_negEfrac","tag fraction with negative energies",200,0,1);
  // Candidates
  TH1D* h_tag_h_Ediff_ = new TH1D("h_tag_h_Ediff","tag h (rechits - cand)/cand",200,-1,8);
  TH1D* h_tag_h0_Ediff_ = new TH1D("h_tag_h0_Ediff","tag h0 (rechits - cand)/cand",200,-1,8);
  TH1D* h_tag_HFhad_Ediff_ = new TH1D("h_tag_HFhad_Ediff","tag HFhad (rechits - cand)/cand",200,-1,8);
  TH1D* h_tag_egammaHF_Ediff_ = new TH1D("h_tag_egammaHF_Ediff","tag egammaHF (rechits - cand)/cand",200,-1,8);
  TH2D* h_tag_h_EvsEdiff_ = new TH2D("h_tag_h_EvsEdiff","tag h cand E vs (rechits - cand)/cand",200,-1,8,200,0,50);
  TH2D* h_tag_h0_EvsEdiff_ = new TH2D("h_tag_h0_EvsEdiff","tag h0 cand E vs (rechits - cand)/cand",200,-1,8,200,0,50);
  TH1D* h_tag_h_Ediff_EcalE_ = new TH1D("h_tag_h_Ediff_EcalE","tag h (rechits + Ecal - cand)/cand",200,-1,8);
  TH1D* h_tag_h0_Ediff_EcalE_ = new TH1D("h_tag_h0_Ediff_EcalE","tag h0 (rechits + Ecal - cand)/cand",200,-1,8);
  TH1D* h_tag_HFhad_Ediff_EcalE_ = new TH1D("h_tag_HFhad_Ediff_EcalE","tag HFhad (rechits +Ecal - cand)/cand",200,-1,8);
  TH1D* h_tag_egammaHF_Ediff_EcalE_ = new TH1D("h_tag_egammaHF_Ediff_EcalE","tag egammaHF (rechits +Ecal - cand)/cand",200,-1,8);
  TH1D* h_tag_rechitspercandidate_ = new TH1D("h_tag_rechitspercandidate","tag rechits per candidate",200,0,200);
  TH1D* h_tag_h_Ediff_rawHcalE_ = new TH1D("h_tag_h_Ediff_rawHcalE","tag h (rechits - rawHcalE)/rawHcalE",200,-1,8);
  TH1D* h_tag_h_Ediff_cut_ = new TH1D("h_tag_h_Ediff_cut","tag h (rechits - cand)/cand with cuts",200,-1,8);
  // Zero h rechits
  TH1D* h_tag_h_pt_rechits_ = new TH1D("h_tag_h_pt_rechits","tag h pt rechits",200,0,100);
  TH1D* h_tag_h_pt_norechits_EcalE_ = new TH1D("h_tag_h_pt_norechits_EcalE","tag h pt norechits EcalE",200,0,100);
  TH1D* h_tag_h_pt_norechits_noEcalE_ = new TH1D("h_tag_h_pt_norechits_noEcalE","tag h pt norechits noEcalE",200,0,100);
  TH1D* h_tag_h_EcalEfrac_norechits_ = new TH1D("h_tag_h_EcalEfrac_norechits","tag h candEcalE/cand E norechits",200,0,1.5);
  TH1D* h_tag_h_EcalEfrac_norechits_lowPt_ = new TH1D("h_tag_h_EcalEfrac_norechits_lowPt","tag h candEcalE/cand E norechits p_{T} < 1 GeV",200,0,1.5);
  TH1D* h_tag_h_EcalEfrac_norechits_midPt_ = new TH1D("h_tag_h_EcalEfrac_norechits_midPt","tag h candEcalE/cand E norechits 1 < p_{T} < 10 GeV",200,0,1.5);
  TH1D* h_tag_h_EcalEfrac_norechits_highPt_ = new TH1D("h_tag_h_EcalEfrac_norechits_highPt","tag h candEcalE/cand E norechits  p_{T} > 10 GeV",200,0,1.5);
  // High h Ediff
  TH1D* h_tag_h_twrE_highEdiff_ = new TH1D("h_tag_h_twrE_highEdiff","tag h tower energies Ediff > 3",200,0,100);
  TH1D* h_tag_h_twrE_lowEdiff_ = new TH1D("h_tag_h_twrE_lowEdiff","tag h tower energies Ediff < 3",200,0,100);
  TH1D* h_tag_h_ntwrs_highEdiff_ = new TH1D("h_tag_h_ntwrs_highEdiff","tag h tower multiplicity Ediff > 3",400,0,400);
  TH1D* h_tag_h_ntwrs_lowEdiff_ = new TH1D("h_tag_h_ntwrs_lowEdiff","tag h tower multiplicity Ediff < 3",400,0,400);
  TH1D* h_tag_h_ntwrs_goodEdiff_ = new TH1D("h_tag_h_ntwrs_goodEdiff","tag h tower multiplicity Ediff < 1.5",400,0,400);
  TH1D* h_tag_h_dR_highEdiff_ = new TH1D("h_tag_h_dR_highEdiff","tag h #DeltaR(candidate,other candidates) Ediff > 3",200,0,5);
  TH1D* h_tag_h_dR_lowEdiff_ = new TH1D("h_tag_h_dR_lowEdiff","tag h #DeltaR(candidate,other candidates) Ediff < 3",200,0,5);
  TH1D* h_tag_h_ncand_highEdiff_ = new TH1D("h_tag_h_ncand_highEdiff","tag h number of candidates in jet Ediff > 3",200,0,200);
  TH1D* h_tag_h_ncand_lowEdiff_ = new TH1D("h_tag_h_ncand_lowEdiff","tag h number of candidates in jet Ediff < 3",200,0,200);
  // h
  TH1D* h_tag_h0_EcalEfrac_ = new TH1D("h_tag_h0_EcalEfrac","tag h0 candidate EcalEfrac",200,0,1.5);
  // HFhad
  TH1D* h_tag_HFhad_eta_ = new TH1D("h_tag_HFhad_eta","tag HFhad candidate #eta",200,-5.5,5.5);
  TH1D* h_tag_HFhad_EcalE_loweta_ = new TH1D("h_tag_HFhad_EcalE_loweta","tag HFhad candidate EcalE #eta < 3",200,0,200);
  TH1D* h_tag_HFhad_EcalE_mideta_ = new TH1D("h_tag_HFhad_EcalE_mideta","tag HFhad candidate EcalE 3 < #eta < 3.5",200,0,200);
  TH1D* h_tag_HFhad_EcalE_higheta_ = new TH1D("h_tag_HFhad_EcalE_higheta","tag HFhad candidate EcalE #eta > 3.5",200,0,200);
  TH1D* h_tag_HFhad_EcalEfrac_ = new TH1D("h_tag_HFhad_EcalEfrac","tag HFhad candidate EcalEfrac",200,0,1.5);
  TH1D* h_tag_HFhad_emf_ = new TH1D("h_tag_HFhad_emf","tag HFhad candidate EM fraction",200,0,1.5);
  // egammaHF
  TH1D* h_tag_egammaHF_EcalEfrac_ = new TH1D("h_tag_egammaHF_EcalEfrac","tag egammaHF candidate EcalEfrac",200,0,1.5);

  int nEvents = tree->GetEntries();
  cout << "Running over " << nEvents << " events" << endl;
  //nEvents = 5;
  for(int iEvent=0; iEvent<nEvents; iEvent++){
    if(iEvent % 1000 == 0){
      cout << "Processing event " << iEvent << endl;
    }
    tree->GetEntry(iEvent);
    
    //////////////////////////
    // Fill tag histograms
    //////////////////////////

    float tag_jet_rechit_E = 0;
    float tag_jet_rechit_E_cut = 0;
    float tag_jet_hadEcalE = 0;
    float tag_jet_negE_n = 0;
    h_tag_rechitspercandidate_->Fill((double)tpfjet_ntwrs_/(double)tpfjet_had_n_);
    for(int i=0; i<tpfjet_had_n_; i++){
      float cand_rechit_E = 0;
      float cand_rechit_E_cut = 0;
      float cand_ntwrs = 0;
      vector<float> cand_rechit_E_vector;
      tag_jet_hadEcalE += tpfjet_had_EcalE_->at(i);
      for(int j=0; j<tpfjet_ntwrs_; j++){
	if(tpfjet_twr_hadind_->at(j) == i &&  tpfjet_twr_hade_->at(j) > 0.0){
	  tag_jet_rechit_E += tpfjet_twr_hade_->at(j)*tpfjet_twr_frac_->at(j);
	  cand_rechit_E += tpfjet_twr_hade_->at(j)*tpfjet_twr_frac_->at(j);
	  cand_rechit_E_vector.push_back(tpfjet_twr_hade_->at(j)*tpfjet_twr_frac_->at(j));
	  if(tpfjet_twr_frac_->at(j) > 0.05){
	    tag_jet_rechit_E_cut += tpfjet_twr_hade_->at(j)*tpfjet_twr_frac_->at(j);
	    cand_rechit_E_cut += tpfjet_twr_hade_->at(j)*tpfjet_twr_frac_->at(j);
	  }
	  //if(tpfjet_twr_frac_->at(j) > 0.05) cand_ntwrs++;
	  cand_ntwrs++;
	}
	else if(tpfjet_twr_hadind_->at(j) == i){
	  h_tag_jet_negEraw_->Fill(tpfjet_twr_hade_->at(j));
	  h_tag_jet_negEtimesFrac_->Fill(tpfjet_twr_hade_->at(j)*tpfjet_twr_frac_->at(j));
	  tag_jet_negE_n++;
	}
      }

      // Candidate types
      
      if(tpfjet_had_id_->at(i) == 0){
	float Ediff = (cand_rechit_E - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i);
	h_tag_h_Ediff_->Fill(Ediff);
	h_tag_h_EvsEdiff_->Fill(Ediff,tpfjet_had_E_->at(i));
	h_tag_h_Ediff_EcalE_->Fill((cand_rechit_E + tpfjet_had_EcalE_->at(i) - tpfjet_had_E_->at(i))/(tpfjet_had_E_->at(i)));// - tpfjet_had_EcalE_->at(i)));
	h_tag_h_Ediff_cut_->Fill((cand_rechit_E_cut - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
	if(cand_rechit_E == 0){
	  h_tag_h_Ediff_rawHcalE_->Fill(-1);
	}
	else{
	  h_tag_h_Ediff_rawHcalE_->Fill((cand_rechit_E - tpfjet_had_rawHcalE_->at(i))/tpfjet_had_rawHcalE_->at(i));
	}

	if(Ediff > 3){
	  h_tag_h_ntwrs_highEdiff_->Fill(cand_ntwrs);
	  h_tag_h_ncand_highEdiff_->Fill(tpfjet_had_n_);
	  for(unsigned int k=0; k<cand_rechit_E_vector.size(); k++){
	    h_tag_h_twrE_highEdiff_->Fill(cand_rechit_E_vector.at(k));
	  }
	  for(int l=0; l<tpfjet_had_n_; l++){
	    if(l != i){
	      h_tag_h_dR_highEdiff_->Fill(deltaR(tpfjet_had_px_->at(i),tpfjet_had_py_->at(i),tpfjet_had_pz_->at(i),tpfjet_had_px_->at(l),tpfjet_had_py_->at(l),tpfjet_had_pz_->at(l)));
	    }
	  }
	}
	else{
	  h_tag_h_ntwrs_lowEdiff_->Fill(cand_ntwrs);
	  h_tag_h_ncand_lowEdiff_->Fill(tpfjet_had_n_);
	  for(unsigned int k=0; k<cand_rechit_E_vector.size(); k++){
	    h_tag_h_twrE_lowEdiff_->Fill(cand_rechit_E_vector.at(k));
	  }
	  for(int l=0; l<tpfjet_had_n_; l++){
	    if(l != i){
	      h_tag_h_dR_lowEdiff_->Fill(deltaR(tpfjet_had_px_->at(i),tpfjet_had_py_->at(i),tpfjet_had_pz_->at(i),tpfjet_had_px_->at(l),tpfjet_had_py_->at(l),tpfjet_had_pz_->at(l)));
	    }
	  }
	  if(Ediff < 1.5){
	    h_tag_h_ntwrs_goodEdiff_->Fill(cand_ntwrs);
	  }
	}

	// Track
	if(tpfjet_had_candtrackind_->at(i) > -1){
	  float track_pt = sqrt(tpfjet_candtrack_px_->at(tpfjet_had_candtrackind_->at(i))*tpfjet_candtrack_px_->at(tpfjet_had_candtrackind_->at(i)) + tpfjet_candtrack_py_->at(tpfjet_had_candtrackind_->at(i))*tpfjet_candtrack_py_->at(tpfjet_had_candtrackind_->at(i)));
	  if(cand_rechit_E == 0){
	    if(tpfjet_had_EcalE_->at(i) == 0){
	      h_tag_h_pt_norechits_noEcalE_->Fill(track_pt);
	    }
	    else{
	      h_tag_h_pt_norechits_EcalE_->Fill(track_pt);
	    }
	    
	    h_tag_h_EcalEfrac_norechits_->Fill(tpfjet_had_EcalE_->at(i)/tpfjet_had_E_->at(i));
	    if(track_pt < 1){
	      h_tag_h_EcalEfrac_norechits_lowPt_->Fill(tpfjet_had_EcalE_->at(i)/tpfjet_had_E_->at(i));
	    }
	    else if(track_pt > 1 && track_pt < 10){
	      h_tag_h_EcalEfrac_norechits_midPt_->Fill(tpfjet_had_EcalE_->at(i)/tpfjet_had_E_->at(i));
	    }
	    else{
	      h_tag_h_EcalEfrac_norechits_highPt_->Fill(tpfjet_had_EcalE_->at(i)/tpfjet_had_E_->at(i));
	    }
	  }
	  else{
	    h_tag_h_pt_rechits_->Fill(track_pt);
	  }
	}
      }
      else if(tpfjet_had_id_->at(i) == 1){
	h_tag_h0_Ediff_->Fill((cand_rechit_E - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
	h_tag_h0_EvsEdiff_->Fill((cand_rechit_E - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i),tpfjet_had_E_->at(i));
	h_tag_h0_Ediff_EcalE_->Fill((cand_rechit_E + tpfjet_had_EcalE_->at(i) - tpfjet_had_E_->at(i))/(tpfjet_had_E_->at(i)));// - tpfjet_had_EcalE_->at(i)));
	h_tag_h0_EcalEfrac_->Fill(tpfjet_had_EcalE_->at(i)/tpfjet_had_E_->at(i));
      }
      else if(tpfjet_had_id_->at(i) == 2){
	h_tag_HFhad_Ediff_->Fill((cand_rechit_E - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
	h_tag_HFhad_Ediff_EcalE_->Fill((cand_rechit_E + tpfjet_had_EcalE_->at(i) - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
	TLorentzVector candvec(tpfjet_had_px_->at(i),tpfjet_had_py_->at(i),tpfjet_had_pz_->at(i),tpfjet_had_E_->at(i));
	h_tag_HFhad_eta_->Fill(candvec.Eta());
	if(fabs(candvec.Eta()) < 3){
	  h_tag_HFhad_EcalE_loweta_->Fill(tpfjet_had_EcalE_->at(i));
	}
	else if(fabs(candvec.Eta()) < 3.5){
	  h_tag_HFhad_EcalE_mideta_->Fill(tpfjet_had_EcalE_->at(i));
	}
	else{
	  h_tag_HFhad_EcalE_higheta_->Fill(tpfjet_had_EcalE_->at(i));
	}
	h_tag_HFhad_EcalEfrac_->Fill(tpfjet_had_EcalE_->at(i)/cand_rechit_E);
	h_tag_HFhad_emf_->Fill(tpfjet_had_emf_->at(i));
      }
      else if(tpfjet_had_id_->at(i) == 3){
	h_tag_egammaHF_Ediff_->Fill((cand_rechit_E - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
	h_tag_egammaHF_Ediff_EcalE_->Fill((cand_rechit_E + tpfjet_had_EcalE_->at(i) - tpfjet_had_E_->at(i))/tpfjet_had_E_->at(i));
	h_tag_egammaHF_EcalEfrac_->Fill(tpfjet_had_EcalE_->at(i)/cand_rechit_E);
      }
    }
    float tag_jet_E = tag_jet_rechit_E + tag_jet_hadEcalE + tpfjet_unkown_E_ + tpfjet_electron_E_ + tpfjet_muon_E_ + tpfjet_photon_E_;
    float tag_jet_E_cut = tag_jet_rechit_E_cut + tag_jet_hadEcalE + tpfjet_unkown_E_ + tpfjet_electron_E_ + tpfjet_muon_E_ + tpfjet_photon_E_;
    h_tag_jet_Ediff_->Fill((tag_jet_E - tpfjet_E_)/tpfjet_E_);
    h_tag_jet_genEdiff_->Fill((tag_jet_E - tpfjet_genE_)/tpfjet_genE_);
    h_tag_jet_Ediff_cut_->Fill((tag_jet_E_cut - tpfjet_E_)/tpfjet_E_);
    h_tag_jet_negEfrac_->Fill((double)tag_jet_negE_n/(double)tpfjet_ntwrs_);

  }
  
  //////////////////////////
  // Save to file
  //////////////////////////

  TFile* fout = new TFile(output,"RECREATE");
  fout->cd();

  h_tag_jet_Ediff_->Write();
  h_tag_jet_genEdiff_->Write();
  h_tag_jet_Ediff_cut_->Write();
  h_tag_jet_negEraw_->Write();
  h_tag_jet_negEtimesFrac_->Write();
  h_tag_jet_negEfrac_->Write();
  h_tag_h_Ediff_->Write();
  h_tag_h0_Ediff_->Write();
  h_tag_HFhad_Ediff_->Write();
  h_tag_egammaHF_Ediff_->Write();
  h_tag_h_EvsEdiff_->Write();
  h_tag_h0_EvsEdiff_->Write();
  h_tag_h_Ediff_EcalE_->Write();
  h_tag_h0_Ediff_EcalE_->Write();
  h_tag_HFhad_Ediff_EcalE_->Write();
  h_tag_h_Ediff_rawHcalE_->Write();
  h_tag_egammaHF_Ediff_EcalE_->Write();
  h_tag_h_Ediff_cut_->Write();
  h_tag_rechitspercandidate_->Write();
  h_tag_h_pt_rechits_->Write();
  h_tag_h_pt_norechits_EcalE_->Write();
  h_tag_h_pt_norechits_noEcalE_->Write();
  h_tag_h_EcalEfrac_norechits_->Write();
  h_tag_h_EcalEfrac_norechits_lowPt_->Write();
  h_tag_h_EcalEfrac_norechits_midPt_->Write();
  h_tag_h_EcalEfrac_norechits_highPt_->Write();
  h_tag_h_twrE_highEdiff_->Write();
  h_tag_h_twrE_lowEdiff_->Write();
  h_tag_h_ntwrs_highEdiff_->Write();
  h_tag_h_ntwrs_lowEdiff_->Write();
  h_tag_h_ntwrs_goodEdiff_->Write();
  h_tag_h_dR_highEdiff_->Write();
  h_tag_h_dR_lowEdiff_->Write();
  h_tag_h_ncand_highEdiff_->Write();
  h_tag_h_ncand_lowEdiff_->Write();
  h_tag_h0_EcalEfrac_->Write();
  h_tag_HFhad_eta_->Write();
  h_tag_HFhad_EcalE_loweta_->Write();
  h_tag_HFhad_EcalE_mideta_->Write();
  h_tag_HFhad_EcalE_higheta_->Write();
  h_tag_HFhad_EcalEfrac_->Write();
  h_tag_HFhad_emf_->Write();
  h_tag_egammaHF_EcalEfrac_->Write();
  
  fout->Close();
  
  return;
}

