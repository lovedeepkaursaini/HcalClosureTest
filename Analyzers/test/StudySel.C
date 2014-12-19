/**
Author  : Lovedeep Kaur Saini, KSU.
Usage : root -l
[] .L StudySel.C++
[] StudySel k
[] k.Loop()

****************************************************************/
#define StudySel_cxx
#include "StudySel.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "iostream"
#include <cmath>
using namespace std;
void StudySel::Loop()
{
   TFile *oFile  = new TFile("Hist_StudySel.root","RECREATE");
   std::map<std::string, TH1D*> hM;
   hM["Counter"] = new TH1D("Counter","Counter",10,0,10);
   TH2D* dphiVs2ndPt    = new TH2D("dphiVs2ndPt","dphiVs2ndPt",100,0.,4.,200,0.,100.);
   hM["resRJoGJ_J2Pt_Sel"] = new TH1D("hRecoJetOverGenJet_Resol_J2Pt_Sel","RecoJet/GenJet Pt ",100,0,2);
   hM["resRJoGJ_J2Pt_Sel_2ndJet"] = new TH1D("hRecoJetOverGenJet_Resol_J2Pt_Sel_2ndJet","RecoJet/GenJet Pt ",100,0,2);
   hM["resRJoGJ_dphi_Sel"] = new TH1D("hRecoJetOverGenJet_Resol_dPhi_Sel","RecoJet/GenJet Pt ",100,0,2);
   hM["dphi"]    = new TH1D("hdphi","dphi (photon,jet)",100,0.,4.);
   //by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (ientry%10000 == 0 ){
        float perc = float(ientry)*100./float(nentries);
        printf("Processed Events          %i       out of %i      %.0f %%\n",int(ientry), int (nentries), perc);
      }
      
      bool jetPassTightId=false;
      if(fabs(ppfjet_eta)<2.4)
	{
	  if((ppfjet_NeutralHadronFrac < 0.90) &&
	     (ppfjet_NeutralEMFrac< 0.90) &&
	     (ppfjet_nConstituents >1))jetPassTightId=true;
	}
      else if(fabs(ppfjet_eta)>2.4)
	{
	  if( (ppfjet_NeutralHadronFrac < 0.90) &&
	      (ppfjet_NeutralEMFrac< 0.90) &&
	      (ppfjet_nConstituents >1)
	      && (ppfjet_ChargedHadronFrac>0)&&
	      (ppfjet_ChargedMultiplicity>0)&&
	      (ppfjet_ChargedEMFrac<0.99))jetPassTightId=true;
	}
      
      //C U T S
      //Let see how energetic the second jet is.
      double jetGamRatio = (pfjet2_pt*pfjet2_scale)/tagPho_pt;//also called alpha.
      //Selection Condition
      hM["Counter"]->Fill(0.);
      if(!tagPho_idTight) continue;
      hM["Counter"]->Fill(1.);
      if(!jetPassTightId) continue;
      hM["Counter"]->Fill(2.);
      float phi1=tagPho_phi;
      float phi2=ppfjet_phi;
      float dphi=fabs(phi1-phi2);
      const float cPi= 4*atan(1);
      while (dphi>cPi) dphi = fabs(2*cPi - dphi);
      hM["dphi"]->Fill(dphi);
      dphiVs2ndPt->Fill(dphi,pfjet2_pt*pfjet2_scale);
      
      //if(photonTrig_fired->at(0)!=1)continue;     //Selection Condition
      bool dphi_Sel = 0;
      if((dphi > 2.95) &&(jetGamRatio<0.2) 
	 && pf_thirdjet_pt*pf_thirdjet_scale<50.
	 && pfjet2_pt*pfjet2_scale<100) dphi_Sel=1;
      
      if(dphi_Sel==1){
	hM["resRJoGJ_dphi_Sel"]->Fill(ppfjet_pt/ppfjet_genpt);
	hM["Counter"]->Fill(9.);
      }
      
      bool J2Pt_Sel = 0;
      if(pfjet2_pt*pfjet2_scale<30 && pf_thirdjet_pt*pf_thirdjet_scale<30.)J2Pt_Sel=1;
      
      if(J2Pt_Sel==1)
	{
	  hM["resRJoGJ_J2Pt_Sel"]->Fill(ppfjet_pt/ppfjet_genpt);
	  double jetEt= ppfjet_pt;
	  double gamEt = tagPho_pt;
	  Use2ndJetPt(jetEt,gamEt);
	  hM["resRJoGJ_J2Pt_Sel_2ndJet"]->Fill(jetEt/ppfjet_genpt);
	  hM["Counter"]->Fill(4.);
	}
      
   }
   std::map<std::string, TH1D*>::iterator it = hM.begin();
   for(; it!=hM.end();it++)
     (it->second)->Write();
   cout<<"effective sigma: J2Pt_Sel: "<<calc_effSigma(hM["resRJoGJ_J2Pt_Sel"])<<endl;
   cout<<"effective sigma: J2Pt_Sel 2ndJet: "<<calc_effSigma(hM["resRJoGJ_J2Pt_Sel_2ndJet"])<<endl;
   cout<<"effective sigma: dphi_Sel: "<<calc_effSigma(hM["resRJoGJ_dphi_Sel"])<<endl;
   dphiVs2ndPt->Write();
   oFile->Close();
}


//from Andrius
//******************************************//
// author: Chris Seez
//
// This function evaluates the effective
// sigma of a histogram.
//
//******************************************//
Double_t StudySel::calc_effSigma(TH1 * hist)
{
  using std::cout;
  using std::endl;
  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
    return 0.;
  }
  Double_t bwid = hist->GetBinWidth(1);
  if(bwid == 0) {
    cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
    return 0.;
  }
  //Double_t xmax = xaxis->GetXmax(); // not used
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();
  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 10.) {
    cout << "effsigma: Too few entries " << total << endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;
  Double_t rlim=0.683*total;
  Int_t nrms= Int_t(2.0*rms/(bwid)); // Set scan size to +/- rms
  //if(nrms > nb/10) nrms=nb/10; // Could be tuned...
  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=Int_t((ave-xmin)/bwid+1+iscan);
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    // ibm -- curent bin value
    // x, xj, xk -- value of the x-axis
    // jbm, kbm -- bin values
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
	jbm++;
	xj+=bwid;
	bin=hist->GetBinContent(jbm);
	total+=bin;
	if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
	kbm--;
	xk-=bwid;
	bin=hist->GetBinContent(kbm);
	total+=bin;
	if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }
  }
  if(ismin == nrms || ismin == -nrms) {
    cout << endl << endl << "CRAP!!!" << endl;
    cout << ismin << '\t' << nrms << endl;
    ierr=3;
  }
  if(ierr != 0 && ierr != 1) cout << "effsigma: Error of type " << ierr <<endl;
  return widmin;
}


// ----------------------------------------------------------------
// --------
void StudySel::Use2ndJetPt(double& jetEt, double& gamPt){
  
  double fJet2Pt =   pfjet2_pt*pfjet2_scale; 
  double fJet2Eta = pfjet2_eta;
  double fJet2Phi =  pfjet2_phi;
  double fProbePhi= ppfjet_phi;
  double fTagPhi = tagPho_phi;
  
  Double_t jetEx = jetEt*std::cos(fProbePhi);
  Double_t jetEy = jetEt*std::sin(fProbePhi);
  
  jetEx +=   0.5*(fJet2Pt)*std::cos(fJet2Phi);
  jetEy +=   0.5*(fJet2Pt)*std::sin(fJet2Phi);
  
  double gamPx = gamPt*std::cos(fTagPhi);
  double gamPy = gamPt*std::sin(fTagPhi);
  
  gamPx +=   0.5*(fJet2Pt)*std::cos(fJet2Phi);
  gamPy +=   0.5*(fJet2Pt)*std::sin(fJet2Phi); 
  
  jetEt = std::sqrt(jetEx*jetEx+jetEy*jetEy);
  gamPt = std::sqrt(gamPx*gamPx+gamPy*gamPy);
  
}
