#define LvdpSkimmer_cxx
/**
Author  : Lovedeep Kaur Saini, KSU.
Adapted From: Original tool by :Andrius

Purpose: Run over the Ntuples, and perform the likelihood minimization
 to obtain HCAL correction factors.

Usage : root -l
[] .L LvdpSkimmer.C++
[] LvdpSkimmer k
[] k.Loop()

****************************************************************/
#include "LvdpSkimmer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TObject.h>
#include <cmath>
#include "TMinuit.h"
#include <sstream>
#include <fstream>
#include <table.h>
const int NUMTOWERS = 83;
const int MAXIETA = 41;
int L=0;

std::ofstream fout("event.txt");

//A very simple structure to hold an event.
class GamJetEvent{
public:
  double fWeight;
  double fTagEcalE;
  double fTagPt, fTagEta, fTagPhi;
  bool   fUse2ndJet;
  double fJet2Pt, fJet2Eta, fJet2Phi, fJet2Scale;
  double fProbePt, fProbeEta, fProbePhi;
  double fProbePUE;
  std::map<Int_t, Double_t> fProbeHcalE;
  double fProbeEcalE;
  void PrintEvt();
} ;


//A simple class to hold a set of GamJetEvent objects.
class GamJetEventSet : public TObject{
public :
  std::vector<GamJetEvent> events;

  double fResolution;
  double EstimateResolution(int setValue, const TArrayD& respcorr);
  double CalculateGamJetBalance(GamJetEvent& e, const TArrayD& respcorr);
  double GetLikelihoodDistance(const TArrayD& respcorr) const;
  void  doFit(TArrayD& respcorr, TArrayD& respcorre, std::vector<int> fix);
  static void FCN (Int_t &npar, 
		   Double_t *gin, 
		   Double_t &f, 
		   Double_t *par, 
		   Int_t flag);
  ClassDef(GamJetEventSet, 1);
};

void LvdpSkimmer::Loop(){
  //  Loop(0,0);
  // Loop(1,0);
  //Loop(0,1);
  Loop(1,1);
}
void LvdpSkimmer::Loop(bool doPfCorr, bool use2ndJet)
{
  
  if (fChain == 0) return;
  
  table* puTable = new table("PFchs_PV_DATA.txt");
  
  
  GamJetEventSet eventsForFit;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  TH1D* h;
  TH1D* histIEta  = new TH1D("histIEta","histIEta",NUMTOWERS+1,-MAXIETA-0.5,MAXIETA+0.5);
  TH1D* hCounter  = new TH1D("hCounter","hCounter",20,0.,20);
  TH1D* h_2ndleadJetCorrEt=new TH1D("h_2ndleadJetCorrEt","h_2ndleadJetCorrEt;corrected #it{E}_{T}^{j2} [GeV];count;",200,0,200);
  TH1D*  h_3rdleadJetCorrEt=new TH1D("h_3rdleadJetCorrEt","h_3rdleadJetCorrEt;corrected #it{E}_{T}^{j3} [GeV];count",200,0,200);
  // Numbers of alpha, pTGamma and eta bins
  const int nPtBins=12; // numbers of Pt bins
  //const int nPtBins=1; // numbers of Pt bins
  const int nAlphaBins=6; // numbers of alpha bins
  const int nEtaBins=2; // number of eta bins
  // Definition of all bin intervals
  const double ptBins[nPtBins+1]={22,36,60,88,105,148.5,165,176,200,250,300,400,1000000};
  // const double ptBins[nPtBins+1]={400,1000000};
  const double alphaBins[nAlphaBins+1]={0.0,0.075,0.10,0.125,0.150,0.175,0.2};
  const double etaBins[nEtaBins+1]={0.0,1.1,2.4};
  if (fChain == 0) return;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (ientry%1000 == 0 ){
      float perc = float(ientry)*100./float(nentries);
      printf("Processed Events          %i       out of %i      %.0f %%\n",int(ientry), int (nentries), perc);
    }
 
    double evtWt = EventWeight;
    double weight = 1;//evtWt;
    //Now let us see if the jet passes tight id criterion.
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
      
      //Selection Condition
      hCounter->Fill(0.,weight);
      
      //Upper and lower Bounds on all Binning variables
      if((tagPho_pt< ptBins[0] || tagPho_pt > ptBins[nPtBins] || fabs(tagPho_eta) > etaBins[nEtaBins]) ||
         (ppfjet_pt*ppfjet_scale)< ptBins[0] || (ppfjet_pt*ppfjet_scale) > ptBins[nPtBins] || fabs(ppfjet_eta) > etaBins[nEtaBins] || (fabs(tagPho_eta) > 1.4442 && fabs(tagPho_eta)<1.566) )
        continue;
      hCounter->Fill(1,weight);
      
      //if(!(photonTrig_fired->at(6) !=-1 && photonTrig_fired->at(6)!=0)) continue;
      //hCounter->Fill(2,weight);

      if(!(pf_NPV>=1))continue;
      hCounter->Fill(3,weight);
      
      //select tight photon
      if(!tagPho_idTight) continue;
      hCounter->Fill(4,weight);
      
      //select leading tight jet
      if(!jetPassTightId) continue;
      hCounter->Fill(5,weight);
      
      
      //calculate dphi (photon, jet)
      double dphi_Pho1stJet = calc_dphi(tagPho_phi, ppfjet_phi);
      if(dphi_Pho1stJet < 2.95)continue;
      hCounter->Fill(6,weight);
      
      double alpha = (pfjet2_pt*pfjet2_scale)/tagPho_pt;//also called alpha.     
      if(alpha > alphaBins[nAlphaBins])continue;
      hCounter->Fill(7,weight);

      h_2ndleadJetCorrEt->Fill(pfjet2_pt*pfjet2_scale,evtWt);
      h_3rdleadJetCorrEt->Fill(pf_thirdjet_pt*pf_thirdjet_scale,evtWt);
            
      if(fabs(tagPho_eta)>2.4)continue;
      hCounter->Fill(8,weight);
      if((tagPho_pt >= ptBins[0] && tagPho_pt < ptBins[1]) && photonTrig_fired->at(0) !=1) continue;
      hCounter->Fill(11,weight);
      if((tagPho_pt >= ptBins[1] && tagPho_pt < ptBins[2]) && photonTrig_fired->at(1) !=1) continue;
      hCounter->Fill(12,weight);
      if((tagPho_pt >= ptBins[2] && tagPho_pt < ptBins[3]) && photonTrig_fired->at(2) !=1) continue;
      hCounter->Fill(13,weight);
      if((tagPho_pt >= ptBins[3] && tagPho_pt < ptBins[4]) && photonTrig_fired->at(3) !=1) continue;
      hCounter->Fill(14,weight);
      if((tagPho_pt >= ptBins[4] && tagPho_pt < ptBins[5]) && photonTrig_fired->at(4) !=1) continue;
      hCounter->Fill(15,weight);
      if((tagPho_pt >= ptBins[5] && tagPho_pt < ptBins[6]) && photonTrig_fired->at(5) !=1) continue;
      hCounter->Fill(16,weight);
      if((tagPho_pt >= ptBins[6]) && photonTrig_fired->at(6) !=1) continue;
      // if(photonTrig_fired->at(6)!=1)continue;     //Selection Condition
      hCounter->Fill(17,weight);
      
      if(photonTrig_fired->at(6)==1)hCounter->Fill(19,weight);
     
      //Cut on second Jetpt > 10 GeV (9. CUT)
      if(pfjet2_pt*pfjet2_scale>10)      hCounter->Fill(9,weight);
      //      if(pfjet2_pt*pfjet2_scale>50)continue;
      //      hCounter->Fill(18,weight);
      if(pf_thirdjet_pt*pf_thirdjet_scale>10)      hCounter->Fill(10,weight);
      // if(pf_thirdjet_pt*pf_thirdjet_scale>50)continue;
      //  hCounter->Fill(20,weight);
      
      //Now let us obtain jet kinematics
      double jetEn = ppfjet_E;
      double jetEta =  ppfjet_eta;
      double jetPhi = ppfjet_phi;
      double jetEcalE=0;
      
      for (unsigned int i=0; i< ppfjet_had_EcalE->size(); i++) 
	jetEcalE+= ppfjet_had_EcalE->at(i);
      double jetHcalEn_noRecHits = 0;//temp
      for (unsigned int i=0; i<ppfjet_had_id->size(); ++i) {
	int trackIdx= ppfjet_had_candtrackind->at(i);
	if ((ppfjet_had_id->at(i) == 0) && // charged hadron
	    (trackIdx>-1) // has a track
	    && (ppfjet_had_ntwrs->at(i) == 0) // has no recHits
	    ) {
	  jetHcalEn_noRecHits += sqrt( pow(ppfjet_candtrack_px->at(trackIdx),2) +
				       pow(ppfjet_candtrack_py->at(trackIdx),2) +
				     pow(ppfjet_candtrack_pz->at(trackIdx),2) )
	    - ppfjet_had_EcalE->at(i);
	}
    }
      
      double jetEcalEn = jetHcalEn_noRecHits+jetEcalE+ppfjet_unkown_E + ppfjet_electron_E + ppfjet_muon_E + ppfjet_photon_E;
      
      //EnMap : E deposit as function of iEta
      //EnMap For leading jet Jet
      std::map<int, double> hcalE;
      //Loop over tower collection in event,
      //for i-th guy, find energy and store as pair,
      // i-Energy. (Under some selection conditions.
      int numTowers = ppfjet_twr_hade->size();
      for(int i=0; i< numTowers; ++i){
	double iTwrEn = ppfjet_twr_hade->at(i);
	if(iTwrEn<=0) continue;
	int iTwrIEta = ppfjet_twr_ieta->at(i);
	double iTwrFrac = ppfjet_twr_frac->at(i);
	int    iTwrClusId = ppfjet_twr_clusterind->at(i);
	double iTwrClusDr = 0;
	if(iTwrClusId>-1) iTwrClusDr = ppfjet_cluster_dR->at(iTwrClusId);
	//Now an if-else situation.
	double iTwrContrib = iTwrEn*iTwrFrac;
	
	if( (iTwrClusId<0) || (iTwrClusDr<0.5)){
	  if(iTwrContrib > 1e-4){
	    hcalE[iTwrIEta] += iTwrContrib;
	  }//if
	}
      }
      //    histNTowers->Fill(hcalE.size());
      
      
      /*
	At this point I have following information:
	iEta of a tower : total energy deposit       
    */
      std::map<int, double>::iterator it = hcalE.begin();
      for(; it!=hcalE.end(); it++){
	histIEta->Fill(it->first,evtWt);
      }
      

      //Determine PU correction factor for the probe;
      double puEnergy = 0;
      if(doPfCorr)
	puEnergy = puTable->density(ppfjet_eta,pf_NPV)*ppfjet_area*cosh(ppfjet_eta);
      
      
      
      GamJetEvent evt;
      evt.fWeight= evtWt;
      evt.fTagEcalE= tagPho_energy;
      evt.fTagPt = tagPho_pt;
      evt.fTagEta= tagPho_eta;
      evt.fTagPhi= tagPho_phi;
      evt.fProbeEta = ppfjet_eta; 
      evt.fProbePhi= ppfjet_phi;
      evt.fProbeHcalE= hcalE;
      evt.fProbeEcalE= jetEcalEn;
      evt.fProbePUE = puEnergy;
      evt.fJet2Pt =   pfjet2_pt; 
      evt.fJet2Scale = pfjet2_scale;
      evt.fJet2Eta = pfjet2_eta;
      evt.fJet2Phi =  pfjet2_phi;
      evt.fUse2ndJet = use2ndJet;
      eventsForFit.events.push_back(evt);
  }   
  
  
  
  
  //Time to start fit procedure.

  std::cout<<">>-----------FITTER COMES HERE--------------<<"<<std::endl;
  
  //First we find very small bins, and take them out of fit.
  std::vector<int> fix(NUMTOWERS,0);
  double max = histIEta->GetMaximum();
  for(int nbin=1; nbin < histIEta->GetNbinsX(); nbin++)
    if(histIEta->GetBinContent(nbin)/max < 0.40)
      fix[nbin-1] = 1;
  
  /*  for(int i=0; i!=eventsForFit.events.size(); i++){
      GamJetEvent e = eventsForFit.events[i];
      e.PrintEvt();
      }*/
  //Command the machinery to do a fit.
  std::cout<<eventsForFit.events.size()<<std::endl;
  TArrayD respcorr, respcorre;
  eventsForFit.doFit(respcorr, respcorre, fix);
  
  

  //Store results in a histogram
  TString fileName  = "Lvdeep_Outfile";
  if(doPfCorr)fileName+="_puCorr";
  else fileName +="_NoPuCorr";

  if(use2ndJet) fileName+="_2ndJetCorr";
  else fileName+="_No2ndJetCorr";

  fileName +=".root";
  
  TFile *oFile = new TFile(fileName,"RECREATE");  
  h=new TH1D("name","title",NUMTOWERS,-MAXIETA-0.5,MAXIETA+0.5);
  for(int i=1; i<=NUMTOWERS; i++) {
    h->SetBinContent(i, respcorr[i-1]);
     h->SetBinError(i, respcorre[i-1]);
  }
  //Write Histogram and close the file.
  h->Write();
  histIEta->Write();
  hCounter->Write(); 
  h_2ndleadJetCorrEt->Write();
  h_3rdleadJetCorrEt->Write();
  oFile->Close();
}


//Miniuit needs this fucntion.
void GamJetEventSet::FCN(Int_t &npar, Double_t*, 
			 Double_t &f, Double_t *par, 
		    Int_t){
  // get the relevant data
  const GamJetEventSet* data= dynamic_cast<const GamJetEventSet*> (gMinuit->GetObjectFit());
  TArrayD respcorr;
  respcorr.Set(NUMTOWERS, par);
  f = data->GetLikelihoodDistance(respcorr);
  return;
}




/**
Adapted From : David et. al. @ Dijet balance team
*********************************************************/
double GamJetEventSet::GetLikelihoodDistance(const TArrayD& respcorr) const
{
  
  double total = 0;
  
  // loop over each events
  int i=-1;
  for(std::vector<GamJetEvent>::const_iterator it=events.begin(); it!=events.end(); ++it) {
    i++;
    //This iterator represents i-th event.
    
    //PHOTON INFORMATION ACCESSED HERE.
    double gamPt = it->fTagPt;
    
    //JET INFORMATION ACCESSED HERE.
    //... AIM IS TO CALC JET'S TRANSVERSE ET.
    double jetEcal = it->fProbeEcalE;
    double jetHcal = 0;
    double jetHf   = 0;
    double jetPUE = it->fProbePUE;
    
    //Access Ecal and Hcal energy of jet.
    for(std::map<int,double>::const_iterator mapit = it->fProbeHcalE.begin();
	mapit!=it->fProbeHcalE.end(); mapit++){
      int ieta=mapit->first;
      double energy=mapit->second;
      int index=ieta+MAXIETA;
      
      if(std::fabs(ieta)>29)
	jetHf += respcorr[index]*energy;
      else
	jetHcal += respcorr[index]*energy;
    }
    
    double jetEt = (jetEcal+jetHcal+jetHf-jetPUE); //not-there yet
    jetEt = jetEt/std::cosh(it->fProbeEta);//thats it.
    
    if(it->fUse2ndJet){
      double tphi = it->fTagPhi;
      double pjphi = it->fProbePhi;
      double sjphi= it->fJet2Phi;
      float dphiPJSJ=fabs(pjphi-sjphi);
      const float cPi= 4*atan(1);
      while (dphiPJSJ>cPi) dphiPJSJ = fabs(2*cPi - dphiPJSJ);
      float dphiTPSJ=fabs(tphi-sjphi);
      while (dphiTPSJ>cPi) dphiTPSJ = fabs(2*cPi - dphiTPSJ);
      
      
      double dphi_Pho2ndJet = dphiTPSJ;
      double dphi_1stJet2ndJet = dphiPJSJ;
      
      bool jetHemiS = 0;
      bool photonHemiS=0;
      if(dphi_1stJet2ndJet < dphi_Pho2ndJet) jetHemiS=1;
      else photonHemiS=1;
 
      Double_t jetEx = jetEt*std::cos(it->fProbePhi);
      Double_t jetEy = jetEt*std::sin(it->fProbePhi);

      double gamPx = gamPt*std::cos(it->fTagPhi);
      double gamPy = gamPt*std::sin(it->fTagPhi);
      
      if(jetHemiS){
	jetEx +=   (it->fJet2Pt*it->fJet2Scale)*std::cos(it->fJet2Phi);
	jetEy +=   (it->fJet2Pt*it->fJet2Scale)*std::sin(it->fJet2Phi);
      }
      
      else if(photonHemiS){
	gamPx +=   (it->fJet2Pt*it->fJet2Scale)*std::cos(it->fJet2Phi);
	gamPy +=   (it->fJet2Pt*it->fJet2Scale)*std::sin(it->fJet2Phi); 
      }
      jetEt = std::sqrt(jetEx*jetEx+jetEy*jetEy);
      gamPt = std::sqrt(gamPx*gamPx+gamPy*gamPy);
    }
    
    //std::cout<<"Photon pt :"<<gamPt<<std::endl;
    //std::cout<<"Jet Et :"<<jetEt<<std::endl;
    //Now calculate the balance

    double balance = (gamPt-jetEt);
    //    double balance = 2*(gamPt-jetEt)/(gamPt+jetEt);
    double resolution = 1.0/std::sqrt(it->fWeight);
    if(i<20 && 0){
      std::cout<<"============"<<std::endl;
      std::cout<<"Probe Et: "<<jetEt<<std::endl;
      std::cout<<"Tag Et: "<<gamPt<<std::endl;
      std::cout<<"Probe Hcal E: "<<jetHcal<<std::endl;
      std::cout<<"Probe Ecal E: "<<jetEcal<<std::endl;
      std::cout<<"i "<<i<<" diff "<<balance<<" w "<<it->fWeight<<std::endl;
      
    }
    
    //Everyone hates long varnames.
    double B = balance;
    double dB = resolution;
    total += std::pow(B/fResolution,2);
    //  cout<<"Balance??????:  "<<events.size()<<'\t'<<B<<'\t'<<dB<<'\t'<<fResolution<<endl;
    //@di-jet people did is as follows:
    //total += 0.5*(std::log(dB*dB)+B*B/dB/dB); //total likelihood
  }
  return total;
}


void GamJetEvent::PrintEvt(){
  /*
  double fWeight;
  double fTagEcalE;
  double fTagPt, fTagEta, fTagPhi;
  double fProbePt, fProbeEta, fProbePhi;
  std::map<Int_t, Double_t> fProbeHcalE;
  double fProbeEcalE; 
  */
  fout<<endl;
  fout<<fProbeEcalE;
  for(int i=0; i!=NUMTOWERS; i++){
    std::map<Int_t, Double_t>::iterator it = this->fProbeHcalE.find(i-MAXIETA);
    if(it==this->fProbeHcalE.end())fout<< " "<<0.0;
    else fout<< "  " << it->second;
  }
  fout<<" "<<this->fTagPt*std::cosh(this->fTagEta);
  fout<<endl;
  return;
}





void GamJetEventSet::doFit(TArrayD& respcorr, TArrayD& respcorre, std::vector<int> fix){
  //Setup the initial response corrections

  int fPrintLevel=5;
  double fParStep=0.10;
  double fParMin=0.0;
  double fParMax=2.0;//0.0;
  

  const int maxIetaFixed=1;
  
  double array[NUMTOWERS] =   { 1.0, 1.0, 1.0, 1.0, 1.0, //1
				1.0, 1.0, 1.0, 1.0, 1.0, //2
				1.0, 1.0, 1.0, 1.0, 1.0, //3
				1.0, 1.0, 1.0, 1.0, 1.0, //4
				1.0, 1.0, 1.0, 1.0, 1.0, //5
				1.0, 1.0, 1.0, 1.0, 1.0, //6
				1.0, 1.0, 1.0, 1.0, 1.0, //7
				1.0, 1.0, 1.0, 1.0, 1.0, //8
				1.0, 1.0, 1.0, 1.0, 1.0, //9
				1.0, 1.0, 1.0, 1.0, 1.0, //10
				1.0, 1.0, 1.0, 1.0, 1.0, //11
				1.0, 1.0, 1.0, 1.0, 1.0, //12
				1.0, 1.0, 1.0, 1.0, 1.0, //13
				1.0, 1.0, 1.0, 1.0, 1.0, //14
				1.0, 1.0, 1.0, 1.0, 1.0, //15
				1.0, 1.0, 1.0, 1.0, 1.0, //16
				1.0, 1.0, 1.0 };  //17
  TArrayD respCorrInit;
  respCorrInit.Set(NUMTOWERS, array);
  //Now estimate the resolution.
  fResolution = EstimateResolution(1,respCorrInit);
  std::cout<<"RESOLUTION: "<<fResolution<<std::endl;
  
  // set the number of parameters to be the number of towers
  TMinuit *gMinuit=new TMinuit(NUMTOWERS);
  //gMinuit->SetPrintLevel(fPrintLevel);
  gMinuit->SetErrorDef(0.5); // for a log likelihood
  gMinuit->SetFCN(FCN);
  gMinuit->SetObjectFit(this);
  
  for(int i=0; i<respCorrInit.GetSize(); i++) {
    int ieta=-MAXIETA+i;
    std::ostringstream oss;
    oss << "Tower ieta: " << ieta;
    gMinuit->DefineParameter(i, oss.str().c_str(), respCorrInit[i], fParStep, fParMin, fParMax);
    if((ieta>=-maxIetaFixed && ieta<=maxIetaFixed)) gMinuit->FixParameter(i);
    if(fix[i]==1)gMinuit->FixParameter(i);
    //if(!(ieta>=-10 && ieta<=10))gMinuit->FixParameter(i);
  }
  
  gMinuit->Migrad();
  
  TArrayD results(NUMTOWERS);
  TArrayD errors(NUMTOWERS);
  for(int i=0; i<results.GetSize(); i++) {
    Double_t val, error;
    gMinuit->GetParameter(i, val, error);
    results[i]=val;
    errors[i]=error;
  }
  respcorr=results;
  respcorre=errors;
  return;
}


double GamJetEventSet::CalculateGamJetBalance( GamJetEvent& e, const TArrayD& respcorr) {
  //PHOTON INFORMATION ACCESSED HERE.
  int gamPt = e.fTagPt;
  
  //JET INFORMATION ACCESSED HERE.
  //... AIM IS TO CALC JET'S TRANSVERSE ET.
  int jetEcal = e.fProbeEcalE;
  int jetHcal = 0;
  int jetHf   = 0;
  //Access Ecal and Hcal energy of jet.
  for(std::map<int,double>::const_iterator mapit = e.fProbeHcalE.begin();
      mapit!=e.fProbeHcalE.end(); mapit++){
    int ieta=mapit->first;
    double energy=mapit->second;
    int index=ieta+MAXIETA;
    if(std::fabs(ieta)>29)
      jetHf += respcorr[index]*energy;
    else
      jetHcal += respcorr[index]*energy;
  }
    
  double jetEt = (jetEcal+jetHcal+jetHf); //not-there yet
  jetEt = jetEt/std::cosh(e.fProbeEta);//thats it.
    //std::cout<<"Photon pt :"<<gamPt<<std::endl;
    //std::cout<<"Jet Et :"<<jetEt<<std::endl;
    //Now calculate the balance
  double balance = (gamPt-jetEt);///(gamPt+jetEt);
  return balance;
}


/*
  Adapted from : Andrius's package.
*************************************************************/
double GamJetEventSet::EstimateResolution(int setValue, const TArrayD& respcorr) {
  if (setValue) fResolution=0.;
  double sumW=0, sumWX=0, sumWXX=0.;
  int allZeros=1;
  for (unsigned int i=0; i<events.size(); ++i) {
    GamJetEvent e= events[i];
    Double_t w= e.fWeight;
    Double_t x= this->CalculateGamJetBalance(e,respcorr);
    if (allZeros && (x!=Double_t(0))) allZeros=0;
    sumW += w;
    sumWX += w*x;
    sumWXX += w*x*x;
  }
  if (allZeros) {
    std::cout << "EstimateResolution: allZeros detected" << std::endl;
   return 0.;
  }
  Double_t sigmaSqr= sumWXX/sumW - pow(sumWX/sumW,2);
  if (sigmaSqr<=Double_t(0)) {
    std::cout << "EstimateResolution: sigmaSqr<=0" << std::endl;
    return 0.;
  }
  Double_t resolution= sqrt(sigmaSqr);
  if (setValue) fResolution = resolution;
  return resolution;
}


double LvdpSkimmer::calc_dphi(double phi1, double phi2){
  float dphi=fabs(phi1-phi2);
  const float cPi= 4*atan(1);
  while (dphi>cPi) dphi = fabs(2*cPi - dphi);
  return dphi;
}
