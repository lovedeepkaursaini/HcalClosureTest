#ifndef GammaJetFitData_H_
#define GammaJetFitData_H_

//
// GammaJetFitData.h
//
//    description: Object which contains data relevant to gamma+jet response 
//                 corrections
//                 Current implementation also allows dijet minimization
//
//


#include <TObject.h>
#include <TMinuit.h>
#include <TH1D.h>
#include <TString.h>
#include <cmath>
#include <algorithm> // std::find

#include <map>
#include <vector>
#include <iostream>

//
//  numerical constants
//

#ifndef __DIJET_RESP_CORR_DATA_H__
const int MAXIETA = 41;
const int NUMTOWERS = 83;
#endif

// --------------------------------------------------------

//
// global variables and functions
//

extern TMinuit *myMinuit;

extern void gammaJet_FCN(Int_t &npar, Double_t* gin, Double_t &f,
			 Double_t *par, Int_t iflag);

extern TArrayD convert(const std::vector<Double_t> &vec);

std::ostream& operator<<(std::ostream& out, const std::map<Int_t,Double_t> &m);
std::ostream& operator<<(std::ostream& out, const std::vector<Double_t> &vec);
std::ostream& operator<<(std::ostream& out, const TArrayD &arr);

// --------------------------------------------------------

inline
Double_t calc_HcalE(const std::map<Int_t,Double_t> &hcalE) {
  double sum=0.;
  for (std::map<Int_t,Double_t>::const_iterator it=hcalE.begin();
       it!=hcalE.end(); ++it) {
    sum+= it->second;
  }
  return sum;
}

// --------------------------------------------------------

inline
Double_t calc_CorrHcalE(const TArrayD &hcalECorrCf,
		       const std::map<Int_t,Double_t> &hcalE) {
  if (hcalECorrCf.GetSize()!=NUMTOWERS) {
    std::cout << "getCorrSum: hcalECorrCf.size=" << hcalECorrCf.GetSize()
	      << " instead of " << NUMTOWERS << "\n";
    return 0.;
  }
  double sum=0.;
  //std::cout << "calc_CorrHcalE\n";
  for (std::map<Int_t,Double_t>::const_iterator it= hcalE.begin();
       it!=hcalE.end(); it++) {
    Int_t ieta= it->first;
    //std::cout << "ieta=" << ieta << " cf=" << hcalECorrCf[ieta+MAXIETA] << " " << "contr=" << it->second << "\n";
    sum += hcalECorrCf[ieta+MAXIETA]*it->second;
  }
  return sum;
}

// --------------------------------------------------------

//
// class definitions
//

// --------------------------------------------------------

//
// Container of misc items
//

class GammaJetEventAuxInfo_t : public TObject {
protected:
  Int_t fEventNo, fRunNo;
public:

  GammaJetEventAuxInfo_t(Int_t eventNo=-1, Int_t runNo=-1)
    : fEventNo(eventNo), fRunNo(runNo)
  {}

  GammaJetEventAuxInfo_t(const GammaJetEventAuxInfo_t &e) :
    TObject(e),
    fEventNo(e.fEventNo),
    fRunNo(e.fRunNo)
  {}

  Int_t GetEventNo() const { return fEventNo; }
  void SetEventNo(Int_t no) { fEventNo=no; }
  Int_t GetRunNo() const { return fRunNo; }
  void SetRunNo(Int_t no) { fRunNo=no; }


  void Assign(const GammaJetEventAuxInfo_t &e);

  GammaJetEventAuxInfo_t& operator=(const GammaJetEventAuxInfo_t &e) {
    this->Assign(e);
    return *this;
  }

  friend
  std::ostream& operator<<(std::ostream& out, const GammaJetEventAuxInfo_t &e);

 public:
  ClassDef(GammaJetEventAuxInfo_t,1)

};


// -----------------------------------------------------------
//
// Container class for a gamma+jet event
//

class GammaJetEvent_t : public TObject
{
 public:
  GammaJetEvent_t();
  GammaJetEvent_t(const GammaJetEvent_t &e);
  ~GammaJetEvent_t();

  Double_t GetWeight() const { return fWeight; }
  void SetWeight(double val) { fWeight=val; }

  Double_t GetTagEta() const { return fTagEta; }
  Double_t GetTagPhi() const { return fTagPhi; }
  Double_t GetTagEcalE() const { return fTagEcalE; }
  Double_t GetTagHcalE(Int_t ieta) const; // slower: uses iterator!
  const std::map<Int_t,Double_t>& GetTagHcalE() const { return fTagHcalE; }
  Double_t CalcTagHcalEtot() const { return calc_HcalE(fTagHcalE); }
  Double_t GetTagE() const { return fTagE; }
  Double_t GetTagEt() const { return fTagEt; }

  // assumes that all energy is ECAL. For gamma+jet events
  void SetTagEEtaPhi(double valE, double valEta, double valPhi);
  // assumes that the energy is distributed. For jet+jet events
  void SetTagEtaPhiEn(double valEta, double valPhi, double valEcalE,
		      const std::map<Int_t,Double_t> &valHcalE);

  Double_t GetTagETtot() const { return fTagEt; }
  Double_t GetTagEtot() const { return fTagE; }

  Double_t GetTagHcalE(const TArrayD &hcalECorrCf) const
  { return calc_CorrHcalE(hcalECorrCf,fTagHcalE); }

  Double_t GetTagEtot(const TArrayD &hcalECorrCf) const
  { return fTagEcalE + calc_CorrHcalE(hcalECorrCf,fTagHcalE); }

  Double_t GetTagETtot(const TArrayD &hcalECorrCf) const
  { return GetTagEtot(hcalECorrCf)/cosh(fTagEta); }

  Double_t GetProbeEta() const { return fProbeEta; }
  Double_t GetProbePhi() const { return fProbePhi; }
  Double_t GetProbeEcalE() const { return fProbeEcalE; }
  Double_t GetProbeHcalE(Int_t ieta) const; // slower: uses iterator!
  const std::map<Int_t,Double_t>& GetProbeHcalE() const { return fProbeHcalE; }
  void SetProbeEtaPhiEn(double valEta, double valPhi, double valEcalE,
			const std::map<Int_t,Double_t> &valHcalE);

  Double_t CalcProbeHcalEtot(const TArrayD &hcalECorrCf) const
  { return calc_CorrHcalE(hcalECorrCf,fProbeHcalE); }

  Double_t GetProbeETtot() const { return fProbeEt; }
  Double_t GetProbeEtot() const { return fProbeE; }

  Double_t GetProbeEtot(const TArrayD &hcalECorrCf) const
  { return fProbeEcalE + calc_CorrHcalE(hcalECorrCf,fProbeHcalE); }

  Double_t GetProbeETtot(const TArrayD &hcalEcorrCf) const
  { return GetProbeEtot(hcalEcorrCf)/cosh(fProbeEta); }

  const GammaJetEventAuxInfo_t& GetAuxInfo() const { return fAInfo; }
  void SetAuxInfo(const GammaJetEventAuxInfo_t &aux) { fAInfo=aux; }

  //
  // calculations
  //

  double CalcDiffEt() const { return (fTagEt - fProbeEt); }

  // for gamma+jet, the photon energy is not changed
  double CalcDiffEt(const TArrayD &hcalEcorrCf) const
  { return (fTagEt - GetProbeETtot(hcalEcorrCf)); }

  // for jet+jet, the tag and probe energies are changed
  double CalcDiffEt_jj(const TArrayD &hcalEcorrCf) const
  { return (GetTagETtot(hcalEcorrCf) - GetProbeETtot(hcalEcorrCf)); }

  //
  // operators
  //

  void MyPrint(int detail=0, std::ostream& out=std::cout) const;

  friend
  std::ostream& operator<<(std::ostream& out, const GammaJetEvent_t &e) {
    e.MyPrint(0, out);
    return out;
  }

 protected:
  // event weight
  Double_t fWeight;

  // photon info
  Double_t fTagEta, fTagPhi, fTagEcalE;
  std::map<Int_t,Double_t> fTagHcalE;
  Double_t fTagE,fTagEt; // calculated from fTagEta,f fTagEcalE, fTagHcalE

  // jet info
  Double_t fProbeEta, fProbePhi, fProbeEcalE;
  std::map<Int_t,Double_t> fProbeHcalE;
  Double_t fProbeE, fProbeEt;// calculated from fProbeEta, fProbeEcalE, fProbeHcalE

  // auxiliary info
  GammaJetEventAuxInfo_t fAInfo;

  ClassDef(GammaJetEvent_t, 1);
};

// -----------------------------------------------------------

//
// A wrapper class to contain data to minimize
//
struct GammaJetCuts_t {
  protected:
  Double_t fEtDiffMin, fEtDiffMax;
  public:
  GammaJetCuts_t();
  GammaJetCuts_t(const GammaJetCuts_t &cuts);

  void SetEtDiff(double min, double max)
  { fEtDiffMin=min; fEtDiffMax=max; }

  int passCuts(const GammaJetEvent_t &e) const;
};

// -----------------------------------------------------------

//
// A wrapper class to contain data to minimize
//

class GammaJetFitter_t : public TObject {
 public:
  GammaJetFitter_t();
  ~GammaJetFitter_t();

  void push_back(const GammaJetEvent_t &e)
  { fData.push_back(new GammaJetEvent_t(e)); }

  const GammaJetEvent_t& GetAt(Int_t idx) const { return *fData.at(idx); }
  const std::vector<GammaJetEvent_t*>& GetData() const { return fData; }
  Int_t GetSize() const { return Int_t(fData.size()); }
  Int_t GetParameterCount() const { return NUMTOWERS; }

  unsigned int size() const { return fData.size(); }
  template<class idx_t>
    const GammaJetEvent_t* at(idx_t idx) const { return fData.at(idx); }


  // ---------------------------------

  Double_t CalcFitValue(const GammaJetEvent_t *e) const {
    Double_t x=0.;
    switch(fFittingProcedure) {
    case 0: // gamma+jet transverse energy diff
    case 1: // jet+jet transverse energy diff
      x= e->CalcDiffEt();
      break;
    default:
      std::cout << "CalcFitValue is not ready for fittingProcedure="
		<< fFittingProcedure << std::endl;
    }
    return x;
  }

  // ---------------------------------

  Double_t EstimateResolution(int setValue=1);

  // get the total distance between the tag and probe
  Double_t GetChi2(const TArrayD& hcalCorrCf) const;
  Double_t GetChi2_jj(const TArrayD& hcalCorrCf) const;
  // depending on fFittingProcedure, get fit value
  Double_t GetFitValue(const GammaJetEvent_t *e,
		       const TArrayD &hcalCorrCf) const;
  Double_t GetFitValue(const TArrayD &hcalCorrCf) const;

  TH1D* doFit(const char *histoName, const char *histoTitle,
	      std::vector<Double_t> &hcalCorrCf,
	      const std::vector<Int_t> *fixTowers=NULL
	      );

  void PrintFitValueChange(const std::vector<Double_t> &finalHcalCorrCf,
		   const std::vector<Double_t> *initialHCalCorrCf=NULL) const;

  // fitting parameters
  inline void SetPrintLevel(Int_t p) { fPrintLevel=p; }
  inline void SetParStep(Double_t p) { fParStep=p; }
  inline void SetParMin(Double_t p) { fParMin=p; }
  inline void SetParMax(Double_t p) { fParMax=p; }
  inline void SetEcalRes(Double_t p=0.07) { fEcalRes=p; }
  inline void SetHcalRes(Double_t p=1.15) { fHcalRes=p; }
  inline void SetHfRes(Double_t p=1.35) { fHfRes=p; }
  void SetErrorDef(Double_t p=0.5) { fErrorDef=p; }
  
  inline void SetFittingProcedure(Int_t p) { fFittingProcedure=p; }

 protected:
  // actual data
  std::vector<GammaJetEvent_t*> fData;

  // guiding parameters
  Int_t fFittingProcedure;
  Double_t fResolution;

  // auxiliary parameters
  Int_t fPrintLevel;
  Double_t fErrorDef;
  Double_t fParStep, fParMin, fParMax;
  Double_t fEcalRes, fHcalRes, fHfRes;

 public:
  ClassDef(GammaJetFitter_t, 1)
};

// -----------------------------------------------------------
// -----------------------------------------------------------

#endif
