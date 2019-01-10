#ifndef ROOSPINZERO_7DCOMPLEX_WITHACCEP_HVV
#define ROOSPINZERO_7DCOMPLEX_WITHACCEP_HVV

#include "RooSpinZero.h"


class RooSpinZero_7DComplex_withAccep_HVV : public RooSpinZero {

public:

  struct accepParameters{
    RooRealVar* aPhi;
    RooRealVar* bPhi;
    RooRealVar* cPhi;
    RooRealVar* dPhi;
    RooRealVar* ePhi;
    RooRealVar* aPhi1;
    RooRealVar* bPhi1;
    RooRealVar* cPhi1;
    RooRealVar* dPhi1;
    RooRealVar* ePhi1;
    RooRealVar* aH1;
    RooRealVar* bH1;
    RooRealVar* cH1;
    RooRealVar* dH1;
    RooRealVar* eH1;
    RooRealVar* aH2;
    RooRealVar* bH2;
    RooRealVar* cH2;
    RooRealVar* dH2;
    RooRealVar* eH2;
    RooRealVar* aHs;
    RooRealVar* bHs;
    RooRealVar* cHs;
    RooRealVar* dHs;
    RooRealVar* eHs;
    RooRealVar* aM1;
    RooRealVar* bM1;
    RooRealVar* cM1;
    RooRealVar* dM1;
    RooRealVar* aM2;
    RooRealVar* bM2;
    RooRealVar* cM2;
    RooRealVar* dM2;
  };

  RooSpinZero_7DComplex_withAccep_HVV();
  RooSpinZero_7DComplex_withAccep_HVV(
    const char *name, const char *title,
    modelMeasurables const& _measurables,
    modelParameters const& _parameters,
    modelCouplings const& _couplings,
    accepParameters const& _accepParams,
    RooSpin::VdecayType _Vdecay1=RooSpin::kVdecayType_Zll, RooSpin::VdecayType _Vdecay2=RooSpin::kVdecayType_Zll,
    TVar::VerbosityLevel verbosity_=TVar::ERROR
  );
  RooSpinZero_7DComplex_withAccep_HVV(const RooSpinZero_7DComplex_withAccep_HVV& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooSpinZero_7DComplex_withAccep_HVV(*this, newname); }
  inline virtual ~RooSpinZero_7DComplex_withAccep_HVV(){}

  Double_t evaluate() const;
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

  void setZZ4fOrdering(Bool_t flag=true);

protected:

  // acceptance parameters
  RooRealProxy aPhi;
  RooRealProxy bPhi;
  RooRealProxy cPhi;
  RooRealProxy dPhi;
  RooRealProxy ePhi;
  RooRealProxy aPhi1;
  RooRealProxy bPhi1;
  RooRealProxy cPhi1;
  RooRealProxy dPhi1;
  RooRealProxy ePhi1;
  RooRealProxy aH1;
  RooRealProxy bH1;
  RooRealProxy cH1;
  RooRealProxy dH1;
  RooRealProxy eH1;
  RooRealProxy aH2;
  RooRealProxy bH2;
  RooRealProxy cH2;
  RooRealProxy dH2;
  RooRealProxy eH2;
  RooRealProxy aHs;
  RooRealProxy bHs;
  RooRealProxy cHs;
  RooRealProxy dHs;
  RooRealProxy eHs;

  RooRealProxy aM1;
  RooRealProxy bM1;
  RooRealProxy cM1;
  RooRealProxy dM1;

  RooRealProxy aM2;
  RooRealProxy bM2;
  RooRealProxy cM2;
  RooRealProxy dM2;

  Bool_t ZZ4fOrdering;

  void evaluatePolarizationTerms(
    Double_t& A00term, Double_t& Appterm, Double_t& Ammterm,
    Double_t& A00ppterm, Double_t& A00mmterm, Double_t& Appmmterm,
    const Int_t code,
    int VGammaVpmode1=0, int VGammaVpmode2=0
  ) const;

};

#endif
