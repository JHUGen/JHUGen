#ifndef ROOSPINTWO_7DCOMPLEX_PPHVV
#define ROOSPINTWO_7DCOMPLEX_PPHVV

#include "RooSpinTwo.h"


class RooSpinTwo_7DComplex_ppHVV : public RooSpinTwo {
public:

  RooSpinTwo_7DComplex_ppHVV(){}
  RooSpinTwo_7DComplex_ppHVV(
    const char *name, const char *title,
    modelMeasurables _measurables,
    modelParameters _parameters,
    modelCouplings _couplings,
    RooSpin::VdecayType _Vdecay1=RooSpin::kVdecayType_Zll, RooSpin::VdecayType _Vdecay2=RooSpin::kVdecayType_Zll
    );
  RooSpinTwo_7DComplex_ppHVV(const RooSpinTwo_7DComplex_ppHVV& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooSpinTwo_7DComplex_ppHVV(*this, newname); }
  inline virtual ~RooSpinTwo_7DComplex_ppHVV(){}

  Double_t evaluate() const;
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

  void setZZ4fOrdering(Bool_t flag=true);

protected:


  Double_t evaluateHSFactor(Int_t di, Int_t dj, Int_t code) const;
  Double_t evaluateH1Factor(Int_t i1, Int_t j1, Int_t helicity, Int_t code) const;
  Double_t evaluateH2Factor(Int_t i2, Int_t j2, Int_t helicity, Int_t code) const;
  Double_t evaluatePhi1PhiFactor(Int_t i1, Int_t i2, Int_t j1, Int_t j2, Int_t code, Double_t extraPhase1, Double_t extraPhase2) const;
  void evaluatePolarizationTerms(std::vector<Double_t>& Axxyyterm, const Int_t code, bool isGammaV1=false, bool isGammaV2=false) const;

  Bool_t ZZ4fOrdering;

};

#endif
